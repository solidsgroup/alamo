//
// This header defines an Allen-Cahn phase field with the same interface and
// driving-force structure as :code:`Model::PhaseField::AllenCahn`, but with a
// mobility that reproduces the eta-evolution behavior of the legacy
// :code:`Model::Propellant::FullFeedback` surrogate combustion model
// (``get_L``): an Arrhenius temperature dependence, with an optional linear
// pressure dependence, and a hard cutoff below a bound temperature.
//
// .. math::
//
//    L(T,P) = L_0 \exp(-E/T) \times
//             \begin{cases} P/P_{ref} & \text{pressure\_dependence} \\
//                           1 & \text{otherwise} \end{cases}
//
// with :math:`L(T,P) = 0` when :math:`T \le T_{bound}`. This mobility
// replaces the constant mobility of :code:`AllenCahn` in the same
// chemical-potential expression, so the interface (:code:`lambda`/:code:`kappa`
// or :code:`sigma`/:code:`epsilon`) and :code:`driving_force` keys behave
// identically to :code:`AllenCahn`.
//
// This model is intended to be used one instance per condensed material
// (e.g. one for AP, one for HTPB) via separate ``phase_change`` mechanisms;
// a cell that is partially occupied by multiple rigid solids is already
// handled upstream by :code:`Model::Mechanism::PhaseChange`, which scales the
// chemical potential evaluated on the shared solid front by this material's
// fraction of that front before returning the source term.
//

#ifndef MODEL_PHASEFIELD_FULLFEEDBACK_H
#define MODEL_PHASEFIELD_FULLFEEDBACK_H

#include "AMReX_Array4.H"
#include "AMReX_Extension.H"
#include "IO/ParmParse.H"
#include "Numeric/Stencil.H"
#include "Set/Set.H"
#include "Unit/Unit.H"
#include "Util/Util.H"

namespace Model
{
namespace PhaseField
{
class FullFeedback
{
public:
    static constexpr const char* name = "fullfeedback";

    static void Parse(FullFeedback& value, IO::ParmParse& pp)
    {
        const bool configured =
            pp.contains("pre_exponential") ||
            pp.contains("lambda") || pp.contains("kappa") ||
            pp.contains("sigma") || pp.contains("epsilon") ||
            pp.contains("driving_force");
        if (!configured) return;

        // Pre-exponential (base) mobility
        pp.query_required("pre_exponential", value.pre_exponential,
                        1.0 / Unit::Pressure() / Unit::Time());

        // Activation temperature (activation energy already divided by the
        // universal gas constant, so this carries units of temperature)
        pp.query_default("activation_temperature", value.activation_temperature,
                        "0.0_K", Unit::Temperature());

        // Whether the mobility additionally scales linearly with pressure
        pp.query_default("pressure_dependence", value.pressure_dependence, false);

        // Reference pressure used to nondimensionalize pressure when
        // pressure_dependence is enabled
        pp.query_default("reference_pressure", value.reference_pressure,
                        "1.0_MPa", Unit::Pressure());

        // Temperature at/below which the mobility is frozen to zero
        pp.query_default("temperature_bound", value.temperature_bound,
                        "0.0_K", Unit::Temperature());

        std::pair<std::string, Set::Scalar> coefficients[2];
        pp.query_exactly<2>({"lambda", "kappa", "sigma", "epsilon"}, coefficients,
                            {Unit::Pressure(), Unit::Energy() / Unit::Length(),
                            Unit::Energy() / Unit::Area(), Unit::Length()});

        if (coefficients[0].first == "lambda" && coefficients[1].first == "kappa")
        {
            value.lambda = coefficients[0].second;
            value.kappa = coefficients[1].second;
        }
        else if (coefficients[0].first == "sigma" && coefficients[1].first == "epsilon")
        {
            const Set::Scalar sigma = coefficients[0].second;
            const Set::Scalar epsilon = coefficients[1].second;
            value.lambda = 4.0 * sigma / epsilon;
            value.kappa = 0.5 * sigma * epsilon;
        }
        else
        {
            Util::Exception(INFO,
                "Full-feedback phase field requires either (lambda,kappa) or (sigma,epsilon)");
        }

        if (value.pre_exponential < 0.0)
            Util::Exception(INFO, "pre_exponential must be nonnegative");
        if (value.activation_temperature < 0.0)
            Util::Exception(INFO, "activation_temperature must be nonnegative");
        if (value.reference_pressure <= 0.0)
            Util::Exception(INFO, "reference_pressure must be positive");
        if (value.temperature_bound < 0.0)
            Util::Exception(INFO, "temperature_bound must be nonnegative");
        if (value.lambda < 0.0 || value.kappa < 0.0)
            Util::Exception(INFO, "Full-feedback energy coefficients must be nonnegative");

        pp.query_default("driving_force", value.driving_force,
                        "0.0", Unit::Pressure());
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar Mobility(const Set::Scalar temperature,
                        const Set::Scalar pressure) const
    {
        if (temperature <= temperature_bound) return 0.0;
        Set::Scalar mobility = pre_exponential *
            exp(-activation_temperature / temperature);
        if (pressure_dependence) mobility *= pressure / reference_pressure;
        return mobility;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar operator()(const amrex::Array4<const Set::Scalar>& field,
                            const int component, const Set::Scalar scale,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM],
                            const Set::Scalar temperature,
                            const Set::Scalar pressure) const
    {
        const Set::Scalar value = scale * field(i,j,k,component);
        const Set::Scalar value2 = value * value;
        const Set::Scalar chemical_potential =
            lambda * (2.0 * value - 6.0 * value2 + 4.0 * value2 * value) -
            kappa * scale * Numeric::Laplacian(field, i, j, k, component, DX) +
            6.0 * driving_force * value * (1.0 - value);
        return -Mobility(temperature, pressure) * chemical_potential;
    }

    Set::Scalar StabilityRate(const Set::Scalar dx,
                            const Set::Scalar temperature,
                            const Set::Scalar pressure) const
    {
        return Mobility(temperature, pressure) *
            (2.0 * lambda + 6.0 * Util::Abs(driving_force) +
            4.0 * AMREX_SPACEDIM * kappa / (dx * dx));
    }

private:
    Set::Scalar pre_exponential = 0.0;
    Set::Scalar activation_temperature = 0.0;
    bool pressure_dependence = false;
    Set::Scalar reference_pressure = NAN;
    Set::Scalar temperature_bound = 0.0;
    Set::Scalar lambda = 0.0;
    Set::Scalar kappa = 0.0;
    Set::Scalar driving_force = 0.0;
};
}
}

#endif
