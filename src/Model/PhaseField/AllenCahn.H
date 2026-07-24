//
// This header defines a standard Allen-Cahn evolution operator for an order
// parameter eta with free energy
//
// .. math::
//
//    F[\eta] = \int_\Omega \left[\lambda\eta^2(1-\eta)^2
//              + \frac{\kappa}{2}|\nabla\eta|^2
//              + \Delta f\eta^2(3-2\eta)\right]d\mathbf{x}.
//
// The resulting nonconservative source is
//
// .. math::
//
//    \dot{\eta} = -L\left[\lambda(2\eta-6\eta^2+4\eta^3)
//                    - \kappa\Delta\eta
//                    + 6\Delta f\eta(1-\eta)\right].
//

#ifndef MODEL_PHASEFIELD_ALLENCAHN_H
#define MODEL_PHASEFIELD_ALLENCAHN_H

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
class AllenCahn
{
public:
    static constexpr const char* name = "allencahn";

    static void Parse(AllenCahn& value, IO::ParmParse& pp)
    {
        const bool configured =
            pp.contains("mobility") ||
            pp.contains("lambda") || pp.contains("kappa") ||
            pp.contains("sigma") || pp.contains("epsilon") ||
            pp.contains("driving_force");
        if (!configured) return;

        pp.query_required("mobility", value.mobility,
                        1.0 / Unit::Pressure() / Unit::Time());

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
                "Allen-Cahn requires either (lambda,kappa) or (sigma,epsilon)");
        }

        if (value.mobility < 0.0)
            Util::Exception(INFO, "mobility must be nonnegative");
        if (value.lambda < 0.0 || value.kappa < 0.0)
            Util::Exception(INFO, "Allen-Cahn energy coefficients must be nonnegative");

        pp.query_default("driving_force", value.driving_force,
                        "0.0", Unit::Pressure());
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar operator()(const amrex::Array4<const Set::Scalar>& field,
                            const int component, const Set::Scalar scale,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM],
                            const Set::Scalar /*temperature*/,
                            const Set::Scalar /*pressure*/) const
    {
        const Set::Scalar value = scale * field(i,j,k,component);
        const Set::Scalar value2 = value * value;
        const Set::Scalar chemical_potential =
            lambda * (2.0 * value - 6.0 * value2 + 4.0 * value2 * value) -
            kappa * scale * Numeric::Laplacian(field, i, j, k, component, DX) +
            6.0 * driving_force * value * (1.0 - value);
        return -mobility * chemical_potential;
    }

    Set::Scalar StabilityRate(const Set::Scalar dx,
                            const Set::Scalar /*temperature*/,
                            const Set::Scalar /*pressure*/) const
    {
        return mobility *
            (2.0 * lambda + 6.0 * Util::Abs(driving_force) +
            4.0 * AMREX_SPACEDIM * kappa / (dx * dx));
    }

private:
    Set::Scalar mobility = 0.0;
    Set::Scalar lambda = 0.0;
    Set::Scalar kappa = 0.0;
    Set::Scalar driving_force = 0.0;
};
}
}

#endif
