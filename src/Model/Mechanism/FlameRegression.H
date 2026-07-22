//
// ``Model::Mechanism::FlameRegression`` applies Flame's common propellant
// phase-field evolution to a multicomponent rigid solid. Solid mass removed by
// regression is transferred one-to-one to corresponding gas species. The
// paired mass sources sum to zero pointwise. The gas velocity source uses
// ``m_dot (u0-u_g)``, with ``u0`` obtained from the solid/gas mass-flux jump,
// so the nonconservative velocity equation retains the conservative momentum
// balance used by Flame.
//

#ifndef MODEL_MECHANISM_FLAMEREGRESSION_H
#define MODEL_MECHANISM_FLAMEREGRESSION_H

#include <array>
#include <cmath>
#include <vector>

#include "Model/Mechanism/Mechanism.H"
#include "Numeric/Stencil.H"
#include "Set/Set.H"

namespace Model
{
namespace Mechanism
{
class FlameRegression
{
public:
    static constexpr const char* name = "flame_regression";
    static constexpr int max_materials = 16;

    static void Parse(FlameRegression& value, IO::ParmParse& pp,
                    const std::vector<std::string>& species_names,
                    const int ngas_species,
                    const std::vector<int>& rigid_solid_species,
                    const std::vector<Set::Scalar>& reference_density,
                    const std::vector<double>& molecular_weight,
                    const Set::Scalar universal_gas_constant)
    {
        std::vector<std::string> inputs;
        std::vector<std::string> outputs;
        std::vector<Set::Scalar> pre_exponential;
        std::vector<Set::Scalar> activation_temperature;
        std::vector<Set::Scalar> ejection_temperature;
        std::vector<double> pressure_exponent;
        pp.queryarr_required("in", inputs);
        pp.queryarr_required("out", outputs);
        pp.queryarr_required("pre_exponential", pre_exponential, Unit::Velocity());
        pp.queryarr_required(
            "activation_temperature", activation_temperature, Unit::Temperature());
        pp.queryarr_default(
            "pressure_exponent", pressure_exponent, std::vector<double>(inputs.size(), 0.0));
        value.use_ejection_temperature = pp.contains("ejection_temperature");
        if (value.use_ejection_temperature)
            pp.queryarr_required(
                "ejection_temperature", ejection_temperature,
                Unit::Temperature());
        value.use_prescribed_speed = pp.contains("prescribed_speed");
        if (value.use_prescribed_speed)
            pp.query_required(
                "prescribed_speed", value.prescribed_speed, Unit::Velocity());

        if (inputs.empty() || inputs.size() > max_materials)
            Util::Exception(INFO, pp.full("in"), " must contain between 1 and ",
                            max_materials, " rigid solid species");
        if (outputs.size() != inputs.size() ||
            pre_exponential.size() != inputs.size() ||
            activation_temperature.size() != inputs.size() ||
            pressure_exponent.size() != inputs.size() ||
            (value.use_ejection_temperature &&
             ejection_temperature.size() != inputs.size()))
            Util::Exception(INFO,
                            "Flame regression requires one output and one set of rate "
                            "parameters per input species");
        if (inputs.size() != rigid_solid_species.size())
            Util::Exception(INFO,
                            "Flame regression must include every configured rigid solid species");

        value.nmaterials = static_cast<int>(inputs.size());
        for (int n = 0; n < value.nmaterials; ++n)
        {
            value.inputs[n] = FindSpecies(inputs[n], species_names, pp.full("in"));
            value.outputs[n] = FindSpecies(outputs[n], species_names, pp.full("out"));

            bool input_is_rigid = false;
            for (const int species : rigid_solid_species)
                if (value.inputs[n] == species) input_is_rigid = true;
            if (!input_is_rigid)
                Util::Exception(INFO, inputs[n], " must identify a rigid solid species");
            if (value.outputs[n] >= ngas_species ||
                value.outputs[n] >= static_cast<int>(molecular_weight.size()))
                Util::Exception(INFO, outputs[n], " must identify a gas species");
            for (int m = 0; m < n; ++m)
            {
                if (value.inputs[n] == value.inputs[m])
                    Util::Exception(INFO, "Duplicate Flame regression input ", inputs[n]);
                if (value.outputs[n] == value.outputs[m])
                    Util::Exception(INFO, "Duplicate Flame regression output ", outputs[n]);
            }

            value.reference_density[n] = reference_density[value.inputs[n]];
            if (!(value.reference_density[n] > 0.0))
                Util::Exception(INFO, "No positive reference density is defined for ",
                                inputs[n]);
            if (pre_exponential[n] < 0.0 || activation_temperature[n] < 0.0)
                Util::Exception(INFO, "Flame regression rate parameters must be nonnegative");
            value.pre_exponential[n] = pre_exponential[n];
            value.activation_temperature[n] = activation_temperature[n];
            value.pressure_exponent[n] = pressure_exponent[n];
            if (value.use_ejection_temperature)
            {
                if (!(ejection_temperature[n] > 0.0))
                    Util::Exception(
                        INFO, "Flame regression ejection temperatures must be positive");
                value.ejection_temperature[n] = ejection_temperature[n];
            }
            value.product_gas_constant[n] =
                universal_gas_constant / molecular_weight[value.outputs[n]];
        }
        for (const int species : rigid_solid_species)
        {
            bool found = false;
            for (int n = 0; n < value.nmaterials; ++n)
                if (value.inputs[n] == species) found = true;
            if (!found)
                Util::Exception(INFO,
                                "Flame regression must include every configured rigid solid species");
        }

        pp.query_default("pressure_reference", value.pressure_reference,
                        "1.0_MPa", Unit::Pressure());
        if (pp.contains("stability_pressure"))
            pp.query_required(
                "stability_pressure", value.stability_pressure, Unit::Pressure());
        else
            value.stability_pressure = value.pressure_reference;
        pp.query_default("temperature_cutoff", value.temperature_cutoff,
                        "0.0_K", Unit::Temperature());
        pp.query_required("phase_field.epsilon", value.epsilon, Unit::Length());
        pp.query_required(
            "phase_field.lambda", value.lambda, Unit::Energy() / Unit::Area());
        pp.query_required(
            "phase_field.kappa", value.kappa, Unit::Energy() / Unit::Area());
        pp.query_required("phase_field.w0", value.w0, Unit::Less());
        pp.query_required("phase_field.w1", value.w1, Unit::Less());
        pp.query_required("phase_field.w12", value.w12, Unit::Less());

        if (!(value.pressure_reference > 0.0) || !(value.stability_pressure > 0.0) ||
            value.temperature_cutoff < 0.0 || !(value.epsilon > 0.0) ||
            value.lambda < 0.0 || value.kappa < 0.0 ||
            (value.use_prescribed_speed && value.prescribed_speed < 0.0))
            Util::Exception(INFO, "Invalid Flame regression model parameter");

        value.a2 = -5.0 * value.w1 + 16.0 * value.w12 - 11.0 * value.w0;
        value.a3 = 14.0 * value.w1 - 32.0 * value.w12 + 18.0 * value.w0;
        value.a4 = -8.0 * value.w1 + 16.0 * value.w12 - 8.0 * value.w0;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar VolumeSource(const State& state,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        const Rates rates = ComputeRates(state, i, j, k, DX);
        Set::Scalar source = 0.0;
        for (int n = 0; n < nmaterials; ++n)
        {
            const Set::Scalar product_specific_volume =
                product_gas_constant[n] *
                EjectionTemperature(state, n) / state.pressure;
            source += rates.mass[n] *
                (1.0 / reference_density[n] - product_specific_volume);
        }
        return source;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Vector MomentumSource(const State& state,
                            const Set::Vector& gas_velocity,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        const Rates rates = ComputeRates(state, i, j, k, DX);
        const Set::Vector interface_velocity =
            ComputeInterfaceVelocity(state, rates, i, j, k, DX);
        Set::Scalar mass_rate = 0.0;
        for (int n = 0; n < nmaterials; ++n)
            mass_rate -= rates.mass[n];
        return mass_rate * (interface_velocity - gas_velocity);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar GasMassSource(const State& state,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        const Rates rates = ComputeRates(state, i, j, k, DX);
        Set::Scalar mass_rate = 0.0;
        for (int n = 0; n < nmaterials; ++n)
            mass_rate -= rates.mass[n];
        return mass_rate;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar GasSpeciesSource(const State& state, const int species,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        const Rates rates = ComputeRates(state, i, j, k, DX);
        for (int n = 0; n < nmaterials; ++n)
            if (outputs[n] == species) return -rates.mass[n];
        return 0.0;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar GasInjectionTemperature(
                            const State& state, const int species,
                            const int /*i*/, const int /*j*/, const int /*k*/,
                            const Set::Scalar /*DX*/[AMREX_SPACEDIM]) const
    {
        for (int n = 0; n < nmaterials; ++n)
            if (outputs[n] == species)
                return EjectionTemperature(state, n);
        return state.temperature;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Vector InterfaceVelocity(const State& state,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        const Rates rates = ComputeRates(state, i, j, k, DX);
        return ComputeInterfaceVelocity(state, rates, i, j, k, DX);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar RegressionSpeed(const State& state,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        const Rates rates = ComputeRates(state, i, j, k, DX);
        if (!(rates.eta_dot < 0.0)) return 0.0;
        const Set::Scalar grad_eta_norm =
            Numeric::Gradient(state.rigid_eta, i, j, k, 0, DX).lpNorm<2>();
        return grad_eta_norm > 0.0 ? -rates.eta_dot / grad_eta_norm : 0.0;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Apply(const amrex::Array4<Set::Scalar>& component_density_rhs,
                const State& state,
                const int i, const int j, const int k,
                const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        const Rates rates = ComputeRates(state, i, j, k, DX);
        for (int n = 0; n < nmaterials; ++n)
        {
            component_density_rhs(i,j,k,inputs[n]) += rates.mass[n];
            component_density_rhs(i,j,k,outputs[n]) -= rates.mass[n];
        }
    }

    Set::Scalar StabilityRate(const Set::Scalar dx,
                            const Set::Scalar temperature) const
    {
        if (use_prescribed_speed) return prescribed_speed / dx;
        if (!(temperature > temperature_cutoff)) return 0.0;

        Set::Scalar mobility = 0.0;
        for (int n = 0; n < nmaterials; ++n)
        {
            const Set::Scalar rate = pre_exponential[n] *
                std::pow(stability_pressure / pressure_reference, pressure_exponent[n]) *
                std::exp(-activation_temperature[n] / temperature);
            mobility = Util::Max(mobility, rate);
        }
        const Set::Scalar local_rate = lambda / epsilon *
            (2.0 * Util::Abs(a2) + 6.0 * Util::Abs(a3) +
            12.0 * Util::Abs(a4));
        const Set::Scalar gradient_rate =
            4.0 * AMREX_SPACEDIM * epsilon * kappa / (dx * dx);
        return mobility * (local_rate + gradient_rate);
    }

private:
    struct Rates
    {
        std::array<Set::Scalar, max_materials> mass{};
        std::array<Set::Scalar, max_materials> fraction{};
        Set::Scalar eta_dot = 0.0;
    };

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Vector ComputeInterfaceVelocity(
                            const State& state, const Rates& rates,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        if (!(rates.eta_dot < 0.0)) return Set::Vector::Zero();

        const Set::Vector grad_eta =
            Numeric::Gradient(state.rigid_eta, i, j, k, 0, DX);
        const Set::Scalar grad_eta_norm = grad_eta.lpNorm<2>();
        if (!(grad_eta_norm > 0.0)) return Set::Vector::Zero();

        Set::Scalar solid_density = 0.0;
        for (int n = 0; n < nmaterials; ++n)
            solid_density += rates.fraction[n] * reference_density[n];
        if (!(solid_density > 0.0) ||
            !(state.temperature > 0.0) || !(state.pressure > 0.0))
            return Set::Vector::Zero();

        Set::Scalar product_specific_volume = 0.0;
        for (int n = 0; n < nmaterials; ++n)
        {
            const Set::Scalar product_mass_fraction =
                rates.fraction[n] * reference_density[n] / solid_density;
            product_specific_volume += product_mass_fraction *
                product_gas_constant[n] * EjectionTemperature(state, n) /
                state.pressure;
        }
        if (!(product_specific_volume > 0.0)) return Set::Vector::Zero();

        const Set::Scalar regression_speed = -rates.eta_dot / grad_eta_norm;
        const Set::Scalar ejection_speed =
            regression_speed * solid_density * product_specific_volume;
        return (-grad_eta / grad_eta_norm) * ejection_speed;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar EjectionTemperature(const State& state, const int material) const
    {
        return use_ejection_temperature ?
            ejection_temperature[material] : state.temperature;
    }

    static int FindSpecies(const std::string& name,
                        const std::vector<std::string>& species_names,
                        const std::string& parameter)
    {
        for (int n = 0; n < static_cast<int>(species_names.size()); ++n)
            if (species_names[n] == name) return n;
        Util::Exception(INFO, "Unknown species identifier ", name, " in ", parameter);
        return -1;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Rates ComputeRates(const State& state,
                    const int i, const int j, const int k,
                    const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        Rates rates;
        if (!use_prescribed_speed &&
            (!(state.temperature > temperature_cutoff) || !(state.pressure > 0.0)))
            return rates;

        const Set::Scalar eta = Util::Clamp(state.rigid_eta(i,j,k), 0.0, 1.0);
        Set::Scalar material_eta = 0.0;
        for (int n = 0; n < nmaterials; ++n)
        {
            rates.fraction[n] = Util::Max(
                state.component_density(i,j,k,inputs[n]) / reference_density[n], 0.0);
            material_eta += rates.fraction[n];
        }
        if (!(eta > 0.0) || !(material_eta > 0.0)) return rates;
        for (int n = 0; n < nmaterials; ++n)
            rates.fraction[n] /= material_eta;

        if (use_prescribed_speed)
        {
            const Set::Scalar grad_eta_norm =
                Numeric::Gradient(state.rigid_eta, i, j, k, 0, DX).lpNorm<2>();
            rates.eta_dot = -prescribed_speed * grad_eta_norm;
        }
        else
        {
            Set::Scalar mobility = 0.0;
            for (int n = 0; n < nmaterials; ++n)
                mobility += rates.fraction[n] * pre_exponential[n] *
                    std::pow(state.pressure / pressure_reference, pressure_exponent[n]) *
                    std::exp(-activation_temperature[n] / state.temperature);
            const Set::Scalar eta2 = eta * eta;
            const Set::Scalar dw =
                2.0 * a2 * eta + 3.0 * a3 * eta2 + 4.0 * a4 * eta2 * eta;
            const Set::Scalar laplacian =
                Numeric::Laplacian(state.rigid_eta, i, j, k, 0, DX);
            Set::Scalar chemical_potential =
                (lambda / epsilon) * dw - epsilon * kappa * laplacian;
            chemical_potential = Util::Max(chemical_potential, 0.0);
            rates.eta_dot = -mobility * chemical_potential;
        }

        for (int n = 0; n < nmaterials; ++n)
            rates.mass[n] =
                reference_density[n] * rates.fraction[n] * rates.eta_dot;
        return rates;
    }

    int nmaterials = 0;
    std::array<int, max_materials> inputs{};
    std::array<int, max_materials> outputs{};
    std::array<Set::Scalar, max_materials> reference_density{};
    std::array<Set::Scalar, max_materials> product_gas_constant{};
    std::array<Set::Scalar, max_materials> pre_exponential{};
    std::array<Set::Scalar, max_materials> activation_temperature{};
    std::array<Set::Scalar, max_materials> ejection_temperature{};
    std::array<Set::Scalar, max_materials> pressure_exponent{};
    Set::Scalar pressure_reference = NAN;
    Set::Scalar stability_pressure = NAN;
    Set::Scalar temperature_cutoff = 0.0;
    Set::Scalar epsilon = NAN;
    Set::Scalar lambda = NAN;
    Set::Scalar kappa = NAN;
    Set::Scalar w0 = NAN;
    Set::Scalar w1 = NAN;
    Set::Scalar w12 = NAN;
    Set::Scalar a2 = NAN;
    Set::Scalar a3 = NAN;
    Set::Scalar a4 = NAN;
    Set::Scalar prescribed_speed = 0.0;
    bool use_prescribed_speed = false;
    bool use_ejection_temperature = false;
};
}
}

#endif
