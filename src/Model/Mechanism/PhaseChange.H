//
// ``Model::Mechanism::PhaseChange`` couples a phase-field rate for one
// condensed input species to mass-conserving production of one or more gas
// species. Product mass fractions also determine the ideal-gas dilation sent
// to the low-Mach projection. An optional Arrhenius factor permits thermally
// activated regression without coupling this generic mechanism to a specific
// propellant model.
//

#ifndef MODEL_MECHANISM_PHASECHANGE_H
#define MODEL_MECHANISM_PHASECHANGE_H

#include <array>
#include <cmath>
#include <vector>

#include "Model/Mechanism/Mechanism.H"
#include "Model/PhaseField/AllenCahn.H"
#include "Model/PhaseField/PhaseField.H"
#include "Numeric/Stencil.H"
#include "Set/Set.H"

namespace Model
{
namespace Mechanism
{
class PhaseChange
{
public:
    static constexpr const char* name = "phase_change";
    static constexpr int max_products = 16;

    static void Parse(PhaseChange& value, IO::ParmParse& pp,
                    const std::vector<std::string>& species_names,
                    const int ngas_species,
                    const std::vector<int>& rigid_solid_species,
                    const std::vector<Set::Scalar>& reference_density,
                    const std::vector<double>& molecular_weight,
                    const Set::Scalar universal_gas_constant)
    {
        std::string input;
        std::vector<std::string> outputs;
        std::vector<Set::Scalar> mass_fractions;
        pp.query_required("in", input);
        pp.queryarr_required("out", outputs);

        for (int n = 0; n < static_cast<int>(species_names.size()); ++n)
            if (species_names[n] == input) value.input = n;
        if (value.input < 0)
            Util::Exception(INFO, "Unknown species identifier ", input,
                            " in ", pp.full("in"));
        if (value.input < ngas_species)
            Util::Exception(INFO, pp.full("in"), " must identify a condensed species");
        for (int n = 0; n < static_cast<int>(rigid_solid_species.size()); ++n)
            if (value.input == rigid_solid_species[n])
            {
                value.input_is_rigid = true;
                value.rigid_component = n;
            }
        value.reference_density = reference_density[value.input];
        if (!(value.reference_density > 0.0))
            Util::Exception(INFO, "No positive reference density is defined for ", input);

        if (outputs.empty() || outputs.size() > max_products)
            Util::Exception(INFO, pp.full("out"), " must contain between 1 and ",
                            max_products, " gas species");
        if (pp.contains("out.mass_fractions"))
            pp.queryarr_required("out.mass_fractions", mass_fractions);
        else if (outputs.size() == 1)
            mass_fractions.push_back(1.0);
        else
            Util::Exception(INFO, pp.full("out.mass_fractions"),
                            " is required for multiple products");
        if (mass_fractions.size() != outputs.size())
            Util::Exception(INFO, pp.full("out.mass_fractions"),
                            " must contain one value per output species");

        value.nproducts = static_cast<int>(outputs.size());
        Set::Scalar sum = 0.0;
        for (int n = 0; n < value.nproducts; ++n)
        {
            value.products[n] = -1;
            for (int m = 0; m < static_cast<int>(species_names.size()); ++m)
                if (species_names[m] == outputs[n]) value.products[n] = m;
            if (value.products[n] < 0)
                Util::Exception(INFO, "Unknown species identifier ", outputs[n],
                                " in ", pp.full("out"));
            if (value.products[n] >= ngas_species ||
                value.products[n] >= static_cast<int>(molecular_weight.size()))
                Util::Exception(INFO, outputs[n], " must identify a gas species");
            if (mass_fractions[n] < 0.0)
                Util::Exception(INFO, pp.full("out.mass_fractions"),
                                " must be nonnegative");
            value.mass_fractions[n] = mass_fractions[n];
            value.product_gas_constant +=
                universal_gas_constant * mass_fractions[n] /
                molecular_weight[value.products[n]];
            sum += mass_fractions[n];
        }
        if (Util::Abs(sum - 1.0) > 1.0e-12)
            Util::Exception(INFO, pp.full("out.mass_fractions"), " must sum to one");

        pp.query_default("rate_multiplier", value.rate_multiplier, 1.0);
        if (value.rate_multiplier < 0.0)
            Util::Exception(INFO, "rate_multiplier must be nonnegative");
        pp.query_default("activation_temperature", value.activation_temperature,
                        "0.0_K", Unit::Temperature());
        pp.query_default("temperature_cutoff", value.temperature_cutoff,
                        "0.0_K", Unit::Temperature());
        if (value.activation_temperature < 0.0 || value.temperature_cutoff < 0.0)
            Util::Exception(INFO, "phase-change temperatures must be nonnegative");
        value.use_prescribed_speed = pp.contains("prescribed_speed");
        if (value.use_prescribed_speed)
        {
            pp.query_required("prescribed_speed", value.prescribed_speed,
                            Unit::Velocity());
            if (value.prescribed_speed < 0.0)
                Util::Exception(INFO,
                    "phase-change prescribed_speed must be nonnegative");
        }
        pp.select<Model::PhaseField::AllenCahn>("phase_field", value.model);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar VolumeSource(const State& state,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        const Set::Scalar source = MassSource(state, i, j, k, DX);
        const Set::Scalar product_specific_volume =
            product_gas_constant * state.temperature / state.pressure;
        return source * (1.0 / reference_density - product_specific_volume);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Apply(const amrex::Array4<Set::Scalar>& component_density_rhs,
                const State& state,
                const int i, const int j, const int k,
                const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        const Set::Scalar source = MassSource(state, i, j, k, DX);
        component_density_rhs(i,j,k,input) += source;
        for (int n = 0; n < nproducts; ++n)
            component_density_rhs(i,j,k,products[n]) -= mass_fractions[n] * source;
    }

    Set::Scalar StabilityRate(const Set::Scalar dx,
                            const Set::Scalar temperature) const
    {
        if (!(temperature > temperature_cutoff)) return 0.0;
        if (use_prescribed_speed) return prescribed_speed / dx;
        const Set::Scalar temperature_factor = activation_temperature > 0.0 ?
            std::exp(-activation_temperature / temperature) : 1.0;
        return rate_multiplier * temperature_factor *
            (input_is_rigid ? model.LocalStabilityRate() :
                              model.StabilityRate(dx));
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar LocalRate(
        const State& state, const int i, const int j, const int k) const
    {
        if (!input_is_rigid || use_prescribed_speed ||
            !(state.temperature > temperature_cutoff))
            return 0.0;

        const Set::Scalar eta =
            Util::Clamp(state.rigid_eta(i,j,k), 0.0, 1.0);
        if (!(eta > 0.0)) return 0.0;
        const Set::Scalar input_eta = Util::Clamp(
            state.rigid_species_eta(i,j,k,rigid_component), 0.0, eta);
        const Set::Scalar composition =
            Model::PhaseField::H(input_eta / eta);
        Set::Scalar rate = composition *
            model.Local(state.rigid_eta, 0, 1.0, i, j, k);
        if (activation_temperature > 0.0)
            rate *= std::exp(
                -activation_temperature / state.temperature);
        return rate_multiplier * rate;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar LocalStabilityRate(
        const State& state, const int i, const int j, const int k) const
    {
        if (!input_is_rigid || use_prescribed_speed ||
            !(state.temperature > temperature_cutoff))
            return 0.0;

        const Set::Scalar eta =
            Util::Clamp(state.rigid_eta(i,j,k), 0.0, 1.0);
        if (!(eta > 0.0) || !(eta < 1.0)) return 0.0;
        const Set::Scalar input_eta = Util::Clamp(
            state.rigid_species_eta(i,j,k,rigid_component), 0.0, eta);
        const Set::Scalar composition =
            Model::PhaseField::H(input_eta / eta);
        const Set::Scalar temperature_factor =
            activation_temperature > 0.0 ?
                std::exp(-activation_temperature / state.temperature) : 1.0;
        return composition * rate_multiplier * temperature_factor *
            model.LocalStabilityRate();
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar GradientCoefficient(
        const State& state, const int i, const int j, const int k) const
    {
        if (!input_is_rigid || use_prescribed_speed ||
            !(state.temperature > temperature_cutoff))
            return 0.0;

        const Set::Scalar eta =
            Util::Clamp(state.rigid_eta(i,j,k), 0.0, 1.0);
        if (!(eta > 0.0)) return 0.0;
        const Set::Scalar input_eta = Util::Clamp(
            state.rigid_species_eta(i,j,k,rigid_component), 0.0, eta);
        const Set::Scalar composition =
            Model::PhaseField::H(input_eta / eta);
        const Set::Scalar temperature_factor =
            activation_temperature > 0.0 ?
                std::exp(-activation_temperature / state.temperature) : 1.0;
        return composition * rate_multiplier * temperature_factor *
            model.GradientCoefficient();
    }

    int RigidComponent() const
    {
        return input_is_rigid && !use_prescribed_speed ?
            rigid_component : -1;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar ApplyImplicitChange(
        const amrex::Array4<Set::Scalar>& component_density,
        const State& state, const Set::Scalar eta_change,
        const int i, const int j, const int k) const
    {
        const Set::Scalar mass_change = reference_density * eta_change;
        component_density(i,j,k,input) += mass_change;
        for (int n = 0; n < nproducts; ++n)
            component_density(i,j,k,products[n]) -=
                mass_fractions[n] * mass_change;

        const Set::Scalar product_specific_volume =
            product_gas_constant * state.temperature / state.pressure;
        return mass_change *
            (1.0 / reference_density - product_specific_volume);
    }

private:
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar MassSource(const State& state,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        if (!(state.temperature > temperature_cutoff)) return 0.0;

        Set::Scalar eta_dot = 0.0;
        if (input_is_rigid)
        {
            const Set::Scalar eta =
                Util::Clamp(state.rigid_eta(i,j,k), 0.0, 1.0);
            if (eta > 0.0 && use_prescribed_speed)
            {
                const Set::Scalar input_eta = Util::Clamp(
                    state.rigid_species_eta(i,j,k,rigid_component),
                    0.0, eta);
                eta_dot = -input_eta / eta * prescribed_speed *
                    Numeric::Gradient(
                        state.rigid_eta, i, j, k, 0, DX).norm();
            }
        }
        else
        {
            if (use_prescribed_speed)
                eta_dot = -prescribed_speed / reference_density *
                    Numeric::Gradient(
                        state.component_density, i, j, k, input, DX).norm();
            else
                eta_dot = model(state.component_density, input,
                                1.0 / reference_density, i, j, k, DX);
        }
        if (activation_temperature > 0.0)
            eta_dot *= std::exp(-activation_temperature / state.temperature);
        return rate_multiplier * reference_density * eta_dot;
    }

    Model::PhaseField::PhaseField<Model::PhaseField::AllenCahn> model;
    int input = -1;
    int rigid_component = -1;
    int nproducts = 0;
    std::array<int, max_products> products{};
    std::array<Set::Scalar, max_products> mass_fractions{};
    Set::Scalar reference_density = NAN;
    Set::Scalar product_gas_constant = 0.0;
    Set::Scalar rate_multiplier = 1.0;
    Set::Scalar activation_temperature = 0.0;
    Set::Scalar temperature_cutoff = 0.0;
    Set::Scalar prescribed_speed = 0.0;
    bool input_is_rigid = false;
    bool use_prescribed_speed = false;
};
}
}

#endif
