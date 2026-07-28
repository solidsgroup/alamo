//
// Conservative reaction between a condensed primary reactant and an
// arbitrary pool of secondary reactants at their diffuse interface.  The
// primary phase field supplies interfacial area, while local secondary
// availability limits the prescribed primary mass flux.  Products may be
// any mixture of gas and condensed species.
//

#ifndef MODEL_MECHANISM_INTERPHASEREACTION_H
#define MODEL_MECHANISM_INTERPHASEREACTION_H

#include <array>
#include <cmath>
#include <vector>

#include "Model/Chemistry/Chemistry.H"
#include "Model/Mechanism/Mechanism.H"
#include "Numeric/Stencil.H"
#include "Unit/Unit.H"

namespace Model
{
namespace Mechanism
{
class InterphaseReaction
{
public:
    static constexpr const char* name = "interphase_reaction";
    static constexpr int max_species = Model::Chemistry::MAX_SPECIES;

    static void Parse(InterphaseReaction& value, IO::ParmParse& pp,
                    const std::vector<std::string>& species_names,
                    const int ngas_species,
                    const std::vector<int>& /*rigid_solid_species*/,
                    const std::vector<Set::Scalar>& reference_density,
                    const std::vector<double>& molecular_weight,
                    const Set::Scalar universal_gas_constant)
    {
        auto find_species = [&](const std::string& name)
        {
            for (int n = 0; n < static_cast<int>(species_names.size()); ++n)
                if (species_names[n] == name) return n;
            Util::Exception(INFO, "Unknown species identifier ", name,
                            " in ", pp.full("type"));
            return -1;
        };

        std::string primary_name;
        std::vector<std::string> secondary_names, product_names;
        std::vector<Set::Scalar> secondary_available_fractions;
        std::vector<Set::Scalar> product_mass_fractions;
        std::vector<Set::Scalar> product_heat_release;
        pp.query_required("reactants.primary", primary_name);
        pp.queryarr_required("reactants.secondary", secondary_names);
        pp.queryarr_required("reactants.secondary.available_mass_fractions",
                            secondary_available_fractions);
        pp.queryarr_required("products", product_names);
        pp.queryarr_required("products.mass_fractions",
                            product_mass_fractions);

        value.nspecies = static_cast<int>(species_names.size());
        value.ngas_species = ngas_species;
        value.primary_reactant = find_species(primary_name);
        if (value.primary_reactant < ngas_species)
            Util::Exception(INFO, pp.full("reactants.primary"),
                " must identify a condensed species whose phase field "
                "defines the reaction interface");
        if (secondary_names.empty() ||
            secondary_names.size() > max_species ||
            secondary_available_fractions.size() != secondary_names.size())
            Util::Exception(INFO,
                pp.full("reactants.secondary.available_mass_fractions"),
                " must contain one value per secondary reactant");
        if (product_names.empty() || product_names.size() > max_species ||
            product_mass_fractions.size() != product_names.size())
            Util::Exception(INFO, pp.full("products.mass_fractions"),
                " must contain one value per product");

        for (int n = 0; n < value.nspecies; ++n)
        {
            if (n < ngas_species)
            {
                if (n >= static_cast<int>(molecular_weight.size()) ||
                    !(molecular_weight[n] > 0.0))
                    Util::Exception(INFO, species_names[n],
                                    " requires a positive molecular weight");
                value.gas_constant[n] =
                    universal_gas_constant / molecular_weight[n];
            }
            else
            {
                if (n >= static_cast<int>(reference_density.size()) ||
                    !(reference_density[n] > 0.0))
                    Util::Exception(INFO, species_names[n],
                                    " requires a positive reference density");
                value.inverse_reference_density[n] =
                    1.0 / reference_density[n];
            }
        }
        value.primary_inverse_reference_density =
            value.inverse_reference_density[value.primary_reactant];

        Set::Scalar available_fraction_sum = 0.0;
        value.nsecondary = static_cast<int>(secondary_names.size());
        for (int n = 0; n < value.nsecondary; ++n)
        {
            value.secondary_reactants[n] = find_species(secondary_names[n]);
            if (value.secondary_reactants[n] == value.primary_reactant)
                Util::Exception(INFO, primary_name,
                    " cannot be both the primary and a secondary reactant");
            for (int m = 0; m < n; ++m)
                if (value.secondary_reactants[n] ==
                    value.secondary_reactants[m])
                    Util::Exception(INFO, "Duplicate secondary reactant ",
                                    secondary_names[n]);
            if (secondary_available_fractions[n] < 0.0 ||
                secondary_available_fractions[n] > 1.0)
                Util::Exception(INFO,
                    "secondary-reactant available mass fractions must lie "
                    "in [0,1]");
            value.secondary_available_fraction[n] =
                secondary_available_fractions[n];
            value.is_secondary[value.secondary_reactants[n]] = 1;
            value.has_gas_secondary = value.has_gas_secondary ||
                value.secondary_reactants[n] < ngas_species;
            value.uses_gas = value.uses_gas ||
                value.secondary_reactants[n] < ngas_species;
            available_fraction_sum += secondary_available_fractions[n];
        }
        if (!(available_fraction_sum > 0.0))
            Util::Exception(INFO,
                "At least one secondary-reactant available mass fraction "
                "must be positive");

        Set::Scalar product_fraction_sum = 0.0;
        value.nproducts = static_cast<int>(product_names.size());
        for (int n = 0; n < value.nproducts; ++n)
        {
            value.products[n] = find_species(product_names[n]);
            for (int m = 0; m < n; ++m)
                if (value.products[n] == value.products[m])
                    Util::Exception(INFO, "Duplicate product ",
                                    product_names[n]);
            if (product_mass_fractions[n] < 0.0)
                Util::Exception(INFO,
                    "product mass fractions must be nonnegative");
            value.product_mass_fraction[n] = product_mass_fractions[n];
            value.uses_gas = value.uses_gas ||
                value.products[n] < ngas_species;
            product_fraction_sum += product_mass_fractions[n];
        }
        if (Util::Abs(product_fraction_sum - 1.0) > 1.0e-12)
            Util::Exception(INFO, pp.full("products.mass_fractions"),
                            " must sum to one");

        if (pp.contains("products.heat_release"))
            pp.queryarr_required("products.heat_release",
                product_heat_release,
                Unit::Energy() / Unit::Mass());
        else
            product_heat_release.assign(value.nproducts, 0.0);
        if (product_heat_release.size() != product_names.size())
            Util::Exception(INFO, pp.full("products.heat_release"),
                            " must contain one value per product");
        for (int n = 0; n < value.nproducts; ++n)
        {
            if (!std::isfinite(product_heat_release[n]))
                Util::Exception(INFO,
                    "product heat-release values must be finite");
            value.product_heat_release[n] = product_heat_release[n];
        }

        pp.query_required("reactants.secondary.mass_ratio",
                        value.secondary_mass_ratio);
        pp.query_required("mass_flux", value.primary_mass_flux,
                        Unit::Mass() / Unit::Area() / Unit::Time());
        pp.query_default("activation_temperature",
                        value.activation_temperature,
                        "0.0_K", Unit::Temperature());
        pp.query_default("temperature_cutoff", value.temperature_cutoff,
                        "0.0_K", Unit::Temperature());
        if (!(value.secondary_mass_ratio > 0.0) ||
            value.primary_mass_flux < 0.0 ||
            value.activation_temperature < 0.0 ||
            value.temperature_cutoff < 0.0)
            Util::Exception(INFO,
                "interphase_reaction requires a positive secondary mass "
                "ratio and nonnegative rate and temperature parameters");
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Apply(const amrex::Array4<Set::Scalar>& rhs,
                const State& state, const int i, const int j, const int k,
                const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        const Set::Scalar primary_rate =
            ComputePrimaryRate(state, i, j, k, DX);
        if (!(primary_rate > 0.0)) return;

        Set::Scalar available_secondary = 0.0;
        for (int n = 0; n < nsecondary; ++n)
            available_secondary += secondary_available_fraction[n] *
                Util::Max(state.component_density(
                    i,j,k,secondary_reactants[n]), 0.0);

        rhs(i,j,k,primary_reactant) -= primary_rate;
        for (int n = 0; n < nsecondary; ++n)
            rhs(i,j,k,secondary_reactants[n]) -=
                secondary_mass_ratio * primary_rate *
                secondary_available_fraction[n] *
                Util::Max(state.component_density(
                    i,j,k,secondary_reactants[n]), 0.0) /
                available_secondary;

        const Set::Scalar product_rate =
            (1.0 + secondary_mass_ratio) * primary_rate;
        for (int n = 0; n < nproducts; ++n)
            rhs(i,j,k,products[n]) +=
                product_mass_fraction[n] * product_rate;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar VolumeSource(const State& state,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        const Set::Scalar primary_rate =
            ComputePrimaryRate(state, i, j, k, DX);
        if (!(primary_rate > 0.0) ||
            (uses_gas && !(state.pressure > 0.0))) return 0.0;

        Set::Scalar available_secondary = 0.0;
        for (int n = 0; n < nsecondary; ++n)
            available_secondary += secondary_available_fraction[n] *
                Util::Max(state.component_density(
                    i,j,k,secondary_reactants[n]), 0.0);

        Set::Scalar volume_rate = -primary_rate *
            primary_inverse_reference_density;
        for (int n = 0; n < nsecondary; ++n)
        {
            const int species = secondary_reactants[n];
            const Set::Scalar specific_volume = species < ngas_species ?
                gas_constant[species] * state.temperature / state.pressure :
                inverse_reference_density[species];
            volume_rate -= secondary_mass_ratio * primary_rate *
                secondary_available_fraction[n] *
                Util::Max(state.component_density(i,j,k,species), 0.0) /
                available_secondary * specific_volume;
        }

        const Set::Scalar product_rate =
            (1.0 + secondary_mass_ratio) * primary_rate;
        for (int n = 0; n < nproducts; ++n)
        {
            const int species = products[n];
            const Set::Scalar specific_volume = species < ngas_species ?
                gas_constant[species] * state.temperature / state.pressure :
                inverse_reference_density[species];
            volume_rate += product_mass_fraction[n] * product_rate *
                specific_volume;
        }
        return volume_rate;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar HeatSource(const State& state,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        const Set::Scalar product_rate =
            (1.0 + secondary_mass_ratio) *
            ComputePrimaryRate(state, i, j, k, DX);
        Set::Scalar heat_source = 0.0;
        for (int n = 0; n < nproducts; ++n)
            heat_source += product_mass_fraction[n] * product_rate *
                product_heat_release[n];
        return heat_source;
    }

    Set::Scalar StabilityRate(const Set::Scalar dx,
                            const Set::Scalar temperature) const
    {
        if (!(temperature > temperature_cutoff) || !(dx > 0.0)) return 0.0;
        const Set::Scalar temperature_factor =
            activation_temperature > 0.0 ?
                std::exp(-activation_temperature / temperature) : 1.0;
        return primary_mass_flux * primary_inverse_reference_density /
            dx * temperature_factor;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar LocalRate(
        const State& /*state*/, const int /*i*/, const int /*j*/,
        const int /*k*/) const
    {
        return 0.0;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar GradientCoefficient(
        const State& /*state*/, const int /*i*/, const int /*j*/,
        const int /*k*/) const
    {
        return 0.0;
    }

    int RigidComponent() const { return -1; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar ApplyImplicitChange(
        const amrex::Array4<Set::Scalar>& /*component_density*/,
        const State& /*state*/, const Set::Scalar /*eta_change*/,
        Set::Scalar& /*integrated_heat*/,
        const int /*i*/, const int /*j*/, const int /*k*/) const
    {
        return 0.0;
    }

private:
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar ComputePrimaryRate(
        const State& state, const int i, const int j, const int k,
        const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        if (!(state.temperature > temperature_cutoff) ||
            !(primary_mass_flux > 0.0)) return 0.0;

        const Set::Scalar primary_eta = Util::Clamp(
            state.component_density(i,j,k,primary_reactant) *
                primary_inverse_reference_density, 0.0, 1.0);
        if (!(primary_eta > 0.0)) return 0.0;

        Set::Scalar condensed_eta = 0.0;
        Set::Scalar contacting_eta = 0.0;
        for (int n = ngas_species; n < nspecies; ++n)
        {
            const Set::Scalar eta = Util::Max(
                state.component_density(i,j,k,n), 0.0) *
                inverse_reference_density[n];
            condensed_eta += eta;
            if (is_secondary[n]) contacting_eta += eta;
        }
        const Set::Scalar gas_eta =
            Util::Clamp(1.0 - condensed_eta, 0.0, 1.0);
        if (has_gas_secondary) contacting_eta += gas_eta;
        const Set::Scalar other_eta =
            Util::Max(gas_eta + condensed_eta - primary_eta, 0.0);
        if (!(contacting_eta > 0.0) || !(other_eta > 0.0)) return 0.0;
        const Set::Scalar exposure =
            Util::Clamp(contacting_eta / other_eta, 0.0, 1.0);

        Set::Scalar available_secondary = 0.0;
        for (int n = 0; n < nsecondary; ++n)
            available_secondary += secondary_available_fraction[n] *
                Util::Max(state.component_density(
                    i,j,k,secondary_reactants[n]), 0.0);
        if (!(available_secondary > 0.0)) return 0.0;

        Set::Scalar primary_rate = primary_mass_flux *
            Numeric::Gradient(state.component_density, i, j, k,
                primary_reactant, DX).norm() *
            primary_inverse_reference_density * exposure;
        if (activation_temperature > 0.0)
            primary_rate *=
                std::exp(-activation_temperature / state.temperature);
        if (state.dt > 0.0)
        {
            primary_rate = Util::Min(primary_rate,
                Util::Max(state.component_density(
                    i,j,k,primary_reactant), 0.0) / state.dt);
            primary_rate = Util::Min(primary_rate,
                available_secondary /
                    (secondary_mass_ratio * state.dt));
        }
        return Util::Max(primary_rate, 0.0);
    }

    int nspecies = 0;
    int ngas_species = 0;
    int primary_reactant = -1;
    int nsecondary = 0;
    int nproducts = 0;
    std::array<int, max_species> secondary_reactants{};
    std::array<int, max_species> products{};
    std::array<int, max_species> is_secondary{};
    std::array<Set::Scalar, max_species> secondary_available_fraction{};
    std::array<Set::Scalar, max_species> product_mass_fraction{};
    std::array<Set::Scalar, max_species> product_heat_release{};
    std::array<Set::Scalar, max_species> inverse_reference_density{};
    std::array<Set::Scalar, max_species> gas_constant{};
    Set::Scalar primary_inverse_reference_density = 0.0;
    Set::Scalar secondary_mass_ratio = 0.0;
    Set::Scalar primary_mass_flux = 0.0;
    Set::Scalar activation_temperature = 0.0;
    Set::Scalar temperature_cutoff = 0.0;
    bool has_gas_secondary = false;
    bool uses_gas = false;
};
}
}

#endif
