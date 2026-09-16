//
// Prescribed unresolved heat flux on a diffuse condensed--gas interface.
//
// The volumetric thermal source is
//
//   qdot = q'' A(T) |eta_g grad(eta_c) - eta_c grad(eta_g)|,
//   A(T) = [1 + tanh((T - T_a) / Delta T_a)] / 2,
//
// so its volume integral converges to q'' times the condensed--gas area without
// selecting a sharp eta contour.  Positive heat_flux adds thermal energy.
// This model intentionally changes no species: it is appropriate only when
// the energy-producing process is outside the resolved chemistry.  Use an
// interphase_reaction when reactant and product masses are represented.
//

#ifndef MODEL_MECHANISM_INTERFACIALHEATSOURCE_H
#define MODEL_MECHANISM_INTERFACIALHEATSOURCE_H

#include <cmath>
#include <string>
#include <vector>

#include "Model/Chemistry/Chemistry.H"
#include "Model/Mechanism/Mechanism.H"
#include "Numeric/Stencil.H"
#include "Unit/Unit.H"

namespace Model
{
namespace Mechanism
{
class InterfacialHeatSource : public MechanismModel
{
public:
    static constexpr const char* name = "interfacial_heat_source";

    static void Parse(InterfacialHeatSource& value, IO::ParmParse& pp,
                    const std::vector<std::string>& species_names,
                    const int ngas_species,
                    const std::vector<int>& rigid_solid_species,
                    const std::vector<int>& liquid_species,
                    const std::vector<Set::Scalar>& reference_density,
                    const std::vector<double>& /*molecular_weight*/,
                    const Set::Scalar /*universal_gas_constant*/)
    {
        // Keep the liquid-only input for existing cases. The condensed key
        // additionally permits prescribed heating of an exposed rigid solid.
        const bool use_condensed_species = pp.contains("condensed");
        if (use_condensed_species && pp.contains("liquid"))
            Util::Exception(INFO, "Specify either condensed or liquid for "
                "an interfacial heat source");
        const std::string species_key = use_condensed_species ? "condensed" : "liquid";
        std::string heated_species_name;
        pp.query_required(species_key, heated_species_name);
        value.heated_species_index = -1;
        for (int n = 0; n < static_cast<int>(species_names.size()); ++n)
            if (species_names[n] == heated_species_name)
                value.heated_species_index = n;
        if (value.heated_species_index < 0)
            Util::Exception(INFO, "Unknown condensed species ", heated_species_name,
                            " in ", pp.full(species_key));

        bool registered_condensed_species = false;
        for (const int species : liquid_species)
            registered_condensed_species = registered_condensed_species ||
                species == value.heated_species_index;
        if (use_condensed_species)
            for (const int species : rigid_solid_species)
                registered_condensed_species = registered_condensed_species ||
                    species == value.heated_species_index;
        if (!registered_condensed_species)
            Util::Exception(INFO, pp.full(species_key),
                use_condensed_species ? " must identify a liquid or rigid solid" :
                                        " must identify a liquid");
        if (value.heated_species_index < ngas_species ||
            value.heated_species_index >=
                static_cast<int>(reference_density.size()) ||
            !(reference_density[value.heated_species_index] > 0.0))
            Util::Exception(INFO, heated_species_name,
                            " requires a positive reference density");

        value.ngas_species = ngas_species;
        value.nspecies = static_cast<int>(species_names.size());
        value.heated_inverse_reference_density =
            1.0 / reference_density[value.heated_species_index];
        for (int n = ngas_species; n < value.nspecies; ++n)
        {
            if (n >= static_cast<int>(reference_density.size()) ||
                !(reference_density[n] > 0.0))
                Util::Exception(INFO, species_names[n],
                                " requires a positive reference density");
            value.inverse_reference_density[n] = 1.0 / reference_density[n];
        }

        pp.query_required("heat_flux", value.heat_flux,
                          Unit::Power() / Unit::Area());
        pp.query_required("activation_temperature",
                          value.activation_temperature,
                          Unit::Temperature());
        pp.query_required("activation_width", value.activation_width,
                          Unit::Temperature());
        if (!(value.heat_flux >= 0.0) || !std::isfinite(value.heat_flux) ||
            !(value.activation_temperature > 0.0) ||
            !std::isfinite(value.activation_temperature) ||
            !(value.activation_width > 0.0) ||
            !std::isfinite(value.activation_width))
            Util::Exception(INFO, pp.full("heat_flux"),
                ", activation_temperature, and activation_width require a "
                "finite nonnegative flux and positive temperatures");
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar HeatSource(const State& state,
                           const int i, const int j, const int k,
                           const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        if (!(heat_flux > 0.0) || !(state.temperature > 0.0)) return 0.0;
        const Set::Scalar heated_volume_fraction =
            state.component_density(i,j,k,heated_species_index) *
                heated_inverse_reference_density;
        if (!(heated_volume_fraction > 0.0)) return 0.0;

        Set::Scalar condensed_eta = 0.0;
        Set::Vector condensed_gradient = Set::Vector::Zero();
        for (int n = ngas_species; n < nspecies; ++n)
        {
            condensed_eta += state.component_density(i,j,k,n) *
                inverse_reference_density[n];
            condensed_gradient += Numeric::Gradient(
                state.component_density, i, j, k, n, DX) *
                inverse_reference_density[n];
        }
        const Set::Scalar gas_eta = 1.0 - condensed_eta;
        if (!(gas_eta > 0.0)) return 0.0;
        const Set::Vector heated_volume_gradient = Numeric::Gradient(
            state.component_density, i, j, k, heated_species_index, DX) *
            heated_inverse_reference_density;
        const Set::Scalar activation = 0.5 *
            (1.0 + std::tanh((state.temperature - activation_temperature) /
                             activation_width));
        return heat_flux * activation *
            (gas_eta * heated_volume_gradient +
             heated_volume_fraction * condensed_gradient).norm();
    }

private:
    int nspecies = 0;
    int ngas_species = 0;
    int heated_species_index = -1;
    Set::Scalar heated_inverse_reference_density = 0.0;
    Set::Scalar heat_flux = 0.0;
    Set::Scalar activation_temperature = 0.0;
    Set::Scalar activation_width = 0.0;
    Model::Chemistry::SpeciesArray inverse_reference_density{};
};
}
}

#endif
