//
// Prescribed unresolved heat flux on a diffuse liquid--gas interface.
//
// The volumetric thermal source is
//
//   qdot = q'' A(T) |eta_g grad(eta_l) - eta_l grad(eta_g)|,
//   A(T) = [1 + tanh((T - T_a) / Delta T_a)] / 2,
//
// so its volume integral converges to q'' times the liquid--gas area without
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
class InterfacialHeatSource
{
public:
    static constexpr const char* name = "interfacial_heat_source";

    static void Parse(InterfacialHeatSource& value, IO::ParmParse& pp,
                    const std::vector<std::string>& species_names,
                    const int ngas_species,
                    const std::vector<int>& /*rigid_solid_species*/,
                    const std::vector<int>& liquid_species,
                    const std::vector<Set::Scalar>& reference_density,
                    const std::vector<double>& /*molecular_weight*/,
                    const Set::Scalar /*universal_gas_constant*/)
    {
        std::string liquid_name;
        pp.query_required("liquid", liquid_name);
        value.liquid_species_index = -1;
        for (int n = 0; n < static_cast<int>(species_names.size()); ++n)
            if (species_names[n] == liquid_name)
                value.liquid_species_index = n;
        if (value.liquid_species_index < 0)
            Util::Exception(INFO, "Unknown liquid species ", liquid_name,
                            " in ", pp.full("liquid"));

        bool registered_liquid = false;
        for (const int species : liquid_species)
            registered_liquid = registered_liquid ||
                species == value.liquid_species_index;
        if (!registered_liquid)
            Util::Exception(INFO, pp.full("liquid"),
                            " must identify a species with mechanics=liquid");
        if (value.liquid_species_index < ngas_species ||
            value.liquid_species_index >=
                static_cast<int>(reference_density.size()) ||
            !(reference_density[value.liquid_species_index] > 0.0))
            Util::Exception(INFO, liquid_name,
                            " requires a positive reference density");

        value.ngas_species = ngas_species;
        value.nspecies = static_cast<int>(species_names.size());
        value.liquid_inverse_reference_density =
            1.0 / reference_density[value.liquid_species_index];
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
    void Apply(const amrex::Array4<Set::Scalar>&,
               const State&, const int, const int, const int,
               const Set::Scalar[AMREX_SPACEDIM]) const
    {}

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar VolumeSource(const State&, const int, const int, const int,
                            const Set::Scalar[AMREX_SPACEDIM]) const
    {
        return 0.0;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar HeatSource(const State& state,
                           const int i, const int j, const int k,
                           const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        if (!(heat_flux > 0.0) || !(state.temperature > 0.0)) return 0.0;
        const Set::Scalar liquid_eta = Util::Clamp(
            state.component_density(i,j,k,liquid_species_index) *
                liquid_inverse_reference_density,
            0.0, 1.0);
        if (!(liquid_eta > 0.0)) return 0.0;

        Set::Scalar condensed_eta = 0.0;
        Set::Vector condensed_gradient = Set::Vector::Zero();
        for (int n = ngas_species; n < nspecies; ++n)
        {
            condensed_eta += Util::Max(
                state.component_density(i,j,k,n), 0.0) *
                inverse_reference_density[n];
            condensed_gradient += Numeric::Gradient(
                state.component_density, i, j, k, n, DX) *
                inverse_reference_density[n];
        }
        const Set::Scalar gas_eta = Util::Clamp(
            1.0 - condensed_eta, 0.0, 1.0);
        const Set::Vector liquid_gradient = Numeric::Gradient(
            state.component_density, i, j, k, liquid_species_index, DX) *
            liquid_inverse_reference_density;
        const Set::Scalar activation = 0.5 *
            (1.0 + std::tanh((state.temperature - activation_temperature) /
                             activation_width));
        return heat_flux * activation *
            (gas_eta * liquid_gradient +
             liquid_eta * condensed_gradient).norm();
    }

    Set::Scalar StabilityRate(const Set::Scalar,
                              const Set::Scalar) const
    {
        // Smooth activation is bounded and saturates at the prescribed flux;
        // it introduces no unbounded thermal feedback rate.
        return 0.0;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar LocalRate(const State&, const int, const int, const int) const
    {
        return 0.0;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar LocalStabilityRate(
        const State&, const int, const int, const int) const
    {
        return 0.0;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar GradientCoefficient(
        const State&, const int, const int, const int) const
    {
        return 0.0;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar ApplyImplicitChange(
        const amrex::Array4<Set::Scalar>&, const State&, const Set::Scalar,
        Set::Scalar&, const int, const int, const int) const
    {
        return 0.0;
    }

    int RigidComponent() const { return -1; }
    bool Equilibrium() const { return false; }
    bool Kinetic() const { return false; }
    Set::Scalar EquilibriumTemperature() const { return NAN; }
    Set::Scalar LatentHeat() const { return 0.0; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::GpuTuple<Set::Scalar,Set::Vector> KineticDiffuseSurface(
        const amrex::Array4<const Set::Scalar>&,
        const int, const int, const int,
        const Set::Scalar[AMREX_SPACEDIM]) const
    {
        return {0.0, Set::Vector::Zero()};
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::GpuTuple<Set::Scalar,Set::Scalar,Set::Scalar>
    ApplyEquilibriumChange(
        const amrex::Array4<Set::Scalar>&, const State&, const Set::Scalar,
        const Set::Scalar, const Set::Scalar,
        const int, const int, const int) const
    {
        return {0.0, 0.0, 0.0};
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar ApplyKineticChange(
        const amrex::Array4<Set::Scalar>&, const State&, const Set::Scalar,
        const Set::Scalar, const Set::Scalar, const Set::Scalar,
        Set::Scalar&, Set::Scalar&, const int, const int, const int) const
    {
        return 0.0;
    }

private:
    int nspecies = 0;
    int ngas_species = 0;
    int liquid_species_index = -1;
    Set::Scalar liquid_inverse_reference_density = 0.0;
    Set::Scalar heat_flux = 0.0;
    Set::Scalar activation_temperature = 0.0;
    Set::Scalar activation_width = 0.0;
    Model::Chemistry::SpeciesArray inverse_reference_density{};
};
}
}

#endif
