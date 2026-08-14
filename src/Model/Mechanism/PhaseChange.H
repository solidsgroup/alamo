//
// Conservative transfer from one condensed species to either one condensed
// species or a fixed-composition gas bundle.  Every mechanism explicitly
// selects one physical closure:
//
//   equilibrium_enthalpy     local enthalpy complementarity at a transition T
//   arrhenius_surface_flux   pressure/temperature-dependent kg m^-2 s^-1 flux
//   hertz_knudsen            vapor-pressure-driven liquid/gas mass flux
//
// The kinetic closures multiply their dimensional surface flux by the local,
// reduction-consistent diffuse phase0/gas surface measure
//
//   delta_0g = |eta_g grad(eta_0) - eta_0 grad(eta_g)|.
//
// This measure integrates to one across an isolated binary interface and
// vanishes on phase0/condensed interfaces that are not exposed to gas.  The
// constitutive flux uses the local diffuse-mixture temperature.  No contour,
// normal extrapolation, signed distance, or sharp-interface state is
// reconstructed.  The local mass/enthalpy system is solved implicitly after
// transport, so kinetic stiffness does not impose an explicit global step.
//
// latent_heat is positive for an endothermic forward transition.  The signed
// coupled_enthalpy_change represents an unresolved process tied only to
// forward phase0 consumption: negative is exothermic and positive is
// endothermic.  It is kept separate from latent heat so equilibrium and vapor
// pressure retain their physical thermodynamic meanings.
//

#ifndef MODEL_MECHANISM_PHASECHANGE_H
#define MODEL_MECHANISM_PHASECHANGE_H

#include <array>
#include <cmath>
#include <string>
#include <vector>

#include "Model/Chemistry/Chemistry.H"
#include "Model/Mechanism/Mechanism.H"
#include "Numeric/Stencil.H"
#include "Set/Set.H"
#include "Unit/Unit.H"

namespace Model
{
namespace Mechanism
{
class PhaseChange : public MechanismModel
{
public:
    static constexpr const char* name = "phase_change";
    static constexpr int max_species = Model::Chemistry::MAX_SPECIES;

    static void Parse(PhaseChange& value, IO::ParmParse& pp,
                    const std::vector<std::string>& species_names,
                    const int ngas_species,
                    const std::vector<int>& /*rigid_solid_species*/,
                    const std::vector<int>& liquid_species,
                    const std::vector<Set::Scalar>& reference_density,
                    const std::vector<double>& molecular_weight,
                    const Set::Scalar universal_gas_constant)
    {
        value.nspecies = static_cast<int>(species_names.size());
        value.ngas_species = ngas_species;
        for (int n = 0; n < ngas_species; ++n)
        {
            if (!(molecular_weight[n] > 0.0))
                Util::Exception(INFO, species_names[n],
                    " requires a positive molecular weight");
            value.inverse_molecular_weight[n] = 1.0 / molecular_weight[n];
        }
        for (int n = ngas_species; n < value.nspecies; ++n)
            if (reference_density[n] > 0.0)
                value.inverse_reference_density[n] =
                    1.0 / reference_density[n];

        std::string phase0_name;
        std::vector<std::string> phase1_names;
        std::vector<Set::Scalar> phase1_mass_fractions;
        pp.query_required("phase0", phase0_name);
        pp.queryarr_required("phase1", phase1_names);

        for (int n = 0; n < value.nspecies; ++n)
            if (species_names[n] == phase0_name) value.phase0 = n;
        if (value.phase0 < 0)
            Util::Exception(INFO, "Unknown species identifier ", phase0_name,
                            " in ", pp.full("phase0"));
        if (value.phase0 < ngas_species)
            Util::Exception(INFO, pp.full("phase0"),
                " must identify a condensed species");
        value.phase0_reference_density = reference_density[value.phase0];
        if (!(value.phase0_reference_density > 0.0))
            Util::Exception(INFO, phase0_name,
                " requires a positive reference density");
        for (const int species : liquid_species)
            value.phase0_is_liquid = value.phase0_is_liquid ||
                species == value.phase0;

        if (phase1_names.empty() || phase1_names.size() > max_species)
            Util::Exception(INFO, pp.full("phase1"),
                " must contain between 1 and ", max_species, " species");
        value.nphase1_species = static_cast<int>(phase1_names.size());
        for (int n = 0; n < value.nphase1_species; ++n)
        {
            value.phase1_species[n] = -1;
            for (int m = 0; m < value.nspecies; ++m)
                if (species_names[m] == phase1_names[n])
                    value.phase1_species[n] = m;
            if (value.phase1_species[n] < 0)
                Util::Exception(INFO, "Unknown species identifier ",
                    phase1_names[n], " in ", pp.full("phase1"));
            if (value.phase1_species[n] == value.phase0)
                Util::Exception(INFO, phase0_name,
                    " cannot appear in both phase0 and phase1");
        }

        value.phase1_is_gas = value.phase1_species[0] < ngas_species;
        if (!value.phase1_is_gas && value.nphase1_species != 1)
            Util::Exception(INFO,
                "A condensed phase1 must contain exactly one species");
        for (int n = 0; n < value.nphase1_species; ++n)
            if ((value.phase1_species[n] < ngas_species) !=
                value.phase1_is_gas)
                Util::Exception(INFO,
                    "phase1 cannot mix gas and condensed species");

        if (value.phase1_is_gas)
        {
            if (pp.contains("phase1.mass_fractions"))
                pp.queryarr_required("phase1.mass_fractions",
                                     phase1_mass_fractions);
            else if (value.nphase1_species == 1)
                phase1_mass_fractions.push_back(1.0);
            else
                Util::Exception(INFO, pp.full("phase1.mass_fractions"),
                    " is required for a multi-species gas phase");
            if (phase1_mass_fractions.size() != phase1_names.size())
                Util::Exception(INFO, pp.full("phase1.mass_fractions"),
                    " must contain one value per gas species");

            Set::Scalar sum = 0.0;
            for (int n = 0; n < value.nphase1_species; ++n)
            {
                const int species = value.phase1_species[n];
                if (phase1_mass_fractions[n] < 0.0)
                    Util::Exception(INFO, pp.full("phase1.mass_fractions"),
                        " must be nonnegative");
                if (!(molecular_weight[species] > 0.0))
                    Util::Exception(INFO, phase1_names[n],
                        " requires a positive molecular weight");
                value.phase1_mass_fractions[n] = phase1_mass_fractions[n];
                value.phase1_gas_constant += universal_gas_constant *
                    phase1_mass_fractions[n] / molecular_weight[species];
                sum += phase1_mass_fractions[n];
            }
            if (Util::Abs(sum - 1.0) > 1.0e-12)
                Util::Exception(INFO, pp.full("phase1.mass_fractions"),
                    " must sum to one");
        }
        else
        {
            if (pp.contains("phase1.mass_fractions"))
                Util::Exception(INFO, pp.full("phase1.mass_fractions"),
                    " is only valid for a gas phase1");
            value.phase1_mass_fractions[0] = 1.0;
            value.phase1_inverse_reference_density =
                1.0 / reference_density[value.phase1_species[0]];
        }

        pp.query_default("reversible", value.reversible, false);
        pp.query_default("latent_heat", value.latent_heat,
                         "0.0_J/kg", Unit::Energy() / Unit::Mass());
        pp.query_default("coupled_enthalpy_change",
                         value.coupled_enthalpy_change,
                         "0.0_J/kg", Unit::Energy() / Unit::Mass());
        if (!(value.latent_heat >= 0.0) ||
            !std::isfinite(value.coupled_enthalpy_change))
            Util::Exception(INFO,
                "latent_heat must be nonnegative and "
                "coupled_enthalpy_change must be finite");

        std::string kinetics;
        pp.query_required("kinetics", kinetics);
        if (kinetics == "equilibrium_enthalpy")
        {
            value.kinetics = Kinetics::EquilibriumEnthalpy;
            pp.query_required("transition_temperature",
                              value.transition_temperature,
                              Unit::Temperature());
            if (!(value.transition_temperature > 0.0) ||
                !(value.latent_heat > 0.0))
                Util::Exception(INFO,
                    "equilibrium_enthalpy requires a positive "
                    "transition_temperature and latent_heat");
        }
        else if (kinetics == "arrhenius_surface_flux")
        {
            value.kinetics = Kinetics::ArrheniusSurfaceFlux;
            pp.query_required("reference_mass_flux",
                value.reference_mass_flux,
                Unit::Mass() / Unit::Area() / Unit::Time());
            pp.query_required("reference_pressure", value.reference_pressure,
                              Unit::Pressure());
            pp.query_required("pressure_exponent", value.pressure_exponent);
            pp.query_required("activation_temperature",
                              value.activation_temperature,
                              Unit::Temperature());
            pp.query_required("reference_temperature",
                              value.reference_temperature,
                              Unit::Temperature());
            if (!(value.reference_mass_flux >= 0.0) ||
                !(value.reference_pressure > 0.0) ||
                !std::isfinite(value.pressure_exponent) ||
                !(value.activation_temperature >= 0.0) ||
                !(value.reference_temperature > 0.0))
                Util::Exception(INFO,
                    "Invalid arrhenius_surface_flux parameter");
        }
        else if (kinetics == "hertz_knudsen")
        {
            value.kinetics = Kinetics::HertzKnudsen;
            if (!value.phase0_is_liquid || !value.phase1_is_gas ||
                value.nphase1_species != 1)
                Util::Exception(INFO,
                    "hertz_knudsen requires a liquid phase0 and one gas "
                    "phase1 species");
            pp.query_required("accommodation_coefficient",
                              value.accommodation_coefficient);
            pp.query_required("saturation_reference_temperature",
                              value.saturation_reference_temperature,
                              Unit::Temperature());
            pp.query_required("saturation_reference_pressure",
                              value.saturation_reference_pressure,
                              Unit::Pressure());
            if (!(value.accommodation_coefficient > 0.0) ||
                value.accommodation_coefficient > 1.0 ||
                !(value.saturation_reference_temperature > 0.0) ||
                !(value.saturation_reference_pressure > 0.0) ||
                !(value.latent_heat > 0.0))
                Util::Exception(INFO,
                    "Invalid hertz_knudsen accommodation, saturation "
                    "reference, or latent heat");
        }
        else
            Util::Exception(INFO, pp.full("kinetics"),
                " must be equilibrium_enthalpy, arrhenius_surface_flux, "
                "or hertz_knudsen");

        if (value.Kinetic())
        {
            pp.query_default("nonlinear.max_iterations",
                value.nonlinear_max_iterations, 40);
            pp.query_default("nonlinear.relative_tolerance",
                value.nonlinear_relative_tolerance, 1.0e-10);
            pp.query_default("nonlinear.absolute_mass_tolerance",
                value.nonlinear_absolute_mass_tolerance,
                "1.0e-12_kg/m^3", Unit::Density());
            if (value.nonlinear_max_iterations < 1 ||
                !(value.nonlinear_relative_tolerance > 0.0) ||
                !(value.nonlinear_absolute_mass_tolerance > 0.0))
                Util::Exception(INFO,
                    "Kinetic phase-change nonlinear controls must be positive");
        }
    }

    bool Equilibrium() const
    {
        return kinetics == Kinetics::EquilibriumEnthalpy;
    }
    bool Kinetic() const
    {
        return kinetics == Kinetics::ArrheniusSurfaceFlux ||
               kinetics == Kinetics::HertzKnudsen;
    }
    Set::Scalar EquilibriumTemperature() const
    {
        return transition_temperature;
    }
    Set::Scalar LatentHeat() const { return latent_heat; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar KineticDiffuseSurfaceMeasure(
        const amrex::Array4<const Set::Scalar>& component_density,
        const int i, const int j, const int k,
        const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        if (!Kinetic()) return 0.0;

        const Set::Scalar eta0 =
            component_density(i,j,k,phase0) / phase0_reference_density;
        if (!(eta0 > 0.0))
            return 0.0;

        const Set::Vector gradient0 = Numeric::Gradient(
            component_density, i, j, k, phase0, DX) /
            phase0_reference_density;
        Set::Scalar condensed_fraction = 0.0;
        Set::Vector condensed_gradient = Set::Vector::Zero();
        for (int n = ngas_species; n < nspecies; ++n)
            if (inverse_reference_density[n] > 0.0)
            {
                condensed_fraction += component_density(i,j,k,n) *
                    inverse_reference_density[n];
                condensed_gradient += Numeric::Gradient(
                    component_density, i, j, k, n, DX) *
                    inverse_reference_density[n];
            }
        const Set::Scalar gas_fraction = 1.0 - condensed_fraction;
        if (!(gas_fraction > 0.0)) return 0.0;
        const Set::Vector surface_vector =
            gas_fraction * gradient0 + eta0 * condensed_gradient;
        return surface_vector.norm();
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::GpuTuple<Set::Scalar,Set::Scalar,Set::Scalar>
    ApplyEquilibriumChange(
        const amrex::Array4<Set::Scalar>& component_density,
        const State& state, Set::Scalar requested_mass_change,
        const Set::Scalar, const Set::Scalar,
        const int i, const int j, const int k) const
    {
        if (!Equilibrium()) return {0.0, 0.0, 0.0};
        Set::Scalar mass_change = requested_mass_change;
        if (!reversible) mass_change = Util::Min(mass_change, 0.0);
        mass_change = Util::Max(mass_change,
            -Util::Max(component_density(i,j,k,phase0), 0.0));
        if (mass_change > 0.0)
        {
            Set::Scalar available = 1.0e300;
            for (int n = 0; n < nphase1_species; ++n)
                if (phase1_mass_fractions[n] > 0.0)
                    available = Util::Min(available,
                        Util::Max(component_density(
                            i,j,k,phase1_species[n]), 0.0) /
                            phase1_mass_fractions[n]);
            mass_change = Util::Min(mass_change, available);
        }
        if (mass_change == 0.0) return {0.0, 0.0, 0.0};

        Set::Scalar phase1_specific_volume =
            phase1_inverse_reference_density;
        if (phase1_is_gas)
        {
            if (!(state.pressure > 0.0)) return {0.0, 0.0, 0.0};
            phase1_specific_volume = phase1_gas_constant *
                transition_temperature / state.pressure;
        }
        component_density(i,j,k,phase0) += mass_change;
        for (int n = 0; n < nphase1_species; ++n)
            component_density(i,j,k,phase1_species[n]) -=
                phase1_mass_fractions[n] * mass_change;

        const Set::Scalar volume_change = mass_change *
            (1.0 / phase0_reference_density - phase1_specific_volume);
        const Set::Scalar coupled_heat = coupled_enthalpy_change *
            Util::Min(mass_change, 0.0);
        return {mass_change, volume_change, coupled_heat};
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar ApplyKineticChange(
        const amrex::Array4<Set::Scalar>& component_density,
        const State& state, const Set::Scalar heat_capacity,
        const Set::Scalar surface_measure,
        Set::Scalar& temperature, Set::Scalar& integrated_heat,
        const int i, const int j, const int k) const
    {
        if (!Kinetic() || !(state.dt > 0.0) ||
            !(heat_capacity > 0.0) || !(temperature > 0.0) ||
            !(surface_measure > 0.0))
            return 0.0;

        const Set::Scalar old_temperature = temperature;
        const Set::Scalar phase0_mass = component_density(i,j,k,phase0);
        if (!(phase0_mass > 0.0)) return 0.0;
        Set::Scalar reverse_mass = 0.0;
        if (reversible)
        {
            reverse_mass = 1.0e300;
            for (int n = 0; n < nphase1_species; ++n)
                if (phase1_mass_fractions[n] > 0.0)
                    reverse_mass = Util::Min(reverse_mass,
                        component_density(i,j,k,phase1_species[n]) /
                            phase1_mass_fractions[n]);
        }

        const auto residual = [=] AMREX_GPU_DEVICE(Set::Scalar mass_change)
        {
            const Set::Scalar heat = latent_heat * mass_change +
                coupled_enthalpy_change * Util::Min(mass_change, 0.0);
            const Set::Scalar kinetic_temperature =
                old_temperature + heat / heat_capacity;
            if (!(kinetic_temperature > 0.0)) return mass_change;

            Set::Scalar mass_flux = 0.0;
            if (kinetics == Kinetics::ArrheniusSurfaceFlux)
            {
                if (!(state.pressure > 0.0)) return mass_change;
                const Set::Scalar pressure_ratio =
                    state.pressure / reference_pressure;
                const Set::Scalar exponent = -activation_temperature *
                    (1.0 / kinetic_temperature -
                     1.0 / reference_temperature);
                mass_flux = reference_mass_flux *
                    std::pow(pressure_ratio, pressure_exponent) *
                    std::exp(exponent);
            }
            else
            {
                const Set::Scalar gas_specific_constant =
                    phase1_gas_constant;
                const Set::Scalar saturation_exponent =
                    -latent_heat / gas_specific_constant *
                    (1.0 / kinetic_temperature -
                     1.0 / saturation_reference_temperature);
                const Set::Scalar saturation_pressure =
                    saturation_reference_pressure *
                    std::exp(saturation_exponent);
                Set::Scalar gas_molar_density = 0.0;
                for (int n = 0; n < ngas_species; ++n)
                    gas_molar_density +=
                        (component_density(i,j,k,n) -
                         (n == phase1_species[0] ? mass_change : 0.0)) *
                        inverse_molecular_weight[n];
                const Set::Scalar vapor_molar_density =
                    (component_density(i,j,k,phase1_species[0]) -
                     mass_change) *
                    inverse_molecular_weight[phase1_species[0]];
                const Set::Scalar vapor_partial_pressure =
                    gas_molar_density > 0.0 ?
                    state.pressure * vapor_molar_density /
                        gas_molar_density : 0.0;
                mass_flux = 2.0 * accommodation_coefficient /
                    (2.0 - accommodation_coefficient) *
                    (saturation_pressure - vapor_partial_pressure) /
                    std::sqrt(2.0 * Set::Constant::Pi *
                              gas_specific_constant * kinetic_temperature);
                if (!reversible) mass_flux = Util::Max(mass_flux, 0.0);
            }
            return mass_change +
                mass_flux * surface_measure * state.dt;
        };

        Set::Scalar lower = -phase0_mass;
        Set::Scalar upper = reversible ? reverse_mass : 0.0;
        const Set::Scalar zero_residual = residual(0.0);
        if (zero_residual > 0.0)
        {
            upper = 0.0;
            if (residual(lower) > 0.0)
                upper = lower;
        }
        else if (zero_residual < 0.0)
        {
            lower = 0.0;
            if (residual(upper) < 0.0)
                lower = upper;
        }
        else
            return 0.0;

        for (int iteration = 0;
             iteration < nonlinear_max_iterations && lower < upper;
             ++iteration)
        {
            const Set::Scalar midpoint = 0.5 * (lower + upper);
            if (residual(midpoint) > 0.0)
                upper = midpoint;
            else
                lower = midpoint;
            const Set::Scalar tolerance = nonlinear_absolute_mass_tolerance +
                nonlinear_relative_tolerance *
                    Util::Max(Util::Abs(lower), Util::Abs(upper));
            if (upper - lower <= tolerance) break;
        }
        const Set::Scalar mass_change = 0.5 * (lower + upper);
        if (mass_change == 0.0) return 0.0;

        component_density(i,j,k,phase0) += mass_change;
        for (int n = 0; n < nphase1_species; ++n)
            component_density(i,j,k,phase1_species[n]) -=
                phase1_mass_fractions[n] * mass_change;
        integrated_heat = latent_heat * mass_change +
            coupled_enthalpy_change * Util::Min(mass_change, 0.0);
        temperature = old_temperature + integrated_heat / heat_capacity;

        Set::Scalar phase1_specific_volume =
            phase1_inverse_reference_density;
        if (phase1_is_gas)
            phase1_specific_volume = phase1_gas_constant *
                temperature / state.pressure;
        return mass_change *
            (1.0 / phase0_reference_density - phase1_specific_volume);
    }

private:
    enum class Kinetics
    {
        EquilibriumEnthalpy,
        ArrheniusSurfaceFlux,
        HertzKnudsen
    };

    Kinetics kinetics = Kinetics::EquilibriumEnthalpy;
    int nspecies = 0;
    int ngas_species = 0;
    int phase0 = -1;
    int nphase1_species = 0;
    std::array<int,max_species> phase1_species{};
    std::array<Set::Scalar,max_species> phase1_mass_fractions{};
    std::array<Set::Scalar,max_species> inverse_reference_density{};
    std::array<Set::Scalar,max_species> inverse_molecular_weight{};
    Set::Scalar phase0_reference_density = NAN;
    Set::Scalar phase1_inverse_reference_density = 0.0;
    Set::Scalar phase1_gas_constant = 0.0;
    Set::Scalar latent_heat = 0.0;
    Set::Scalar coupled_enthalpy_change = 0.0;
    Set::Scalar transition_temperature = NAN;
    Set::Scalar reference_mass_flux = NAN;
    Set::Scalar reference_pressure = NAN;
    Set::Scalar pressure_exponent = NAN;
    Set::Scalar activation_temperature = NAN;
    Set::Scalar reference_temperature = NAN;
    Set::Scalar accommodation_coefficient = NAN;
    Set::Scalar saturation_reference_temperature = NAN;
    Set::Scalar saturation_reference_pressure = NAN;
    Set::Scalar nonlinear_relative_tolerance = 1.0e-10;
    Set::Scalar nonlinear_absolute_mass_tolerance = 1.0e-12;
    int nonlinear_max_iterations = 40;
    bool phase0_is_liquid = false;
    bool phase1_is_gas = false;
    bool reversible = false;
};
}
}

#endif
