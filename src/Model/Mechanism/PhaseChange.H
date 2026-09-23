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
// reduction-consistent diffuse phase0/phase1 surface measure
//
//   delta_01 = |eta_1 grad(eta_0) - eta_0 grad(eta_1)|.
//
// Irreversible solid-to-gas regression evaluates the normal derivative with
// the integrator's input-selected advection operator and limiter. The local
// inward normal is frozen over that reconstruction; recoil retains the
// centered geometric pair gradient. This measure integrates to one across a
// resolved monotone planar binary interface and
// vanishes on interfaces that do not contain both selected phases.  The
// constitutive flux defaults to the local diffuse-mixture temperature.
// Irreversible solid-to-gas Arrhenius kinetics may instead use temperature
// reconstructed at the total-condensed eta=0.5 contour. The reconstruction
// is frozen during each phase-change update; the local enthalpy increment
// is still included in the implicit kinetic temperature. This split treatment
// requires timestep convergence, rather than a fully implicit surface solve.
// A kinetic transition sends the occupied-volume change from that same mass
// and heat update to the projection and supplies the associated recoil.
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
#include "Numeric/Advect/Advect.H"
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

    // Chen's steady planar reference balance (q is heat INTO the surface,
    // Q is signed condensed heat release):
    //   r = A exp(-Ta/Ts),  cp (Ts-T0) = Q + q/(rho*r).
    // These reference calculations do not replace the transient energy
    // equation or the reconstructed kinetic temperature during a run.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static bool ChenSurfaceReference(
        const Set::Scalar rho, const Set::Scalar cp, const Set::Scalar Q,
        const Set::Scalar A, const Set::Scalar Ta, const Set::Scalar T0,
        const Set::Scalar q, Set::Scalar& Ts, Set::Scalar& r)
    {
        Ts = r = NAN;
        if (!(rho > 0.0) || !std::isfinite(rho) ||
            !(cp > 0.0) || !std::isfinite(cp) || !std::isfinite(Q) ||
            !(A > 0.0) || !std::isfinite(A) ||
            !(Ta >= 0.0) || !std::isfinite(Ta) ||
            !(T0 > 0.0) || !std::isfinite(T0) ||
            !(q > 0.0) || !std::isfinite(q)) return false;
        const Set::Scalar equilibrium_temperature = T0 + Q / cp;
        if (!std::isfinite(equilibrium_temperature)) return false;
        Set::Scalar lower = Util::Max(0.0, equilibrium_temperature);
        Set::Scalar upper = Util::Max(1000.0, 2.0 * lower);
        const Set::Scalar log_scale =
            std::log(rho) + std::log(A) + std::log(cp) - std::log(q);
        // The log residual is monotone for Ta >= 0. It avoids Arrhenius
        // underflow at the lower bracket and mass-flux overflow at the upper.
        bool bracketed = false;
        for (int iteration = 0; iteration < 256; ++iteration)
        {
            if (!std::isfinite(upper)) return false;
            if (log_scale - Ta / upper +
                std::log(upper - equilibrium_temperature) >= 0.0)
            {
                bracketed = true;
                break;
            }
            upper *= 2.0;
        }
        if (!bracketed) return false;
        // With Ta=0 and a negative equilibrium temperature, a positive
        // root need not exist even for positive incident flux.
        if (Ta == 0.0 && equilibrium_temperature < 0.0 &&
            log_scale + std::log(-equilibrium_temperature) >= 0.0)
            return false;
        for (int iteration = 0; iteration < 128; ++iteration)
        {
            const Set::Scalar trial = lower + 0.5 * (upper - lower);
            if (trial == lower || trial == upper) break;
            if (log_scale - Ta / trial +
                std::log(trial - equilibrium_temperature) > 0.0)
                upper = trial;
            else lower = trial;
        }
        Ts = lower + 0.5 * (upper - lower);
        r = std::exp(std::log(A) - Ta / Ts);
        return Ts > 0.0 && std::isfinite(Ts) && r > 0.0 && std::isfinite(r);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    static bool ChenReferencePrefactor(
        const Set::Scalar rho, const Set::Scalar cp, const Set::Scalar Q,
        const Set::Scalar Ta, const Set::Scalar T0, const Set::Scalar q,
        const Set::Scalar r, Set::Scalar& Ts, Set::Scalar& A)
    {
        Ts = A = NAN;
        if (!(rho > 0.0) || !std::isfinite(rho) ||
            !(cp > 0.0) || !std::isfinite(cp) || !std::isfinite(Q) ||
            !(Ta >= 0.0) || !std::isfinite(Ta) ||
            !(T0 > 0.0) || !std::isfinite(T0) ||
            !(q > 0.0) || !std::isfinite(q) ||
            !(r > 0.0) || !std::isfinite(r)) return false;
        Ts = T0 + (Q + q / rho / r) / cp;
        if (!(Ts > 0.0) || !std::isfinite(Ts)) return false;
        A = std::exp(std::log(r) + Ta / Ts);
        return A > 0.0 && std::isfinite(A);
    }

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
        value.gas_constant = universal_gas_constant;
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

        std::string thermal_closure;
        pp.query_default("thermal_closure", thermal_closure, "local_heat");
        if (thermal_closure != "local_heat" && thermal_closure != "reference_enthalpy")
            Util::Exception(INFO, "Unknown phase-change thermal_closure: ", thermal_closure);
        value.reference_enthalpy = thermal_closure == "reference_enthalpy";
        pp.query_default("enthalpy_reference_temperature", value.enthalpy_reference_temperature,
            "300_K", Unit::Temperature());
        if (!(value.enthalpy_reference_temperature > 0.0) ||
            !std::isfinite(value.enthalpy_reference_temperature))
            Util::Exception(INFO, "enthalpy_reference_temperature must be positive and finite");

        std::string kinetics;
        pp.query_required("kinetics", kinetics);
        if (value.reference_enthalpy &&
            (kinetics != "arrhenius_surface_flux" || value.phase0_is_liquid ||
             !value.phase1_is_gas || value.reversible))
            Util::Exception(INFO, "reference_enthalpy requires irreversible solid-to-gas Arrhenius kinetics");
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
            pp.query_default("reconstruct_surface_temperature",
                value.reconstruct_surface_temperature, false);
            if (value.reconstruct_surface_temperature &&
                (value.phase0_is_liquid || !value.phase1_is_gas || value.reversible))
                Util::Exception(INFO,
                    "reconstruct_surface_temperature requires irreversible solid-to-gas kinetics");
            pp.query_required("reference_pressure", value.reference_pressure,
                              Unit::Pressure());
            pp.query_required("pressure_exponent", value.pressure_exponent);
            pp.query_required("activation_temperature",
                              value.activation_temperature,
                              Unit::Temperature());
            std::string reconstruction = "contour";
            pp.query_default("surface_temperature_reconstruction", reconstruction, "contour");
            if (reconstruction == "contour")
                value.surface_temperature_options.mode = SurfaceTemperatureOptions::Contour;
            else if (reconstruction == "solid_mass_weighted")
                value.surface_temperature_options.mode = SurfaceTemperatureOptions::SolidMass;
            else if (reconstruction == "conversion_weighted")
                value.surface_temperature_options.mode = SurfaceTemperatureOptions::Conversion;
            else if (reconstruction == "arrhenius_equivalent")
                value.surface_temperature_options.mode = SurfaceTemperatureOptions::ArrheniusEquivalent;
            else Util::Exception(INFO, "Unknown surface_temperature_reconstruction: ", reconstruction);
            if (reconstruction != "contour" && !value.reconstruct_surface_temperature)
                Util::Exception(INFO, "Weighted surface temperature requires reconstruct_surface_temperature=true");
            pp.query_default("surface_temperature_band_cutoff",
                value.surface_temperature_options.band_cutoff, 0.01);
            if (!(value.surface_temperature_options.band_cutoff > 1.0e-10 &&
                  value.surface_temperature_options.band_cutoff < 0.5))
                Util::Exception(INFO, "surface_temperature_band_cutoff must lie in (1e-10,0.5)");
            pp.query_default("surface_temperature_contour_fraction",
                value.surface_temperature_options.contour_fraction, 0.5);
            const Set::Scalar contour = value.surface_temperature_options.contour_fraction;
            if (!(contour > 1.0e-10 && contour < 1.0-1.0e-10))
                Util::Exception(INFO, "surface_temperature_contour_fraction must lie in (1e-10,1-1e-10)");
            if (contour != 0.5 && !value.reconstruct_surface_temperature)
                Util::Exception(INFO, "Nondefault contour requires reconstruct_surface_temperature=true");
            if (reconstruction != "contour" &&
                !(contour > value.surface_temperature_options.band_cutoff &&
                  contour < 1.0-value.surface_temperature_options.band_cutoff))
                Util::Exception(INFO, "Weighted surface contour must lie inside its diffuse band");
            value.surface_temperature_options.activation_temperature = value.activation_temperature;
            value.surface_temperature_options.solid_component = value.phase0;
            if (pp.contains("pre_exponential_speed"))
                // A fixes the absolute Arrhenius law; this temperature only
                // normalizes its stored reference mass flux.
                pp.query_default("reference_temperature",
                                 value.reference_temperature, "1000_K",
                                 Unit::Temperature());
            else
                // A supplied reference mass flux has meaning only together
                // with its explicitly specified reference temperature.
                pp.query_required("reference_temperature",
                                  value.reference_temperature, Unit::Temperature());
            if (pp.contains("pre_exponential_speed"))
            {
                if (pp.contains("reference_mass_flux"))
                    Util::Exception(INFO, "Specify either pre_exponential_speed "
                        "or reference_mass_flux, not both");
                Set::Scalar pre_exponential_speed;
                pp.query_required("pre_exponential_speed",
                                  pre_exponential_speed, Unit::Velocity());
                if (!(pre_exponential_speed >= 0.0) ||
                    !std::isfinite(pre_exponential_speed))
                    Util::Exception(INFO,
                        "pre_exponential_speed must be finite and nonnegative");
                // r = A exp(-E/(R T)) is a recession speed. Convert it once
                // to the dimensional mass flux at the reference temperature;
                // the diffuse surface measure supplies the spatial support.
                value.reference_mass_flux = value.phase0_reference_density *
                    pre_exponential_speed * std::exp(
                        -value.activation_temperature / value.reference_temperature);
            }
            else
                pp.query_required("reference_mass_flux",
                    value.reference_mass_flux,
                    Unit::Mass() / Unit::Area() / Unit::Time());
            if (!(value.reference_mass_flux >= 0.0) ||
                !std::isfinite(value.reference_mass_flux) ||
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

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool Equilibrium() const
    {
        return kinetics == Kinetics::EquilibriumEnthalpy;
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
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
    bool ReconstructSurfaceTemperature() const
    {
        return reconstruct_surface_temperature;
    }
    SurfaceTemperatureOptions SurfaceTemperatureClosure() const
    {
        return surface_temperature_options;
    }
    template<class Advection>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::GpuTuple<Set::Scalar,Set::Vector>
    KineticDiffuseSurfaceGeometry(
        const amrex::Array4<const Set::Scalar>& component_density,
        const int i, const int j, const int k,
        const Set::Scalar DX[AMREX_SPACEDIM], const Advection& advection) const
    {
        if (!Kinetic())
            return {0.0, Set::Vector::Zero()};

        const Set::Scalar eta0 =
            component_density(i,j,k,phase0) / phase0_reference_density;
        if (!(eta0 > 0.0))
            return {0.0, Set::Vector::Zero()};

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
        Set::Scalar eta1 = 0.0;
        Set::Vector gradient1 = Set::Vector::Zero();
        if (phase1_is_gas)
        {
            eta1 = 1.0 - condensed_fraction;
            gradient1 = -condensed_gradient;
        }
        else
        {
            const int phase1 = phase1_species[0];
            eta1 = component_density(i,j,k,phase1) *
                inverse_reference_density[phase1];
            gradient1 = Numeric::Gradient(
                component_density, i, j, k, phase1, DX) *
                inverse_reference_density[phase1];
        }
        const bool irreversible_solid_regression =
            kinetics == Kinetics::ArrheniusSurfaceFlux &&
            !phase0_is_liquid && phase1_is_gas && !reversible;
        if (!(eta1 > 0.0) && !irreversible_solid_regression)
            return {0.0, Set::Vector::Zero()};
        eta1 = Util::Max(eta1, 0.0);
        const Set::Vector phase_pair_gradient =
            eta1 * gradient0 - eta0 * gradient1;
        Set::Scalar surface_measure = phase_pair_gradient.norm();
        if (irreversible_solid_regression && surface_measure > 0.0)
        {
            // With the central fractions held fixed, grad(pair_fraction) is
            // eta_g grad(eta_s) - eta_s grad(eta_g). Regression is transport
            // along its inward normal. Reuse the configured reconstruction
            // instead of introducing a separate Hamilton-Jacobi stencil.
            const Set::Vector inward_normal = phase_pair_gradient / surface_measure;
            const int regressing_species = phase0;
            const Set::Scalar regressing_density = phase0_reference_density;
            const int first_condensed_species = ngas_species;
            const int total_species = nspecies;
            const auto inverse_densities = inverse_reference_density;
            const auto pair_fraction = [=] AMREX_GPU_HOST_DEVICE(
                int sample_i, int sample_j, int sample_k, int) -> Set::Scalar
            {
                Set::Scalar gas_fraction = 1.0;
                for (int species = first_condensed_species;
                     species < total_species; ++species)
                    gas_fraction -= component_density(
                        sample_i, sample_j, sample_k, species) * inverse_densities[species];
                return eta1 * component_density(
                    sample_i, sample_j, sample_k, regressing_species) / regressing_density -
                    eta0 * Util::Max(gas_fraction, 0.0);
            };
            const auto normal_velocity = [=] AMREX_GPU_HOST_DEVICE(
                int, int, int, int direction) -> Set::Scalar
            {
                return inward_normal(direction);
            };
            // The advection operator returns -n.grad(pair_fraction). Its
            // negative is the surface measure multiplying the implicit mass
            // flux. No velocity-divergence/curvature source belongs here.
            surface_measure = -advection(pair_fraction, normal_velocity,
                i, j, k, 0, DX, {Numeric::Advect::Form::Advective});
        }
        if (!(surface_measure > 0.0))
            return {0.0, Set::Vector::Zero()};
        return {surface_measure, phase_pair_gradient};
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

    bool UsesReferenceEnthalpy() const { return reference_enthalpy; }

    void ConfigurePhaseHeatCapacities(const std::vector<Set::Scalar>& cp)
    {
        if (!reference_enthalpy) return;
        if (cp.size() != static_cast<std::size_t>(nspecies) || !(cp[phase0] > 0.0))
            Util::Exception(INFO, "reference_enthalpy requires material heat capacities");
        reactant_minus_product_cp = cp[phase0];
        for (int n = 0; n < nphase1_species; ++n)
        {
            if (!(cp[phase1_species[n]] > 0.0))
                Util::Exception(INFO, "reference_enthalpy requires positive product heat capacities");
            reactant_minus_product_cp -= phase1_mass_fractions[n] * cp[phase1_species[n]];
        }
        if (!std::isfinite(reactant_minus_product_cp))
            Util::Exception(INFO, "Non-finite phase heat-capacity difference");
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::GpuTuple<Set::Scalar,Set::Scalar,Set::Scalar,Set::Scalar>
    ApplyKineticChange(
        const amrex::Array4<Set::Scalar>& component_density,
        const State& state, const Set::Scalar heat_capacity,
        const Set::Scalar surface_measure,
        Set::Scalar& temperature, Set::Scalar& integrated_heat,
        const int i, const int j, const int k) const
    {
        if (!Kinetic() || !(state.dt > 0.0) ||
            !(heat_capacity > 0.0) || !(temperature > 0.0) ||
            !(surface_measure > 0.0))
            return {0.0, 0.0, 0.0, 0.0};

        const Set::Scalar old_temperature = temperature;
        const Set::Scalar kinetic_base_temperature =
            reconstruct_surface_temperature ? state.surface_temperature : old_temperature;
        if (!(kinetic_base_temperature > 0.0) || !std::isfinite(kinetic_base_temperature))
            Util::Abort(INFO, "Missing positive reconstructed surface temperature");
        const Set::Scalar phase0_mass = component_density(i,j,k,phase0);
        if (!(phase0_mass > 0.0)) return {0.0, 0.0, 0.0, 0.0};
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

        if (reference_enthalpy &&
            (!std::isfinite(reactant_minus_product_cp) || state.caloric.Enabled()))
            Util::Abort(INFO, "reference_enthalpy requires configured constant phase heat capacities");
        const auto thermal_increment = [=] AMREX_GPU_DEVICE(
            Set::Scalar mass_change, Set::Scalar heat)
        {
            return reference_enthalpy ? ReferenceEnthalpyDeltaT(old_temperature, heat,
                heat_capacity, mass_change, reactant_minus_product_cp,
                enthalpy_reference_temperature) :
                state.caloric.DeltaT(old_temperature, heat, heat_capacity);
        };
        const auto kinetic_mass_flux =
            [=] AMREX_GPU_DEVICE(Set::Scalar mass_change)
        {
            const Set::Scalar heat = latent_heat * mass_change +
                coupled_enthalpy_change * Util::Min(mass_change, 0.0);
            const Set::Scalar delta_temperature =
                thermal_increment(mass_change, heat);
            const Set::Scalar kinetic_temperature =
                kinetic_base_temperature + delta_temperature;
            if (!(kinetic_temperature > 0.0) ||
                !(old_temperature + delta_temperature > 0.0)) return 0.0;

            Set::Scalar mass_flux = 0.0;
            if (kinetics == Kinetics::ArrheniusSurfaceFlux)
            {
                if (!(state.pressure > 0.0)) return 0.0;
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
            return mass_flux;
        };
        const auto residual = [=] AMREX_GPU_DEVICE(Set::Scalar mass_change)
        {
            return mass_change + kinetic_mass_flux(mass_change) *
                surface_measure * state.dt;
        };

        Set::Scalar lower = 0.0;
        Set::Scalar upper = 0.0;
        const Set::Scalar zero_residual = residual(0.0);
        if (zero_residual > 0.0)
        {
            // Endothermic forward transfer reduces its own kinetic rate, so
            // its zero-state increment is the natural local root bracket.
            lower = Util::Max(-phase0_mass, -zero_residual);
            if (residual(lower) > 0.0)
            {
                lower = -phase0_mass;
                if (residual(lower) > 0.0) upper = lower;
            }
        }
        else if (zero_residual < 0.0)
        {
            upper = Util::Min(reverse_mass, -zero_residual);
            if (residual(upper) < 0.0)
            {
                upper = reverse_mass;
                if (residual(upper) < 0.0) lower = upper;
            }
        }
        else
            return {0.0, 0.0, 0.0, 0.0};

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
        if (mass_change == 0.0) return {0.0, 0.0, 0.0, 0.0};

        Set::Scalar old_gas_molar_density = 0.0;
        for (int n = 0; n < ngas_species; ++n)
            old_gas_molar_density += component_density(i,j,k,n) *
                inverse_molecular_weight[n];
        component_density(i,j,k,phase0) += mass_change;
        for (int n = 0; n < nphase1_species; ++n)
            component_density(i,j,k,phase1_species[n]) -=
                phase1_mass_fractions[n] * mass_change;
        integrated_heat = latent_heat * mass_change +
            coupled_enthalpy_change * Util::Min(mass_change, 0.0);
        temperature = old_temperature + thermal_increment(mass_change, integrated_heat);

        Set::Scalar phase1_specific_volume =
            phase1_inverse_reference_density;
        if (phase1_is_gas)
            phase1_specific_volume = phase1_gas_constant *
                temperature / state.pressure;
        const Set::Scalar specific_volume_jump =
            phase1_specific_volume - 1.0 / phase0_reference_density;
        // Products use the final physical temperature above; pre-existing gas
        // also expands/contracts with the heat update. This is the exact EOS
        // volume difference, including the mass/temperature cross term once.
        const Set::Scalar thermal_volume_change = old_gas_molar_density *
            gas_constant * (temperature - old_temperature) / state.pressure;
        const Set::Scalar surface_mass_flux = -mass_change /
            (surface_measure * state.dt);
        return {mass_change, -mass_change * specific_volume_jump + thermal_volume_change,
            surface_mass_flux, surface_mass_flux * specific_volume_jump};
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar RecoilPressure(const Set::Scalar mass_flux,
        const Set::Scalar temperature, const Set::Scalar pressure) const
    {
        Set::Scalar phase1_specific_volume =
            phase1_inverse_reference_density;
        if (phase1_is_gas)
        {
            if (!(temperature > 0.0) || !(pressure > 0.0)) return 0.0;
            phase1_specific_volume =
                phase1_gas_constant * temperature / pressure;
        }
        return mass_flux * mass_flux *
            (phase1_specific_volume - 1.0 / phase0_reference_density);
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
    Set::Scalar gas_constant = 0.0;
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
    bool reconstruct_surface_temperature = false;
    bool reference_enthalpy = false;
    Set::Scalar enthalpy_reference_temperature = 300.0;
    Set::Scalar reactant_minus_product_cp = NAN;
    SurfaceTemperatureOptions surface_temperature_options{};
};
}
}

#endif
