#include "LowMach.H"

#include <cctype>
#include <cstring>
#include "AMReX_MultiFabUtil.H"
#include "AMReX_SPACE.H"
#include "AMReX_TimeIntegrator.H"
#include "Model/Chemistry/Chemistry.H"
#include "Model/PhaseField/PhaseField.H"
#include "Numeric/Advect/Advect.H"
#include "Numeric/Stencil.H"

#include "Model/Gas/Thermo/Thermo.H"
#include "Model/Gas/Thermo/CpConstant.H"
#include "Model/Gas/Transport/Transport.H"
#include "Model/Gas/Transport/Mixture_Averaged.H"
#include "Model/Gas/EOS/EOS.H"
#include "Model/Gas/EOS/CPG.H"

namespace Integrator
{
namespace
{
AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
Set::Scalar MinimumImageDisplacement(
    Set::Scalar displacement, Set::Scalar period, int periodic)
{
    if (!periodic || !(period > 0.0)) return displacement;
    return displacement - period *
        std::floor(displacement / period + 0.5);
}

// Lightweight velocity accessor used directly by the advection stencil.  It
// avoids allocating and filling an auxiliary vector MultiFab for an affine
// rigid velocity that is cheaper to evaluate than to load from memory.
struct RigidBodyVelocity
{
    Set::Vector translation = Set::Vector::Zero();
    Set::Vector rotation = Set::Vector::Zero();
    Set::Vector center = Set::Vector::Zero();
    amrex::GpuArray<Set::Scalar,AMREX_SPACEDIM> prob_lo{};
    amrex::GpuArray<Set::Scalar,AMREX_SPACEDIM> dx{};
    amrex::GpuArray<Set::Scalar,AMREX_SPACEDIM> period{};
    amrex::GpuArray<int,AMREX_SPACEDIM> periodic{};

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar operator()(int i, int j, int k, int d) const
    {
        (void)k;
        Set::Scalar value = translation(d);
#if AMREX_SPACEDIM == 2
        const Set::Scalar x = MinimumImageDisplacement(
            prob_lo[0] + (i + 0.5) * dx[0] - center(0),
            period[0], periodic[0]);
        const Set::Scalar y = MinimumImageDisplacement(
            prob_lo[1] + (j + 0.5) * dx[1] - center(1),
            period[1], periodic[1]);
        value += d == 0 ? -rotation(0) * y : rotation(0) * x;
#elif AMREX_SPACEDIM == 3
        const Set::Scalar x = MinimumImageDisplacement(
            prob_lo[0] + (i + 0.5) * dx[0] - center(0),
            period[0], periodic[0]);
        const Set::Scalar y = MinimumImageDisplacement(
            prob_lo[1] + (j + 0.5) * dx[1] - center(1),
            period[1], periodic[1]);
        const Set::Scalar z = MinimumImageDisplacement(
            prob_lo[2] + (k + 0.5) * dx[2] - center(2),
            period[2], periodic[2]);
        if (d == 0) value += rotation(1) * z - rotation(2) * y;
        if (d == 1) value += rotation(2) * x - rotation(0) * z;
        if (d == 2) value += rotation(0) * y - rotation(1) * x;
#else
        (void)i;
        (void)j;
        (void)k;
#endif
        return value;
    }
};
}

LowMach::LowMach(IO::ParmParse& pp) : LowMach()
{
    pp_queryclass(*this);
}

void
LowMach::Parse(LowMach& value, IO::ParmParse& pp)
{
#if !AMREX_DEVICE_COMPILE
    BL_PROFILE("Integrator::LowMach::Parse");

    // Plotfiles retain the current RK state, not the previous-stage buffers.
    value.synchronize_restart_state = !value.restart_file_cell.empty();

    pp.query_required("cfl", value.cfl);
    pp.query_default(
        "phase_field.cfl", value.phase_field_cfl, value.cfl);
    if (!(value.phase_field_cfl > 0.0))
        Util::Exception(INFO, "phase_field.cfl must be positive");
    pp.query_default("cfl_v", value.cfl_v, 1.0e100);
    pp.query_default("small", value.small, 1.0e-12);
    pp.query_default("density_floor", value.density_floor,
        "1.0e-12_kg/m^3", Unit::Density());
    pp.query_default("pressure_floor", value.pressure_floor,
        "1.0e-12_Pa", Unit::Pressure());
    pp.query_default("pressure_scale", value.pressure_scale, 1.0);
    pp.query_default("projection.enabled", value.projection_enabled, true);
    pp.query_default("diagnostics.interval", value.diagnostics_interval, 0);
    pp.query_default("diagnostics.extended_fields", value.diagnostics_extended_fields, false);
    if (value.projection_enabled)
    {
        pp.query_default("projection.update_pressure", value.projection_update_pressure, true);
        pp.queryclass("projection", value.pressure_poisson);
    }
    pp.query_default("include_viscosity", value.include_viscosity, true);
    pp.query_default("implicit_viscosity", value.implicit_momentum_diffusion, false);
    pp.query_default("include_conduction", value.include_conduction, true);
    pp.query_default("advect_temperature", value.advect_temperature, true);
    if (value.implicit_momentum_diffusion && !value.include_viscosity)
        Util::Exception(INFO,
            "implicit_viscosity requires include_viscosity=1");

    pp.queryclass<Model::Gas::Gas>("gas", value.gas);
    value.gas_device_data =
        value.gas.GetDeviceData(value.gas_device_storage);
    value.ngas_species = value.gas.nspecies;
    if (value.ngas_species > Model::Chemistry::MAX_SPECIES)
        Util::Exception(INFO, "LowMach chemistry supports at most ",
                        Model::Chemistry::MAX_SPECIES, " gas species");

    pp.queryarr_required("species.names", value.species_names);
    value.nspecies = value.species_names.size();
    if (value.nspecies < value.ngas_species)
        Util::Exception(INFO, "species.names must contain one identifier for every gas species");
    value.component_density_ic.resize(value.nspecies, nullptr);
    value.reference_density.assign(value.nspecies, NAN);
    value.condensed_specific_heat.fill(NAN);
    value.condensed_thermal_conductivity.fill(NAN);
    value.condensed_inverse_reference_density.fill(NAN);
    value.condensed_temperature_override.fill(-1.0);
    value.condensed_dynamic_viscosity.fill(NAN);
    value.liquid_inverse_reference_density.fill(0.0);
    if (value.nspecies > Model::Chemistry::MAX_SPECIES)
        Util::Exception(INFO, "LowMach supports at most ",
                        Model::Chemistry::MAX_SPECIES, " total species");
    for (int n = 0; n < value.nspecies; ++n)
    {
        const std::string& name = value.species_names[n];
        if (name.empty() ||
            !(std::isalpha(static_cast<unsigned char>(name[0])) || name[0] == '_'))
            Util::Exception(INFO, name, " is not a valid species identifier");
        for (const char c : name)
            if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_'))
                Util::Exception(INFO, name, " is not a valid species identifier");
        for (int m = 0; m < n; ++m)
            if (name == value.species_names[m])
                Util::Exception(INFO, "Duplicate species identifier ", name);

        std::string mechanics;
        pp.query_required(name + ".mechanics", mechanics);
        if (mechanics == "fluid")
        {
            if (n >= value.ngas_species)
                Util::Exception(INFO, "Fluid species ", name,
                                " has no corresponding entry in gas.mw");
        }
        else if (mechanics == "deformable_solid")
        {
            if (value.deformable_solid_species >= 0)
                Util::Exception(INFO,
                                "LowMach currently supports one deformable solid species");
            value.deformable_solid_species = n;
        }
        else if (mechanics == "rigid_solid")
            value.rigid_solid_species.push_back(n);
        else if (mechanics == "liquid")
        {
            value.liquid_species.push_back(n);
            pp.query_required(name + ".reference_density",
                            value.reference_density[n], Unit::Density());
            pp.query_required(name + ".dynamic_viscosity",
                            value.condensed_dynamic_viscosity[n],
                            Unit::Pressure() * Unit::Time());
            if (!(value.reference_density[n] > 0.0) ||
                !(value.condensed_dynamic_viscosity[n] > 0.0))
                Util::Exception(INFO, name,
                    " liquid density and dynamic viscosity must be positive");
            value.interfacial_reference_density_min = Util::Min(
                value.interfacial_reference_density_min,
                value.reference_density[n]);
        }
        else
            Util::Exception(INFO, mechanics,
                            " is not a valid mechanics type for species ", name);
        if (n < value.ngas_species && mechanics != "fluid")
            Util::Exception(INFO, "The first ", value.ngas_species,
                            " species must be fluid species described by the gas model");

        if (mechanics != "fluid")
        {
            pp.query_default(name + ".temperature_override",
                value.condensed_temperature_override[n], "-1.0",
                Unit::Temperature());
            if (value.condensed_temperature_override[n] == 0.0)
                Util::Exception(INFO, name,
                    ".temperature_override must be negative or positive");

            const bool has_specific_heat =
                pp.contains(name + ".specific_heat");
            const bool has_thermal_conductivity =
                pp.contains(name + ".thermal_conductivity");
            if (has_specific_heat != has_thermal_conductivity)
                Util::Exception(INFO, name,
                    " requires both specific_heat and thermal_conductivity");
            if (has_specific_heat)
            {
                pp.query_required(name + ".specific_heat",
                    value.condensed_specific_heat[n],
                    Unit::SpecificHeatCapacity());
                pp.query_required(name + ".thermal_conductivity",
                    value.condensed_thermal_conductivity[n],
                    Unit::ThermalConductivity());
                if (!(value.condensed_specific_heat[n] > 0.0) ||
                    !(value.condensed_thermal_conductivity[n] > 0.0))
                    Util::Exception(INFO, name,
                        " thermal properties must be positive");
                value.condensed_thermal_transport = true;
            }
        }

        if (pp.contains(name + ".density.ic.type"))
            pp.select<IC::Constant,IC::Expression,IC::PNG,IC::PSRead>(
                name + ".density.ic", value.component_density_ic[n],
                pp.forward_args(value.geom, Unit::Density()));
    }

    std::string chemistry_timestep_mode;
    pp.query_default(
        "chemistry.timestep.mode", chemistry_timestep_mode, "off");
    if (chemistry_timestep_mode == "off")
        value.chemistry_timestep_mode = ChemistryTimestepMode::Off;
    else if (chemistry_timestep_mode == "report")
        value.chemistry_timestep_mode = ChemistryTimestepMode::Report;
    else if (chemistry_timestep_mode == "limit")
        value.chemistry_timestep_mode = ChemistryTimestepMode::Limit;
    else
        Util::Exception(INFO, "chemistry.timestep.mode must be off, "
            "report, or limit");
    pp.query_default("chemistry.timestep.max_fractional_change",
        value.chemistry_max_fractional_change, 0.1);
    pp.query_default("chemistry.timestep.reactant_mass_fraction_floor",
        value.chemistry_reactant_mass_fraction_floor, 1.0e-4);
    if (!(value.chemistry_max_fractional_change > 0.0) ||
        value.chemistry_max_fractional_change > 1.0)
        Util::Exception(INFO,
            "chemistry.timestep.max_fractional_change must be in (0,1]");
    if (!(value.chemistry_reactant_mass_fraction_floor > 0.0) ||
        value.chemistry_reactant_mass_fraction_floor >= 1.0)
        Util::Exception(INFO,
            "chemistry.timestep.reactant_mass_fraction_floor must be in "
            "(0,1)");

    // Consume the timestep-control keys before the nested chemistry parser
    // performs its strict unused-input check.
    value.chemistry.Define(value.ngas_species);
    pp.queryclass("chemistry", value.chemistry);
    if (value.chemistry.Split() && !value.projection_enabled)
        Util::Exception(INFO, "Locally integrated LowMach chemistry requires projection.enabled=1");
    if (value.chemistry_timestep_mode != ChemistryTimestepMode::Off &&
        !value.chemistry.Reactive())
        Util::Exception(INFO, "chemistry timestep reporting requires a "
            "reactive chemistry model");
    if (value.chemistry_timestep_mode != ChemistryTimestepMode::Off &&
        !value.dynamictimestep.on)
        Util::Exception(INFO, "chemistry timestep reporting requires "
            "dynamictimestep.on=1");

    value.implicit_thermal_diffusion = value.include_conduction;
    value.implicit_species_diffusion = value.ngas_species > 1 &&
        value.gas.transport.CommonDiffusivity();
    if (value.deformable_solid_species >= 0)
    {
        pp.queryclass("reference_map", value.reference_map_reconstruction);

        std::string solid_model_type;
        pp.query_required("solid.model.type", solid_model_type);
        if (solid_model_type != "finite.neohookean")
            Util::Exception(INFO, solid_model_type,
                            " is not a valid deformable solid model for LowMach");
        pp.query_default("solid.model.eta_threshold", value.finite_solid_eta_threshold, 0.5);
        pp.query_default("solid.model.J_floor",       value.finite_solid_J_floor, 1.0e-6);
        pp.query_default("solid.model.viscosity", value.finite_solid_viscosity,
                        "0.0", Unit::Pressure() * Unit::Time());
        pp.query_default("solid.model.interface_viscosity", value.finite_solid_interface_viscosity,
                        "0.0", Unit::Pressure() * Unit::Time());
        pp.query_default("solid.model.deviatoric_stress_divergence_sign",
                        value.finite_solid_deviatoric_stress_divergence_sign, 0.0);
        pp.query_required("solid.model.reference_density",
                        value.finite_solid_reference_density, Unit::Density());
        if (value.finite_solid_reference_density <= 0.0)
            Util::Exception(INFO, "solid.model.reference_density must be positive");
        value.reference_density[value.deformable_solid_species] =
            value.finite_solid_reference_density;
        pp.queryclass<Model::Solid::Finite::NeoHookean>("solid.model.finite.neohookean", value.finite_solid_model);
    }
    if (!value.rigid_solid_species.empty())
    {
        if (!value.projection_enabled)
            Util::Exception(INFO, "Rigid solid mechanics requires projection.enabled=1");
        const int nrigid = static_cast<int>(value.rigid_solid_species.size());
        value.rigid_relaxation_time.resize(nrigid);
        value.rigid_velocity.resize(nrigid, Set::Vector::Zero());
        value.rigid_anderson_coupling.resize(nrigid, true);
        value.rigid_coupling_max_iterations.resize(nrigid, 3);
        value.rigid_coupling_relative_tolerance.resize(nrigid, 1.0e-3);
        value.rigid_coupling_absolute_tolerance.resize(nrigid, 1.0e-6);
        for (int m = 0;
             m < nrigid; ++m)
        {
            const int n = value.rigid_solid_species[m];
            const std::string prefix = value.species_names[n] + ".rigid";
            pp.query_required(value.species_names[n] + ".reference_density",
                            value.reference_density[n], Unit::Density());
            if (value.reference_density[n] <= 0.0)
                Util::Exception(INFO, value.species_names[n],
                                ".reference_density must be positive");

            pp.query_required(prefix + ".relaxation_time",
                            value.rigid_relaxation_time[m], Unit::Time());
            if (!(value.rigid_relaxation_time[m] > 0.0))
                Util::Exception(INFO, prefix,
                    ".relaxation_time must be positive");

            std::string motion;
            pp.query_default(prefix + ".motion", motion, "fixed");
            if (motion == "fixed")
            {
                pp.queryarr_required(prefix + ".velocity",
                                    value.rigid_velocity[m], Unit::Velocity());
                value.fixed_rigid_solid_species.push_back(n);
            }
            else if (motion == "free")
            {
                value.free_rigid_solid_species.push_back(n);
                value.free_rigid_solid_components.push_back(m);

                std::string coupling;
                pp.query_default(prefix + ".coupling", coupling, "anderson");
                if (coupling == "anderson")
                    value.rigid_anderson_coupling[m] = true;
                else if (coupling == "picard")
                    value.rigid_anderson_coupling[m] = false;
                else
                    Util::Exception(INFO, prefix,
                        ".coupling must be anderson or picard");

                pp.query_default(prefix + ".max_iterations",
                                value.rigid_coupling_max_iterations[m], 3);
                pp.query_default(prefix + ".relative_tolerance",
                                value.rigid_coupling_relative_tolerance[m],
                                1.0e-3);
                pp.query_default(prefix + ".absolute_tolerance",
                                value.rigid_coupling_absolute_tolerance[m],
                                "1.0e-6_m/s", Unit::Velocity());
                if (value.rigid_coupling_max_iterations[m] < 1)
                    Util::Exception(INFO, prefix,
                        ".max_iterations must be at least one");
                if (value.rigid_coupling_relative_tolerance[m] < 0.0 ||
                    value.rigid_coupling_absolute_tolerance[m] < 0.0 ||
                    !(value.rigid_coupling_relative_tolerance[m] > 0.0 ||
                      value.rigid_coupling_absolute_tolerance[m] > 0.0))
                    Util::Exception(INFO, prefix,
                        " coupling tolerances must be nonnegative and at "
                        "least one must be positive");
            }
            else
                Util::Exception(INFO, prefix,
                    ".motion must be fixed or free");
        }
        const int nfree = static_cast<int>(
            value.free_rigid_solid_species.size());
        value.last_free_rigid_coupling_iterations.assign(nfree, 0);
        value.last_free_rigid_coupling_residual.assign(nfree, 0.0);
        value.last_free_rigid_coupling_converged.assign(nfree, false);
        value.free_rigid_body_time.assign(nfree, NAN);
        value.free_rigid_bodies.resize(nfree);
        value.previous_free_rigid_projection_increments.resize(nfree);
    }
    // The same diffuse thickness defines capillary profiles, conservative
    // interface transport, and the physical support of kinetic interfacial
    // sources.  Parse it independently of the liquid model so gas--solid
    // phase change uses a physical length rather than a grid-cell stencil.
    Set::Scalar interface_thickness = 0.0;
    if (pp.contains("interface.thickness"))
    {
        pp.query_required("interface.thickness",
                          interface_thickness, Unit::Length());
        if (!(interface_thickness > 0.0))
            Util::Exception(INFO,
                "interface.thickness must be positive");
        value.interfacial_thickness = interface_thickness;
    }

    if (!value.liquid_species.empty())
    {
        pp.query_default("interface.enabled",
                         value.interfacial_forces_enabled, true);
        bool include_liquid_solid_interfaces = false;
        pp.query_default("interface.liquid_solid.enabled",
                         include_liquid_solid_interfaces, false);
        value.liquid_solid_interface_enabled =
            value.interfacial_forces_enabled &&
            include_liquid_solid_interfaces;

        std::string interfacial_transport;
        pp.query_default("interface.transport.type",
                         interfacial_transport, "advected");
        if (interfacial_transport == "advected")
            value.interfacial_transport =
                InterfacialTransport::Advected;
        else if (interfacial_transport == "conservative_allen_cahn")
        {
            value.interfacial_transport =
                InterfacialTransport::ConservativeAllenCahn;
            pp.query_required("interface.transport.diffusivity",
                value.conservative_allen_cahn_diffusivity,
                Unit::Area() / Unit::Time());
            if (!(value.conservative_allen_cahn_diffusivity > 0.0))
                Util::Exception(INFO,
                    "interface.transport.diffusivity must be positive");
        }
        else
            Util::Exception(INFO, "interface.transport.type must be "
                "advected or conservative_allen_cahn");

        if (value.interfacial_forces_enabled ||
            value.interfacial_transport ==
                InterfacialTransport::ConservativeAllenCahn)
        {
            if (!(interface_thickness > 0.0))
                Util::Exception(INFO, "interface.thickness is required by "
                    "capillary forces and conservative interface transport");
        }

        // Every solid is represented explicitly in the capillary simplex so
        // it is never folded into the aggregate gas phase.  Only interfaces
        // involving at least one liquid carry an interfacial energy below;
        // solid-solid and solid-gas surface energies are intentionally absent.
        for (int species = value.ngas_species;
             species < value.nspecies; ++species)
        {
            bool is_solid = species == value.deformable_solid_species;
            for (const int n : value.rigid_solid_species)
                is_solid = is_solid || species == n;
            if (!is_solid) continue;
            value.interfacial_solid_species.push_back(species);
            if (value.liquid_solid_interface_enabled)
                value.interfacial_reference_density_min = Util::Min(
                    value.interfacial_reference_density_min,
                    value.reference_density[species]);
        }

        std::vector<std::string> phase_names;
        for (const int n : value.liquid_species)
            phase_names.push_back(value.species_names[n]);
        for (const int n : value.interfacial_solid_species)
            phase_names.push_back(value.species_names[n]);
        phase_names.push_back("gas");
        const int nliquid = static_cast<int>(value.liquid_species.size());
        const int nphase = static_cast<int>(phase_names.size());
        std::vector<Set::Scalar> surface_tension(nphase * nphase, 0.0);
        std::vector<Set::Scalar> solid_interface_regularization(
            nphase * nphase, 0.0);
        std::vector<Set::Scalar> solid_surface_energy_difference(
            nphase * nphase, 0.0);
        for (int a = 0; a < nliquid; ++a)
            for (int b = a + 1; b < nphase; ++b)
            {
                const bool solid_pair =
                    b >= nliquid && b < nphase - 1;
                if (solid_pair) continue;
                const std::string key = "interface.surface_tension." +
                    phase_names[a] + "_" + phase_names[b];
                if (!value.interfacial_forces_enabled)
                {
                    pp.ignore(key);
                    continue;
                }
                Set::Scalar sigma = NAN;
                pp.query_required(key, sigma,
                    Unit::Energy() / Unit::Area());
                if (!(sigma > 0.0))
                    Util::Exception(INFO, "liquid surface tension for ",
                        phase_names[a], "_", phase_names[b],
                        " must be positive");
                surface_tension[a * nphase + b] = sigma;
                surface_tension[b * nphase + a] = sigma;
            }
        for (int a = 0; a < nliquid; ++a)
            for (int b = nliquid; b < nphase - 1; ++b)
            {
                const std::string pair =
                    phase_names[a] + "_" + phase_names[b];
                const std::string regularization_key =
                    "interface.liquid_solid.regularization." + pair;
                const std::string contact_angle_key =
                    "interface.liquid_solid.contact_angle." + pair;
                const std::string surface_difference_key =
                    "interface.liquid_solid.surface_energy_difference." + pair;
                if (!value.liquid_solid_interface_enabled)
                {
                    pp.ignore(regularization_key);
                    pp.ignore(contact_angle_key);
                    pp.ignore(surface_difference_key);
                    continue;
                }

                Set::Scalar regularization = NAN;
                pp.query_required(regularization_key, regularization,
                    Unit::Energy() / Unit::Area());
                if (!(regularization > 0.0))
                    Util::Exception(INFO,
                        "liquid-solid interface regularization for ", pair,
                        " must be positive");

                const bool has_contact_angle =
                    pp.contains(contact_angle_key);
                const bool has_surface_difference =
                    pp.contains(surface_difference_key);
                if (has_contact_angle == has_surface_difference)
                    Util::Exception(INFO, "liquid-solid pair ", pair,
                        " requires exactly one of ", contact_angle_key,
                        " or ", surface_difference_key);

                Set::Scalar surface_difference = NAN;
                if (has_contact_angle)
                {
                    Set::Scalar contact_angle = NAN;
                    pp.query_required(contact_angle_key, contact_angle,
                                      Unit::Angle());
                    if (!std::isfinite(contact_angle) ||
                        contact_angle < 0.0 ||
                        contact_angle > Set::Constant::Pi)
                        Util::Exception(INFO, "contact angle for ", pair,
                            " must be between 0 and 180 degrees");
                    const int gas = nphase - 1;
                    surface_difference =
                        Model::PhaseField::MultiphaseInterface::
                        SurfaceEnergyDifferenceFromContactAngle(
                            surface_tension[a * nphase + gas],
                            contact_angle);
                }
                else
                {
                    pp.query_required(surface_difference_key,
                        surface_difference,
                        Unit::Energy() / Unit::Area());
                    if (!std::isfinite(surface_difference))
                        Util::Exception(INFO,
                            "solid surface energy difference for ", pair,
                            " must be finite");
                }

                solid_interface_regularization[a * nphase + b] =
                    regularization;
                solid_interface_regularization[b * nphase + a] =
                    regularization;
                solid_surface_energy_difference[a * nphase + b] =
                    surface_difference;
                solid_surface_energy_difference[b * nphase + a] =
                    surface_difference;
            }
        if (value.interfacial_forces_enabled)
            value.interfacial_model.Define(
                phase_names, nliquid,
                value.liquid_solid_interface_enabled,
                interface_thickness, surface_tension,
                solid_interface_regularization,
                solid_surface_energy_difference);
        for (const int n : value.liquid_species)
            value.liquid_inverse_reference_density[n] =
                1.0 / value.reference_density[n];
    }

    for (int n = value.ngas_species; n < value.nspecies; ++n)
    {
        if (!(value.reference_density[n] > 0.0))
            Util::Exception(INFO, value.species_names[n],
                " requires a positive reference density");
        value.condensed_inverse_reference_density[n] =
            1.0 / value.reference_density[n];
        if (value.condensed_thermal_transport)
        {
            if (!(value.condensed_specific_heat[n] > 0.0) ||
                !(value.condensed_thermal_conductivity[n] > 0.0))
                Util::Exception(INFO, value.species_names[n],
                    " requires specific_heat and thermal_conductivity");
        }
    }

    std::vector<std::string> mechanism_names;
    pp.queryarr_default("mechanisms.names", mechanism_names, {});
    value.mechanisms.resize(mechanism_names.size());
    for (int n = 0; n < static_cast<int>(mechanism_names.size()); ++n)
    {
        const std::string& id = mechanism_names[n];
        if (id.empty() ||
            !(std::isalpha(static_cast<unsigned char>(id[0])) || id[0] == '_'))
            Util::Exception(INFO, id, " is not a valid mechanism identifier");
        for (const char c : id)
            if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_'))
                Util::Exception(INFO, id, " is not a valid mechanism identifier");
        for (int m = 0; m < n; ++m)
            if (id == mechanism_names[m])
                Util::Exception(INFO, "Duplicate mechanism identifier ", id);
        pp.select<Model::Mechanism::PhaseChange,
                Model::Mechanism::InterphaseReaction,
                Model::Mechanism::InterfacialHeatSource>(
            id, value.mechanisms[n],
            pp.forward_args(value.species_names, value.ngas_species,
                            value.rigid_solid_species,
                            value.liquid_species,
                            value.reference_density, value.gas.MW,
                            value.gas.Rg));
        value.has_equilibrium_phase_change =
            value.has_equilibrium_phase_change ||
            value.mechanisms[n].Equilibrium();
        value.has_kinetic_phase_change =
            value.has_kinetic_phase_change ||
            value.mechanisms[n].Kinetic();
    }
    value.has_split_phase_change = value.has_equilibrium_phase_change ||
        value.has_kinetic_phase_change;
    if (value.has_kinetic_phase_change &&
        !(value.interfacial_thickness > 0.0))
        Util::Exception(INFO, "interface.thickness is required by kinetic "
            "diffuse-interface phase change");
    if (value.implicit_momentum_diffusion || value.implicit_thermal_diffusion ||
        value.implicit_species_diffusion)
        pp.queryclass("diffusion", value.diffusion);

    pp.select<Numeric::Advect::MUSCL,
            Numeric::Advect::Upwind,
            Numeric::Advect::Centered,
            Numeric::Advect::QUICK,
            Numeric::Advect::WENO5>("advection",value.advect);
    if (value.advect.PhiLocation() != Set::HC::Cell ||
        value.advect.VelocityLocation() != Set::HC::Cell)
        Util::Exception(INFO, "LowMach currently requires cell-centered phi and velocity advection data");

    pp.query_default("velocity_refinement_criterion",
        value.velocity_refinement_criterion, "1.0e100_m/s",
        Unit::Velocity());
    pp.query_default("pressure_refinement_criterion",
        value.pressure_refinement_criterion, "1.0e100_Pa",
        Unit::Pressure());
    pp.query_default("temperature_refinement_criterion",
        value.temperature_refinement_criterion, "1.0e100_K",
        Unit::Temperature());
    pp.query_default("reaction_refinement_criterion",
                    value.reaction_refinement_criterion,
                    "1.0e100_1/s", 1.0 / Unit::Time());
    if (value.deformable_solid_species >= 0 ||
        !value.rigid_solid_species.empty() || !value.liquid_species.empty())
        pp.query_default("eta_refinement_criterion", value.eta_refinement_criterion, 1.0e100);
    pp.query_default("amr.reinitialize_condensed_composition",
                    value.reinitialize_condensed_composition, false);
    if (value.reinitialize_condensed_composition)
    {
        pp.query_default("amr.reinitialize_condensed_composition_eta_min",
                        value.reinitialize_condensed_composition_eta_min, 0.99);
        if (value.reinitialize_condensed_composition_eta_min < 0.0 ||
            value.reinitialize_condensed_composition_eta_min > 1.0)
            Util::Exception(INFO,
                "amr.reinitialize_condensed_composition_eta_min must lie in [0,1]");
        if (value.deformable_solid_species >= 0 ||
            static_cast<int>(value.rigid_solid_species.size()) !=
                value.nspecies - value.ngas_species)
            Util::Exception(INFO,
                "amr.reinitialize_condensed_composition currently requires "
                "all condensed species to use rigid_solid mechanics");
        for (int n = value.ngas_species; n < value.nspecies; ++n)
            if (value.component_density_ic[n] == nullptr)
                Util::Exception(INFO,
                    "amr.reinitialize_condensed_composition requires a density IC for ",
                    value.species_names[n]);
    }
    pp.queryarr_default("g", value.g, "0.0_m/s^2 0.0_m/s^2 0.0_m/s^2",
        Unit::Length() / Unit::Time() / Unit::Time());

    int kinetic_diffuse_moment_nghost = 0;
    if (value.has_kinetic_phase_change)
    {
        // For eta = 0.5[1-tanh(2x/ell)], the eta=0.1--0.9 band has
        // half-width atanh(0.8)*ell/2.  Allocate enough ghosts to represent
        // that same physical support on every configured AMR level.
        constexpr Set::Scalar support_per_thickness =
            0.54930614433405484570;
        const Set::Scalar support = support_per_thickness *
            value.interfacial_thickness;
        for (int lev = 0; lev <= value.maxLevel(); ++lev)
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
                kinetic_diffuse_moment_nghost = Util::Max(
                    kinetic_diffuse_moment_nghost,
                    static_cast<int>(std::ceil(
                        support / value.geom[lev].CellSize(d))) + 1);
    }

    int nghost = value.advect.NGhost();
    if (nghost < 3) nghost = 3;
    pp.select_default<BC::Constant,BC::Expression>(
        "velocity.bc", value.velocity_bc,
        pp.forward_args(AMREX_SPACEDIM,
            Unit::Length() / Unit::Time()));
    pp.select_default<BC::Constant,BC::Expression>(
        "temperature.bc", value.temperature_bc,
        pp.forward_args(1, Unit::Temperature()));
    pp.select_default<BC::Constant::ZeroNeumann,BC::Constant,BC::Expression>(
        "component_density.bc", value.component_density_bc,
        pp.forward_args(value.nspecies, Unit::Density()));
    pp.select_default<BC::Constant,BC::Expression>(
        "pressure.bc", value.pressure_bc,
        pp.forward_args(1, Unit::Pressure()));

    pp.select_default<IC::Constant,IC::Expression>(
        "velocity.ic", value.velocity_ic,
        pp.forward_args(value.geom, Unit::Length() / Unit::Time()));
    pp.select_default<IC::Constant,IC::Expression>(
        "temperature.ic", value.temperature_ic,
        pp.forward_args(value.geom, Unit::Temperature()));
    if (pp.contains("heat_source.ic.type"))
        pp.select<IC::Constant,IC::Expression>(
            "heat_source.ic", value.heat_source_ic,
            pp.forward_args(value.geom,Unit::Power() / Unit::Volume()));
    pp.select_default<IC::Constant,IC::Expression>(
        "pressure.ic", value.pressure_ic,
        pp.forward_args(value.geom, Unit::Pressure()));
    if (value.deformable_solid_species >= 0)
    {
        pp.select_default<BC::Constant::ZeroNeumann,BC::Constant,BC::Expression>(
            "xi.bc", value.xi_bc,
            pp.forward_args(AMREX_SPACEDIM, Unit::Length()));
        pp.select_default<IC::Expression::X,IC::Constant,IC::Expression>(
            "xi.ic", value.xi_ic,
            pp.forward_args(value.geom, Unit::Length()));
    }

    std::vector<std::string> species_suffix(value.nspecies);
    for (int n = 0; n < value.nspecies; ++n)
        species_suffix[n] = "_" + value.species_names[n];
    std::vector<std::string> gas_species_suffix(
        species_suffix.begin(), species_suffix.begin() + value.ngas_species);
    std::vector<std::string> vector_suffix = {"x", "y", "z"};
    vector_suffix.resize(AMREX_SPACEDIM);

    value.AddField<Set::Scalar,Set::HC::Cell>(value.velocity_mf,        value.velocity_bc,      AMREX_SPACEDIM, nghost, "velocity",          true,  true, vector_suffix);
    value.AddField<Set::Scalar,Set::HC::Cell>(value.velocity_old_mf,    value.velocity_bc,      AMREX_SPACEDIM, nghost, "velocity_old",      false, true, vector_suffix);
    value.AddField<Set::Scalar,Set::HC::Cell>(value.temperature_mf,       value.temperature_bc,   1,              nghost, "temperature",       true,  true);
    value.AddField<Set::Scalar,Set::HC::Cell>(value.temperature_old_mf,   value.temperature_bc,   1,              nghost, "temperature_old",   false, true);
    if (value.heat_source_ic)
        value.AddField<Set::Scalar,Set::HC::Cell>(
            value.heat_source_mf, &value.bc_nothing, 1, 0,
            "heat_source", false, false);
    value.AddField<Set::Scalar,Set::HC::Cell>(value.component_density_mf,     value.component_density_bc, value.nspecies, nghost, "component_density",     true,  true, species_suffix);
    value.AddField<Set::Scalar,Set::HC::Cell>(value.component_density_old_mf, value.component_density_bc, value.nspecies, nghost, "component_density_old", false, true, species_suffix);
    if (value.deformable_solid_species >= 0)
    {
        value.AddField<Set::Scalar,Set::HC::Cell>(value.eta_mf, &value.bc_nothing, 1, nghost, "eta", true, false);
        value.AddField<Set::Scalar,Set::HC::Cell>(value.xi_mf,     value.xi_bc, AMREX_SPACEDIM, nghost, "xi",     true,  true, vector_suffix);
        value.AddField<Set::Scalar,Set::HC::Cell>(value.xi_old_mf, value.xi_bc, AMREX_SPACEDIM, nghost, "xi_old", false, true, vector_suffix);
    }
    if (!value.rigid_solid_species.empty())
    {
        std::vector<std::string> rigid_species_suffix;
        for (const int n : value.rigid_solid_species)
            rigid_species_suffix.push_back(species_suffix[n]);
        value.AddField<Set::Scalar,Set::HC::Cell>(value.rigid_eta_mf, &value.bc_nothing, 1, 1, "rigid_eta", true, false);
        value.AddField<Set::Scalar,Set::HC::Cell>(
            value.rigid_species_eta_mf, &value.bc_nothing,
            value.rigid_solid_species.size(), 1, "rigid_species_eta",
            true, false, rigid_species_suffix);
    }
    if (!value.fixed_rigid_solid_species.empty())
        value.AddField<Set::Scalar,Set::HC::Cell>(value.fixed_rigid_eta_mf, &value.bc_nothing, 1, 1, "fixed_rigid_eta", true, false);
    if (!value.liquid_species.empty())
    {
        std::vector<std::string> interfacial_phase_suffix;
        std::vector<std::string> liquid_species_suffix;
        for (const int n : value.liquid_species)
        {
            interfacial_phase_suffix.push_back(
                "_liquid_" + value.species_names[n]);
            liquid_species_suffix.push_back(species_suffix[n]);
        }
        for (const int n : value.interfacial_solid_species)
            interfacial_phase_suffix.push_back(
                "_solid_" + value.species_names[n]);
        interfacial_phase_suffix.push_back("_gas");
        const int nphase = static_cast<int>(value.liquid_species.size() +
            value.interfacial_solid_species.size() + 1);
        value.AddField<Set::Scalar,Set::HC::Cell>(
            value.interfacial_volume_fraction_mf,
            &value.bc_nothing, nphase, nghost,
            "interfacial_volume_fraction_internal", false, false,
            interfacial_phase_suffix);
        // Mechanical liquid volume fractions are reconstructed directly from
        // partial densities, just like rigid_species_eta.  The capillary
        // simplex below may normalize its private copy, but that must not
        // change the diagnostic/mechanical liquid occupancy.
        value.AddField<Set::Scalar,Set::HC::Cell>(
            value.liquid_species_eta_mf, &value.bc_nothing,
            value.liquid_species.size(), nghost, "liquid_species_eta",
            true, false, liquid_species_suffix);
        // Solid components of the capillary simplex duplicate eta_mf or
        // rigid_species_eta_mf.  Liquid components are written by
        // liquid_species_eta_mf.  Only gas remains to be copied out of the
        // normalized capillary simplex.
        value.AddField<Set::Scalar,Set::HC::Cell>(
            value.gas_volume_fraction_mf,
            &value.bc_nothing, 1, 0, "gas_eta", true, false);
        if (value.interfacial_forces_enabled)
            value.AddField<Set::Scalar,Set::HC::Cell>(
                value.interfacial_chemical_potential_mf, &value.bc_nothing,
                nphase, 1, "chemical_potential", true, false,
                interfacial_phase_suffix);
    }

    value.AddField<Set::Scalar,Set::HC::Cell>(value.density_mf,             &value.bc_nothing, 1,              1,      "density",             true,  false);
    value.AddField<Set::Scalar,Set::HC::Cell>(value.pressure_mf,            value.pressure_bc, 1,              nghost, "pressure",            true,  true);
    value.AddField<Set::Scalar,Set::HC::Cell>(value.pressure_correction_mf, &value.bc_nothing, 1,              nghost, "pressure_correction", true, false);
    if (value.chemistry.Split())
        value.AddField<Set::Scalar,Set::HC::Cell>(value.chemistry_dilatation_mf,
            &value.bc_nothing, 1, 0, "chemistry_dilatation", false, false);
    if (value.has_split_phase_change)
    {
        value.AddField<Set::Scalar,Set::HC::Cell>(
            value.phase_change_dilatation_mf, &value.bc_nothing, 1, 0,
            "phase_change_dilatation", false, false);
        value.AddField<Set::Scalar,Set::HC::Cell>(
            value.phase_change_heat_mf, &value.bc_nothing, 1, 0,
            "phase_change_heat", false, false);
    }
    if (value.has_kinetic_phase_change)
    {
        // Compact diffuse moments are reused by each configured kinetic
        // mechanism.  Ghost depth is derived from the physical interface
        // thickness and the finest configured AMR spacing.
        value.AddField<Set::Scalar,Set::HC::Cell>(
            value.kinetic_diffuse_moment_mf, &value.bc_nothing,
            4 + AMREX_SPACEDIM,
            kinetic_diffuse_moment_nghost,
            "kinetic_diffuse_moment", false, false);
    }
    if (value.implicit_thermal_diffusion || value.implicit_species_diffusion)
        value.AddField<Set::Scalar,Set::HC::Cell>(value.diffusion_dilatation_mf,
            &value.bc_nothing, 1, 0, "diffusion_dilatation", false, false);
    if (value.deformable_solid_species >= 0)
        value.AddField<Set::Matrix,Set::HC::Cell>(value.solid_deviatoric_stress_mf, nullptr, 1, 1, "solid_deviatoric_stress", true, false);
    if (value.diagnostics_extended_fields)
    {
        value.AddField<Set::Scalar,Set::HC::Cell>(value.mass_fraction_mf, &value.bc_nothing, value.ngas_species, 1, "mass_fraction", true, false, gas_species_suffix);
        value.AddField<Set::Scalar,Set::HC::Cell>(value.mole_fraction_mf, &value.bc_nothing, value.ngas_species, 1, "mole_fraction", true, false, gas_species_suffix);
        value.AddField<Set::Scalar,Set::HC::Cell>(value.momentum_mf, &value.bc_nothing, AMREX_SPACEDIM, 1, "momentum", true, false, vector_suffix);
        value.AddField<Set::Scalar,Set::HC::Cell>(value.energy_mf, &value.bc_nothing, 1, 1, "energy", true, false);
        value.AddField<Set::Scalar,Set::HC::Cell>(value.vorticity_mf, &value.bc_nothing, 1, 1, "vorticity", true, false);
        value.AddField<Set::Scalar,Set::HC::Cell>(value.viscosity_mf, &value.bc_nothing, 1, 1, "viscosity", true, false);
        value.AddField<Set::Scalar,Set::HC::Cell>(value.thermal_conductivity_coeff_mf, &value.bc_nothing, 1, 1, "thermal_conductivity_coeff", true, false);
        value.AddField<Set::Scalar,Set::HC::Cell>(value.diffusion_coeff_mf, &value.bc_nothing, value.ngas_species, 1, "diffusion_coeff", true, false);
        value.AddField<Set::Scalar,Set::HC::Cell>(value.wdot_mf, &value.bc_nothing, value.ngas_species, 1, "wdot", true, false);
        value.AddField<Set::Scalar,Set::HC::Cell>(value.qdot_mf, &value.bc_nothing, 1, 1, "qdot", true, false);
        if (value.deformable_solid_species >= 0)
            value.AddField<Set::Matrix,Set::HC::Cell>(value.deformation_gradient_mf, nullptr, 1, 1, "F", true, false);
    }

    const bool gross_model_chemistry =
        std::strcmp(
            value.chemistry.model_name(),
            Model::Chemistry::GrossModel::name) == 0;
#ifdef ALAMO_GPU
    if (value.chemistry.Reactive() && !gross_model_chemistry)
        Util::Exception(
            INFO, "CUDA LowMach currently supports frozen and gross_model chemistry");
    const Chemistry* host_chemistry = nullptr;
    const Model::Gas::Gas* host_gas = nullptr;
#else
    const Chemistry* host_chemistry = &value.chemistry;
    const Model::Gas::Gas* host_gas = &value.gas;
#endif
    const Set::Scalar* specific_heat =
        value.condensed_specific_heat.data();
    const Set::Scalar* thermal_conductivity =
        value.condensed_thermal_conductivity.data();
    const Set::Scalar* inverse_reference_density =
        value.condensed_inverse_reference_density.data();
    const Set::Scalar* liquid_inverse_reference_density =
        value.liquid_inverse_reference_density.data();
    const Set::Scalar* dynamic_viscosity =
        value.condensed_dynamic_viscosity.data();
#ifdef ALAMO_GPU
    std::array<Set::Scalar, 5 * Model::Chemistry::MAX_SPECIES>
        thermal_host{};
    for (int n = 0; n < Model::Chemistry::MAX_SPECIES; ++n)
    {
        thermal_host[n] = value.condensed_specific_heat[n];
        thermal_host[Model::Chemistry::MAX_SPECIES + n] =
            value.condensed_thermal_conductivity[n];
        thermal_host[2 * Model::Chemistry::MAX_SPECIES + n] =
            value.condensed_inverse_reference_density[n];
        thermal_host[3 * Model::Chemistry::MAX_SPECIES + n] =
            value.liquid_inverse_reference_density[n];
        thermal_host[4 * Model::Chemistry::MAX_SPECIES + n] =
            value.condensed_dynamic_viscosity[n];
    }
    value.thermal_device_storage.resize(thermal_host.size());
    amrex::Gpu::copy(
        amrex::Gpu::hostToDevice, thermal_host.begin(), thermal_host.end(),
        value.thermal_device_storage.begin());
    specific_heat = value.thermal_device_storage.data();
    thermal_conductivity =
        specific_heat + Model::Chemistry::MAX_SPECIES;
    inverse_reference_density =
        specific_heat + 2 * Model::Chemistry::MAX_SPECIES;
    liquid_inverse_reference_density =
        specific_heat + 3 * Model::Chemistry::MAX_SPECIES;
    dynamic_viscosity =
        specific_heat + 4 * Model::Chemistry::MAX_SPECIES;
#endif
    value.thermal_data = {
        value.gas_device_data, value.nspecies, value.ngas_species,
        value.density_floor, value.pressure_reference,
        value.condensed_thermal_transport,
        specific_heat, thermal_conductivity, inverse_reference_density,
        liquid_inverse_reference_density, dynamic_viscosity};
    value.chemistry_device_data = {
        gross_model_chemistry,
        value.chemistry.Get<Model::Chemistry::GrossModel>(),
        value.chemistry.Solver(), host_chemistry, host_gas};
    value.thermochemical_data = {
        value.thermal_data, value.chemistry_device_data,
        value.implicit_species_diffusion, value.include_conduction,
        value.implicit_thermal_diffusion};

    bool allow_unused;
    pp.query_default("allow_unused", allow_unused, false);
    if (!allow_unused && pp.AnyUnusedInputs(true, false))
    {
        Util::Warning(INFO, "The following inputs were specified but not used:");
        pp.AllUnusedInputs();
        Util::Exception(INFO, "Aborting. Specify 'allow_unused=True` to ignore this error.");
    }
#else
    Util::IgnoreUnused(value, pp);
#endif
}

void
LowMach::UpdateSolidStress(int lev,
                            const amrex::MultiFab& u_mf,
                            const amrex::MultiFab& eta_mf,
                            const amrex::MultiFab& xi_mf,
                            bool write_diagnostics)
{
    BL_PROFILE("Integrator::LowMach::UpdateSolidStress");

    solid_deviatoric_stress_mf[lev]->setVal(Set::Matrix::Zero());
    if (write_diagnostics)
        deformation_gradient_mf[lev]->setVal(Set::Matrix::Zero());

    const auto dx = geom[lev].CellSizeArray();
    amrex::Box domain = geom[lev].Domain();
    amrex::GpuArray<int,AMREX_SPACEDIM> periodic{};
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
        periodic[d] = geom[lev].isPeriodic(d);
    const Model::Solid::Finite::NeoHookean solid_model = finite_solid_model;
    const Set::Scalar eta_threshold = Util::Clamp(finite_solid_eta_threshold, 0.0, 1.0);
    const Set::Scalar stress_eta_min =
        Util::Clamp(reference_map_reconstruction.EtaExtension(), 0.0, 1.0);
    const Set::Scalar det_floor = Util::Max(small, finite_solid_J_floor);
    const Set::Scalar solid_viscosity = finite_solid_viscosity;
    const Set::Scalar solid_interface_viscosity = finite_solid_interface_viscosity;

    amrex::MultiFab xi_stress_mf;
    const amrex::MultiFab* xi_for_stress_mf = &xi_mf;
    if (reference_map_reconstruction.SmoothsStressMap())
    {
        xi_stress_mf.define(xi_mf.boxArray(), xi_mf.DistributionMap(), AMREX_SPACEDIM, xi_mf.nGrow());
        amrex::MultiFab::Copy(xi_stress_mf, xi_mf, 0, 0, AMREX_SPACEDIM, xi_mf.nGrow());
        reference_map_reconstruction.SmoothForStress(geom[lev], eta_mf, xi_stress_mf);
        xi_for_stress_mf = &xi_stress_mf;
    }

    for (amrex::MFIter mfi(*xi_for_stress_mf, true); mfi.isValid(); ++mfi)
    {
        amrex::Box bx = mfi.growntilebox(1);
        bx &= domain;
        Set::Patch<const Set::Scalar> u = u_mf.array(mfi);
        Set::Patch<const Set::Scalar> eta = eta_mf.array(mfi);
        Set::Patch<const Set::Scalar> xi = xi_for_stress_mf->array(mfi);
        Set::Patch<Set::Matrix> solid_deviatoric_stress = solid_deviatoric_stress_mf.Patch(lev,mfi);
        Set::Patch<Set::Matrix> F_field;
        if (write_diagnostics) F_field = deformation_gradient_mf.Patch(lev,mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            auto sten = Numeric::GetStencil(i, j, k, domain, periodic);
            Set::Matrix F = Set::Matrix::Identity();
            Set::Matrix grad_u =
                Numeric::Gradient(u, i, j, k, dx.data(), sten);

            Set::Matrix solid_sigma = Set::Matrix::Zero();

            Set::Scalar eta_val = Util::Clamp(eta(i,j,k), 0.0, 1.0);
            Set::Scalar solid_weight = 0.0;
            if (eta_val > eta_threshold)
            {
                Set::Scalar denom = Util::Max(1.0 - eta_threshold, 1.0e-12);
                solid_weight = Util::SmootherStep((eta_val - eta_threshold) / denom);
            }

            if (eta_val > stress_eta_min)
            {
                Set::Matrix grad_xi =
                    Numeric::Gradient(xi, i, j, k, dx.data(), sten);
                Set::Scalar det_grad_xi = grad_xi.determinant();
                bool valid = Util::Abs(det_grad_xi) > det_floor;

                if (valid) F = grad_xi.inverse();

                Set::Scalar J = F.determinant();
                if (valid && Util::Abs(J) > det_floor)
                    solid_sigma = (solid_model.DW(F) * F.transpose()) / J;

                if (solid_viscosity != 0.0)
                    solid_sigma += solid_viscosity * (grad_u + grad_u.transpose());
            }

            Set::Matrix solid_sigma_dev = solid_sigma -
                                        (solid_sigma.trace() / Set::Scalar(AMREX_SPACEDIM)) *
                                            Set::Matrix::Identity();
            Set::Matrix strain_rate_dev = grad_u + grad_u.transpose();
            strain_rate_dev -= (strain_rate_dev.trace() / Set::Scalar(AMREX_SPACEDIM)) *
                                Set::Matrix::Identity();
            Set::Matrix interface_damping_stress =
                solid_interface_viscosity * 4.0 * eta_val * (1.0 - eta_val) * strain_rate_dev;
            solid_deviatoric_stress(i,j,k) =
                solid_weight * solid_sigma_dev + interface_damping_stress;
            if (write_diagnostics)
                F_field(i,j,k) = F;
        });
    }

    solid_deviatoric_stress_mf[lev]->FillBoundary(geom[lev].periodicity());
    if (write_diagnostics)
        deformation_gradient_mf[lev]->FillBoundary(geom[lev].periodicity());
}

//
// Set the following calculated variables:
// - eta_mf         (if there ia a deformable solid present, calculated based on densities)
// - rigid_species_eta_mf and rigid_eta_mf
//                  (if rigid solids are present, calculated based on densities)
// - liquid_species_eta_mf
//                  (if liquids are present, calculated based on densities)
// - density        (based on partial densities)
//
// If writing diagonistics, also calculate:
// - mass_fraction_mf
// - mole_fraction_mf
//
void
LowMach::UpdateComponentState(int lev, const amrex::MultiFab& component_density_mf)
{
    density_mf[lev]->setVal(0.0);
    const bool deformable_solid = deformable_solid_species >= 0;
    const bool rigid_solid = !rigid_solid_species.empty();
    const bool fixed_rigid_solid = !fixed_rigid_solid_species.empty();
    const bool liquid = !liquid_species.empty();
    const int ngas = ngas_species;
    const int number_of_species = nspecies;
    const int solid = deformable_solid_species;
    const Set::Scalar solid_reference_density =
        finite_solid_reference_density;
    const bool write_diagnostics = diagnostics_extended_fields;
    const Set::Scalar rho_floor = density_floor;
    const Set::Scalar p_reference = pressure_reference;
    const Set::Scalar* inverse_reference_density =
        amrex::get<8>(thermal_data);
    const auto gas_data = gas_device_data;
    if (deformable_solid) eta_mf[lev]->setVal(0.0);
    if (rigid_solid)
    {
        rigid_eta_mf[lev]->setVal(0.0);
        rigid_species_eta_mf[lev]->setVal(0.0);
    }
    if (fixed_rigid_solid) fixed_rigid_eta_mf[lev]->setVal(0.0);
    if (liquid)
    {
        liquid_species_eta_mf[lev]->setVal(0.0);
        interfacial_volume_fraction_mf[lev]->setVal(0.0);
    }
    if (diagnostics_extended_fields)
    {
        mass_fraction_mf[lev]->setVal(0.0);
        mole_fraction_mf[lev]->setVal(0.0);
    }

    for (amrex::MFIter mfi(component_density_mf, true); mfi.isValid(); ++mfi)
    {
        amrex::Box bx = mfi.growntilebox(1);
        Set::Patch<const Set::Scalar> component_density = component_density_mf.array(mfi);
        Set::Patch<const Set::Scalar> T = temperature_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> rho = density_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> eta = eta_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> mass_fraction = mass_fraction_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> mole_fraction = mole_fraction_mf.Patch(lev,mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Scalar gas_partial_density = 0.0;
            Set::Scalar total_partial_density = 0.0;
            for (int n = 0; n < ngas; ++n)
                gas_partial_density += component_density(i,j,k,n);
            for (int n = 0; n < number_of_species; ++n)
                total_partial_density += component_density(i,j,k,n);

            // Partial densities are the conserved species state. Mechanical
            // volume fractions are reconstructed from them when needed.
            rho(i,j,k) = total_partial_density;
            if (deformable_solid)
                eta(i,j,k) =
                    component_density(i,j,k,solid) /
                    solid_reference_density;

            if (write_diagnostics)
            {
                Set::Scalar gas_volume_fraction = 0.0;
                if (gas_partial_density > rho_floor && T(i,j,k) > 0.0 &&
                    p_reference > 0.0)
                    gas_volume_fraction = Util::Clamp(
                        gas_partial_density *
                            Model::Gas::Gas::GasConstant(
                                gas_data, component_density, i, j, k) *
                            T(i,j,k) / p_reference,
                        0.0, 1.0);
                Set::Scalar condensed_volume_fraction = 0.0;
                for (int n = ngas; n < number_of_species; ++n)
                    condensed_volume_fraction += Util::Max(
                        component_density(i,j,k,n), 0.0) *
                        inverse_reference_density[n];
                gas_volume_fraction = Util::Min(gas_volume_fraction,
                    1.0 - Model::PhaseField::H(condensed_volume_fraction));
                // Diagnostic fractions include phase occupancy; gas-model
                // calculations continue to use the conditional composition X.
                Set::Scalar moles = 0.0;
                for (int n = 0; n < ngas; ++n)
                    moles += component_density(i,j,k,n) /
                        Model::Gas::Gas::MolecularWeight(gas_data, n);
                for (int n = 0; n < ngas; ++n)
                {
                    mass_fraction(i,j,k,n) =
                        total_partial_density > rho_floor ?
                        component_density(i,j,k,n) / total_partial_density : 0.0;
                    mole_fraction(i,j,k,n) = moles > 0.0 ?
                        gas_volume_fraction * component_density(i,j,k,n) /
                            Model::Gas::Gas::MolecularWeight(gas_data, n) /
                            moles : 0.0;
                }
            }
        });
    }

    density_mf[lev]->FillBoundary(geom[lev].periodicity());
    if (deformable_solid)
        eta_mf[lev]->FillBoundary(geom[lev].periodicity());
    if (rigid_solid)
    {
        for (int m = 0; m < static_cast<int>(rigid_solid_species.size()); ++m)
        {
            const int n = rigid_solid_species[m];
            amrex::MultiFab::Copy(
                *rigid_species_eta_mf[lev], component_density_mf,
                n, m, 1, 1);
            rigid_species_eta_mf[lev]->mult(
                1.0 / reference_density[n], m, 1, 1);
            amrex::MultiFab::Saxpy(
                *rigid_eta_mf[lev], 1.0, *rigid_species_eta_mf[lev],
                m, 0, 1, 1);
        }
        rigid_species_eta_mf[lev]->FillBoundary(geom[lev].periodicity());
        rigid_eta_mf[lev]->FillBoundary(geom[lev].periodicity());
    }
    if (fixed_rigid_solid)
    {
        for (const int n : fixed_rigid_solid_species)
            amrex::MultiFab::Saxpy(*fixed_rigid_eta_mf[lev], 1.0 / reference_density[n],
                                    component_density_mf, n, 0, 1, 1);
        fixed_rigid_eta_mf[lev]->FillBoundary(geom[lev].periodicity());
    }
    if (liquid)
    {
        const int nliquid = static_cast<int>(liquid_species.size());
        const int nsolid =
            static_cast<int>(interfacial_solid_species.size());
        const int nmaterial = nliquid + nsolid;
        for (int m = 0; m < nliquid; ++m)
        {
            const int n = liquid_species[m];
            amrex::MultiFab::Copy(
                *liquid_species_eta_mf[lev], component_density_mf,
                n, m, 1, liquid_species_eta_mf[lev]->nGrow());
            liquid_species_eta_mf[lev]->mult(
                1.0 / reference_density[n], m, 1,
                liquid_species_eta_mf[lev]->nGrow());
        }
        liquid_species_eta_mf[lev]->FillBoundary(
            geom[lev].periodicity());
        amrex::MultiFab::Copy(
            *interfacial_volume_fraction_mf[lev],
            *liquid_species_eta_mf[lev], 0, 0, nliquid,
            interfacial_volume_fraction_mf[lev]->nGrow());
        for (int m = 0; m < nsolid; ++m)
        {
            const int n = interfacial_solid_species[m];
            const int phase = nliquid + m;
            amrex::MultiFab::Copy(
                *interfacial_volume_fraction_mf[lev], component_density_mf,
                n, phase, 1,
                interfacial_volume_fraction_mf[lev]->nGrow());
            interfacial_volume_fraction_mf[lev]->mult(
                1.0 / reference_density[n], phase, 1,
                interfacial_volume_fraction_mf[lev]->nGrow());
        }
        for (amrex::MFIter mfi(*interfacial_volume_fraction_mf[lev], false);
             mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.fabbox();
            Set::Patch<Set::Scalar> phase =
                interfacial_volume_fraction_mf[lev]->array(mfi);
            amrex::ParallelFor(
                bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar material_sum = 0.0;
                    for (int n = 0; n < nmaterial; ++n)
                    {
                        phase(i,j,k,n) = Util::Clamp(
                            phase(i,j,k,n), 0.0, 1.0);
                        material_sum += phase(i,j,k,n);
                    }
                    if (material_sum > 1.0)
                    {
                        const Set::Scalar scale = 1.0 / material_sum;
                        for (int n = 0; n < nmaterial; ++n)
                            phase(i,j,k,n) *= scale;
                        phase(i,j,k,nmaterial) = 0.0;
                    }
                    else
                        phase(i,j,k,nmaterial) = 1.0 - material_sum;
                });
        }
        interfacial_volume_fraction_mf[lev]->FillBoundary(
            geom[lev].periodicity());
    }
    if (diagnostics_extended_fields)
    {
        mass_fraction_mf[lev]->FillBoundary(geom[lev].periodicity());
        mole_fraction_mf[lev]->FillBoundary(geom[lev].periodicity());
    }
}

void
LowMach::UpdateInterfacialChemicalPotential(int lev)
{
    BL_PROFILE("Integrator::LowMach::UpdateInterfacialChemicalPotential");
    if (!interfacial_forces_enabled) return;

    const int nphase = interfacial_model.NumberOfPhases();
    const Set::Scalar epsilon =
        interfacial_model.InterfaceThickness();
    const Set::Scalar* DX = geom[lev].CellSize();
    const amrex::Box domain = geom[lev].Domain();
    amrex::GpuArray<int,AMREX_SPACEDIM> periodic{};
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
        periodic[d] = geom[lev].isPeriodic(d);
    amrex::MultiFab& chemical_potential =
        *interfacial_chemical_potential_mf[lev];
    chemical_potential.setVal(0.0, 0, nphase,
                              chemical_potential.nGrow());

    // Pair kernels avoid placing an N-by-N coefficient matrix in every GPU
    // launch.  Cost scales only with the configured material pairs, while the
    // common two- and three-liquid cases remain small.
    for (int a = 0; a < nphase; ++a)
        for (int b = a + 1; b < nphase; ++b)
        {
            const Set::Scalar pair_regularization =
                interfacial_model.PairRegularization(a,b);
            const Set::Scalar surface_correction =
                interfacial_model.SolidSurfaceCorrection(a,b);
            if (!interfacial_model.HasInterfacialEnergy(a,b)) continue;
            for (amrex::MFIter mfi(chemical_potential, false);
                 mfi.isValid(); ++mfi)
            {
                const amrex::Box bx = mfi.fabbox() & domain;
                Set::Patch<const Set::Scalar> phase =
                    interfacial_volume_fraction_mf.Patch(lev,mfi);
                Set::Patch<Set::Scalar> mu =
                    chemical_potential.array(mfi);
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        const auto stencil = Numeric::GetStencil(
                            i, j, k, domain, periodic);
                        mu(i,j,k,a) +=
                            Model::PhaseField::MultiphaseInterface::
                                PairChemicalPotential(
                                    phase, a, b, pair_regularization, epsilon,
                                    i, j, k, DX, stencil);
                        mu(i,j,k,b) +=
                            Model::PhaseField::MultiphaseInterface::
                                PairChemicalPotential(
                                    phase, b, a, pair_regularization, epsilon,
                                    i, j, k, DX, stencil);
                        if (surface_correction != 0.0)
                        {
                            mu(i,j,k,a) +=
                                Model::PhaseField::MultiphaseInterface::
                                    LiquidSurfaceCorrectionChemicalPotential(
                                        phase, a, b, surface_correction,
                                        epsilon,
                                        i, j, k, DX, stencil);
                            mu(i,j,k,b) +=
                                Model::PhaseField::MultiphaseInterface::
                                    SolidSurfaceCorrectionChemicalPotential(
                                        phase, a, b, surface_correction,
                                        epsilon,
                                        i, j, k, DX, stencil);
                        }
                    });
            }
        }
    chemical_potential.FillBoundary(geom[lev].periodicity());
}

void
LowMach::ApplyInterfacialTransport(Set::Scalar time, Set::Scalar dt)
{
    BL_PROFILE("Integrator::LowMach::ApplyInterfacialTransport");
    if (interfacial_transport == InterfacialTransport::Advected ||
        liquid_species.empty() || !(dt > 0.0)) return;

    const int nlev = finest_level + 1;
    const int nliquid = static_cast<int>(liquid_species.size());
    Set::Scalar fastest_rate = 0.0;
    for (int lev = 0; lev < nlev; ++lev)
    {
        const auto dx = geom[lev].CellSizeArray();
        Set::Scalar dxmin = dx[0];
        for (int d = 1; d < AMREX_SPACEDIM; ++d)
            dxmin = Util::Min(dxmin, dx[d]);
        fastest_rate = Util::Max(fastest_rate,
            conservative_allen_cahn_diffusivity / (dxmin * dxmin));
    }
    int substeps = static_cast<int>(
        std::ceil(dt * fastest_rate / 0.02));
    if (substeps < 1) substeps = 1;
    const Set::Scalar subdt = dt / substeps;

    for (int step = 0; step < substeps; ++step)
    {
        for (int lev = 0; lev < nlev; ++lev)
            UpdateComponentState(lev, *component_density_mf[lev]);
        for (int lev = 0; lev < nlev; ++lev)
        {
            const Set::Scalar* DX = geom[lev].CellSize();
            const amrex::Box domain = geom[lev].Domain();
            amrex::GpuArray<int,AMREX_SPACEDIM> periodic{};
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
                periodic[d] = geom[lev].isPeriodic(d);
            for (int liquid = 0; liquid < nliquid; ++liquid)
            {
                amrex::MultiFab compression_flux(
                    component_density_mf[lev]->boxArray(),
                    component_density_mf[lev]->DistributionMap(),
                    AMREX_SPACEDIM, 1);
                for (amrex::MFIter mfi(compression_flux,
                        amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const amrex::Box bx = mfi.growntilebox(1) & domain;
                    Set::Patch<const Set::Scalar> volume_fraction =
                        interfacial_volume_fraction_mf.Patch(lev,mfi);
                    Set::Patch<Set::Scalar> flux =
                        compression_flux.array(mfi);
                    amrex::ParallelFor(
                        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                        {
                            const auto stencil = Numeric::GetStencil(
                                i, j, k, domain, periodic);
                            const Set::Scalar eta = Util::Clamp(
                                volume_fraction(i,j,k,liquid), 0.0, 1.0);
                            const Set::Vector gradient = Numeric::Gradient(
                                volume_fraction, i, j, k, liquid,
                                DX, stencil);
                            const Set::Scalar target =
                                4.0 * eta * (1.0 - eta) /
                                interfacial_thickness;
                            const Set::Scalar norm = gradient.norm();
                            for (int d = 0; d < AMREX_SPACEDIM; ++d)
                                flux(i,j,k,d) = norm > 1.0e-14 /
                                    interfacial_thickness ?
                                    target * gradient(d) / norm : 0.0;
                        });
                }
                compression_flux.FillBoundary(geom[lev].periodicity());
                const int species = liquid_species[liquid];
                const Set::Scalar reference_density =
                    this->reference_density[species];
                for (amrex::MFIter mfi(*component_density_mf[lev],
                        amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const amrex::Box& bx = mfi.tilebox();
                    Set::Patch<Set::Scalar> component_density =
                        component_density_mf.Patch(lev,mfi);
                    Set::Patch<const Set::Scalar> volume_fraction =
                        interfacial_volume_fraction_mf.Patch(lev,mfi);
                    Set::Patch<const Set::Scalar> flux =
                        compression_flux.const_array(mfi);
                    amrex::ParallelFor(
                        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                        {
                            const auto stencil = Numeric::GetStencil(
                                i, j, k, domain, periodic);
                            Set::Scalar divergence = 0.0;
                            for (int d = 0; d < AMREX_SPACEDIM; ++d)
                                divergence += Numeric::Gradient(
                                    flux, i, j, k, d, DX, stencil)(d);
                            const Set::Scalar laplacian =
                                Numeric::Laplacian(volume_fraction,
                                    i, j, k, liquid, DX, stencil);
                            component_density(i,j,k,species) +=
                                subdt * reference_density *
                                conservative_allen_cahn_diffusivity *
                                (laplacian - divergence);
                        });
                }
            }
        }

        for (int lev = nlev - 2; lev >= 0; --lev)
            amrex::average_down(*component_density_mf[lev + 1],
                *component_density_mf[lev], geom[lev + 1], geom[lev],
                0, nspecies, refRatio(lev));
        for (int lev = 0; lev < nlev; ++lev)
        {
            component_density_bc->FillBoundary(
                *component_density_mf[lev], 0, nspecies,
                time + (step + 1) * subdt, 0);
            component_density_mf[lev]->FillBoundary(
                geom[lev].periodicity());
        }
    }
    for (int lev = 0; lev < nlev; ++lev)
        UpdateComponentState(lev, *component_density_mf[lev]);
}

//
// Calculate derived quantities (only if diagnostics_extended_fields is enabled)
// - energy
// - vorticity
// - momentum
// - mixture transport properties
// - instantaneous chemistry source terms
//
void
LowMach::UpdateDerivedDiagnostics(int lev, const amrex::MultiFab& u_mf, const amrex::MultiFab& T_mf)
{
    if (!diagnostics_extended_fields) return;

    const auto dx = geom[lev].CellSizeArray();
    amrex::Box domain = geom[lev].Domain();
    amrex::GpuArray<int,AMREX_SPACEDIM> periodic{};
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
        periodic[d] = geom[lev].isPeriodic(d);
    const int ngas = ngas_species;
    const int number_of_species = nspecies;
    const Set::Scalar rho_floor = density_floor;
    const Set::Scalar p_reference = pressure_reference;
    const Set::Scalar* inverse_reference_density =
        amrex::get<8>(thermal_data);
    const auto gas_data = gas_device_data;
    const auto thermal = thermal_data;
    const auto chemistry_data = chemistry_device_data;
    momentum_mf[lev]->setVal(0.0);
    energy_mf[lev]->setVal(0.0);
    vorticity_mf[lev]->setVal(0.0);
    viscosity_mf[lev]->setVal(0.0);
    thermal_conductivity_coeff_mf[lev]->setVal(0.0);
    diffusion_coeff_mf[lev]->setVal(0.0);
    wdot_mf[lev]->setVal(0.0);
    qdot_mf[lev]->setVal(0.0);

    for (amrex::MFIter mfi(u_mf, true); mfi.isValid(); ++mfi)
    {
        amrex::Box bx = mfi.growntilebox(1);
        bx &= domain;
        Set::Patch<const Set::Scalar> u = u_mf.array(mfi);
        Set::Patch<const Set::Scalar> T = T_mf.array(mfi);
        Set::Patch<const Set::Scalar> component_density = component_density_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> rho = density_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> M = momentum_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> E = energy_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> omega = vorticity_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> mu = viscosity_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> kappa = thermal_conductivity_coeff_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> diffusion = diffusion_coeff_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> wdot = wdot_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> qdot = qdot_mf.Patch(lev,mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            const Set::Scalar density = rho(i,j,k);
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
                M(i,j,k,d) = density * u(i,j,k,d);
            E(i,j,k) = Model::Gas::Gas::ComputeEnergy(
                gas_data, density, M(i,j,k,0), M(i,j,k,1),
                T(i,j,k), component_density, i, j, k);

            auto sten = Numeric::GetStencil(i, j, k, domain, periodic);
            Set::Matrix grad_u =
                Numeric::Gradient(u, i, j, k, dx.data(), sten);
            omega(i,j,k) = grad_u(1,0) - grad_u(0,1);

            mu(i,j,k) = ComputeViscosity(
                component_density, T(i,j,k), i, j, k, thermal);
            auto [thermal_gas_volume_fraction, gas_heat_capacity,
                heat_capacity, conductivity, cp] = ComputeThermalState(
                    component_density, T(i,j,k), i, j, k, thermal);
            (void)thermal_gas_volume_fraction;
            (void)gas_heat_capacity;
            (void)heat_capacity;
            (void)cp;
            kappa(i,j,k) = conductivity;
            for (int n = 0; n < ngas; ++n)
                diffusion(i,j,k,n) =
                    Model::Gas::Gas::DiffusionCoefficient(
                        gas_data, T(i,j,k), p_reference,
                        component_density, i, j, k, n);

            Model::Chemistry::SpeciesArray rhoY{};
            Set::Scalar gas_density = 0.0;
            for (int n = 0; n < ngas; ++n)
            {
                rhoY[n] = Util::Max(component_density(i,j,k,n), 0.0);
                gas_density += rhoY[n];
            }
            Set::Scalar raw_gas_volume_fraction = 0.0;
            if (gas_density > rho_floor && T(i,j,k) > 0.0 &&
                p_reference > 0.0)
                raw_gas_volume_fraction = Util::Max(
                    gas_density * Model::Gas::Gas::GasConstant(
                        gas_data, component_density, i, j, k) *
                        T(i,j,k) / p_reference, 0.0);
            Set::Scalar condensed_volume_fraction = 0.0;
            for (int n = ngas; n < number_of_species; ++n)
                condensed_volume_fraction += Util::Max(
                    component_density(i,j,k,n), 0.0) *
                    inverse_reference_density[n];
            const Set::Scalar gas_volume_fraction = Util::Min(
                raw_gas_volume_fraction,
                1.0 - Model::PhaseField::H(condensed_volume_fraction));
            if (gas_volume_fraction > 0.0 &&
                raw_gas_volume_fraction > 0.0)
            {
                for (int n = 0; n < ngas; ++n)
                    rhoY[n] /= raw_gas_volume_fraction;
                Model::Chemistry::Source reaction{};
#ifdef ALAMO_GPU
                const auto& [gross_model, model, solver, host_chemistry,
                             host_gas] = chemistry_data;
                if (gross_model)
                    reaction = model.ComputeChemistrySources(
                        p_reference, T(i,j,k), rhoY, 0.0, nullptr);
#else
                const auto& [gross_model, model, solver, host_chemistry,
                             host_gas] = chemistry_data;
                reaction = host_chemistry->ComputeChemistrySources(
                    p_reference, T(i,j,k), rhoY, 0.0, host_gas);
#endif
                for (int n = 0; n < ngas; ++n)
                    wdot(i,j,k,n) =
                        gas_volume_fraction * reaction.first[n];
                qdot(i,j,k) = gas_volume_fraction * reaction.second;
            }
        });
    }

    momentum_mf[lev]->FillBoundary(geom[lev].periodicity());
    energy_mf[lev]->FillBoundary(geom[lev].periodicity());
    vorticity_mf[lev]->FillBoundary(geom[lev].periodicity());
    viscosity_mf[lev]->FillBoundary(geom[lev].periodicity());
    thermal_conductivity_coeff_mf[lev]->FillBoundary(geom[lev].periodicity());
    diffusion_coeff_mf[lev]->FillBoundary(geom[lev].periodicity());
    wdot_mf[lev]->FillBoundary(geom[lev].periodicity());
    qdot_mf[lev]->FillBoundary(geom[lev].periodicity());
}

void
LowMach::AdvanceChemistry(int lev, amrex::MultiFab& T_mf,
                            amrex::MultiFab& component_density_mf, Set::Scalar dt)
{
    if (!chemistry.Split() || !(dt > 0.0)) return;
    const int ngas = ngas_species;
    const int number_of_species = nspecies;
    const Set::Scalar rho_floor = density_floor;
    const Set::Scalar p_reference = pressure_reference;
    const Set::Scalar* inverse_reference_density =
        amrex::get<8>(thermal_data);
    const auto gas_data = gas_device_data;
    const auto thermal = thermal_data;
    const auto chemistry_data = chemistry_device_data;

    for (amrex::MFIter mfi(T_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        Set::Patch<Set::Scalar> T = T_mf.array(mfi);
        Set::Patch<Set::Scalar> component_density = component_density_mf.array(mfi);
        Set::Patch<Set::Scalar> integrated_dilatation = chemistry_dilatation_mf.Patch(lev,mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Model::Chemistry::SpeciesArray rhoY{};
            Set::Scalar gas_density = 0.0;
            for (int n = 0; n < ngas; ++n)
            {
                rhoY[n] = Util::Max(component_density(i,j,k,n), 0.0);
                gas_density += rhoY[n];
            }
            if (!(gas_density > rho_floor) || !(T(i,j,k) > 0.0)) return;

            // Keep the EOS volume unbounded here.  If the transported state is
            // temporarily overfilled, clipping it before time weighting would
            // multiply the chemistry heat release by the overfill ratio.
            const Set::Scalar raw_gas_volume_fraction = Util::Max(
                gas_density * Model::Gas::Gas::GasConstant(
                    gas_data, component_density, i, j, k) *
                    T(i,j,k) / p_reference, 0.0);
            Set::Scalar condensed_volume_fraction = 0.0;
            for (int n = ngas; n < number_of_species; ++n)
                condensed_volume_fraction += Util::Max(
                    component_density(i,j,k,n), 0.0) *
                    inverse_reference_density[n];
            const Set::Scalar gas_accessibility =
                1.0 - Model::PhaseField::H(condensed_volume_fraction);
            const Set::Scalar reacting_volume_fraction = Util::Min(
                raw_gas_volume_fraction, gas_accessibility);
            const Set::Scalar chemistry_weight =
                raw_gas_volume_fraction > 0.0 ?
                    reacting_volume_fraction /
                        raw_gas_volume_fraction : 0.0;
            if (!(chemistry_weight > 0.0)) return;

            Set::Scalar temperature = T(i,j,k);
            auto [thermal_gas_volume_fraction, gas_heat_capacity, heat_capacity,
                conductivity, cp] = ComputeThermalState(
                    component_density, temperature, i, j, k, thermal);
            (void)thermal_gas_volume_fraction;
            (void)gas_heat_capacity;
            (void)conductivity;
            const Set::Scalar mixture_density = cp > 0.0 ?
                heat_capacity / cp : gas_density;

            const auto& [gross_model, model, solver, host_chemistry, host_gas] =
                chemistry_data;
#ifdef ALAMO_GPU
            auto result = model.Advance(
                dt * chemistry_weight, p_reference, mixture_density,
                rhoY, temperature, gas_data, solver);
#else
            auto result = host_chemistry->Advance(
                dt * chemistry_weight, p_reference, mixture_density,
                ngas, rhoY, temperature, host_gas);
#endif
            if (!result.converged)
                Util::Abort(INFO, "Local chemistry integration failed at ",
                    i, ",", j, ",", k, " dt=", result.failed_dt,
                    " residual=", result.residual_norm);

            for (int n = 0; n < ngas; ++n)
                component_density(i,j,k,n) = rhoY[n];
            T(i,j,k) = temperature;
            integrated_dilatation(i,j,k) += result.volume_fraction_change;
        });
    }
}

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
amrex::GpuTuple<Set::Scalar, Set::Scalar, Set::Scalar, Set::Scalar, Set::Scalar>
LowMach::ComputeThermalState(
    Set::Patch<const Set::Scalar> component_density,
    Set::Scalar temperature, int i, int j, int k,
    const ThermalData& data)
{
    const auto& [gas_data, nspecies, ngas_species, density_floor,
                 pressure_reference, condensed_thermal_transport,
                 condensed_specific_heat, condensed_thermal_conductivity,
                 condensed_inverse_reference_density,
                 liquid_inverse_reference_density,
                 condensed_dynamic_viscosity] = data;
    (void)liquid_inverse_reference_density;
    (void)condensed_dynamic_viscosity;
    Set::Scalar density = 0.0;
    Set::Scalar gas_density = 0.0;
    for (int n = 0; n < nspecies; ++n)
        density += Util::Max(component_density(i,j,k,n), 0.0);
    for (int n = 0; n < ngas_species; ++n)
        gas_density += Util::Max(component_density(i,j,k,n), 0.0);

    const Set::Scalar cp = Model::Gas::Gas::CpMass(
        gas_data, temperature, component_density, i, j, k);
    const Set::Scalar gas_conductivity =
        Model::Gas::Gas::ThermalConductivity(
            gas_data, temperature, component_density, i, j, k);
    Set::Scalar gas_volume_fraction = 0.0;
    if (gas_density > density_floor && temperature > 0.0 &&
        pressure_reference > 0.0)
        gas_volume_fraction = Util::Clamp(
            gas_density * Model::Gas::Gas::GasConstant(
                gas_data, component_density, i, j, k) * temperature /
                pressure_reference, 0.0, 1.0);

    if (!condensed_thermal_transport)
        return {gas_volume_fraction, density * cp, density * cp,
                gas_conductivity, cp};

    Set::Scalar condensed_volume_fraction = 0.0;
    for (int n = ngas_species; n < nspecies; ++n)
        condensed_volume_fraction += Util::Max(
            component_density(i,j,k,n), 0.0) *
            condensed_inverse_reference_density[n];
    const Set::Scalar solid_fraction =
        Model::PhaseField::H(condensed_volume_fraction);
    gas_volume_fraction = 1.0 - solid_fraction;

    Set::Scalar intrinsic_gas_density = 0.0;
    const Set::Scalar gas_constant = Model::Gas::Gas::GasConstant(
        gas_data, component_density, i, j, k);
    if (gas_density > density_floor && temperature > 0.0 &&
        pressure_reference > 0.0 && gas_constant > 0.0)
        intrinsic_gas_density =
            pressure_reference / (gas_constant * temperature);
    const Set::Scalar gas_heat_capacity =
        gas_volume_fraction * intrinsic_gas_density * cp;
    Set::Scalar heat_capacity = gas_heat_capacity;
    Set::Scalar conductivity = gas_volume_fraction * gas_conductivity;
    const Set::Scalar solid_scale = condensed_volume_fraction > 0.0 ?
        solid_fraction / condensed_volume_fraction : 0.0;
    for (int n = ngas_species; n < nspecies; ++n)
    {
        const Set::Scalar partial_density =
            Util::Max(component_density(i,j,k,n), 0.0);
        heat_capacity += solid_scale * partial_density *
            condensed_specific_heat[n];
        conductivity += solid_scale * partial_density *
            condensed_inverse_reference_density[n] *
            condensed_thermal_conductivity[n];
    }
    return {gas_volume_fraction, gas_heat_capacity, heat_capacity,
            conductivity, cp};
}

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
Set::Scalar
LowMach::ComputeViscosity(
    Set::Patch<const Set::Scalar> component_density,
    Set::Scalar temperature, int i, int j, int k,
    const ThermalData& data)
{
    const auto& [gas_data, nspecies, ngas_species, density_floor,
                 pressure_reference, condensed_thermal_transport,
                 condensed_specific_heat, condensed_thermal_conductivity,
                 condensed_inverse_reference_density,
                 liquid_inverse_reference_density,
                 condensed_dynamic_viscosity] = data;
    (void)density_floor;
    (void)pressure_reference;
    (void)condensed_thermal_transport;
    (void)condensed_specific_heat;
    (void)condensed_thermal_conductivity;
    (void)condensed_inverse_reference_density;
    const Set::Scalar gas_viscosity =
        Model::Gas::Gas::DynamicViscosity(
            gas_data, temperature, component_density, i, j, k);
    Set::Scalar liquid_volume = 0.0;
    for (int n = ngas_species; n < nspecies; ++n)
        liquid_volume += Util::Max(component_density(i,j,k,n), 0.0) *
            liquid_inverse_reference_density[n];
    if (!(liquid_volume > 0.0)) return gas_viscosity;

    const Set::Scalar liquid_fraction = Model::PhaseField::H(liquid_volume);
    const Set::Scalar scale = liquid_fraction / liquid_volume;
    Set::Scalar viscosity = (1.0 - liquid_fraction) * gas_viscosity;
    for (int n = ngas_species; n < nspecies; ++n)
        if (liquid_inverse_reference_density[n] > 0.0)
            viscosity += scale *
                Util::Max(component_density(i,j,k,n), 0.0) *
                liquid_inverse_reference_density[n] *
                condensed_dynamic_viscosity[n];
    return viscosity;
}

void
LowMach::ApplyImplicitPhaseChange(Set::Scalar time, Set::Scalar dt)
{
    BL_PROFILE("Integrator::LowMach::ApplyImplicitPhaseChange");
    if (!has_split_phase_change || !(dt > 0.0))
        return;

    const int nlev = finest_level + 1;
    for (int lev = 0; lev < nlev; ++lev)
    {
        phase_change_dilatation_mf[lev]->setVal(0.0);
        phase_change_heat_mf[lev]->setVal(0.0);
    }

    if (has_equilibrium_phase_change)
        ApplyEquilibriumPhaseChange(time, dt);
    if (has_kinetic_phase_change)
        ApplyKineticPhaseChange(time, dt);
}

void
LowMach::ApplyEquilibriumPhaseChange(Set::Scalar time, Set::Scalar dt)
{
    BL_PROFILE("Integrator::LowMach::ApplyEquilibriumPhaseChange");
    const int nlev = finest_level + 1;
    const auto thermal = thermal_data;
    const Set::Scalar p_reference = pressure_reference;

    // Mechanisms are applied in their configured order.  This is important
    // for solid -> liquid -> gas chains: newly melted material is immediately
    // visible to the following equilibrium transition without introducing a
    // kinetic interface-speed parameter.
    for (const auto& configured_mechanism : mechanisms)
    {
        const auto mechanism = configured_mechanism;
        if (!mechanism.Equilibrium()) continue;
        const Set::Scalar transition_temperature =
            mechanism.EquilibriumTemperature();
        const Set::Scalar latent_heat = mechanism.LatentHeat();

        for (int lev = 0; lev < nlev; ++lev)
            for (amrex::MFIter mfi(*component_density_mf[lev],
                                    amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                Set::Patch<Set::Scalar> component_density =
                    component_density_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> component_density_state =
                    component_density_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> rigid_eta;
                Set::Patch<const Set::Scalar> rigid_species_eta;
                if (!rigid_solid_species.empty())
                {
                    rigid_eta = rigid_eta_mf.Patch(lev,mfi);
                    rigid_species_eta = rigid_species_eta_mf.Patch(lev,mfi);
                }
                Set::Patch<Set::Scalar> temperature =
                    temperature_mf.Patch(lev,mfi);
                Set::Patch<Set::Scalar> integrated_dilatation =
                    phase_change_dilatation_mf.Patch(lev,mfi);

                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        auto [gas_volume_fraction, gas_heat_capacity,
                            heat_capacity, conductivity, cp] =
                            ComputeThermalState(component_density_state,
                                temperature(i,j,k), i, j, k, thermal);
                        (void)gas_heat_capacity;
                        (void)conductivity;
                        (void)cp;
                        if (!(heat_capacity > 0.0)) return;

                        const Set::Scalar old_temperature =
                            temperature(i,j,k);
                        const Set::Scalar sensible_energy = heat_capacity *
                            (old_temperature - transition_temperature);
                        const Set::Scalar requested_mass_change =
                            -sensible_energy / latent_heat;
                        if (requested_mass_change == 0.0) return;

                        const Model::Mechanism::State state = {
                            component_density_state, rigid_eta,
                            rigid_species_eta, old_temperature,
                            p_reference, dt};
                        auto [mass_change, volume_change, coupled_heat] =
                            mechanism.ApplyEquilibriumChange(
                                component_density, state,
                                requested_mass_change, heat_capacity,
                                gas_volume_fraction, i, j, k);
                        if (mass_change == 0.0) return;

                        auto [new_gas_volume_fraction,
                            new_gas_heat_capacity, new_heat_capacity,
                            new_conductivity, new_cp] =
                            ComputeThermalState(component_density_state,
                                old_temperature, i, j, k, thermal);
                        (void)new_gas_volume_fraction;
                        (void)new_gas_heat_capacity;
                        (void)new_conductivity;
                        (void)new_cp;
                        if (new_heat_capacity > 0.0)
                        {
                            const Set::Scalar remaining_energy =
                                sensible_energy +
                                latent_heat * mass_change + coupled_heat;
                            temperature(i,j,k) = Util::Max(
                                transition_temperature +
                                    remaining_energy / new_heat_capacity,
                                1.0e-12);
                        }
                        integrated_dilatation(i,j,k) += volume_change;
                    });
            }

        // The finest representation owns covered regions.  Average it down
        // before the next transition so a chained mechanism sees one
        // consistent composite state on every level.
        for (int lev = nlev - 2; lev >= 0; --lev)
        {
            amrex::average_down(*component_density_mf[lev + 1],
                *component_density_mf[lev], geom[lev + 1], geom[lev],
                0, nspecies, refRatio(lev));
            amrex::average_down(*temperature_mf[lev + 1],
                *temperature_mf[lev], geom[lev + 1], geom[lev],
                0, 1, refRatio(lev));
            amrex::average_down(*phase_change_dilatation_mf[lev + 1],
                *phase_change_dilatation_mf[lev], geom[lev + 1], geom[lev],
                0, 1, refRatio(lev));
        }
        for (int lev = 0; lev < nlev; ++lev)
        {
            component_density_bc->FillBoundary(
                *component_density_mf[lev], 0, nspecies, time, 0);
            component_density_mf[lev]->FillBoundary(
                geom[lev].periodicity());
            temperature_bc->FillBoundary(
                *temperature_mf[lev], 0, 1, time, 0);
            temperature_mf[lev]->FillBoundary(geom[lev].periodicity());
            UpdateComponentState(lev, *component_density_mf[lev]);
        }
    }
}

void
LowMach::ApplyKineticPhaseChange(Set::Scalar time, Set::Scalar dt)
{
    BL_PROFILE("Integrator::LowMach::ApplyKineticPhaseChange");
    const int nlev = finest_level + 1;
    for (int lev = 0; lev < nlev; ++lev)
    {
        component_density_bc->FillBoundary(
            *component_density_mf[lev], 0, nspecies, time, 0);
        component_density_mf[lev]->FillBoundary(geom[lev].periodicity());
        temperature_bc->FillBoundary(
            *temperature_mf[lev], 0, 1, time, 0);
        temperature_mf[lev]->FillBoundary(geom[lev].periodicity());
        UpdateComponentState(lev, *component_density_mf[lev]);
    }

    const auto thermal = thermal_data;
    const Set::Scalar p_reference = pressure_reference;
    for (const auto& configured_mechanism : mechanisms)
    {
        const auto mechanism = configured_mechanism;
        if (!mechanism.Kinetic()) continue;

        // Form the local diffuse surface moments on every level before any
        // level is updated.  This keeps the AMR hierarchy synchronized while
        // FillPatch supplies coarse/fine and periodic ghost values.
        for (int lev = 0; lev < nlev; ++lev)
        {
            const Set::Scalar* DX = geom[lev].CellSize();
            kinetic_diffuse_moment_mf[lev]->setVal(0.0);

            for (amrex::MFIter mfi(*component_density_mf[lev]);
                 mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                Set::Patch<const Set::Scalar> component_density_state =
                    component_density_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> temperature_state =
                    temperature_mf.Patch(lev,mfi);
                Set::Patch<Set::Scalar> moment =
                    kinetic_diffuse_moment_mf.Patch(lev,mfi);
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        auto [surface_measure, surface_normal] =
                            mechanism.KineticDiffuseSurface(
                                component_density_state,
                                i, j, k, DX);
                        if (!(surface_measure > 0.0) ||
                            !(temperature_state(i,j,k) > 0.0))
                            return;

                        auto [gas_volume_fraction, gas_heat_capacity,
                            heat_capacity, conductivity, cp] =
                            ComputeThermalState(component_density_state,
                                temperature_state(i,j,k), i, j, k, thermal);
                        (void)gas_volume_fraction;
                        (void)gas_heat_capacity;
                        (void)conductivity;
                        (void)cp;
                        if (!(heat_capacity > 0.0)) return;
                        moment(i,j,k,0) = surface_measure;
                        moment(i,j,k,1) = surface_measure *
                            temperature_state(i,j,k);
                        moment(i,j,k,2) = surface_measure *
                            surface_measure / heat_capacity;
                        moment(i,j,k,3) = surface_measure;
                        for (int d = 0; d < AMREX_SPACEDIM; ++d)
                            moment(i,j,k,4+d) = surface_normal(d);
                    });
            }
        }

        if (nlev == 1)
            kinetic_diffuse_moment_mf[0]->FillBoundary(
                geom[0].periodicity());
        else
            for (int lev = 0; lev < nlev; ++lev)
                FillPatch(lev, time, kinetic_diffuse_moment_mf,
                    *kinetic_diffuse_moment_mf[lev], bc_nothing, 0);

        for (int lev = 0; lev < nlev; ++lev)
        {
            const Set::Scalar* DX = geom[lev].CellSize();
            constexpr Set::Scalar support_per_thickness =
                0.54930614433405484570;
            const Set::Scalar support = support_per_thickness *
                interfacial_thickness;

            for (amrex::MFIter mfi(*component_density_mf[lev],
                                    amrex::TilingIfNotGPU());
                 mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                Set::Patch<Set::Scalar> component_density =
                    component_density_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> component_density_state =
                    component_density_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> rigid_eta;
                Set::Patch<const Set::Scalar> rigid_species_eta;
                if (!rigid_solid_species.empty())
                {
                    rigid_eta = rigid_eta_mf.Patch(lev,mfi);
                    rigid_species_eta = rigid_species_eta_mf.Patch(lev,mfi);
                }
                Set::Patch<Set::Scalar> temperature =
                    temperature_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> moment =
                    kinetic_diffuse_moment_mf.Patch(lev,mfi);
                Set::Patch<Set::Scalar> integrated_dilatation =
                    phase_change_dilatation_mf.Patch(lev,mfi);
                Set::Patch<Set::Scalar> integrated_heat =
                    phase_change_heat_mf.Patch(lev,mfi);
                amrex::ParallelFor(
                    bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        const Set::Scalar surface_measure = moment(i,j,k,3);
                        if (!(surface_measure > 0.0)) return;

                        Set::Scalar diffuse_moment[3] = {0.0, 0.0, 0.0};
                        // Five-point Gauss--Legendre quadrature on the
                        // physical normal interval [-support,support].  The
                        // quadrature order is fixed, but its locations scale
                        // with interface.thickness and are unrelated to a
                        // cell stencil width.
                        for (int sample = 0; sample < 5; ++sample)
                        {
                            Set::Scalar point = 0.0;
                            Set::Scalar quadrature_weight =
                                0.56888888888888888889;
                            if (sample == 1 || sample == 2)
                            {
                                point = sample == 1 ?
                                    -0.53846931010568309104 :
                                     0.53846931010568309104;
                                quadrature_weight =
                                    0.47862867049936646804;
                            }
                            else if (sample == 3 || sample == 4)
                            {
                                point = sample == 3 ?
                                    -0.90617984593866399280 :
                                     0.90617984593866399280;
                                quadrature_weight =
                                    0.23692688505618908751;
                            }
                            const Set::Scalar distance = point * support;
                            const Set::Scalar interval_weight =
                                quadrature_weight * support;

                            const Set::Scalar shift_x = distance *
                                moment(i,j,k,4) / DX[0];
                            const int i0 = i + static_cast<int>(
                                std::floor(shift_x));
                            const Set::Scalar fraction_x =
                                shift_x - std::floor(shift_x);
#if AMREX_SPACEDIM >= 2
                            const Set::Scalar shift_y = distance *
                                moment(i,j,k,5) / DX[1];
                            const int j0 = j + static_cast<int>(
                                std::floor(shift_y));
                            const Set::Scalar fraction_y =
                                shift_y - std::floor(shift_y);
#endif
#if AMREX_SPACEDIM == 3
                            const Set::Scalar shift_z = distance *
                                moment(i,j,k,6) / DX[2];
                            const int k0 = k + static_cast<int>(
                                std::floor(shift_z));
                            const Set::Scalar fraction_z =
                                shift_z - std::floor(shift_z);
#endif

#if AMREX_SPACEDIM == 1
                            for (int ii = 0; ii <= 1; ++ii)
                            {
                                const Set::Scalar weight = interval_weight *
                                    (ii ? fraction_x : 1.0 - fraction_x);
                                for (int n = 0; n < 3; ++n)
                                    diffuse_moment[n] += weight *
                                        moment(i0+ii,j,k,n);
                            }
#elif AMREX_SPACEDIM == 2
                            for (int jj = 0; jj <= 1; ++jj)
                                for (int ii = 0; ii <= 1; ++ii)
                                {
                                    const Set::Scalar weight =
                                        interval_weight *
                                        (ii ? fraction_x :
                                              1.0 - fraction_x) *
                                        (jj ? fraction_y :
                                              1.0 - fraction_y);
                                    for (int n = 0; n < 3; ++n)
                                        diffuse_moment[n] += weight *
                                            moment(i0+ii,j0+jj,k,n);
                                }
#else
                            for (int kk = 0; kk <= 1; ++kk)
                                for (int jj = 0; jj <= 1; ++jj)
                                    for (int ii = 0; ii <= 1; ++ii)
                                    {
                                        const Set::Scalar weight =
                                            interval_weight *
                                            (ii ? fraction_x :
                                                  1.0 - fraction_x) *
                                            (jj ? fraction_y :
                                                  1.0 - fraction_y) *
                                            (kk ? fraction_z :
                                                  1.0 - fraction_z);
                                        for (int n = 0; n < 3; ++n)
                                            diffuse_moment[n] += weight *
                                                moment(i0+ii,j0+jj,k0+kk,n);
                                    }
#endif
                        }

                        if (!(diffuse_moment[0] > 0.0)) return;
                        const Set::Scalar diffuse_temperature =
                            diffuse_moment[1] / diffuse_moment[0];
                        const Set::Scalar diffuse_thermal_response =
                            diffuse_moment[2] / diffuse_moment[0];

                        auto [gas_volume_fraction, gas_heat_capacity,
                            heat_capacity, conductivity, cp] =
                            ComputeThermalState(component_density_state,
                                temperature(i,j,k), i, j, k, thermal);
                        (void)gas_volume_fraction;
                        (void)gas_heat_capacity;
                        (void)conductivity;
                        (void)cp;
                        Set::Scalar local_heat = 0.0;
                        const Model::Mechanism::State state = {
                            component_density_state, rigid_eta,
                            rigid_species_eta, temperature(i,j,k),
                            p_reference, dt};
                        integrated_dilatation(i,j,k) +=
                            mechanism.ApplyKineticChange(
                                component_density, state, heat_capacity,
                                surface_measure, diffuse_temperature,
                                diffuse_thermal_response,
                                temperature(i,j,k), local_heat,
                                i, j, k);
                        integrated_heat(i,j,k) += local_heat;
                    });
            }
        }

        for (int lev = nlev - 2; lev >= 0; --lev)
        {
            amrex::average_down(*component_density_mf[lev + 1],
                *component_density_mf[lev], geom[lev + 1], geom[lev],
                0, nspecies, refRatio(lev));
            amrex::average_down(*temperature_mf[lev + 1],
                *temperature_mf[lev], geom[lev + 1], geom[lev],
                0, 1, refRatio(lev));
            amrex::average_down(*phase_change_dilatation_mf[lev + 1],
                *phase_change_dilatation_mf[lev], geom[lev + 1], geom[lev],
                0, 1, refRatio(lev));
            amrex::average_down(*phase_change_heat_mf[lev + 1],
                *phase_change_heat_mf[lev], geom[lev + 1], geom[lev],
                0, 1, refRatio(lev));
        }
        for (int lev = 0; lev < nlev; ++lev)
        {
            component_density_bc->FillBoundary(
                *component_density_mf[lev], 0, nspecies, time, 0);
            component_density_mf[lev]->FillBoundary(
                geom[lev].periodicity());
            temperature_bc->FillBoundary(
                *temperature_mf[lev], 0, 1, time, 0);
            temperature_mf[lev]->FillBoundary(geom[lev].periodicity());
            UpdateComponentState(lev, *component_density_mf[lev]);
        }
    }

}

void
LowMach::ApplyImplicitDiffusion(Set::Scalar time, Set::Scalar dt)
{
    BL_PROFILE("Integrator::LowMach::ApplyImplicitDiffusion");
    if ((!implicit_momentum_diffusion && !implicit_thermal_diffusion &&
        !implicit_species_diffusion) || !(dt > 0.0)) return;

    const int nlev = finest_level + 1;
    const bool thermochemical_diffusion =
        implicit_thermal_diffusion || implicit_species_diffusion;
    const int ngas = ngas_species;
    const int number_of_species = nspecies;
    const Set::Scalar rho_floor = density_floor;
    const Set::Scalar p_reference = pressure_reference;
    const Set::Scalar* inverse_reference_density =
        amrex::get<8>(thermal_data);
    const auto gas_data = gas_device_data;
    const auto thermal = thermal_data;
    diffusion.SetLayout(geom, refRatio(), temperature_mf, nlev, 1);
    diffusion.FillBoundary(temperature_mf, *temperature_bc, time, 1);
    diffusion.FillBoundary(
        component_density_mf, *component_density_bc, time, nspecies);

    if (implicit_momentum_diffusion)
        diffusion.FillBoundary(
            velocity_mf, *velocity_bc, time, AMREX_SPACEDIM);

    //
    // This updates the dilatation field due to diffusion (no integration)
    //
    if (thermochemical_diffusion)
    {
        for (int lev = 0; lev < nlev; ++lev)
        {
            diffusion_dilatation_mf[lev]->setVal(0.0);

            for (amrex::MFIter mfi(*diffusion_dilatation_mf[lev],
                                    amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                Set::Patch<const Set::Scalar> component_density = component_density_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> T = temperature_mf.Patch(lev,mfi);
                Set::Patch<Set::Scalar> diffusion_dilatation = diffusion_dilatation_mf.Patch(lev,mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar gas_density = 0.0;
                    for (int n = 0; n < ngas; ++n)
                        gas_density += component_density(i,j,k,n);
                    diffusion_dilatation(i,j,k) = -gas_density *
                        Model::Gas::Gas::GasConstant(
                            gas_data, component_density, i, j, k) *
                        T(i,j,k) / p_reference;
                });
            }
        }
    }

    //
    // This implicitly solves the species diffusion equation over an
    // interval dt
    //
    if (implicit_species_diffusion)
    {
        diffusion.SetLayout(
            geom, refRatio(), component_density_mf, nlev, ngas_species);
        for (int lev = 0; lev < nlev; ++lev)
        {
            diffusion.State(lev, ngas_species).setVal(0.0);
            diffusion.Mass(lev,ngas_species).setVal(0.0);
            diffusion.Mobility(lev, ngas_species).setVal(0.0);

            for (amrex::MFIter mfi(diffusion.State(lev, ngas_species), amrex::TilingIfNotGPU());
                mfi.isValid(); ++mfi)
            {
                const amrex::Box& grown_box = mfi.growntilebox(1);
                const amrex::Box& valid_box = mfi.tilebox();
                Set::Patch<const Set::Scalar> component_density = component_density_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> T = temperature_mf.Patch(lev,mfi);
                Set::Patch<Set::Scalar> Y = diffusion.State(lev, ngas_species).array(mfi);
                Set::Patch<Set::Scalar> a = diffusion.Mass(lev,ngas_species).array(mfi);
                Set::Patch<Set::Scalar> b = diffusion.Mobility(lev, ngas_species).array(mfi);

                amrex::ParallelFor(
                    grown_box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        Set::Scalar gas_density = 0.0;
                        for (int n = 0; n < ngas; ++n)
                            gas_density += component_density(i,j,k,n);
                        for (int n = 0; n < ngas; ++n)
                            Y(i,j,k,n) = gas_density > rho_floor ?
                                component_density(i,j,k,n) / gas_density :
                                (n == 0 ? 1.0 : 0.0);
                    });
                amrex::ParallelFor(
                    valid_box, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        Set::Scalar gas_density = 0.0;
                        for (int n = 0; n < ngas; ++n)
                            gas_density += component_density(i,j,k,n);
                        const Set::Scalar raw_gas_volume_fraction = Util::Max(
                            gas_density * Model::Gas::Gas::GasConstant(
                                gas_data, component_density, i, j, k) *
                                T(i,j,k) / p_reference, 0.0);
                        Set::Scalar condensed_volume_fraction = 0.0;
                        for (int n = ngas; n < number_of_species; ++n)
                            condensed_volume_fraction += Util::Max(
                                component_density(i,j,k,n), 0.0) *
                                inverse_reference_density[n];
                        const Set::Scalar gas_accessibility =
                            1.0 - Model::PhaseField::H(
                                condensed_volume_fraction);
                        const Set::Scalar diffusion_weight =
                            raw_gas_volume_fraction > 0.0 ? Util::Min(
                                1.0, gas_accessibility /
                                    raw_gas_volume_fraction) : 0.0;
                        a(i,j,k) = Util::Max(gas_density, rho_floor);
                        b(i,j,k) = diffusion_weight * gas_density *
                            Model::Gas::Gas::DiffusionCoefficient(
                                gas_data, T(i,j,k), p_reference,
                                component_density, i, j, k, 0);
                    });
            }
        }

        amrex::Vector<amrex::BCRec> species_boundary_conditions(ngas_species);
        for (int n = 0; n < ngas_species; ++n)
            species_boundary_conditions[n] = component_density_bc->GetBCRec(n);
        diffusion.Solve(
            time, dt, species_boundary_conditions, ngas_species);
        for (int lev = 0; lev < nlev; ++lev)
        {
            const amrex::MultiFab& state = diffusion.State(lev, ngas_species);
            for (amrex::MFIter mfi(*component_density_mf[lev],
                                    amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                Set::Patch<Set::Scalar> component_density =
                    component_density_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> Y = state.array(mfi);
                const int ngas = ngas_species;

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar gas_density = 0.0;
                    for (int n = 0; n < ngas; ++n)
                        gas_density += component_density(i,j,k,n);
                    for (int n = 0; n < ngas; ++n)
                        component_density(i,j,k,n) = gas_density * Y(i,j,k,n);
                });
            }
        }
        diffusion.Synchronize(component_density_mf, ngas_species);
        diffusion.FillBoundary(
            component_density_mf, *component_density_bc, time, nspecies);
    }

    //
    // This implicitly solves the thermal transport equation over the
    // interval dt
    //
    if (implicit_thermal_diffusion)
    {
        diffusion.SetLayout(geom, refRatio(), temperature_mf, nlev, 1);
        for (int lev = 0; lev < nlev; ++lev)
        {
            amrex::MultiFab::Copy(
                diffusion.State(lev, 1), *temperature_mf[lev], 0, 0, 1, diffusion.State(lev,1).nGrow());

            for (amrex::MFIter mfi(diffusion.Mass(lev, 1), amrex::TilingIfNotGPU());
                mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                Set::Patch<const Set::Scalar> component_density = component_density_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> T                 = temperature_mf.Patch(lev,mfi);
                Set::Patch<Set::Scalar> a = diffusion.Mass(lev,1).array(mfi);
                Set::Patch<Set::Scalar> b = diffusion.Mobility(lev,1).array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    auto [gas_volume_fraction, gas_heat_capacity,
                        heat_capacity, conductivity, cp] = ComputeThermalState(
                            component_density, T(i,j,k), i, j, k, thermal);
                    (void)gas_volume_fraction;
                    (void)gas_heat_capacity;
                    a(i,j,k) = Util::Max(heat_capacity, rho_floor * cp);
                    b(i,j,k) = conductivity;
                });
            }
        }

        diffusion.Solve(time, dt, temperature_bc->GetBCRec(), 1);
        for (int lev = 0; lev < nlev; ++lev)
            amrex::MultiFab::Copy(*temperature_mf[lev], diffusion.State(lev, 1),
                                    0, 0, 1, 0);
        diffusion.Synchronize(temperature_mf, 1);
        diffusion.FillBoundary(temperature_mf, *temperature_bc, time, 1);
    }

    if (implicit_momentum_diffusion)
    {
        diffusion.SetLayout(
            geom, refRatio(), velocity_mf, nlev, AMREX_SPACEDIM);
        for (int lev = 0; lev < nlev; ++lev)
        {
            amrex::MultiFab& state = diffusion.State(lev, AMREX_SPACEDIM);
            amrex::MultiFab& mass = diffusion.Mass(lev, AMREX_SPACEDIM);
            amrex::MultiFab& mobility =
                diffusion.Mobility(lev, AMREX_SPACEDIM);
            amrex::MultiFab::Copy(
                state, *velocity_mf[lev], 0, 0, AMREX_SPACEDIM, state.nGrow());

            for (amrex::MFIter mfi(mass, amrex::TilingIfNotGPU());
                mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                Set::Patch<const Set::Scalar> component_density =
                    component_density_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> T = temperature_mf.Patch(lev,mfi);
                Set::Patch<Set::Scalar> a = mass.array(mfi);
                Set::Patch<Set::Scalar> b = mobility.array(mfi);
                const Set::Scalar rho_floor = density_floor;
                const int number_of_species = nspecies;

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar density = 0.0;
                    for (int n = 0; n < number_of_species; ++n)
                        density += component_density(i,j,k,n);
                    a(i,j,k) = Util::Max(density, rho_floor);
                    b(i,j,k) = ComputeViscosity(
                        component_density, T(i,j,k), i, j, k, thermal);
                });
            }
        }
        amrex::Vector<amrex::BCRec> velocity_boundary_conditions(AMREX_SPACEDIM);
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
            velocity_boundary_conditions[d] = velocity_bc->GetBCRec(d);
        diffusion.Solve(
            time, dt, velocity_boundary_conditions, AMREX_SPACEDIM);
        for (int lev = 0; lev < nlev; ++lev)
            amrex::MultiFab::Copy(*velocity_mf[lev],
                diffusion.State(lev, AMREX_SPACEDIM), 0, 0,
                AMREX_SPACEDIM, 0);
        diffusion.Synchronize(velocity_mf, AMREX_SPACEDIM);
        diffusion.FillBoundary(
            velocity_mf, *velocity_bc, time, AMREX_SPACEDIM);
    }

    //
    // This updates the diffusion dilatation (no integration)
    //
    if (thermochemical_diffusion)
    {
        for (int lev = 0; lev < nlev; ++lev)
        {
            for (amrex::MFIter mfi(*diffusion_dilatation_mf[lev],
                                    amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                Set::Patch<const Set::Scalar> component_density = component_density_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> T = temperature_mf.Patch(lev,mfi);
                Set::Patch<Set::Scalar> volume_change = diffusion_dilatation_mf.Patch(lev,mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar gas_density = 0.0;
                    for (int n = 0; n < ngas; ++n)
                        gas_density += component_density(i,j,k,n);
                    volume_change(i,j,k) += gas_density *
                        Model::Gas::Gas::GasConstant(
                            gas_data, component_density, i, j, k) *
                        T(i,j,k) / p_reference;
                });
            }
        }
        diffusion.Synchronize(diffusion_dilatation_mf, 1);
    }
}

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
amrex::GpuTuple<Model::Chemistry::SpeciesArray, Set::Scalar /*temperature*/, Set::Scalar /*dilatation*/>
LowMach::ComputeThermochemicalSource(
    Set::Patch<const Set::Scalar> component_density,
    Set::Patch<const Set::Scalar> T,
    int i, int j, int k, const Set::Scalar* DX, Set::Scalar dt,
    const ThermochemicalData& data, bool include_reaction,
    bool project_negative_gas_density)
{
    const auto& [thermal, chemistry_data, implicit_species_diffusion,
                 include_conduction, implicit_thermal_diffusion] = data;
    const auto& [gas_data, nspecies, ngas_species, density_floor,
                 pressure_reference, condensed_thermal_transport,
                 condensed_specific_heat, condensed_thermal_conductivity,
                 condensed_inverse_reference_density,
                 liquid_inverse_reference_density,
                 condensed_dynamic_viscosity] = thermal;
    (void)liquid_inverse_reference_density;
    (void)condensed_dynamic_viscosity;

    // return values
    Model::Chemistry::SpeciesArray species{};
    Set::Scalar temperature = 0.0;
    Set::Scalar dilatation = 0.0;

    Model::Chemistry::SpeciesArray rhoY{};
    Set::Scalar gas_density = 0.0;
    Set::Scalar molar_density = 0.0;

    // calculate density (and molar density) of gas species only
    for (int n = 0; n < ngas_species; ++n)
    {
        rhoY[n] = project_negative_gas_density ?
            Util::Max(component_density(i,j,k,n), 0.0) :
            component_density(i,j,k,n);
        gas_density += rhoY[n];
        molar_density += rhoY[n] /
            Model::Gas::Gas::MolecularWeight(gas_data, n);
    }

    // The raw EOS volume converts extrinsic partial densities to intrinsic
    // gas densities.  Only the occupied/reacting volume is bounded by the
    // space available outside condensed phases.
    Set::Scalar raw_gas_volume_fraction = 0.0;
    if (gas_density > density_floor && T(i,j,k) > 0.0 && pressure_reference > 0.0)
        raw_gas_volume_fraction = Util::Max(
            gas_density * Model::Gas::Gas::GasConstant(
                gas_data, component_density, i, j, k) *
                T(i,j,k) / pressure_reference,
            0.0);
    Set::Scalar condensed_volume_fraction = 0.0;
    for (int n = ngas_species; n < nspecies; ++n)
        condensed_volume_fraction += Util::Max(
            component_density(i,j,k,n), 0.0) *
            condensed_inverse_reference_density[n];
    const Set::Scalar gas_accessibility =
        1.0 - Model::PhaseField::H(condensed_volume_fraction);
    const Set::Scalar gas_volume_fraction =
        Util::Min(raw_gas_volume_fraction, gas_accessibility);
    const Set::Scalar reacting_volume_fraction = gas_volume_fraction;
    const Set::Scalar chemistry_weight = raw_gas_volume_fraction > 0.0 ?
        reacting_volume_fraction / raw_gas_volume_fraction : 0.0;

    // Calculate relative density with respect to the gas volume fraction
    Model::Chemistry::SpeciesArray intrinsic_rhoY{};
    if (raw_gas_volume_fraction > 0.0)
        for (int n = 0; n < ngas_species; ++n)
            intrinsic_rhoY[n] = rhoY[n] / raw_gas_volume_fraction;

    Model::Chemistry::Source reaction{};
    if (include_reaction && reacting_volume_fraction > 0.0)
    {
        const auto& [gross_model, model, solver, host_chemistry, host_gas] =
            chemistry_data;
        (void)solver;
#ifdef ALAMO_GPU
        if (gross_model)
            reaction = model.ComputeChemistrySources(
                pressure_reference, T(i,j,k), intrinsic_rhoY,
                dt * chemistry_weight, nullptr);
#else
        reaction = host_chemistry->ComputeChemistrySources(
            pressure_reference, T(i,j,k), intrinsic_rhoY,
            dt * chemistry_weight, host_gas);
#endif
    }

    for (int n = 0; n < ngas_species; ++n)
        species[n] = reacting_volume_fraction * reaction.first[n];

    Set::Scalar enthalpy_diffusion = 0.0;
    if (ngas_species > 1 && !implicit_species_diffusion)
    {
        auto gas_density_at = [=] AMREX_GPU_HOST_DEVICE(
            int ii, int jj, int kk)
        {
            Set::Scalar value = 0.0;
            for (int n = 0; n < ngas_species; ++n)
                value += component_density(ii,jj,kk,n);
            return value;
        };
        auto mass_fraction_at = [=] AMREX_GPU_HOST_DEVICE(
            int ii, int jj, int kk, int n)
        {
            const Set::Scalar value = gas_density_at(ii,jj,kk);
            return value > density_floor ?
                component_density(ii,jj,kk,n) / value :
                (n == 0 ? 1.0 : 0.0);
        };
        auto transport_coefficient_at = [=] AMREX_GPU_HOST_DEVICE(
            int ii, int jj, int kk, int n)
        {
            const Set::Scalar local_gas_density =
                gas_density_at(ii,jj,kk);
            const Set::Scalar local_raw_gas_volume_fraction = Util::Max(
                local_gas_density * Model::Gas::Gas::GasConstant(
                    gas_data, component_density, ii, jj, kk) *
                    T(ii,jj,kk) / pressure_reference, 0.0);
            Set::Scalar local_condensed_volume_fraction = 0.0;
            for (int m = ngas_species; m < nspecies; ++m)
                local_condensed_volume_fraction += Util::Max(
                    component_density(ii,jj,kk,m), 0.0) *
                    condensed_inverse_reference_density[m];
            const Set::Scalar local_accessibility =
                1.0 - Model::PhaseField::H(
                    local_condensed_volume_fraction);
            const Set::Scalar weight = local_raw_gas_volume_fraction > 0.0 ?
                Util::Min(1.0, local_accessibility /
                    local_raw_gas_volume_fraction) : 0.0;
            return weight * local_gas_density *
                Model::Gas::Gas::DiffusionCoefficient(
                    gas_data, T(ii,jj,kk), pressure_reference,
                    component_density, ii, jj, kk, n);
        };
        const bool has_condensed_species = nspecies > ngas_species;
        auto face_coefficient = [=] AMREX_GPU_HOST_DEVICE(
            Set::Scalar a, Set::Scalar b)
        {
            if (!has_condensed_species) return 0.5 * (a + b);
            return a > 0.0 && b > 0.0 ? 2.0 * a * b / (a + b) : 0.0;
        };
        auto species_flux = [=] AMREX_GPU_HOST_DEVICE(
            int ia, int ja, int ka, int ib, int jb, int kb,
            int d, int n)
        {
            const Set::Scalar rhoDa =
                transport_coefficient_at(ia,ja,ka,n);
            const Set::Scalar rhoDb =
                transport_coefficient_at(ib,jb,kb,n);
            const Set::Scalar raw_flux = face_coefficient(rhoDa, rhoDb) *
                (mass_fraction_at(ib,jb,kb,n) - mass_fraction_at(ia,ja,ka,n)) /
                DX[d];

            Set::Scalar total_flux = 0.0;
            for (int m = 0; m < ngas_species; ++m)
            {
                const Set::Scalar rhoDma =
                    transport_coefficient_at(ia,ja,ka,m);
                const Set::Scalar rhoDmb =
                    transport_coefficient_at(ib,jb,kb,m);
                total_flux += face_coefficient(rhoDma, rhoDmb) *
                    (mass_fraction_at(ib,jb,kb,m) -
                    mass_fraction_at(ia,ja,ka,m)) / DX[d];
            }
            const Set::Scalar face_fraction = 0.5 *
                (mass_fraction_at(ia,ja,ka,n) + mass_fraction_at(ib,jb,kb,n));
            return raw_flux - face_fraction * total_flux;
        };
        auto conservative_species_flux = [=] AMREX_GPU_HOST_DEVICE(
            int ia, int ja, int ka, int ib, int jb, int kb,
            int d, int n)
        {
            if (n < ngas_species - 1)
                return species_flux(ia, ja, ka, ib, jb, kb, d, n);
            Set::Scalar flux = 0.0;
            for (int m = 0; m < ngas_species - 1; ++m)
                flux -= species_flux(ia, ja, ka, ib, jb, kb, d, m);
            return flux;
        };

        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            const int di = d == 0;
            const int dj = d == 1;
            const int dk = d == 2;
            for (int n = 0; n < ngas_species; ++n)
            {
                const Set::Scalar flux_hi = conservative_species_flux(
                    i, j, k, i+di, j+dj, k+dk, d, n);
                const Set::Scalar flux_lo = conservative_species_flux(
                    i-di, j-dj, k-dk, i, j, k, d, n);
                species[n] += (flux_hi - flux_lo) / DX[d];

                const Set::Scalar h =
                    Model::Gas::Gas::EnthalpyMassSpecies(
                        gas_data, T(i,j,k), n);
                const Set::Scalar h_hi = 0.5 * (h +
                    Model::Gas::Gas::EnthalpyMassSpecies(
                        gas_data, T(i+di,j+dj,k+dk), n));
                const Set::Scalar h_lo = 0.5 * (
                    Model::Gas::Gas::EnthalpyMassSpecies(
                        gas_data, T(i-di,j-dj,k-dk), n) + h);
                enthalpy_diffusion +=
                    ((h_hi - h) * flux_hi + (h - h_lo) * flux_lo) / DX[d];
            }
        }
    }

    auto [thermal_gas_volume_fraction, gas_heat_capacity, heat_capacity,
        conductivity, cp] = ComputeThermalState(
            component_density, T(i,j,k), i, j, k, thermal);
    (void)thermal_gas_volume_fraction;
    (void)gas_heat_capacity;
    (void)cp;
    if (include_conduction && !implicit_thermal_diffusion &&
        heat_capacity > 0.0)
    {
        temperature += conductivity / heat_capacity *
            Numeric::Laplacian(T, i, j, k, 0, DX);
    }
    if (heat_capacity > 0.0)
        temperature +=
            (reacting_volume_fraction * reaction.second + enthalpy_diffusion) /
            heat_capacity;

    if (T(i,j,k) > 0.0)
        dilatation += gas_volume_fraction * temperature / T(i,j,k);
    if (molar_density > 0.0)
    {
        Set::Scalar molar_source = 0.0;
        for (int n = 0; n < ngas_species; ++n)
            molar_source += species[n] /
                Model::Gas::Gas::MolecularWeight(gas_data, n);
        dilatation += gas_volume_fraction * molar_source / molar_density;
    }
    return {species,temperature,dilatation};
}

//
// Compute the mass-weighted rigid velocity associated with each freely moving
// condensed species.  Coarse cells covered by a finer AMR level are excluded
// so that every body's moments describe the composite hierarchy exactly once.
//
void
LowMach::UpdateFreeRigidBodyStates()
{
    BL_PROFILE("Integrator::LowMach::UpdateFreeRigidBodyStates");
    if (free_rigid_solid_species.empty()) return;

    for (int free_body_index = 0;
         free_body_index < static_cast<int>(free_rigid_solid_species.size());
         ++free_body_index)
    {
        // Keep the preceding center as an unwrapping reference.  A compact body
        // crossing a periodic boundary must not suddenly acquire a domain-sized
        // radius of gyration or a center near the middle of the box.
        const FreeRigidBodyState reference_body =
            free_rigid_bodies[free_body_index];
        free_rigid_bodies[free_body_index] = FreeRigidBodyState{};
        FreeRigidBodyState& body = free_rigid_bodies[free_body_index];
        const int free_component =
            free_rigid_solid_components[free_body_index];

        // M, Mx[3], P[3], raw symmetric second moment[6], raw angular momentum[3].
        std::array<Set::Scalar, 16> moment{};
        const int component_count = nspecies;

        Set::Vector unwrapping_center = reference_body.center;
        bool have_unwrapping_center = reference_body.valid;
        if (!have_unwrapping_center)
        {
            // A body may already straddle a periodic boundary in its initial
            // condition.  Obtain a topology-aware first center from circular
            // moments; subsequent calls use the cheaper preceding-center path.
            std::array<Set::Scalar, 6> circular_moment{};
            for (int lev = 0; lev <= finest_level; ++lev)
            {
                std::unique_ptr<amrex::iMultiFab> uncovered;
                if (lev < finest_level)
                    uncovered = std::make_unique<amrex::iMultiFab>(
                        amrex::makeFineMask(*component_density_mf[lev],
                            component_density_mf[lev + 1]->boxArray(),
                            refRatio(lev), geom[lev].periodicity(), 1, 0));

                const auto prob_lo = geom[lev].ProbLoArray();
                const auto dx = geom[lev].CellSizeArray();
                amrex::GpuArray<Set::Scalar,AMREX_SPACEDIM> period{};
                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    period[d] = geom[lev].ProbLength(d);
                Set::Scalar cell_volume = 1.0;
                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    cell_volume *= dx[d];

                for (amrex::MFIter mfi(*component_density_mf[lev],
                        amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const amrex::Box& bx = mfi.tilebox();
                    Set::Patch<const Set::Scalar> component_density =
                        component_density_mf.Patch(lev,mfi);
                    Set::Patch<const Set::Scalar> rigid_species_eta =
                        rigid_species_eta_mf.Patch(lev,mfi);
                    amrex::Array4<const int> uncovered_patch;
                    const bool has_fine_coverage = uncovered != nullptr;
                    if (has_fine_coverage)
                        uncovered_patch = uncovered->const_array(mfi);

                    using Sum = amrex::ReduceOpSum;
                    amrex::ReduceOps<Sum,Sum,Sum,Sum,Sum,Sum> reduce_op;
                    amrex::ReduceData<Set::Scalar,Set::Scalar,Set::Scalar,
                        Set::Scalar,Set::Scalar,Set::Scalar> reduce_data(reduce_op);
                    using ReduceTuple = typename decltype(reduce_data)::Type;
                    reduce_op.eval(bx, reduce_data,
                        [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple
                    {
                        if (has_fine_coverage && uncovered_patch(i,j,k) == 0)
                            return {0.0,0.0,0.0,0.0,0.0,0.0};
                        Set::Scalar mixture_density = 0.0;
                        for (int n = 0; n < component_count; ++n)
                            mixture_density += Util::Max(
                                component_density(i,j,k,n), 0.0);
                        const Set::Scalar indicator = Model::PhaseField::H(
                            Util::Clamp(rigid_species_eta(
                                i,j,k,free_component), 0.0, 1.0));
                        const Set::Scalar dm =
                            mixture_density * indicator * cell_volume;
                        const Set::Scalar ax = 2.0 * Set::Constant::Pi *
                            (prob_lo[0] + (i + 0.5) * dx[0] - prob_lo[0]) /
                            period[0];
                        Set::Scalar ay = 0.0, az = 0.0;
#if AMREX_SPACEDIM > 1
                        ay = 2.0 * Set::Constant::Pi *
                            (prob_lo[1] + (j + 0.5) * dx[1] - prob_lo[1]) /
                            period[1];
#endif
#if AMREX_SPACEDIM > 2
                        az = 2.0 * Set::Constant::Pi *
                            (prob_lo[2] + (k + 0.5) * dx[2] - prob_lo[2]) /
                            period[2];
#endif
                        return {dm*std::cos(ax), dm*std::sin(ax),
                                dm*std::cos(ay), dm*std::sin(ay),
                                dm*std::cos(az), dm*std::sin(az)};
                    });
                    ReduceTuple value = reduce_data.value();
                    circular_moment[0] += amrex::get<0>(value);
                    circular_moment[1] += amrex::get<1>(value);
                    circular_moment[2] += amrex::get<2>(value);
                    circular_moment[3] += amrex::get<3>(value);
                    circular_moment[4] += amrex::get<4>(value);
                    circular_moment[5] += amrex::get<5>(value);
                }
            }
            amrex::ParallelDescriptor::ReduceRealSum(
                circular_moment.data(), circular_moment.size());
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
                if (geom[0].isPeriodic(d))
                {
                    Set::Scalar angle = std::atan2(
                        circular_moment[2*d + 1], circular_moment[2*d]);
                    if (angle < 0.0) angle += 2.0 * Set::Constant::Pi;
                    unwrapping_center(d) = geom[0].ProbLo(d) +
                        geom[0].ProbLength(d) * angle /
                            (2.0 * Set::Constant::Pi);
                }
            have_unwrapping_center = true;
        }

        for (int lev = 0; lev <= finest_level; ++lev)
        {
            std::unique_ptr<amrex::iMultiFab> uncovered;
            if (lev < finest_level)
                uncovered = std::make_unique<amrex::iMultiFab>(
                    amrex::makeFineMask(*component_density_mf[lev],
                        component_density_mf[lev + 1]->boxArray(), refRatio(lev),
                        geom[lev].periodicity(), 1, 0));

            const auto prob_lo = geom[lev].ProbLoArray();
            const auto dx = geom[lev].CellSizeArray();
            amrex::GpuArray<Set::Scalar,AMREX_SPACEDIM> period{};
            amrex::GpuArray<int,AMREX_SPACEDIM> periodic{};
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
            {
                period[d] = geom[lev].ProbLength(d);
                periodic[d] = geom[lev].isPeriodic(d);
            }
            const bool unwrap_about_reference = have_unwrapping_center;
            const Set::Vector reference_center = unwrapping_center;
            Set::Scalar cell_volume = 1.0;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) cell_volume *= dx[d];

            for (amrex::MFIter mfi(*component_density_mf[lev],
                    amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                Set::Patch<const Set::Scalar> component_density =
                    component_density_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> rigid_species_eta =
                    rigid_species_eta_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> u = velocity_mf.Patch(lev,mfi);
                amrex::Array4<const int> uncovered_patch;
                const bool has_fine_coverage = uncovered != nullptr;
                if (has_fine_coverage) uncovered_patch = uncovered->const_array(mfi);

                using Sum = amrex::ReduceOpSum;
                amrex::ReduceOps<Sum,Sum,Sum,Sum,Sum,Sum,Sum,Sum,
                                Sum,Sum,Sum,Sum,Sum,Sum,Sum,Sum> reduce_op;
                amrex::ReduceData<Set::Scalar,Set::Scalar,Set::Scalar,Set::Scalar,
                                Set::Scalar,Set::Scalar,Set::Scalar,Set::Scalar,
                                Set::Scalar,Set::Scalar,Set::Scalar,Set::Scalar,
                                Set::Scalar,Set::Scalar,Set::Scalar,Set::Scalar>
                    reduce_data(reduce_op);
                using ReduceTuple = typename decltype(reduce_data)::Type;

                reduce_op.eval(bx, reduce_data,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple
                {
                    if (has_fine_coverage && uncovered_patch(i,j,k) == 0)
                        return {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,
                                0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};

                    Set::Scalar mixture_density = 0.0;
                    for (int n = 0; n < component_count; ++n)
                        mixture_density += Util::Max(component_density(i,j,k,n), 0.0);
                    const Set::Scalar indicator = Model::PhaseField::H(
                        Util::Clamp(rigid_species_eta(
                            i,j,k,free_component), 0.0, 1.0));
                    const Set::Scalar dm = mixture_density * indicator * cell_volume;
                    Set::Scalar x = prob_lo[0] + (i + 0.5) * dx[0];
                    Set::Scalar y = 0.0, z = 0.0;
                    Set::Scalar ux = u(i,j,k,0), uy = 0.0, uz = 0.0;
#if AMREX_SPACEDIM > 1
                    y = prob_lo[1] + (j + 0.5) * dx[1];
                    uy = u(i,j,k,1);
#endif
#if AMREX_SPACEDIM > 2
                    z = prob_lo[2] + (k + 0.5) * dx[2];
                    uz = u(i,j,k,2);
#endif
                    if (unwrap_about_reference)
                    {
                        if (periodic[0])
                            x = reference_center(0) + MinimumImageDisplacement(
                                x - reference_center(0), period[0], periodic[0]);
#if AMREX_SPACEDIM > 1
                        if (periodic[1])
                            y = reference_center(1) + MinimumImageDisplacement(
                                y - reference_center(1), period[1], periodic[1]);
#endif
#if AMREX_SPACEDIM > 2
                        if (periodic[2])
                            z = reference_center(2) + MinimumImageDisplacement(
                                z - reference_center(2), period[2], periodic[2]);
#endif
                    }
                    return {
                        dm,
                        dm*x, dm*y, dm*z,
                        dm*ux, dm*uy, dm*uz,
                        dm*x*x, dm*x*y, dm*x*z,
                        dm*y*y, dm*y*z, dm*z*z,
                        dm*(y*uz-z*uy), dm*(z*ux-x*uz), dm*(x*uy-y*ux)};
                });

                ReduceTuple value = reduce_data.value();
                moment[0]  += amrex::get<0>(value);
                moment[1]  += amrex::get<1>(value);
                moment[2]  += amrex::get<2>(value);
                moment[3]  += amrex::get<3>(value);
                moment[4]  += amrex::get<4>(value);
                moment[5]  += amrex::get<5>(value);
                moment[6]  += amrex::get<6>(value);
                moment[7]  += amrex::get<7>(value);
                moment[8]  += amrex::get<8>(value);
                moment[9]  += amrex::get<9>(value);
                moment[10] += amrex::get<10>(value);
                moment[11] += amrex::get<11>(value);
                moment[12] += amrex::get<12>(value);
                moment[13] += amrex::get<13>(value);
                moment[14] += amrex::get<14>(value);
                moment[15] += amrex::get<15>(value);
            }
        }
        amrex::ParallelDescriptor::ReduceRealSum(moment.data(), moment.size());

        const Set::Scalar mass = moment[0];
        if (!(mass > density_floor) || !std::isfinite(mass)) continue;

        body.mass = mass;
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            body.center(d) = moment[1 + d] / mass;
            body.velocity(d) = moment[4 + d] / mass;
        }

        Set::Matrix central_second_moment = Set::Matrix::Zero();
        central_second_moment(0,0) = moment[7] -
            mass * body.center(0) * body.center(0);
#if AMREX_SPACEDIM > 1
        central_second_moment(0,1) = central_second_moment(1,0) = moment[8] -
            mass * body.center(0) * body.center(1);
        central_second_moment(1,1) = moment[10] -
            mass * body.center(1) * body.center(1);
#endif
#if AMREX_SPACEDIM > 2
        central_second_moment(0,2) = central_second_moment(2,0) = moment[9] -
            mass * body.center(0) * body.center(2);
        central_second_moment(1,2) = central_second_moment(2,1) = moment[11] -
            mass * body.center(1) * body.center(2);
        central_second_moment(2,2) = moment[12] -
            mass * body.center(2) * body.center(2);
#endif
        body.inertia = central_second_moment.trace() *
            Set::Matrix::Identity() - central_second_moment;
        body.radius_of_gyration = std::sqrt(Util::Max(
            central_second_moment.trace() / mass, 0.0));

#if AMREX_SPACEDIM == 2
        const Set::Scalar angular_momentum = moment[15] -
            (body.center(0) * moment[5] - body.center(1) * moment[4]);
        const Set::Scalar scalar_inertia = central_second_moment.trace();
        if (scalar_inertia > density_floor)
            body.angular_velocity(0) =
                angular_momentum / scalar_inertia;
#elif AMREX_SPACEDIM == 3
        Set::Vector angular_momentum;
        angular_momentum(0) = moment[13] -
            (body.center(1) * moment[6] - body.center(2) * moment[5]);
        angular_momentum(1) = moment[14] -
            (body.center(2) * moment[4] - body.center(0) * moment[6]);
        angular_momentum(2) = moment[15] -
            (body.center(0) * moment[5] - body.center(1) * moment[4]);
        if (Util::Abs(body.inertia.determinant()) > density_floor)
            body.angular_velocity = body.inertia.inverse() * angular_momentum;
#endif
        body.valid = true;
    }
}


void
LowMach::ProjectVelocity(Set::Scalar time, Set::Scalar dt)
{
    BL_PROFILE("Integrator::LowMach::ProjectVelocity");
    if (!projection_enabled || !(dt > 0.0)) return;

    const int nlev = finest_level + 1;
    const bool deformable_solid = deformable_solid_species >= 0;
    const bool rigid_solid = !rigid_solid_species.empty();
    const int nfree = static_cast<int>(free_rigid_solid_species.size());
    const bool free_rigid_solid = nfree > 0;
    const bool liquid = !liquid_species.empty();
    const bool capillary = interfacial_forces_enabled;
    const bool mixed_phase = deformable_solid || rigid_solid || liquid;
    const bool split_chemistry = chemistry.Split();
    const bool split_diffusion =
        implicit_thermal_diffusion || implicit_species_diffusion;
    const bool split_phase_change = has_split_phase_change;
    const int nrigid = static_cast<int>(rigid_solid_species.size());
    int projection_iterations = 1;
    for (const int component : free_rigid_solid_components)
        projection_iterations = Util::Max(projection_iterations,
            rigid_coupling_max_iterations[component]);
    Model::Chemistry::SpeciesArray rigid_inverse_relaxation_time{};
    std::array<Model::Chemistry::SpeciesArray, AMREX_SPACEDIM>
        rigid_prescribed_velocity{};
    for (int m = 0; m < nrigid; ++m)
    {
        rigid_inverse_relaxation_time[m] =
            1.0 / rigid_relaxation_time[m];
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
            rigid_prescribed_velocity[d][m] =
                rigid_velocity[m](d);
    }
    // When every rigid phase has the same relaxation time, the Brinkman
    // mobility depends only on aggregate rigid occupancy.  If all phases are
    // also fixed at the same velocity, the complete penalty update can use the
    // aggregate field and avoid a per-cell loop over rigid species.  The
    // common AP/HTPB configuration follows this path and recovers the original
    // aggregate-solid work and behavior.
    bool shared_rigid_relaxation = rigid_solid;
    bool shared_fixed_rigid_penalty = rigid_solid && !free_rigid_solid;
    Set::Scalar shared_rigid_inverse_relaxation_time = 0.0;
    Set::Vector shared_rigid_velocity = Set::Vector::Zero();
    if (rigid_solid)
    {
        shared_rigid_inverse_relaxation_time =
            rigid_inverse_relaxation_time[0];
        shared_rigid_velocity = rigid_velocity[0];
        for (int m = 1; m < nrigid; ++m)
        {
            shared_rigid_relaxation = shared_rigid_relaxation &&
                rigid_relaxation_time[m] == rigid_relaxation_time[0];
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
                shared_fixed_rigid_penalty =
                    shared_fixed_rigid_penalty &&
                    rigid_velocity[m](d) == rigid_velocity[0](d);
        }
        shared_fixed_rigid_penalty = shared_fixed_rigid_penalty &&
            shared_rigid_relaxation;
    }
    if (!(pressure_reference == pressure_reference))
        pressure_reference = pressure_mf[0]->sum(0, false) /
                            static_cast<Set::Scalar>(geom[0].Domain().numPts());
    amrex::get<4>(thermal_data) = pressure_reference;
    amrex::get<0>(thermochemical_data) = thermal_data;
    const auto thermal = thermal_data;
    const auto gas_data = gas_device_data;
    const auto thermochemical = thermochemical_data;
    const int number_of_species = nspecies;
    const Set::Scalar* liquid_inverse_density =
        amrex::get<9>(thermal_data);
    const amrex::BCRec pressure_boundary = pressure_bc->GetBCRec();
    amrex::GpuArray<int, AMREX_SPACEDIM> pressure_outlet_lo{};
    amrex::GpuArray<int, AMREX_SPACEDIM> pressure_outlet_hi{};
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
    {
        pressure_outlet_lo[d] = !geom[0].isPeriodic(d) &&
            BC::BCUtil::IsDirichlet(pressure_boundary.lo(d));
        pressure_outlet_hi[d] = !geom[0].isPeriodic(d) &&
            BC::BCUtil::IsDirichlet(pressure_boundary.hi(d));
    }
    pressure_poisson.SetLayout(geom, refRatio(), velocity_mf, nlev);
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> projection_source(nlev);
    amrex::Vector<std::unique_ptr<amrex::MultiFab>> provisional_velocity(nlev);

    for (int lev = 0; lev < nlev; ++lev)
    {
        velocity_bc->define(geom[lev]);
        velocity_bc->FillBoundary(*velocity_mf[lev], 0, AMREX_SPACEDIM, time, 0);
        velocity_mf[lev]->FillBoundary(geom[lev].periodicity());
        UpdateComponentState(lev, *component_density_mf[lev]);
    }

    // Preserve the momentum already accumulated by advection and diffusion,
    // then predict only the missing pressure/constraint impulse from the
    // preceding step.  This is substantially more accurate than extrapolating
    // old final velocities because the fresh provisional state is retained.
    std::vector<FreeRigidBodyState> provisional_body_states;
    std::vector<FreeRigidBodyState> target_bodies;
    if (free_rigid_solid)
    {
        UpdateFreeRigidBodyStates();
        provisional_body_states = free_rigid_bodies;
        target_bodies = provisional_body_states;
        for (int b = 0; b < nfree; ++b)
        {
            FreeRigidBodyState& target_body = target_bodies[b];
            const FreeRigidBodyState& previous_increment =
                previous_free_rigid_projection_increments[b];
            if (!target_body.valid || !previous_increment.valid ||
                !(previous_free_rigid_dt > 0.0))
                continue;
            const Set::Scalar step_ratio = dt / previous_free_rigid_dt;
            const Set::Scalar mass_scale = Util::Max(
                Util::Max(target_body.mass, previous_increment.mass),
                density_floor);
            const Set::Scalar radius_scale = Util::Max(
                Util::Max(target_body.radius_of_gyration,
                    previous_increment.radius_of_gyration),
                1.0e-12);
            const bool compatible_history =
                step_ratio >= 0.25 && step_ratio <= 4.0 &&
                Util::Abs(target_body.mass -
                    previous_increment.mass) <=
                    0.25 * mass_scale &&
                Util::Abs(target_body.radius_of_gyration -
                    previous_increment.radius_of_gyration) <=
                    0.25 * radius_scale;
            if (compatible_history)
            {
                target_body.velocity += step_ratio *
                    previous_increment.velocity;
                target_body.angular_velocity += step_ratio *
                    previous_increment.angular_velocity;
            }
        }

        if (projection_iterations > 1)
        {
            for (int lev = 0; lev < nlev; ++lev)
            {
                provisional_velocity[lev] = std::make_unique<amrex::MultiFab>(
                    velocity_mf[lev]->boxArray(),
                    velocity_mf[lev]->DistributionMap(), AMREX_SPACEDIM,
                    velocity_mf[lev]->nGrow());
                amrex::MultiFab::Copy(*provisional_velocity[lev],
                    *velocity_mf[lev], 0, 0, AMREX_SPACEDIM,
                    velocity_mf[lev]->nGrow());
            }
        }
    }

    for (int lev = 0; lev < nlev; ++lev)
    {
        const auto dx = geom[lev].CellSizeArray();
        const amrex::Dim3 domain_lo = amrex::lbound(geom[lev].Domain());
        const amrex::Dim3 domain_hi = amrex::ubound(geom[lev].Domain());
        amrex::MultiFab& beta_mf = pressure_poisson.Coefficient(lev);
        const Set::Scalar rho_floor = density_floor;

        pressure_poisson.RHS(lev).setVal(0.0);
        for (const auto& configured_mechanism : mechanisms)
        {
            const auto mechanism = configured_mechanism;
            const Set::Scalar p_reference = pressure_reference;
            for (amrex::MFIter mfi(pressure_poisson.RHS(lev), amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                Set::Patch<const Set::Scalar> T = temperature_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> component_density = component_density_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> rigid_eta = rigid_eta_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> rigid_species_eta =
                    rigid_species_eta_mf.Patch(lev,mfi);
                Set::Patch<Set::Scalar> rhs = pressure_poisson.RHS(lev).array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    const Model::Mechanism::State state = {
                        component_density, rigid_eta, rigid_species_eta,
                        T(i,j,k), p_reference, dt};
                    rhs(i,j,k) += mechanism.VolumeSource(
                        state, i, j, k, dx.data());
                    const Set::Scalar heat_source = mechanism.HeatSource(
                        state, i, j, k, dx.data());
                    if (heat_source != 0.0 && T(i,j,k) > 0.0)
                    {
                        auto [gas_volume_fraction, gas_heat_capacity,
                            heat_capacity, conductivity, cp] =
                            ComputeThermalState(component_density,
                                T(i,j,k), i, j, k, thermal);
                        (void)gas_heat_capacity;
                        (void)conductivity;
                        (void)cp;
                        if (heat_capacity > 0.0)
                            rhs(i,j,k) += gas_volume_fraction * heat_source /
                                (heat_capacity * T(i,j,k));
                    }
                });
            }
        }

        for (amrex::MFIter mfi(pressure_poisson.RHS(lev), amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            Set::Patch<const Set::Scalar> T = temperature_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> component_density = component_density_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> chemistry_dilatation =
                chemistry_dilatation_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> phase_change_dilatation;
            if (split_phase_change)
                phase_change_dilatation =
                    phase_change_dilatation_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> diffusion_dilatation;
            if (split_diffusion)
                diffusion_dilatation = diffusion_dilatation_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> eta = eta_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> rigid_eta = rigid_eta_mf.Patch(lev,mfi);
            Set::Patch<Set::Scalar> rhs = pressure_poisson.RHS(lev).array(mfi);
            const int ngas = ngas_species;
            const Set::Scalar p_reference = pressure_reference;
            const Set::Scalar inverse_dt = 1.0 / dt;

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                auto [species,temperature,dilatation] = ComputeThermochemicalSource(
                    component_density, T, i, j, k, dx.data(), dt,
                    thermochemical, !split_chemistry);

                rhs(i,j,k) += dilatation;
                if (split_chemistry)
                    rhs(i,j,k) += chemistry_dilatation(i,j,k) / dt;
                if (split_phase_change)
                    rhs(i,j,k) += phase_change_dilatation(i,j,k) / dt;
                if (split_diffusion)
                    rhs(i,j,k) += diffusion_dilatation(i,j,k) / dt;

                if (mixed_phase)
                {
                    Set::Scalar gas_density = 0.0;
                    for (int n = 0; n < ngas; ++n)
                        gas_density += component_density(i,j,k,n);
                    Set::Scalar volume_fraction = gas_density *
                        Model::Gas::Gas::GasConstant(
                            gas_data, component_density, i, j, k) *
                        T(i,j,k) / p_reference;
                    if (deformable_solid) volume_fraction += eta(i,j,k);
                    if (rigid_solid) volume_fraction += rigid_eta(i,j,k);
                    if (liquid)
                        for (int n = ngas; n < number_of_species; ++n)
                            volume_fraction += Util::Max(
                                component_density(i,j,k,n), 0.0) *
                                liquid_inverse_density[n];

                    // Correct only the unexplained part of the mixture-volume
                    // defect.  Split operators have already placed their
                    // exact integrated volume changes in both the state and
                    // the projection source above; counting those changes a
                    // second time produces an impulsive velocity whenever a
                    // dense phase melts or evaporates.
                    Set::Scalar unexplained_volume_defect =
                        volume_fraction - 1.0;
                    if (split_chemistry)
                        unexplained_volume_defect -=
                            chemistry_dilatation(i,j,k);
                    if (split_phase_change)
                        unexplained_volume_defect -=
                            phase_change_dilatation(i,j,k);
                    if (split_diffusion)
                        unexplained_volume_defect -=
                            diffusion_dilatation(i,j,k);
                    // A prescribed-pressure boundary is an open boundary for
                    // the mixture constraint.  A diffuse phase clipped there
                    // has an algebraic volume defect that must leave with the
                    // phase; forcing it to zero locally in one shrinking time
                    // step creates an unbounded outlet velocity.  Fade only
                    // this drift-control term where the diffuse interface can
                    // overlap the outlet, then restore it smoothly over the
                    // next physical interface thickness.  The
                    // thermochemical and phase-change divergence sources
                    // above remain active at the outlet.
                    Set::Scalar outlet_weight = 1.0;
                    if (interfacial_thickness > 0.0)
                    {
                        const int index[AMREX_SPACEDIM] =
                            {AMREX_D_DECL(i,j,k)};
                        const int lo[AMREX_SPACEDIM] =
                            {AMREX_D_DECL(domain_lo.x,domain_lo.y,domain_lo.z)};
                        const int hi[AMREX_SPACEDIM] =
                            {AMREX_D_DECL(domain_hi.x,domain_hi.y,domain_hi.z)};
                        for (int d = 0; d < AMREX_SPACEDIM; ++d)
                        {
                            // Two interface thicknesses represent the
                            // fully restored correction.  Only a physical
                            // pressure outlet can reduce this distance;
                            // periodic directions and closed boundaries must
                            // leave the weight unchanged.
                            Set::Scalar distance =
                                2.0 * interfacial_thickness;
                            if (pressure_outlet_lo[d])
                                distance = Util::Min(distance,
                                    (index[d] - lo[d] + 0.5) * dx[d]);
                            if (pressure_outlet_hi[d])
                                distance = Util::Min(distance,
                                    (hi[d] - index[d] + 0.5) * dx[d]);
                            const Set::Scalar q = Util::Clamp(
                                distance / interfacial_thickness - 1.0,
                                0.0, 1.0);
                            outlet_weight *= q * q * (3.0 - 2.0 * q);
                        }
                    }
                    rhs(i,j,k) += outlet_weight *
                        unexplained_volume_defect * inverse_dt /
                        Util::Max(volume_fraction, 0.1);
                }
            });
        }

        for (amrex::MFIter mfi(pressure_poisson.RHS(lev), amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            Set::Patch<const Set::Scalar> rho = density_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> rigid_eta;
            Set::Patch<const Set::Scalar> rigid_species_eta;
            if (rigid_solid)
            {
                rigid_eta = rigid_eta_mf.Patch(lev,mfi);
                rigid_species_eta = rigid_species_eta_mf.Patch(lev,mfi);
            }
            Set::Patch<Set::Scalar> beta = beta_mf.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar penalty_rate = 0.0;
                const Set::Scalar rigid_composition_sum = rigid_solid ?
                    Util::Max(rigid_eta(i,j,k), 0.0) : 0.0;
                const Set::Scalar aggregate_eta =
                    Util::Clamp(rigid_composition_sum, 0.0, 1.0);
                const Set::Scalar aggregate_weight =
                    Model::PhaseField::H(aggregate_eta);
                if (shared_rigid_relaxation)
                    penalty_rate = aggregate_weight *
                        shared_rigid_inverse_relaxation_time;
                else if (rigid_composition_sum > 0.0)
                {
                    for (int m = 0; m < nrigid; ++m)
                    {
                        const Set::Scalar composition = Util::Clamp(
                            rigid_species_eta(i,j,k,m) /
                                rigid_composition_sum,
                            0.0, 1.0);
                        penalty_rate += aggregate_weight * composition *
                            rigid_inverse_relaxation_time[m];
                    }
                }
                const Set::Scalar mobility =
                    1.0 / (1.0 + dt * penalty_rate);
                beta(i,j,k) = mobility / Util::Max(rho(i,j,k), rho_floor);
            });
        }
        if (free_rigid_solid)
        {
            projection_source[lev] = std::make_unique<amrex::MultiFab>(
                pressure_poisson.RHS(lev).boxArray(),
                pressure_poisson.RHS(lev).DistributionMap(), 1, 0);
            amrex::MultiFab::Copy(*projection_source[lev],
                pressure_poisson.RHS(lev), 0, 0, 1, 0);
        }
    }

    // Capillarity and pressure use the same face gradient and face mobility.
    // The variational mu_i grad(eta_i) force is the discrete counterpart of
    // the configured multiphase interfacial free energy.  It is well-defined
    // where
    // |grad eta| is small and remains balanced with the face-centered
    // projection.
    pressure_poisson.PrepareCoefficients(time);
    using OwnedFaceField =
        amrex::Array<std::unique_ptr<amrex::MultiFab>, AMREX_SPACEDIM>;
    amrex::Vector<OwnedFaceField> capillary_face_acceleration(nlev);
    if (capillary)
    {
        for (int lev = 0; lev < nlev; ++lev)
        {
            UpdateInterfacialChemicalPotential(lev);
            const Set::Scalar* DX = geom[lev].CellSize();
            const amrex::Box domain = geom[lev].Domain();
            const int nphase = interfacial_model.NumberOfPhases();

            for (int d = 0; d < AMREX_SPACEDIM; ++d)
            {
                amrex::BoxArray face_grids = velocity_mf[lev]->boxArray();
                face_grids.surroundingNodes(d);
                capillary_face_acceleration[lev][d] =
                    std::make_unique<amrex::MultiFab>(
                        face_grids, velocity_mf[lev]->DistributionMap(), 1, 0);
                amrex::MultiFab& capillary_acceleration =
                    *capillary_face_acceleration[lev][d];
                const amrex::MultiFab& beta =
                    pressure_poisson.FaceCoefficient(lev, d);
                const int di = d == 0;
                const int dj = d == 1;
                const int dk = d == 2;
                const int domain_face_lo = domain.smallEnd(d);
                const int domain_face_hi = domain.bigEnd(d) + 1;
                const bool periodic = geom[lev].isPeriodic(d);

                for (amrex::MFIter mfi(capillary_acceleration,
                        amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const amrex::Box& bx = mfi.tilebox();
                    Set::Patch<const Set::Scalar> eta =
                        interfacial_volume_fraction_mf.Patch(lev,mfi);
                    Set::Patch<const Set::Scalar> mu =
                        interfacial_chemical_potential_mf.Patch(lev,mfi);
                    Set::Patch<const Set::Scalar> face_beta =
                        beta.const_array(mfi);
                    Set::Patch<Set::Scalar> capillary_face_force =
                        capillary_acceleration.array(mfi);

                    amrex::ParallelFor(
                        bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                        {
                            const int face_index =
                                d == 0 ? i : (d == 1 ? j : k);
                            if (!periodic && (face_index == domain_face_lo ||
                                             face_index == domain_face_hi))
                            {
                                capillary_face_force(i,j,k) = 0.0;
                                return;
                            }

                            Set::Scalar force = 0.0;
                            for (int phase = 0; phase < nphase; ++phase)
                            {
                                const Set::Scalar mu_face = 0.5 *
                                    (mu(i-di,j-dj,k-dk,phase) +
                                     mu(i,j,k,phase));
                                force += mu_face *
                                    (eta(i,j,k,phase) -
                                     eta(i-di,j-dj,k-dk,phase)) / DX[d];
                            }
                            capillary_face_force(i,j,k) =
                                face_beta(i,j,k) * force;
                        });
                }
            }
        }
    }

    auto ApplyRigidPenalty = [&](const std::vector<FreeRigidBodyState>& bodies)
    {
        Model::Chemistry::SpeciesArray free_rigid_flag{};
        Model::Chemistry::SpeciesArray active_free_rigid{};
        std::array<Model::Chemistry::SpeciesArray, AMREX_SPACEDIM>
            free_center{}, free_translation{}, free_rotation{};
        bool any_active_free_rigid = false;
        for (int b = 0; b < nfree; ++b)
        {
            const int m = free_rigid_solid_components[b];
            free_rigid_flag[m] = 1.0;
            active_free_rigid[m] = bodies[b].valid ? 1.0 : 0.0;
            any_active_free_rigid = any_active_free_rigid || bodies[b].valid;
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
            {
                free_center[d][m] = bodies[b].center(d);
                free_translation[d][m] = bodies[b].velocity(d);
                free_rotation[d][m] = bodies[b].angular_velocity(d);
            }
        }
        for (int lev = 0; lev < nlev; ++lev)
        {
            const auto prob_lo = geom[lev].ProbLoArray();
            const auto dx = geom[lev].CellSizeArray();
            amrex::GpuArray<Set::Scalar,AMREX_SPACEDIM> period{};
            amrex::GpuArray<int,AMREX_SPACEDIM> periodic{};
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
            {
                period[d] = geom[lev].ProbLength(d);
                periodic[d] = geom[lev].isPeriodic(d);
            }
            for (amrex::MFIter mfi(*velocity_mf[lev],
                    amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                Set::Patch<Set::Scalar> u = velocity_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> rigid_eta =
                    rigid_eta_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> rigid_species_eta =
                    rigid_species_eta_mf.Patch(lev,mfi);

                amrex::ParallelFor(bx,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    const Set::Scalar rigid_composition_sum = Util::Max(
                        rigid_eta(i,j,k), 0.0);
                    const Set::Scalar aggregate_eta = Util::Clamp(
                        rigid_composition_sum, 0.0, 1.0);
                    const Set::Scalar aggregate_weight =
                        Model::PhaseField::H(aggregate_eta);
                    if (shared_fixed_rigid_penalty)
                    {
                        const Set::Scalar penalty_rate = aggregate_weight *
                            shared_rigid_inverse_relaxation_time;
                        const Set::Scalar mobility = 1.0 /
                            (1.0 + dt * penalty_rate);
                        for (int d = 0; d < AMREX_SPACEDIM; ++d)
                            u(i,j,k,d) = mobility * u(i,j,k,d) +
                                (1.0 - mobility) *
                                    shared_rigid_velocity(d);
                        return;
                    }

                    Set::Scalar total_penalty_rate = 0.0;
                    Set::Scalar weighted_target[AMREX_SPACEDIM] = {0.0};
                    if (!(rigid_composition_sum > 0.0)) return;
                    for (int m = 0; m < nrigid; ++m)
                    {
                        const Set::Scalar composition = Util::Clamp(
                            rigid_species_eta(i,j,k,m) /
                                rigid_composition_sum,
                            0.0, 1.0);
                        const Set::Scalar weight =
                            aggregate_weight * composition;
                        const bool is_free = free_rigid_flag[m] > 0.5;
                        if (!(weight > 0.0) ||
                            (is_free && !(active_free_rigid[m] > 0.5)))
                            continue;
                        const Set::Scalar penalty_rate = weight *
                            rigid_inverse_relaxation_time[m];
                        total_penalty_rate += penalty_rate;
                        Set::Scalar target[AMREX_SPACEDIM];
                        for (int d = 0; d < AMREX_SPACEDIM; ++d)
                            target[d] = is_free ?
                                free_translation[d][m] :
                                rigid_prescribed_velocity[d][m];
                        if (is_free)
                        {
#if AMREX_SPACEDIM == 2
                            const Set::Scalar x = MinimumImageDisplacement(
                                prob_lo[0] + (i + 0.5) * dx[0] -
                                    free_center[0][m],
                                period[0], periodic[0]);
                            const Set::Scalar y = MinimumImageDisplacement(
                                prob_lo[1] + (j + 0.5) * dx[1] -
                                    free_center[1][m],
                                period[1], periodic[1]);
                            target[0] -= free_rotation[0][m] * y;
                            target[1] += free_rotation[0][m] * x;
#elif AMREX_SPACEDIM == 3
                            const Set::Scalar x = MinimumImageDisplacement(
                                prob_lo[0] + (i + 0.5) * dx[0] -
                                    free_center[0][m],
                                period[0], periodic[0]);
                            const Set::Scalar y = MinimumImageDisplacement(
                                prob_lo[1] + (j + 0.5) * dx[1] -
                                    free_center[1][m],
                                period[1], periodic[1]);
                            const Set::Scalar z = MinimumImageDisplacement(
                                prob_lo[2] + (k + 0.5) * dx[2] -
                                    free_center[2][m],
                                period[2], periodic[2]);
                            target[0] += free_rotation[1][m] * z -
                                free_rotation[2][m] * y;
                            target[1] += free_rotation[2][m] * x -
                                free_rotation[0][m] * z;
                            target[2] += free_rotation[0][m] * y -
                                free_rotation[1][m] * x;
#endif
                        }
                        for (int d = 0; d < AMREX_SPACEDIM; ++d)
                            weighted_target[d] += penalty_rate * target[d];
                    }
                    const Set::Scalar mobility = 1.0 /
                        (1.0 + dt * total_penalty_rate);
                    for (int d = 0; d < AMREX_SPACEDIM; ++d)
                        u(i,j,k,d) = mobility * (u(i,j,k,d) + dt *
                            weighted_target[d]);
                });
            }

            velocity_bc->define(geom[lev]);
            velocity_bc->FillBoundary(*velocity_mf[lev], 0,
                AMREX_SPACEDIM, time, 0);
            velocity_mf[lev]->FillBoundary(geom[lev].periodicity());
        }

        // A composite projection requires covered coarse cells to contain the
        // same penalized predictor represented by the fine level.
        if (any_active_free_rigid)
        {
            for (int lev = nlev - 2; lev >= 0; --lev)
                amrex::average_down(*velocity_mf[lev + 1], *velocity_mf[lev],
                    geom[lev + 1], geom[lev], 0, AMREX_SPACEDIM,
                    refRatio(lev));
            for (int lev = 0; lev < nlev; ++lev)
            {
                velocity_bc->define(geom[lev]);
                velocity_bc->FillBoundary(*velocity_mf[lev], 0,
                    AMREX_SPACEDIM, time, 0);
                velocity_mf[lev]->FillBoundary(geom[lev].periodicity());
            }
        }
    };

    auto RigidBodyDifferenceMagnitude = [](
        const FreeRigidBodyState& a, const FreeRigidBodyState& b,
        Set::Scalar radius)
    {
        if (!a.valid || !b.valid) return Set::Scalar(0.0);
        Set::Scalar magnitude_squared = 0.0;
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            const Set::Scalar velocity = a.velocity(d) - b.velocity(d);
            const Set::Scalar rotation =
                radius * (a.angular_velocity(d) - b.angular_velocity(d));
            magnitude_squared += velocity * velocity + rotation * rotation;
        }
        return std::sqrt(magnitude_squared);
    };

    // The fixed-point residual is attenuated by the implicit Brinkman
    // mobility.  Undo that attenuation when deciding convergence, and allow
    // Aitken acceleration to span the corresponding inverse-mobility scale.
    std::vector<Set::Scalar> coupling_mobility(nfree, 1.0);
    std::vector<Set::Scalar> maximum_anderson_weight(nfree, 1.0);
    std::vector<FreeRigidBodyState> previous_coupling_output(nfree);
    std::vector<Set::Vector> previous_residual_velocity(
        nfree, Set::Vector::Zero());
    std::vector<Set::Vector> previous_residual_rotation(
        nfree, Set::Vector::Zero());
    std::vector<bool> have_previous_coupling_residual(nfree, false);
    for (int b = 0; b < nfree; ++b)
    {
        const int m = free_rigid_solid_components[b];
        coupling_mobility[b] =
            1.0 / (1.0 + dt / rigid_relaxation_time[m]);
        maximum_anderson_weight[b] = 1.0 / coupling_mobility[b];
        last_free_rigid_coupling_iterations[b] = 0;
        last_free_rigid_coupling_residual[b] = 0.0;
        last_free_rigid_coupling_converged[b] = false;
    }
    for (int iteration = 0; iteration < projection_iterations; ++iteration)
    {
        // Coupling corrections are alternative solutions of the same step.
        // Always restart from u* rather than projecting an already corrected
        // velocity a second time.
        if (iteration > 0)
            for (int lev = 0; lev < nlev; ++lev)
                amrex::MultiFab::Copy(*velocity_mf[lev],
                    *provisional_velocity[lev], 0, 0, AMREX_SPACEDIM,
                    velocity_mf[lev]->nGrow());

        if (rigid_solid) ApplyRigidPenalty(target_bodies);

        for (int lev = 0; lev < nlev; ++lev)
        {
            if (free_rigid_solid)
                amrex::MultiFab::Copy(pressure_poisson.RHS(lev),
                    *projection_source[lev], 0, 0, 1, 0);
            Operator::PressurePoisson::FaceField capillary_acceleration;
            const Operator::PressurePoisson::FaceField*
                capillary_acceleration_ptr = nullptr;
            if (capillary)
            {
                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    capillary_acceleration[d] =
                        capillary_face_acceleration[lev][d].get();
                capillary_acceleration_ptr = &capillary_acceleration;
            }
            pressure_poisson.PrepareRHS(
                lev, *velocity_mf[lev], dt, capillary_acceleration_ptr);
        }
        pressure_poisson.Solve(time, pressure_bc->GetBCRec());

        for (int lev = 0; lev < nlev; ++lev)
        {
            pressure_poisson.ApplyCorrection(lev, *velocity_mf[lev], dt);
            velocity_bc->define(geom[lev]);
            velocity_bc->FillBoundary(*velocity_mf[lev], 0,
                AMREX_SPACEDIM, time, 0);
            velocity_mf[lev]->FillBoundary(geom[lev].periodicity());
        }
        for (int lev = nlev - 2; lev >= 0; --lev)
            amrex::average_down(*velocity_mf[lev + 1], *velocity_mf[lev],
                geom[lev + 1], geom[lev], 0, AMREX_SPACEDIM,
                refRatio(lev));

        if (!free_rigid_solid) break;

        UpdateFreeRigidBodyStates();
        bool all_bodies_finished = true;
        for (int b = 0; b < nfree; ++b)
        {
            const int m = free_rigid_solid_components[b];
            const FreeRigidBodyState& body = free_rigid_bodies[b];
            FreeRigidBodyState& target_body = target_bodies[b];
            const bool exhausted =
                iteration + 1 >= rigid_coupling_max_iterations[m];
            last_free_rigid_coupling_iterations[b] = iteration + 1;
            if (!body.valid || !target_body.valid)
            {
                last_free_rigid_coupling_converged[b] =
                    !body.valid && !target_body.valid;
                all_bodies_finished = all_bodies_finished &&
                    (last_free_rigid_coupling_converged[b] || exhausted);
                continue;
            }

            const Set::Scalar residual_radius = Util::Max(
                Util::Max(body.radius_of_gyration,
                    target_body.radius_of_gyration), 1.0e-12);
            Set::Vector residual_velocity = Set::Vector::Zero();
            Set::Vector residual_rotation = Set::Vector::Zero();
            Set::Scalar residual_squared = 0.0;
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
            {
                residual_velocity(d) =
                    body.velocity(d) - target_body.velocity(d);
                residual_rotation(d) = body.angular_velocity(d) -
                    target_body.angular_velocity(d);
                residual_squared +=
                    residual_velocity(d) * residual_velocity(d) +
                    residual_radius * residual_radius *
                    residual_rotation(d) * residual_rotation(d);
            }
            const Set::Scalar residual = std::sqrt(residual_squared);
            const Set::Scalar estimated_error =
                residual / coupling_mobility[b];
            const Set::Scalar correction_scale = Util::Max(
                RigidBodyDifferenceMagnitude(body,
                    provisional_body_states[b], residual_radius),
                RigidBodyDifferenceMagnitude(target_body,
                    provisional_body_states[b], residual_radius));
            const Set::Scalar absolute_tolerance =
                rigid_coupling_absolute_tolerance[m];
            last_free_rigid_coupling_residual[b] = estimated_error /
                Util::Max(correction_scale, absolute_tolerance);
            last_free_rigid_coupling_converged[b] = estimated_error <=
                absolute_tolerance + rigid_coupling_relative_tolerance[m] *
                    correction_scale;
            all_bodies_finished = all_bodies_finished &&
                (last_free_rigid_coupling_converged[b] || exhausted);

            if (!last_free_rigid_coupling_converged[b] && !exhausted)
            {
                if (rigid_anderson_coupling[m] &&
                    have_previous_coupling_residual[b])
                {
                    Set::Scalar numerator = 0.0;
                    Set::Scalar denominator = 0.0;
                    for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    {
                        const Set::Scalar delta_velocity =
                            residual_velocity(d) -
                            previous_residual_velocity[b](d);
                        const Set::Scalar delta_rotation = residual_radius *
                            (residual_rotation(d) -
                            previous_residual_rotation[b](d));
                        numerator += previous_residual_velocity[b](d) *
                            delta_velocity + residual_radius *
                            previous_residual_rotation[b](d) *
                            delta_rotation;
                        denominator += delta_velocity * delta_velocity +
                            delta_rotation * delta_rotation;
                    }
                    if (denominator > 1.0e-24)
                    {
                        const Set::Scalar weight = Util::Clamp(
                            -numerator / denominator, 0.0,
                            maximum_anderson_weight[b]);
                        target_body = previous_coupling_output[b];
                        target_body.velocity += weight *
                            (body.velocity -
                            previous_coupling_output[b].velocity);
                        target_body.angular_velocity += weight *
                            (body.angular_velocity -
                            previous_coupling_output[b].angular_velocity);
                    }
                    else
                        target_body = body;
                }
                else
                    target_body = body;
            }

            previous_coupling_output[b] = body;
            previous_residual_velocity[b] = residual_velocity;
            previous_residual_rotation[b] = residual_rotation;
            have_previous_coupling_residual[b] = true;
        }
        if (all_bodies_finished) break;
    }

    // Only the final coupling candidate contributes a pressure correction.
    for (int lev = 0; lev < nlev; ++lev)
    {
        amrex::MultiFab::Copy(*pressure_correction_mf[lev],
            pressure_poisson.Solution(lev), 0, 0, 1, 0);
        amrex::MultiFab& p_mf = *pressure_mf[lev];
        const Set::Scalar p_floor = pressure_floor;
        const Set::Scalar p_scale_inv = 1.0 / pressure_scale;
        const Set::Scalar p_reference = pressure_reference;
        const bool update_pressure = projection_update_pressure;
        for (amrex::MFIter mfi(p_mf, amrex::TilingIfNotGPU());
            mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            Set::Patch<Set::Scalar> p = p_mf.array(mfi);
            Set::Patch<const Set::Scalar> phi =
                pressure_poisson.Solution(lev).array(mfi);
            amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                if (update_pressure)
                    p(i,j,k) = Util::Max(
                        p_reference + p_scale_inv * phi(i,j,k), p_floor);
            });
        }
        pressure_bc->define(geom[lev]);
        pressure_bc->FillBoundary(p_mf, 0, 1, time, 0);
        p_mf.FillBoundary(geom[lev].periodicity());
    }

    if (free_rigid_solid)
    {
        bool have_projection_history = false;
        for (int b = 0; b < nfree; ++b)
        {
            const FreeRigidBodyState& body = free_rigid_bodies[b];
            const FreeRigidBodyState& provisional =
                provisional_body_states[b];
            FreeRigidBodyState& previous_increment =
                previous_free_rigid_projection_increments[b];
            if (last_free_rigid_coupling_converged[b] &&
                body.valid && provisional.valid)
            {
                previous_increment = provisional;
                previous_increment.velocity =
                    body.velocity - provisional.velocity;
                previous_increment.angular_velocity =
                    body.angular_velocity - provisional.angular_velocity;
                have_projection_history = true;
            }
            else
                previous_increment = FreeRigidBodyState{};
            free_rigid_body_time[b] = body.valid ? time : NAN;
        }
        previous_free_rigid_dt = have_projection_history ? dt : NAN;
    }
}

void
LowMach::Initialize(int lev)
{
    BL_PROFILE("Integrator::LowMach::Initialize");
    const bool deformable_solid = deformable_solid_species >= 0;

    velocity_mf[lev]->setVal(0.0, 0, velocity_mf[lev]->nComp(), velocity_mf[lev]->nGrow());
    velocity_old_mf[lev]->setVal(0.0, 0, velocity_old_mf[lev]->nComp(), velocity_old_mf[lev]->nGrow());
    temperature_mf[lev]->setVal(0.0, 0, temperature_mf[lev]->nComp(), temperature_mf[lev]->nGrow());
    temperature_old_mf[lev]->setVal(0.0, 0, temperature_old_mf[lev]->nComp(), temperature_old_mf[lev]->nGrow());
    component_density_mf[lev]->setVal(0.0, 0, component_density_mf[lev]->nComp(), component_density_mf[lev]->nGrow());
    component_density_old_mf[lev]->setVal(0.0, 0, component_density_old_mf[lev]->nComp(), component_density_old_mf[lev]->nGrow());
    if (deformable_solid)
    {
        eta_mf[lev]->setVal(0.0, 0, eta_mf[lev]->nComp(), eta_mf[lev]->nGrow());
        xi_mf[lev]->setVal(0.0, 0, xi_mf[lev]->nComp(), xi_mf[lev]->nGrow());
        xi_old_mf[lev]->setVal(0.0, 0, xi_old_mf[lev]->nComp(), xi_old_mf[lev]->nGrow());
    }
    pressure_mf[lev]->setVal(0.0, 0, pressure_mf[lev]->nComp(), pressure_mf[lev]->nGrow());

    velocity_ic->Initialize(lev, velocity_mf, 0.0);
    velocity_ic->Initialize(lev, velocity_old_mf, 0.0);
    temperature_ic->Initialize(lev, temperature_mf, 0.0);
    temperature_ic->Initialize(lev, temperature_old_mf, 0.0);
    Set::Field<Set::Scalar> species_density(component_density_mf.size());
    for (int n = 0; n < nspecies; ++n)
    {
        if (component_density_ic[n] == nullptr) continue;
        species_density[lev] = std::make_unique<amrex::MultiFab>(
            *component_density_mf[lev], amrex::MakeType::make_alias, n, 1);
        component_density_ic[n]->Initialize(lev, species_density, 0.0);
    }
    pressure_ic->Initialize(lev, pressure_mf, 0.0);

    bool has_temperature_override = false;
    for (int n = ngas_species; n < nspecies; ++n)
        if (condensed_temperature_override[n] >= 0.0)
            has_temperature_override = true;
    if (has_temperature_override)
    {
        const auto inverse_reference_density =
            condensed_inverse_reference_density;
        const auto temperature_override = condensed_temperature_override;
        const int first_condensed = ngas_species;
        const int number_of_species = nspecies;
        const int number_of_gas_species = ngas_species;
        for (amrex::MFIter mfi(*temperature_mf[lev],
                                amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            Set::Patch<Set::Scalar> temperature =
                temperature_mf.Patch(lev,mfi);
            Set::Patch<Set::Scalar> component_density =
                component_density_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> pressure =
                pressure_mf.Patch(lev,mfi);

            amrex::ParallelFor(
                bx, [=,this] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar override_volume = 0.0;
                    Set::Scalar override_temperature = 0.0;
                    Set::Scalar condensed_volume = 0.0;
                    for (int n = first_condensed;
                         n < number_of_species; ++n)
                    {
                        const Set::Scalar volume = Util::Max(
                            component_density(i,j,k,n), 0.0) *
                            inverse_reference_density[n];
                        condensed_volume += volume;
                        if (temperature_override[n] >= 0.0)
                        {
                            override_volume += volume;
                            override_temperature +=
                                volume * temperature_override[n];
                        }
                    }
                    if (!(override_volume > 0.0)) return;

                    const Set::Scalar override_scale =
                        override_volume > 1.0 ?
                        1.0 / override_volume : 1.0;
                    const Set::Scalar weight = Util::Min(
                        override_volume, 1.0);
                    temperature(i,j,k) =
                        (1.0 - weight) * temperature(i,j,k) +
                        override_scale * override_temperature;

                    if (!(pressure(i,j,k) > 0.0) ||
                        !(temperature(i,j,k) > 0.0)) return;
                    Set::Scalar gas_molar_density = 0.0;
                    for (int n = 0; n < number_of_gas_species; ++n)
                        gas_molar_density += Util::Max(
                            component_density(i,j,k,n), 0.0) / gas.MW[n];
                    const Set::Scalar gas_volume = gas_molar_density *
                        Set::Constant::Rg * temperature(i,j,k) /
                        pressure(i,j,k);
                    const Set::Scalar available_volume =
                        Util::Max(1.0 - condensed_volume, 0.0);
                    if (gas_volume > 0.0)
                    {
                        const Set::Scalar scale =
                            available_volume / gas_volume;
                        for (int n = 0; n < number_of_gas_species; ++n)
                            component_density(i,j,k,n) *= scale;
                    }
                });
        }
    }
    if (deformable_solid)
    {
        xi_ic->Initialize(lev, xi_mf, 0.0);
        xi_ic->Initialize(lev, xi_old_mf, 0.0);
    }
    pressure_correction_mf[lev]->setVal(0.0);
    if (has_split_phase_change)
    {
        phase_change_dilatation_mf[lev]->setVal(0.0);
        phase_change_heat_mf[lev]->setVal(0.0);
    }
    if (lev == 0 && !(pressure_reference == pressure_reference))
        pressure_reference = pressure_mf[0]->sum(0, false) /
                            static_cast<Set::Scalar>(geom[0].Domain().numPts());
    if (lev == 0)
    {
        amrex::get<4>(thermal_data) = pressure_reference;
        amrex::get<0>(thermochemical_data) = thermal_data;
    }

    velocity_bc->define(geom[lev]);
    temperature_bc->define(geom[lev]);
    component_density_bc->define(geom[lev]);
    if (deformable_solid)
        xi_bc->define(geom[lev]);
    pressure_bc->define(geom[lev]);

    velocity_bc->FillBoundary(*velocity_mf[lev], 0, AMREX_SPACEDIM, 0.0, 0);
    velocity_mf[lev]->FillBoundary(geom[lev].periodicity());
    temperature_bc->FillBoundary(*temperature_mf[lev], 0, 1, 0.0, 0);
    temperature_mf[lev]->FillBoundary(geom[lev].periodicity());
    component_density_bc->FillBoundary(*component_density_mf[lev], 0, nspecies, 0.0, 0);
    component_density_mf[lev]->FillBoundary(geom[lev].periodicity());
    amrex::MultiFab::Copy(*component_density_old_mf[lev], *component_density_mf[lev],
                        0, 0, nspecies, component_density_mf[lev]->nGrow());
    amrex::MultiFab::Copy(*temperature_old_mf[lev], *temperature_mf[lev],
                        0, 0, 1, temperature_mf[lev]->nGrow());
    UpdateComponentState(lev, *component_density_mf[lev]);

    if (deformable_solid)
    {
        xi_bc->FillBoundary(*xi_mf[lev], 0, AMREX_SPACEDIM, 0.0, 0);
        xi_mf[lev]->FillBoundary(geom[lev].periodicity());
    }

    velocity_bc->FillBoundary(*velocity_old_mf[lev], 0, AMREX_SPACEDIM, 0.0, 0);
    velocity_old_mf[lev]->FillBoundary(geom[lev].periodicity());
    temperature_bc->FillBoundary(*temperature_old_mf[lev], 0, 1, 0.0, 0);
    temperature_old_mf[lev]->FillBoundary(geom[lev].periodicity());
    if (deformable_solid)
    {
        xi_bc->FillBoundary(*xi_old_mf[lev], 0, AMREX_SPACEDIM, 0.0, 0);
        xi_old_mf[lev]->FillBoundary(geom[lev].periodicity());
    }

    pressure_bc->FillBoundary(*pressure_mf[lev], 0, 1, 0.0, 0);
    pressure_mf[lev]->FillBoundary(geom[lev].periodicity());
    if (deformable_solid)
    {
        reference_map_reconstruction(
            geom[lev], *eta_mf[lev], *xi_mf[lev], *xi_bc, 0.0);
        reference_map_reconstruction(
            geom[lev], *eta_mf[lev], *xi_old_mf[lev], *xi_bc, 0.0);
        UpdateSolidStress(lev, *velocity_mf[lev], *eta_mf[lev], *xi_mf[lev], diagnostics_extended_fields);
    }
    UpdateDerivedDiagnostics(lev, *velocity_mf[lev], *temperature_mf[lev]);
}

void
LowMach::Regrid(int lev, Set::Scalar time)
{
    if (!reinitialize_condensed_composition || lev == 0 ||
        ngas_species >= nspecies)
        return;

    const int ncondensed = nspecies - ngas_species;
    amrex::MultiFab initial_density(
        component_density_mf[lev]->boxArray(),
        component_density_mf[lev]->DistributionMap(),
        ncondensed, 0);
    initial_density.setVal(0.0);

    Set::Field<Set::Scalar> species_density(component_density_mf.size());
    for (int n = ngas_species; n < nspecies; ++n)
    {
        if (component_density_ic[n] == nullptr) continue;
        species_density[lev] = std::make_unique<amrex::MultiFab>(
            initial_density, amrex::MakeType::make_alias,
            n - ngas_species, 1);
        component_density_ic[n]->Initialize(lev, species_density, 0.0);
    }

    const Set::Scalar* inverse_reference_density =
        amrex::get<8>(thermal_data);
    const int first_condensed = ngas_species;
    const int total_species = nspecies;
    const Set::Scalar eta_min =
        time == 0.0 ? 0.0 : reinitialize_condensed_composition_eta_min;
    for (amrex::MFIter mfi(*component_density_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        Set::Patch<Set::Scalar> density = component_density_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> density_ic = initial_density.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Scalar eta = 0.0;
            Set::Scalar eta_ic = 0.0;
            for (int n = first_condensed; n < total_species; ++n)
            {
                const Set::Scalar inverse_rho_ref =
                    inverse_reference_density[n];
                eta += Util::Max(density(i,j,k,n), 0.0) * inverse_rho_ref;
                eta_ic += Util::Max(
                    density_ic(i,j,k,n - first_condensed), 0.0) *
                    inverse_rho_ref;
            }

            if (eta > 0.0 && eta >= eta_min && eta_ic > 0.0)
            {
                // Restore the IC composition without moving the evolved
                // condensed/gas interface represented by the total eta.
                const Set::Scalar scale = eta / eta_ic;
                for (int n = first_condensed; n < total_species; ++n)
                    density(i,j,k,n) =
                        scale * density_ic(i,j,k,n - first_condensed);
            }
        });
    }

    component_density_bc->define(geom[lev]);
    component_density_bc->FillBoundary(
        *component_density_mf[lev], 0, nspecies, time, 0);
    component_density_mf[lev]->FillBoundary(geom[lev].periodicity());
    // Composition reconstruction is a remap, so the next RK stage must not
    // begin from the pre-remap buffer.
    amrex::MultiFab::Copy(*component_density_old_mf[lev],
                          *component_density_mf[lev], 0, 0, nspecies,
                          component_density_mf[lev]->nGrow());
    UpdateComponentState(lev, *component_density_mf[lev]);
}

void
LowMach::RHS(int lev, Set::Scalar time, Set::Scalar dt,
            amrex::MultiFab& u_rhs_mf,
            amrex::MultiFab& T_rhs_mf,
            amrex::MultiFab& component_density_rhs_mf,
            amrex::MultiFab* xi_rhs_mf,
            const amrex::MultiFab& u_mf,
            const amrex::MultiFab& T_mf,
            const amrex::MultiFab& component_density_mf,
            const amrex::MultiFab* xi_mf)
{
    const bool deformable_solid = deformable_solid_species >= 0;
    const bool split_chemistry = chemistry.Split();
    const bool external_heat_source = heat_source_ic != nullptr;
    if (deformable_solid && finite_solid_deviatoric_stress_divergence_sign != 0.0)
        UpdateSolidStress(lev, u_mf, *eta_mf[lev], *xi_mf);

    const auto dx = geom[lev].CellSizeArray();
    amrex::Box domain = geom[lev].Domain();
    amrex::GpuArray<int,AMREX_SPACEDIM> periodic{};
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
        periodic[d] = geom[lev].isPeriodic(d);
    const auto gas_data = gas_device_data;
    const auto thermal = thermal_data;
    const auto thermochemical = thermochemical_data;
    const auto advect_scheme = advect;
    const Set::Vector gravity = g;
    const Set::Scalar rho_floor = density_floor;
    const bool explicit_viscosity =
        include_viscosity && !implicit_momentum_diffusion;
    const bool transport_condensed = amrex::get<5>(thermal);
    const Set::Scalar* condensed_cp = amrex::get<6>(thermal);
    const Set::Scalar* condensed_inverse_density = amrex::get<8>(thermal);
    const Set::Scalar* liquid_inverse_density = amrex::get<9>(thermal);
    const Set::Scalar stress_sign =
        finite_solid_deviatoric_stress_divergence_sign;
    if (external_heat_source)
        heat_source_ic->Initialize(lev, heat_source_mf, time);

    // Each freely moving rigid species must be transported by the same rigid
    // translation and rotation used by its momentum penalty.  Transporting
    // one with the mixture velocity lets shear in the diffuse interface peel
    // off dilute solid filaments, which subsequently behave as independent
    // material.  Body velocities are held explicit over one flow step and
    // evaluated about their translated stage centers.
    const int nfree = static_cast<int>(free_rigid_solid_species.size());
    Model::Chemistry::SpeciesArray free_rigid_species_flag{};
    std::vector<bool> advect_free_rigid(nfree, false);
    std::vector<RigidBodyVelocity> free_rigid_velocity(nfree);
    for (int b = 0; b < nfree; ++b)
    {
        const int n = free_rigid_solid_species[b];
        free_rigid_species_flag[n] = 1.0;
        advect_free_rigid[b] = free_rigid_bodies[b].valid &&
            std::isfinite(free_rigid_body_time[b]);
        free_rigid_velocity[b].translation = free_rigid_bodies[b].velocity;
        free_rigid_velocity[b].rotation =
            free_rigid_bodies[b].angular_velocity;
        free_rigid_velocity[b].center = free_rigid_bodies[b].center +
            (time - free_rigid_body_time[b]) *
            free_rigid_velocity[b].translation;
        free_rigid_velocity[b].prob_lo = geom[lev].ProbLoArray();
        free_rigid_velocity[b].dx = geom[lev].CellSizeArray();
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
        {
            free_rigid_velocity[b].period[d] = geom[lev].ProbLength(d);
            free_rigid_velocity[b].periodic[d] =
                geom[lev].isPeriodic(d);
        }
    }

    u_rhs_mf.setVal(0.0, 0, AMREX_SPACEDIM, u_rhs_mf.nGrow());

    //
    // Momentum equation
    //
    for (amrex::MFIter mfi(u_rhs_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        Set::Patch<const Set::Scalar> u = u_mf.array(mfi);
        Set::Patch<const Set::Scalar> T = T_mf.array(mfi);
        Set::Patch<const Set::Scalar> component_density = component_density_mf.array(mfi);
        Set::Patch<const Set::Scalar> rho = density_mf.Patch(lev,mfi);
        Set::Patch<const Set::Matrix> solid_deviatoric_stress =
            solid_deviatoric_stress_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> u_rhs = u_rhs_mf.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            auto sten = Numeric::GetStencil(i, j, k, domain, periodic);
            Set::Scalar density = Util::Max(rho(i,j,k), rho_floor);
            Set::Scalar mu = explicit_viscosity ?
                ComputeViscosity(
                    component_density, T(i,j,k), i, j, k, thermal) : 0.0;

            Set::Vector div_sigma = Set::Vector::Zero();
            if (deformable_solid)
                div_sigma = Numeric::Divergence(
                    solid_deviatoric_stress, i, j, k, dx.data(), sten);

            Set::Vector adv_u = advect_scheme.Vector(
                u, u, i, j, k, 0, dx.data(),
                {Numeric::Advect::Form::Advective}, sten);
            Set::Vector lap_u =
                Numeric::VectorLaplacian(u, i, j, k, 0, dx.data());

            // Pressure is applied once by the nonincremental projection after
            // this predictor; including the stored pressure here feeds the
            // projection solution back into the next time step.
            Set::Vector rhs_vec = adv_u
                                + gravity
                                + (mu / density) * lap_u
                                + (stress_sign / density) * div_sigma;

            // Put the result into the multicomponent field
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
                u_rhs(i,j,k,d) = rhs_vec(d);
        });
    }

    //
    // Mechanism equations (phase field)
    //
    T_rhs_mf.setVal(0.0);
    component_density_rhs_mf.setVal(0.0);
    for (const auto& mechanism : mechanisms)
    {
        const auto mechanism_view = mechanism;
        const Set::Scalar p_reference = pressure_reference;
        for (amrex::MFIter mfi(component_density_rhs_mf, false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            Set::Patch<const Set::Scalar> T = T_mf.array(mfi);
            Set::Patch<const Set::Scalar> component_density = component_density_mf.array(mfi);
            Set::Patch<const Set::Scalar> rigid_eta = rigid_eta_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> rigid_species_eta =
                rigid_species_eta_mf.Patch(lev,mfi);
            Set::Patch<Set::Scalar> component_density_rhs = component_density_rhs_mf.array(mfi);
            Set::Patch<Set::Scalar> T_rhs = T_rhs_mf.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                const Model::Mechanism::State state = {
                    component_density, rigid_eta, rigid_species_eta,
                    T(i,j,k), p_reference, dt};
                mechanism_view.Apply(
                    component_density_rhs, state, i, j, k, dx.data());
                auto [gas_volume_fraction, gas_heat_capacity, heat_capacity,
                    conductivity, cp] = ComputeThermalState(
                        component_density, T(i,j,k), i, j, k, thermal);
                (void)gas_volume_fraction;
                (void)gas_heat_capacity;
                (void)conductivity;
                (void)cp;
                if (heat_capacity > 0.0)
                    T_rhs(i,j,k) += mechanism_view.HeatSource(
                        state, i, j, k, dx.data()) / heat_capacity;
            });
        }
    }

    //
    // Advection equations for other fields including:
    // - temperature
    // - component densities
    // - reference map (xi)
    //
    for (amrex::MFIter mfi(T_mf, false); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        Set::Patch<const Set::Scalar> u = u_mf.array(mfi);
        Set::Patch<const Set::Scalar> T = T_mf.array(mfi);
        Set::Patch<const Set::Scalar> component_density = component_density_mf.array(mfi);
        Set::Patch<const Set::Scalar> heat_source;
        if (external_heat_source)
            heat_source = heat_source_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> xi;
        if (deformable_solid)
            xi = xi_mf->array(mfi);
        Set::Patch<Set::Scalar> T_rhs = T_rhs_mf.array(mfi);
        Set::Patch<Set::Scalar> component_density_rhs = component_density_rhs_mf.array(mfi);
        Set::Patch<Set::Scalar> xi_rhs;
        if (deformable_solid)
            xi_rhs = xi_rhs_mf->array(mfi);
        const int ngas = ngas_species;
        const int component_count = nspecies;
        const int deformable_component = deformable_solid_species;
        const bool advect_T = advect_temperature;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            auto sten = Numeric::GetStencil(i, j, k, domain, periodic);

            auto [species,temperature,dilatation] =
                ComputeThermochemicalSource(
                    component_density, T, i, j, k, dx.data(), dt,
                    thermochemical, !split_chemistry);
            (void)dilatation; // ignore unused

            // Species are finite-volume conserved variables.  Retain their
            // transport contribution separately so temperature advection can
            // enforce the matching conservative gas-enthalpy balance without
            // folding chemistry or phase-change sources into that balance.
            Model::Chemistry::SpeciesArray gas_advection{};
            for (int n = 0; n < ngas; ++n)
                gas_advection[n] = advect_scheme(
                    component_density, u, i, j, k, n, dx.data(),
                    {Numeric::Advect::Form::Conservative}, sten);

            const Set::Scalar mechanism_temperature_source = T_rhs(i,j,k);
            T_rhs(i,j,k) = mechanism_temperature_source;
            if (advect_T)
            {
                auto [gas_volume_fraction, gas_heat_capacity, heat_capacity,
                    conductivity, cp] = ComputeThermalState(
                        component_density, T(i,j,k), i, j, k, thermal);
                (void)gas_volume_fraction;
                (void)gas_heat_capacity;
                (void)conductivity;
                (void)cp;
                if (heat_capacity > 0.0)
                {
                    // Advect volumetric gas enthalpy with the same
                    // finite-volume operator used for partial densities.  A
                    // direct advective update of T is not energy consistent
                    // when a hot, light gas and a cold, dense gas mix
                    // numerically: it combines a volume-weighted temperature
                    // with mass-conservative density and creates sensible
                    // enthalpy.  The chain-rule subtraction below converts
                    // d(rho*h)/dt back to dT/dt at fixed composition.
                    auto gas_enthalpy = [=] AMREX_GPU_HOST_DEVICE(
                        int ii, int jj, int kk, int /*comp*/)
                    {
                        Set::Scalar value = 0.0;
                        for (int n = 0; n < ngas; ++n)
                            value += component_density(ii,jj,kk,n) *
                                Model::Gas::Gas::EnthalpyMassSpecies(
                                    gas_data, T(ii,jj,kk), n);
                        return value;
                    };
                    Set::Scalar gas_enthalpy_rhs = advect_scheme(
                        gas_enthalpy, u, i, j, k, 0, dx.data(),
                        {Numeric::Advect::Form::Conservative}, sten);
                    for (int n = 0; n < ngas; ++n)
                        gas_enthalpy_rhs -=
                            Model::Gas::Gas::EnthalpyMassSpecies(
                                gas_data, T(i,j,k), n) * gas_advection[n];
                    T_rhs(i,j,k) += gas_enthalpy_rhs / heat_capacity;

                    Set::Scalar moving_condensed_heat_capacity = 0.0;
                    if (transport_condensed)
                    {
                        Set::Scalar condensed_volume = 0.0;
                        for (int n = ngas; n < component_count; ++n)
                            condensed_volume += Util::Max(
                                component_density(i,j,k,n), 0.0) *
                                condensed_inverse_density[n];
                        const Set::Scalar solid_fraction =
                            Model::PhaseField::H(condensed_volume);
                        const Set::Scalar solid_scale =
                            condensed_volume > 0.0 ?
                            solid_fraction / condensed_volume : 0.0;
                        for (int n = ngas; n < component_count; ++n)
                        {
                            const bool moving_condensed_phase =
                                n == deformable_component ||
                                free_rigid_species_flag[n] > 0.5 ||
                                liquid_inverse_density[n] > 0.0;
                            if (moving_condensed_phase)
                                moving_condensed_heat_capacity += solid_scale *
                                    Util::Max(
                                        component_density(i,j,k,n), 0.0) *
                                        condensed_cp[n];
                        }
                    }
                    T_rhs(i,j,k) += moving_condensed_heat_capacity /
                        heat_capacity * advect_scheme(
                            T, u, i, j, k, 0, dx.data(),
                            {Numeric::Advect::Form::Advective}, sten);
                }
            }
            T_rhs(i,j,k) += temperature;
            if (external_heat_source)
            {
                auto [gas_volume_fraction, gas_heat_capacity, heat_capacity,
                    conductivity, cp] = ComputeThermalState(
                        component_density, T(i,j,k), i, j, k, thermal);
                (void)gas_volume_fraction;
                (void)gas_heat_capacity;
                (void)conductivity;
                (void)cp;
                if (heat_capacity > 0.0)
                    T_rhs(i,j,k) += heat_source(i,j,k) / heat_capacity;
            }

            for (int n = 0; n < ngas; ++n)
            {
                const Set::Scalar mechanism_source = component_density_rhs(i,j,k,n);
                component_density_rhs(i,j,k,n) =
                    gas_advection[n] + mechanism_source + species[n];
            }
            for (int n = ngas; n < component_count; ++n)
            {
                // Prescribed/fixed rigid phases keep their existing spatial
                // profiles.  Their mechanism source was assembled above;
                // only deformable and freely moving phases are transported.
                const bool liquid_component =
                    liquid_inverse_density[n] > 0.0;
                const bool free_rigid_component =
                    free_rigid_species_flag[n] > 0.5;
                if (free_rigid_component)
                    continue;
                if (n != deformable_component && !liquid_component)
                    continue;
                const Set::Scalar mechanism_source = component_density_rhs(i,j,k,n);
                component_density_rhs(i,j,k,n) =
                    advect_scheme(component_density, u, i, j, k, n,
                        dx.data(),
                        {Numeric::Advect::Form::Conservative}, sten) +
                    mechanism_source;
            }
            if (deformable_solid)
                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    xi_rhs(i,j,k,d) = advect_scheme(
                        xi, u, i, j, k, d, dx.data(),
                        {Numeric::Advect::Form::Advective}, sten);
        });
        for (int b = 0; b < nfree; ++b)
        {
            const int n = free_rigid_solid_species[b];
            const bool use_rigid_velocity = advect_free_rigid[b];
            const RigidBodyVelocity body_velocity = free_rigid_velocity[b];
            amrex::ParallelFor(bx,
                [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                auto sten = Numeric::GetStencil(
                    i, j, k, domain, periodic);
                const Set::Scalar mechanism_source =
                    component_density_rhs(i,j,k,n);
                const Set::Scalar advection = use_rigid_velocity ?
                    advect_scheme(component_density, body_velocity,
                        i, j, k, n, dx.data(),
                        {Numeric::Advect::Form::Conservative}, sten) :
                    advect_scheme(component_density, u,
                        i, j, k, n, dx.data(),
                        {Numeric::Advect::Form::Conservative}, sten);
                component_density_rhs(i,j,k,n) =
                    advection + mechanism_source;
            });
        }
    }

    u_rhs_mf.FillBoundary(geom[lev].periodicity());
}

void
LowMach::Advance(int lev, Set::Scalar time, Set::Scalar dt)
{
    const bool deformable_solid = deformable_solid_species >= 0;
    const bool split_chemistry = chemistry.Split();
    if (split_chemistry && time == t_new[0])
        chemistry_dilatation_mf[lev]->setVal(0.0);
    std::swap(velocity_old_mf[lev], velocity_mf[lev]);
    std::swap(temperature_old_mf[lev], temperature_mf[lev]);
    std::swap(component_density_old_mf[lev], component_density_mf[lev]);
    if (deformable_solid)
        std::swap(xi_old_mf[lev], xi_mf[lev]);

    if (split_chemistry)
    {
        AdvanceChemistry(lev, *temperature_old_mf[lev],
                            *component_density_old_mf[lev], 0.5 * dt);
        temperature_bc->FillBoundary(*temperature_old_mf[lev], 0, 1,
                                    time + 0.5 * dt, 0);
        temperature_old_mf[lev]->FillBoundary(geom[lev].periodicity());
        component_density_bc->FillBoundary(*component_density_old_mf[lev],
                                            0, nspecies, time + 0.5 * dt, 0);
        component_density_old_mf[lev]->FillBoundary(geom[lev].periodicity());
    }

    amrex::Vector<amrex::MultiFab> solution_new;
    solution_new.emplace_back(*velocity_mf[lev].get(), amrex::MakeType::make_alias, 0, AMREX_SPACEDIM);
    solution_new.emplace_back(*temperature_mf[lev].get(), amrex::MakeType::make_alias, 0, 1);
    solution_new.emplace_back(*component_density_mf[lev].get(), amrex::MakeType::make_alias, 0, nspecies);
    if (deformable_solid)
        solution_new.emplace_back(*xi_mf[lev].get(), amrex::MakeType::make_alias, 0, AMREX_SPACEDIM);

    amrex::Vector<amrex::MultiFab> solution_old;
    solution_old.emplace_back(*velocity_old_mf[lev].get(), amrex::MakeType::make_alias, 0, AMREX_SPACEDIM);
    solution_old.emplace_back(*temperature_old_mf[lev].get(), amrex::MakeType::make_alias, 0, 1);
    solution_old.emplace_back(*component_density_old_mf[lev].get(), amrex::MakeType::make_alias, 0, nspecies);
    if (deformable_solid)
        solution_old.emplace_back(*xi_old_mf[lev].get(), amrex::MakeType::make_alias, 0, AMREX_SPACEDIM);

    velocity_bc->define(geom[lev]);
    temperature_bc->define(geom[lev]);
    component_density_bc->define(geom[lev]);
    if (deformable_solid)
        xi_bc->define(geom[lev]);

    amrex::TimeIntegrator timeintegrator(solution_new, time);
    timeintegrator.set_rhs([&](amrex::Vector<amrex::MultiFab>& rhs_mf,
                                amrex::Vector<amrex::MultiFab>& state_mf,
                                const Set::Scalar rhs_time)
    {
        velocity_bc->FillBoundary(state_mf[0], 0, AMREX_SPACEDIM, rhs_time, 0);
        state_mf[0].FillBoundary(geom[lev].periodicity());
        temperature_bc->FillBoundary(state_mf[1], 0, 1, rhs_time, 0);
        state_mf[1].FillBoundary(geom[lev].periodicity());
        component_density_bc->FillBoundary(state_mf[2], 0, nspecies, rhs_time, 0);
        state_mf[2].FillBoundary(geom[lev].periodicity());
        if (deformable_solid)
        {
            xi_bc->FillBoundary(state_mf[3], 0, AMREX_SPACEDIM, rhs_time, 0);
            state_mf[3].FillBoundary(geom[lev].periodicity());
        }

        UpdateComponentState(lev, state_mf[2]);
        if (deformable_solid)
            reference_map_reconstruction(
                geom[lev], *eta_mf[lev], state_mf[3], *xi_bc, rhs_time);
        RHS(lev, rhs_time, dt, rhs_mf[0], rhs_mf[1], rhs_mf[2],
            deformable_solid ? &rhs_mf[3] : nullptr,
            state_mf[0], state_mf[1], state_mf[2],
            deformable_solid ? &state_mf[3] : nullptr);
    });

    timeintegrator.set_post_stage_action([&](amrex::Vector<amrex::MultiFab>& stage_mf, Set::Scalar stage_time)
    {
        velocity_bc->FillBoundary(stage_mf[0], 0, AMREX_SPACEDIM, stage_time, 0);
        stage_mf[0].FillBoundary(geom[lev].periodicity());
        temperature_bc->FillBoundary(stage_mf[1], 0, 1, stage_time, 0);
        stage_mf[1].FillBoundary(geom[lev].periodicity());
        component_density_bc->FillBoundary(stage_mf[2], 0, nspecies, stage_time, 0);
        stage_mf[2].FillBoundary(geom[lev].periodicity());
        if (deformable_solid)
        {
            xi_bc->FillBoundary(stage_mf[3], 0, AMREX_SPACEDIM, stage_time, 0);
            stage_mf[3].FillBoundary(geom[lev].periodicity());
        }

        UpdateComponentState(lev, stage_mf[2]);
        if (deformable_solid)
            reference_map_reconstruction(
                geom[lev], *eta_mf[lev], stage_mf[3], *xi_bc, stage_time);
    });

    timeintegrator.advance(solution_old, solution_new, time, dt);
    if (split_chemistry)
    {
        AdvanceChemistry(lev, *temperature_mf[lev],
                            *component_density_mf[lev], 0.5 * dt);
    }
    Set::Scalar new_time = time + dt;
    velocity_bc->FillBoundary(*velocity_mf[lev], 0, AMREX_SPACEDIM, new_time, 0);
    velocity_mf[lev]->FillBoundary(geom[lev].periodicity());
    temperature_bc->FillBoundary(*temperature_mf[lev], 0, 1, new_time, 0);
    temperature_mf[lev]->FillBoundary(geom[lev].periodicity());
    component_density_bc->FillBoundary(*component_density_mf[lev], 0, nspecies, new_time, 0);
    component_density_mf[lev]->FillBoundary(geom[lev].periodicity());
    if (deformable_solid)
    {
        xi_bc->FillBoundary(*xi_mf[lev], 0, AMREX_SPACEDIM, new_time, 0);
        xi_mf[lev]->FillBoundary(geom[lev].periodicity());
    }

    UpdateComponentState(lev, *component_density_mf[lev]);
    if (deformable_solid)
    {
        reference_map_reconstruction(
            geom[lev], *eta_mf[lev], *xi_mf[lev], *xi_bc, new_time);
        xi_bc->FillBoundary(*xi_mf[lev], 0, AMREX_SPACEDIM, new_time, 0);
        xi_mf[lev]->FillBoundary(geom[lev].periodicity());
        reference_map_reconstruction(
            geom[lev], *eta_mf[lev], *xi_mf[lev], *xi_bc, new_time);
    }
}

//
// Calculate dynamic timestep
//
void
LowMach::TimeStepBegin(Set::Scalar time, int /*iter*/)
{
    if (synchronize_restart_state)
    {
        if (!(pressure_reference == pressure_reference))
        {
            Set::Field<Set::Scalar> reference_pressure(pressure_mf.size());
            reference_pressure[0] = std::make_unique<amrex::MultiFab>(
                pressure_mf[0]->boxArray(), pressure_mf[0]->DistributionMap(),
                1, 0);
            pressure_ic->Initialize(0, reference_pressure, 0.0);
            pressure_reference = reference_pressure[0]->sum(0, false) /
                static_cast<Set::Scalar>(geom[0].Domain().numPts());
        }
        amrex::get<4>(thermal_data) = pressure_reference;
        amrex::get<0>(thermochemical_data) = thermal_data;
        for (int lev = 0; lev <= finest_level; ++lev)
        {
            amrex::MultiFab::Copy(*velocity_old_mf[lev], *velocity_mf[lev],
                                  0, 0, AMREX_SPACEDIM, velocity_mf[lev]->nGrow());
            amrex::MultiFab::Copy(*temperature_old_mf[lev], *temperature_mf[lev],
                                  0, 0, 1, temperature_mf[lev]->nGrow());
            amrex::MultiFab::Copy(*component_density_old_mf[lev],
                                  *component_density_mf[lev], 0, 0, nspecies,
                                  component_density_mf[lev]->nGrow());
            if (deformable_solid_species >= 0)
                amrex::MultiFab::Copy(*xi_old_mf[lev], *xi_mf[lev],
                                      0, 0, AMREX_SPACEDIM, xi_mf[lev]->nGrow());
        }
        synchronize_restart_state = false;
    }

    // Initialize the body kinematics before the first Runge-Kutta stage (and
    // refresh them after a restart).  Later steps already carry the state
    // produced by the preceding projection.
    bool initialize_free_rigid_bodies = false;
    for (int b = 0;
         b < static_cast<int>(free_rigid_solid_species.size()); ++b)
        initialize_free_rigid_bodies = initialize_free_rigid_bodies ||
            !free_rigid_bodies[b].valid ||
            !std::isfinite(free_rigid_body_time[b]);
    if (initialize_free_rigid_bodies)
    {
        for (int lev = 0; lev <= finest_level; ++lev)
            UpdateComponentState(lev, *component_density_mf[lev]);
        UpdateFreeRigidBodyStates();
        for (int b = 0;
             b < static_cast<int>(free_rigid_solid_species.size()); ++b)
            if (free_rigid_bodies[b].valid)
                free_rigid_body_time[b] = time;
    }
    if (!dynamictimestep.on) return;

    const bool deformable_solid = deformable_solid_species >= 0;
    const bool explicit_solid_deviatoric_stress = deformable_solid &&
        finite_solid_deviatoric_stress_divergence_sign != 0.0;
    const auto gas_data = gas_device_data;
    const auto thermal = thermal_data;
    const auto chemistry_data = chemistry_device_data;
    const bool report_chemistry_timescale =
        chemistry_timestep_mode != ChemistryTimestepMode::Off;
    const bool limit_chemistry_timestep =
        chemistry_timestep_mode == ChemistryTimestepMode::Limit;
    chemistry_timescale_min = report_chemistry_timescale ? cfl_v : NAN;
    chemistry_timestep_candidate = report_chemistry_timescale ? cfl_v : NAN;
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        Set::Scalar advmax = 0.0;
        Set::Scalar viscmax = 0.0;
        Set::Scalar elasticmax = 0.0;
        Set::Scalar phasefieldmax = 0.0;
        Set::Scalar capillarymax = 0.0;
        Set::Scalar chemistryratemax = 0.0;
        UpdateComponentState(lev, *component_density_mf[lev]);
        const auto dx = geom[lev].CellSizeArray();
        Set::Scalar dxmin = std::min(dx[0], dx[1]);
        if (interfacial_forces_enabled)
        {
            capillarymax = std::sqrt(
                interfacial_model.MaximumCapillaryEnergy() /
                (interfacial_reference_density_min * dxmin * dxmin * dxmin));
        }
        Set::Scalar temperaturemax = 0.0;

        for (amrex::MFIter mfi(*velocity_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            Set::Patch<const Set::Scalar> u = velocity_mf.Patch(lev,mfi);
            amrex::ReduceOps<amrex::ReduceOpMax> reduce_op;
            amrex::ReduceData<Set::Scalar> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;
            reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple
            {
                // The cell-centered MUSCL operator applies all directional
                // fluxes in one unsplit update.  Its multidimensional CFL
                // rate is therefore the sum of the directional Courant
                // rates, not |u| divided by the smallest cell width.  The
                // latter can underpredict the rate by sqrt(dim) for diagonal
                // flow and is especially visible when a curved interface
                // crosses a periodic seam.
                Set::Scalar advection_rate = 0.0;
                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    advection_rate += Util::Abs(u(i,j,k,d)) / dx[d];
                return {advection_rate};
            });
            ReduceTuple hv = reduce_data.value();
            advmax = std::max(advmax, amrex::get<0>(hv));
        }

        for (amrex::MFIter mfi(*temperature_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            Set::Patch<const Set::Scalar> T = temperature_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> component_density = component_density_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> rho = density_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> eta = eta_mf.Patch(lev,mfi);
            const Set::Scalar rho_floor = density_floor;
            const bool viscous =
                include_viscosity && !implicit_momentum_diffusion;
            const bool conductive =
                include_conduction && !implicit_thermal_diffusion;
            const bool species_diffusive =
                ngas_species > 1 && !implicit_species_diffusion;
            const bool elastic = explicit_solid_deviatoric_stress;
            const int ngas = ngas_species;
            const Set::Scalar p_reference = pressure_reference;
            const Set::Scalar eta_threshold = Util::Clamp(finite_solid_eta_threshold, 0.0, 1.0);
            const Set::Scalar mu_solid = finite_solid_model.mu;
            const Set::Scalar kappa_solid = finite_solid_model.kappa;
            const Set::Scalar solid_viscosity = finite_solid_viscosity;
            const Set::Scalar solid_interface_viscosity = finite_solid_interface_viscosity;
            const bool evaluate_chemistry = report_chemistry_timescale;
            const Set::Scalar reactant_fraction_floor =
                chemistry_reactant_mass_fraction_floor;

            amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpMax,
                            amrex::ReduceOpMax, amrex::ReduceOpMax> reduce_op;
            amrex::ReduceData<Set::Scalar, Set::Scalar, Set::Scalar,
                            Set::Scalar> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;
            reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple
            {
                Set::Scalar nu = 0.0;
                if (viscous)
                {
                    Set::Scalar mu = ComputeViscosity(
                        component_density, T(i,j,k), i, j, k, thermal);
                    nu = mu / Util::Max(rho(i,j,k), rho_floor);
                }
                if (conductive)
                {
                    const Set::Scalar cp = Model::Gas::Gas::CpMass(
                        gas_data, T(i,j,k), component_density, i, j, k);
                    const Set::Scalar kappa =
                        Model::Gas::Gas::ThermalConductivity(
                            gas_data, T(i,j,k),
                            component_density, i, j, k);
                    nu = Util::Max(nu, kappa /
                        (Util::Max(rho(i,j,k), rho_floor) * cp));
                }
                if (species_diffusive)
                    for (int n = 0; n < ngas; ++n)
                        nu = Util::Max(
                            nu, Model::Gas::Gas::DiffusionCoefficient(
                                gas_data, T(i,j,k), p_reference,
                                component_density, i, j, k, n));
                Set::Scalar elastic_rate = 0.0;
                if (elastic && Util::Clamp(eta(i,j,k), 0.0, 1.0) > eta_threshold)
                {
                    Set::Scalar density = Util::Max(rho(i,j,k), rho_floor);
                    Set::Scalar wave_speed = std::sqrt(Util::Max(kappa_solid + (4.0 / 3.0) * mu_solid, mu_solid) / density);
                    elastic_rate = wave_speed / dxmin;
                    nu = Util::Max(nu, (solid_viscosity + solid_interface_viscosity) / density);
                }

                // This is an accuracy timescale for the locally implicit,
                // split chemistry solve, not an explicit stability bound.
                // It limits a forward estimate of fractional temperature
                // rise or major-reactant consumption.  The mass-fraction
                // floor prevents trace reactants from setting the global
                // timestep while retaining the temperature signal from a
                // strongly exothermic trace reaction.
                Set::Scalar chemistry_rate = 0.0;
                if (evaluate_chemistry && T(i,j,k) > 0.0 &&
                    p_reference > 0.0)
                {
                    Model::Chemistry::SpeciesArray rhoY{};
                    Set::Scalar gas_density = 0.0;
                    for (int n = 0; n < ngas; ++n)
                    {
                        rhoY[n] = Util::Max(
                            component_density(i,j,k,n), 0.0);
                        gas_density += rhoY[n];
                    }

                    Set::Scalar raw_gas_volume_fraction = 0.0;
                    if (gas_density > rho_floor)
                        raw_gas_volume_fraction = Util::Max(
                            gas_density * Model::Gas::Gas::GasConstant(
                                gas_data, component_density, i, j, k) *
                                T(i,j,k) / p_reference,
                            0.0);
                    Set::Scalar condensed_volume_fraction = 0.0;
                    const Set::Scalar* inverse_reference_density =
                        amrex::get<8>(thermal);
                    for (int n = ngas; n < amrex::get<1>(thermal); ++n)
                        condensed_volume_fraction += Util::Max(
                            component_density(i,j,k,n), 0.0) *
                            inverse_reference_density[n];
                    const Set::Scalar reacting_volume_fraction = Util::Min(
                        raw_gas_volume_fraction,
                        1.0 - Model::PhaseField::H(
                            condensed_volume_fraction));

                    if (raw_gas_volume_fraction > 0.0 &&
                        reacting_volume_fraction > 0.0)
                    {
                        Model::Chemistry::SpeciesArray intrinsic_rhoY{};
                        for (int n = 0; n < ngas; ++n)
                            intrinsic_rhoY[n] =
                                rhoY[n] / raw_gas_volume_fraction;

                        Model::Chemistry::Source reaction{};
                        const auto& [gross_model, model, solver,
                                     host_chemistry, host_gas] =
                            chemistry_data;
                        (void)solver;
#ifdef ALAMO_GPU
                        if (gross_model)
                            reaction = model.ComputeChemistrySources(
                                p_reference, T(i,j,k), intrinsic_rhoY,
                                0.0, nullptr);
#else
                        reaction = host_chemistry->ComputeChemistrySources(
                            p_reference, T(i,j,k), intrinsic_rhoY,
                            0.0, host_gas);
#endif

                        auto [gas_volume_fraction, gas_heat_capacity,
                            heat_capacity, conductivity, cp] =
                            ComputeThermalState(
                                component_density, T(i,j,k), i, j, k,
                                thermal);
                        (void)gas_volume_fraction;
                        (void)gas_heat_capacity;
                        (void)conductivity;
                        (void)cp;
                        if (heat_capacity > 0.0)
                            chemistry_rate = Util::Abs(
                                reacting_volume_fraction * reaction.second /
                                heat_capacity) / T(i,j,k);

                        const Set::Scalar density_scale =
                            reactant_fraction_floor * gas_density;
                        for (int n = 0; n < ngas; ++n)
                            if (reaction.first[n] < 0.0)
                                chemistry_rate = Util::Max(
                                    chemistry_rate,
                                    reacting_volume_fraction *
                                        (-reaction.first[n]) /
                                        Util::Max(rhoY[n], density_scale));
                    }
                }
                return {nu / (dxmin * dxmin), elastic_rate, T(i,j,k),
                        chemistry_rate};
            });
            ReduceTuple hv = reduce_data.value();
            viscmax = std::max(viscmax, amrex::get<0>(hv));
            elasticmax = std::max(elasticmax, amrex::get<1>(hv));
            temperaturemax = std::max(temperaturemax, amrex::get<2>(hv));
            chemistryratemax = std::max(
                chemistryratemax, amrex::get<3>(hv));
        }
        amrex::ParallelDescriptor::ReduceRealMax(advmax);
        amrex::ParallelDescriptor::ReduceRealMax(viscmax);
        amrex::ParallelDescriptor::ReduceRealMax(elasticmax);
        amrex::ParallelDescriptor::ReduceRealMax(temperaturemax);
        amrex::ParallelDescriptor::ReduceRealMax(chemistryratemax);
        for (const auto& mechanism : mechanisms)
            if (mechanism.RigidComponent() < 0)
                phasefieldmax = std::max(
                    phasefieldmax,
                    mechanism.StabilityRate(dxmin, temperaturemax));

        Set::Scalar adv_dt = advmax > 0.0 ? cfl / advmax : cfl_v;
        Set::Scalar visc_dt =
            viscmax > 0.0 ? 0.5 * cfl / viscmax : cfl_v;
        Set::Scalar elastic_dt =
            elasticmax > 0.0 ? cfl / elasticmax : cfl_v;
        Set::Scalar phasefield_dt = phasefieldmax > 0.0 ?
            phase_field_cfl / phasefieldmax : cfl_v;
        Set::Scalar capillary_dt =
            capillarymax > 0.0 ? cfl / capillarymax : cfl_v;
        Set::Scalar chemistry_timescale = chemistryratemax > 0.0 ?
            1.0 / chemistryratemax : cfl_v;
        Set::Scalar chemistry_dt = report_chemistry_timescale ?
            chemistry_max_fractional_change * chemistry_timescale : cfl_v;
        chemistry_timescale_min = std::min(
            chemistry_timescale_min, chemistry_timescale);
        chemistry_timestep_candidate = std::min(
            chemistry_timestep_candidate, chemistry_dt);
        if (dynamictimestep.verbose &&
            amrex::ParallelDescriptor::IOProcessor())
            amrex::Print() << "LowMach timestep level " << lev
                << " advective " << adv_dt
                << " viscous " << visc_dt
                << " elastic " << elastic_dt
                << " phase_field " << phasefield_dt
                << " capillary " << capillary_dt
                << " chemistry_timescale " << chemistry_timescale
                << " chemistry_candidate " << chemistry_dt
                << " chemistry_limits " << limit_chemistry_timestep
                << "\n";
        DynamicTimestep_SyncTimeStep(
            lev, std::min({adv_dt, visc_dt, elastic_dt, phasefield_dt,
                           capillary_dt,
                           limit_chemistry_timestep ? chemistry_dt : cfl_v}));
    }
    DynamicTimestep_Update();
}


void
LowMach::PrintDiagnostics(Set::Scalar time, int iter)
{
    if (diagnostics_interval <= 0 || iter % diagnostics_interval != 0) return;

    Set::Scalar vmax = 0.0;
    Set::Scalar uxmin = 1.0e300;
    Set::Scalar uxmax = -1.0e300;
    Set::Scalar uymin = 1.0e300;
    Set::Scalar uymax = -1.0e300;
    Set::Scalar divmax = 0.0;
    Set::Scalar div2 = 0.0;
    Set::Scalar ncell = 0.0;
    Set::Scalar divmax_interior = 0.0;
    Set::Scalar div2_interior = 0.0;
    Set::Scalar ncell_interior = 0.0;

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        const auto dx = geom[lev].CellSizeArray();
        amrex::Box domain = geom[lev].Domain();
        amrex::GpuArray<int,AMREX_SPACEDIM> periodic{};
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
            periodic[d] = geom[lev].isPeriodic(d);

        for (amrex::MFIter mfi(*velocity_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = velocity_mf[lev]->ixType() == amrex::IndexType::TheNodeType() ? mfi.tilebox() : mfi.tilebox();
            Set::Patch<const Set::Scalar> u = velocity_mf.Patch(lev,mfi);
            amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpMin, amrex::ReduceOpMax, amrex::ReduceOpMin, amrex::ReduceOpMax> reduce_op;
            amrex::ReduceData<Set::Scalar, Set::Scalar, Set::Scalar, Set::Scalar, Set::Scalar> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;
            reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple
            {
                Set::Scalar speed = std::sqrt(u(i,j,k,0) * u(i,j,k,0) + u(i,j,k,1) * u(i,j,k,1));
                return {speed, u(i,j,k,0), u(i,j,k,0), u(i,j,k,1), u(i,j,k,1)};
            });
            ReduceTuple hv = reduce_data.value();
            vmax = std::max(vmax, amrex::get<0>(hv));
            uxmin = std::min(uxmin, amrex::get<1>(hv));
            uxmax = std::max(uxmax, amrex::get<2>(hv));
            uymin = std::min(uymin, amrex::get<3>(hv));
            uymax = std::max(uymax, amrex::get<4>(hv));
        }

        for (amrex::MFIter mfi(*velocity_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            Set::Patch<const Set::Scalar> u = velocity_mf.Patch(lev,mfi);
            amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpSum, amrex::ReduceOpSum,
                            amrex::ReduceOpMax, amrex::ReduceOpSum, amrex::ReduceOpSum> reduce_op;
            amrex::ReduceData<Set::Scalar, Set::Scalar, Set::Scalar,
                            Set::Scalar, Set::Scalar, Set::Scalar> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;
            reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple
            {
                auto sten = Numeric::GetStencil(
                    i, j, k, domain, periodic);
                Set::Matrix grad_u =
                    Numeric::Gradient(u, i, j, k, dx.data(), sten);
                Set::Scalar div = 0.0;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) div += grad_u(d,d);
                bool interior = i > domain.smallEnd(0) + 2 && i < domain.bigEnd(0) - 2 &&
                                j > domain.smallEnd(1) + 2 && j < domain.bigEnd(1) - 2;
                Set::Scalar abs_div = Util::Abs(div);
                return {abs_div, div * div, 1.0,
                        interior ? abs_div : 0.0,
                        interior ? div * div : 0.0,
                        interior ? 1.0 : 0.0};
            });
            ReduceTuple hv = reduce_data.value();
            divmax = std::max(divmax, amrex::get<0>(hv));
            div2 += amrex::get<1>(hv);
            ncell += amrex::get<2>(hv);
            divmax_interior = std::max(divmax_interior, amrex::get<3>(hv));
            div2_interior += amrex::get<4>(hv);
            ncell_interior += amrex::get<5>(hv);
        }
    }

    amrex::ParallelDescriptor::ReduceRealMax(vmax);
    amrex::ParallelDescriptor::ReduceRealMin(uxmin);
    amrex::ParallelDescriptor::ReduceRealMax(uxmax);
    amrex::ParallelDescriptor::ReduceRealMin(uymin);
    amrex::ParallelDescriptor::ReduceRealMax(uymax);
    amrex::ParallelDescriptor::ReduceRealMax(divmax);
    amrex::ParallelDescriptor::ReduceRealSum(div2);
    amrex::ParallelDescriptor::ReduceRealSum(ncell);
    amrex::ParallelDescriptor::ReduceRealMax(divmax_interior);
    amrex::ParallelDescriptor::ReduceRealSum(div2_interior);
    amrex::ParallelDescriptor::ReduceRealSum(ncell_interior);
    Set::Scalar divrms = ncell > 0.0 ? std::sqrt(div2 / ncell) : 0.0;
    Set::Scalar divrms_interior = ncell_interior > 0.0 ? std::sqrt(div2_interior / ncell_interior) : 0.0;
    if (amrex::ParallelDescriptor::IOProcessor())
        amrex::Print() << "LowMach diagnostics step " << iter
                        << " time " << time
                        << " vmax " << vmax
                        << " uxmin " << uxmin
                        << " uxmax " << uxmax
                        << " uymin " << uymin
                        << " uymax " << uymax
                        << " divmax " << divmax
                        << " divrms " << divrms
                        << " divmax_interior " << divmax_interior
                        << " divrms_interior " << divrms_interior;
    if (amrex::ParallelDescriptor::IOProcessor() &&
        chemistry_timestep_mode != ChemistryTimestepMode::Off)
        amrex::Print() << " chemistry_timescale_min "
                        << chemistry_timescale_min
                        << " chemistry_timestep_candidate "
                        << chemistry_timestep_candidate
                        << " chemistry_timestep_limits "
                        << (chemistry_timestep_mode ==
                            ChemistryTimestepMode::Limit);
    if (amrex::ParallelDescriptor::IOProcessor())
    {
        for (int b = 0;
             b < static_cast<int>(free_rigid_solid_species.size()); ++b)
        {
            const FreeRigidBodyState& body = free_rigid_bodies[b];
            if (!body.valid) continue;
            const std::string prefix = "free_rigid_" +
                species_names[free_rigid_solid_species[b]];
            amrex::Print() << " " << prefix << "_mass " << body.mass
                            << " " << prefix << "_center_x " << body.center(0)
                            << " " << prefix << "_velocity_x " << body.velocity(0);
#if AMREX_SPACEDIM > 1
            amrex::Print() << " " << prefix << "_center_y " << body.center(1)
                            << " " << prefix << "_velocity_y " << body.velocity(1)
                            << " " << prefix << "_omega " << body.angular_velocity(0);
#endif
#if AMREX_SPACEDIM > 2
            amrex::Print() << " " << prefix << "_center_z " << body.center(2)
                            << " " << prefix << "_velocity_z " << body.velocity(2)
                            << " " << prefix << "_omega_y " << body.angular_velocity(1)
                            << " " << prefix << "_omega_z " << body.angular_velocity(2);
#endif
            amrex::Print() << " " << prefix << "_radius_of_gyration "
                            << body.radius_of_gyration
                            << " " << prefix << "_coupling_iterations "
                            << last_free_rigid_coupling_iterations[b]
                            << " " << prefix << "_coupling_residual "
                            << last_free_rigid_coupling_residual[b]
                            << " " << prefix << "_coupling_converged "
                            << last_free_rigid_coupling_converged[b];
        }
    }
    if (amrex::ParallelDescriptor::IOProcessor()) amrex::Print() << "\n";
}

void
LowMach::TimeStepComplete(Set::Scalar time, int iter)
{
    ApplyImplicitDiffusion(time + dt[0], dt[0]);
    ApplyInterfacialTransport(time, dt[0]);
    ApplyImplicitPhaseChange(time + dt[0], dt[0]);
    ProjectVelocity(time + dt[0], dt[0]);
    PrintDiagnostics(time, iter);
}

void
LowMach::PreparePlotFile(Set::Scalar /*time*/, const amrex::Vector<int>& /*iter*/)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        UpdateComponentState(lev, *component_density_mf[lev]);
        if (!liquid_species.empty())
        {
            const int nliquid = static_cast<int>(liquid_species.size());
            const int gas_phase = nliquid +
                static_cast<int>(interfacial_solid_species.size());
            amrex::MultiFab::Copy(
                *gas_volume_fraction_mf[lev],
                *interfacial_volume_fraction_mf[lev],
                gas_phase, 0, 1, 0);
        }
        if (interfacial_forces_enabled)
            UpdateInterfacialChemicalPotential(lev);
        if (deformable_solid_species >= 0)
            UpdateSolidStress(lev, *velocity_mf[lev], *eta_mf[lev], *xi_mf[lev], diagnostics_extended_fields);
        UpdateDerivedDiagnostics(lev, *velocity_mf[lev], *temperature_mf[lev]);
    }
}

void
LowMach::TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/)
{
    const auto dx = geom[lev].CellSizeArray();
    Set::Scalar dr = 0.0;
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
        dr += dx[d] * dx[d];
    dr = std::sqrt(dr);
    const bool deformable_solid = deformable_solid_species >= 0;
    const bool rigid_solid = !rigid_solid_species.empty();
    const bool liquid = !liquid_species.empty();
    const int nliquid = static_cast<int>(liquid_species.size());
    const int ngas = ngas_species;
    const Set::Scalar rho_floor = density_floor;
    const ThermochemicalData reaction_data = {
        thermal_data, chemistry_device_data, true, false, true};
    amrex::GpuArray<int,AMREX_SPACEDIM> periodic{};
    for (int d = 0; d < AMREX_SPACEDIM; ++d)
        periodic[d] = geom[lev].isPeriodic(d);

    for (amrex::MFIter mfi(*temperature_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        Set::Patch<char> tag = tags.array(mfi);
        Set::Patch<const Set::Scalar> u = velocity_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> pressure = pressure_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> T = temperature_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> component_density =
            component_density_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> eta = eta_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> rigid_eta = rigid_eta_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> interfacial_volume_fraction =
            interfacial_volume_fraction_mf.Patch(lev,mfi);
        const Set::Scalar vcrit = velocity_refinement_criterion;
        const Set::Scalar pcrit = pressure_refinement_criterion;
        const Set::Scalar Tcrit = temperature_refinement_criterion;
        const Set::Scalar reaction_crit = reaction_refinement_criterion;
        const Set::Scalar etacrit = eta_refinement_criterion;
        amrex::Box domain = geom[lev].Domain();

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            auto sten = Numeric::GetStencil(i, j, k, domain, periodic);
            Set::Matrix grad_u =
                Numeric::Gradient(u, i, j, k, dx.data(), sten);
            Set::Vector grad_p =
                Numeric::Gradient(
                    pressure, i, j, k, 0, dx.data(), sten);
            Set::Vector grad_T =
                Numeric::Gradient(T, i, j, k, 0, dx.data(), sten);
            Set::Scalar reaction_rate = 0.0;
            if (reaction_crit < 1.0e100)
            {
                Set::Scalar gas_density = 0.0;
                for (int n = 0; n < ngas; ++n)
                    gas_density +=
                        Util::Max(component_density(i,j,k,n), 0.0);
                if (gas_density > rho_floor)
                {
                    // Conservative transport can leave small negative species
                    // undershoots where split chemistry is inactive.  Project
                    // the diagnostic state just as AdvanceChemistry does; AMR
                    // tagging must not turn those undershoots into a chemistry
                    // model abort.
                    const auto source = ComputeThermochemicalSource(
                        component_density, T, i, j, k, dx.data(), 0.0,
                        reaction_data, true, true);
                    const auto& species_source = amrex::get<0>(source);
                    for (int n = 0; n < ngas; ++n)
                        reaction_rate += Util::Abs(species_source[n]);
                    reaction_rate /= gas_density;
                }
            }
            bool refine_eta = false;
            if (deformable_solid)
            {
                Set::Vector grad_eta =
                    Numeric::Gradient(
                        eta, i, j, k, 0, dx.data(), sten);
                const Set::Scalar eta_val = eta(i,j,k);
                refine_eta = grad_eta.lpNorm<2>() * dr * 2.0 > etacrit ||
                            (eta_val > etacrit && eta_val < 1.0 - etacrit);
            }
            if (rigid_solid)
            {
                Set::Vector grad_eta =
                    Numeric::Gradient(
                        rigid_eta, i, j, k, 0, dx.data(), sten);
                const Set::Scalar eta_val = rigid_eta(i,j,k);
                refine_eta = refine_eta || grad_eta.lpNorm<2>() * dr * 2.0 > etacrit ||
                            (eta_val > etacrit && eta_val < 1.0 - etacrit);
            }
            if (liquid)
                for (int n = 0; n < nliquid; ++n)
                {
                    Set::Vector grad_eta = Numeric::Gradient(
                        interfacial_volume_fraction, i, j, k, n,
                        dx.data(), sten);
                    const Set::Scalar eta_val =
                        interfacial_volume_fraction(i,j,k,n);
                    refine_eta = refine_eta ||
                        grad_eta.lpNorm<2>() * dr * 2.0 > etacrit ||
                        (eta_val > etacrit && eta_val < 1.0 - etacrit);
                }
            if (grad_u.norm() * dr > vcrit ||
                grad_p.lpNorm<2>() * dr > pcrit ||
                grad_T.lpNorm<2>() * dr > Tcrit ||
                reaction_rate > reaction_crit ||
                refine_eta)
                tag(i,j,k) = amrex::TagBox::SET;
        });
    }
}

}
