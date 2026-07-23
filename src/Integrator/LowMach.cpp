#include "LowMach.H"

#include <cctype>
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
LowMach::LowMach(IO::ParmParse& pp) : LowMach()
{
    pp_queryclass(*this);
}

void
LowMach::Parse(LowMach& value, IO::ParmParse& pp)
{
    BL_PROFILE("Integrator::LowMach::Parse");

    pp.query_required("cfl", value.cfl);
    pp.query_default("cfl_v", value.cfl_v, 1.0e100);
    pp.query_default("small", value.small, 1.0e-12);
    pp.query_default("density_floor", value.density_floor, value.small);
    pp.query_default("pressure_floor", value.pressure_floor, value.small);
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
        else
            Util::Exception(INFO, mechanics,
                            " is not a valid mechanics type for species ", name);
        if (n < value.ngas_species && mechanics != "fluid")
            Util::Exception(INFO, "The first ", value.ngas_species,
                            " species must be fluid species described by the gas model");

        if (mechanics != "fluid")
        {
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
            pp.select<IC::Constant,IC::Expression>(
                name + ".density.ic", value.component_density_ic[n],
                value.geom, Unit::Density());
    }

    value.chemistry.Define(value.ngas_species);
    pp.queryclass("chemistry", value.chemistry);
    if (value.chemistry.Split() && !value.projection_enabled)
        Util::Exception(INFO, "Locally integrated LowMach chemistry requires projection.enabled=1");
    value.implicit_thermal_diffusion = value.include_conduction;
    value.implicit_species_diffusion = value.ngas_species > 1 &&
        value.gas.transport.CommonDiffusivity();
    if (value.implicit_momentum_diffusion || value.implicit_thermal_diffusion ||
        value.implicit_species_diffusion)
        pp.queryclass("diffusion", value.diffusion);

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
        pp.query_required("rigid.relaxation_time",
                        value.rigid_relaxation_time, Unit::Time());
        pp.queryarr_default("rigid.velocity", value.rigid_velocity, Set::Vector::Zero());
        if (value.rigid_relaxation_time <= 0.0)
            Util::Exception(INFO, "rigid.relaxation_time must be positive");
        for (const int n : value.rigid_solid_species)
        {
            pp.query_required(value.species_names[n] + ".reference_density",
                            value.reference_density[n], Unit::Density());
            if (value.reference_density[n] <= 0.0)
                Util::Exception(INFO, value.species_names[n],
                                ".reference_density must be positive");
        }
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

    const bool has_thermal_bridge_width =
        pp.contains("thermal_bridge.width");
    const bool has_thermal_bridge_peclet =
        pp.contains("thermal_bridge.peclet");
    if (has_thermal_bridge_width != has_thermal_bridge_peclet)
        Util::Exception(INFO,
            "thermal_bridge requires both width and peclet");
    if (has_thermal_bridge_width)
    {
        if (!value.condensed_thermal_transport || !value.include_conduction)
            Util::Exception(INFO,
                "thermal_bridge requires condensed thermal transport");
        pp.query_required("thermal_bridge.width",
            value.thermal_bridge_width, Unit::Length());
        pp.query_required("thermal_bridge.peclet",
            value.thermal_bridge_peclet);
        if (!(value.thermal_bridge_width > 0.0) ||
            !(value.thermal_bridge_peclet > 0.0))
            Util::Exception(INFO,
                "thermal_bridge width and peclet must be positive");
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
        pp.select<Model::Mechanism::PhaseChange>(
            id, value.mechanisms[n], value.species_names, value.ngas_species,
            value.rigid_solid_species, value.reference_density, value.gas.MW,
            value.gas.Rg);
    }

    pp.select<Numeric::Advect::MUSCL,
            Numeric::Advect::Upwind,
            Numeric::Advect::Centered,
            Numeric::Advect::QUICK,
            Numeric::Advect::WENO5>("advection",value.advect);
    if (value.advect.PhiLocation() != Set::HC::Cell ||
        value.advect.VelocityLocation() != Set::HC::Cell)
        Util::Exception(INFO, "LowMach currently requires cell-centered phi and velocity advection data");

    pp.query_default("velocity_refinement_criterion", value.velocity_refinement_criterion, 1.0e100);
    pp.query_default("pressure_refinement_criterion", value.pressure_refinement_criterion, 1.0e100);
    pp.query_default("temperature_refinement_criterion", value.temperature_refinement_criterion, 1.0e100);
    if (value.deformable_solid_species >= 0 || !value.rigid_solid_species.empty())
        pp.query_default("eta_refinement_criterion", value.eta_refinement_criterion, 1.0e100);
    pp.queryarr_default("g", value.g, Set::Vector::Zero());

    int nghost = value.advect.NGhost();
    if (nghost < 3) nghost = 3;
    pp.select_default<BC::Constant,BC::Expression>("velocity.bc", value.velocity_bc, AMREX_SPACEDIM);
    pp.select_default<BC::Constant,BC::Expression>("temperature.bc", value.temperature_bc, 1);
    pp.select_default<BC::Constant::ZeroNeumann,BC::Constant,BC::Expression>("component_density.bc", value.component_density_bc, value.nspecies);
    pp.select_default<BC::Constant,BC::Expression>("pressure.bc", value.pressure_bc, 1);

    pp.select_default<IC::Constant,IC::Expression>("velocity.ic", value.velocity_ic, value.geom);
    pp.select_default<IC::Constant,IC::Expression>("temperature.ic", value.temperature_ic, value.geom);
    pp.select_default<IC::Constant,IC::Expression>("pressure.ic", value.pressure_ic, value.geom);
    if (value.deformable_solid_species >= 0)
    {
        pp.select_default<BC::Constant::ZeroNeumann,BC::Constant,BC::Expression>("xi.bc", value.xi_bc, AMREX_SPACEDIM);
        pp.select_default<IC::Expression::X,IC::Constant,IC::Expression>("xi.ic", value.xi_ic, value.geom);
    }

    std::vector<std::string> species_suffix(value.nspecies);
    for (int n = 0; n < value.nspecies; ++n)
        species_suffix[n] = "_" + value.species_names[n];
    std::vector<std::string> gas_species_suffix(
        species_suffix.begin(), species_suffix.begin() + value.ngas_species);

    value.AddField<Set::Scalar,Set::HC::Cell>(value.velocity_mf,        value.velocity_bc,      AMREX_SPACEDIM, nghost, "velocity",          true,  true, {"x","y"});
    value.AddField<Set::Scalar,Set::HC::Cell>(value.velocity_old_mf,    value.velocity_bc,      AMREX_SPACEDIM, nghost, "velocity_old",      false, true, {"x","y"});
    value.AddField<Set::Scalar,Set::HC::Cell>(value.temperature_mf,       value.temperature_bc,   1,              nghost, "temperature",       true,  true);
    value.AddField<Set::Scalar,Set::HC::Cell>(value.temperature_old_mf,   value.temperature_bc,   1,              nghost, "temperature_old",   false, true);
    value.AddField<Set::Scalar,Set::HC::Cell>(value.component_density_mf,     value.component_density_bc, value.nspecies, nghost, "component_density",     true,  true, species_suffix);
    value.AddField<Set::Scalar,Set::HC::Cell>(value.component_density_old_mf, value.component_density_bc, value.nspecies, nghost, "component_density_old", false, true, species_suffix);
    if (value.deformable_solid_species >= 0)
    {
        value.AddField<Set::Scalar,Set::HC::Cell>(value.eta_mf, &value.bc_nothing, 1, nghost, "eta", true, false);
        value.AddField<Set::Scalar,Set::HC::Cell>(value.xi_mf,     value.xi_bc, AMREX_SPACEDIM, nghost, "xi",     true,  true, {"x","y"});
        value.AddField<Set::Scalar,Set::HC::Cell>(value.xi_old_mf, value.xi_bc, AMREX_SPACEDIM, nghost, "xi_old", false, true, {"x","y"});
    }
    if (!value.rigid_solid_species.empty())
        value.AddField<Set::Scalar,Set::HC::Cell>(value.rigid_eta_mf, &value.bc_nothing, 1, 1, "rigid_eta", true, false);

    value.AddField<Set::Scalar,Set::HC::Cell>(value.density_mf,             &value.bc_nothing, 1,              1,      "density",             true,  false);
    value.AddField<Set::Scalar,Set::HC::Cell>(value.pressure_mf,            value.pressure_bc, 1,              nghost, "pressure",            true,  true);
    value.AddField<Set::Scalar,Set::HC::Cell>(value.pressure_correction_mf, &value.bc_nothing, 1,              nghost, "pressure_correction", true, false);
    if (value.chemistry.Split())
        value.AddField<Set::Scalar,Set::HC::Cell>(value.chemistry_dilatation_mf,
            &value.bc_nothing, 1, 0, "chemistry_dilatation", false, false);
    if (value.implicit_thermal_diffusion || value.implicit_species_diffusion)
        value.AddField<Set::Scalar,Set::HC::Cell>(value.diffusion_dilatation_mf,
            &value.bc_nothing, 1, 0, "diffusion_dilatation", false, false);
    if (value.deformable_solid_species >= 0)
        value.AddField<Set::Matrix,Set::HC::Cell>(value.solid_deviatoric_stress_mf, nullptr, 1, 1, "solid_deviatoric_stress", true, false);
    if (value.diagnostics_extended_fields)
    {
        value.AddField<Set::Scalar,Set::HC::Cell>(value.mass_fraction_mf, &value.bc_nothing, value.ngas_species, 1, "mass_fraction", true, false, gas_species_suffix);
        value.AddField<Set::Scalar,Set::HC::Cell>(value.mole_fraction_mf, &value.bc_nothing, value.ngas_species, 1, "mole_fraction", true, false, gas_species_suffix);
        value.AddField<Set::Scalar,Set::HC::Cell>(value.momentum_mf, &value.bc_nothing, AMREX_SPACEDIM, 1, "momentum", true, false, {"x","y"});
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

    bool allow_unused;
    pp.query_default("allow_unused", allow_unused, false);
    if (!allow_unused && pp.AnyUnusedInputs(true, false))
    {
        Util::Warning(INFO, "The following inputs were specified but not used:");
        pp.AllUnusedInputs();
        Util::Exception(INFO, "Aborting. Specify 'allow_unused=True` to ignore this error.");
    }
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

    const Set::Scalar* DX = geom[lev].CellSize();
    amrex::Box domain = geom[lev].Domain();
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
            auto sten = Numeric::GetStencil(i, j, k, domain);
            Set::Matrix F = Set::Matrix::Identity();
            Set::Matrix grad_u = Numeric::Gradient(u, i, j, k, DX, sten);

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
                Set::Matrix grad_xi = Numeric::Gradient(xi, i, j, k, DX, sten);
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
// - rigid_eta_mf   (if there is a rigid solid present, calculated based on densities)
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
    if (deformable_solid) eta_mf[lev]->setVal(0.0);
    if (rigid_solid) rigid_eta_mf[lev]->setVal(0.0);
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

        amrex::ParallelFor(bx, [=,this] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Scalar gas_partial_density = 0.0;
            Set::Scalar total_partial_density = 0.0;
            for (int n = 0; n < ngas_species; ++n)
                gas_partial_density += component_density(i,j,k,n);
            for (int n = 0; n < nspecies; ++n)
                total_partial_density += component_density(i,j,k,n);

            // Partial densities are the conserved species state. Mechanical
            // volume fractions are reconstructed from them when needed.
            rho(i,j,k) = total_partial_density;
            if (deformable_solid)
                eta(i,j,k) = component_density(i,j,k,deformable_solid_species) / finite_solid_reference_density;

            if (diagnostics_extended_fields)
            {
                const Model::Gas::MoleFraction X = gas.MoleFractions(component_density, i, j, k);
                Set::Scalar gas_volume_fraction = 0.0;
                if (gas_partial_density > density_floor && T(i,j,k) > 0.0 &&
                    pressure_reference > 0.0)
                    gas_volume_fraction = Util::Clamp(
                        gas_partial_density * gas.R(X, i, j, k) * T(i,j,k) /
                            pressure_reference,
                        0.0, 1.0);
                Set::Scalar condensed_volume_fraction = 0.0;
                for (int n = ngas_species; n < nspecies; ++n)
                    condensed_volume_fraction += Util::Max(
                        component_density(i,j,k,n), 0.0) *
                        condensed_inverse_reference_density[n];
                gas_volume_fraction = Util::Min(gas_volume_fraction,
                    1.0 - Model::PhaseField::H(condensed_volume_fraction));
                // Diagnostic fractions include phase occupancy; gas-model
                // calculations continue to use the conditional composition X.
                for (int n = 0; n < ngas_species; ++n)
                {
                    mass_fraction(i,j,k,n) = total_partial_density > density_floor ?
                        component_density(i,j,k,n) / total_partial_density : 0.0;
                    mole_fraction(i,j,k,n) = gas_volume_fraction * X(i,j,k,n);
                }
            }
        });
    }

    density_mf[lev]->FillBoundary(geom[lev].periodicity());
    if (deformable_solid)
        eta_mf[lev]->FillBoundary(geom[lev].periodicity());
    if (rigid_solid)
    {
        for (const int n : rigid_solid_species)
            amrex::MultiFab::Saxpy(*rigid_eta_mf[lev], 1.0 / reference_density[n],
                                    component_density_mf, n, 0, 1, 1);
        rigid_eta_mf[lev]->FillBoundary(geom[lev].periodicity());
    }
    if (diagnostics_extended_fields)
    {
        mass_fraction_mf[lev]->FillBoundary(geom[lev].periodicity());
        mole_fraction_mf[lev]->FillBoundary(geom[lev].periodicity());
    }
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

    const Set::Scalar* DX = geom[lev].CellSize();
    amrex::Box domain = geom[lev].Domain();
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

        amrex::ParallelFor(bx, [=,this] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            const Model::Gas::MoleFraction X = gas.MoleFractions(component_density, i, j, k);
            const Set::Scalar density = rho(i,j,k);
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
                M(i,j,k,d) = density * u(i,j,k,d);
            E(i,j,k) = gas.ComputeE(density, M(i,j,k,0), M(i,j,k,1), T(i,j,k), X, i, j, k);

            auto sten = Numeric::GetStencil(i, j, k, domain);
            Set::Matrix grad_u = Numeric::Gradient(u, i, j, k, DX, sten);
            omega(i,j,k) = grad_u(1,0) - grad_u(0,1);

            mu(i,j,k) = gas.dynamic_viscosity(T(i,j,k), X, i, j, k);
            auto [thermal_gas_volume_fraction, gas_heat_capacity,
                heat_capacity, conductivity, cp] = ComputeThermalState(
                    component_density, T(i,j,k), i, j, k);
            (void)thermal_gas_volume_fraction;
            (void)gas_heat_capacity;
            (void)heat_capacity;
            (void)cp;
            kappa(i,j,k) = conductivity;
            for (int n = 0; n < ngas_species; ++n)
                diffusion(i,j,k,n) = gas.diffusion_coefficient(
                    T(i,j,k), pressure_reference, X, i, j, k, n);

            Model::Chemistry::SpeciesArray rhoY{};
            Set::Scalar gas_density = 0.0;
            for (int n = 0; n < ngas_species; ++n)
            {
                rhoY[n] = Util::Max(component_density(i,j,k,n), 0.0);
                gas_density += rhoY[n];
            }
            Set::Scalar gas_volume_fraction = 0.0;
            if (gas_density > density_floor && T(i,j,k) > 0.0 &&
                pressure_reference > 0.0)
                gas_volume_fraction = Util::Clamp(
                    gas_density * gas.R(X, i, j, k) * T(i,j,k) /
                        pressure_reference, 0.0, 1.0);
            const Set::Scalar intrinsic_gas_volume_fraction =
                gas_volume_fraction;
            Set::Scalar condensed_volume_fraction = 0.0;
            for (int n = ngas_species; n < nspecies; ++n)
                condensed_volume_fraction += Util::Max(
                    component_density(i,j,k,n), 0.0) *
                    condensed_inverse_reference_density[n];
            gas_volume_fraction = Util::Min(gas_volume_fraction,
                1.0 - Model::PhaseField::H(condensed_volume_fraction));
            if (gas_volume_fraction > 0.0 &&
                intrinsic_gas_volume_fraction > 0.0)
            {
                for (int n = 0; n < ngas_species; ++n)
                    rhoY[n] /= intrinsic_gas_volume_fraction;
                const Model::Chemistry::Source reaction =
                    chemistry.ComputeChemistrySources(
                        pressure_reference, T(i,j,k), rhoY, 0.0, &gas);
                for (int n = 0; n < ngas_species; ++n)
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

    for (amrex::MFIter mfi(T_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        Set::Patch<Set::Scalar> T = T_mf.array(mfi);
        Set::Patch<Set::Scalar> component_density = component_density_mf.array(mfi);
        Set::Patch<Set::Scalar> integrated_dilatation = chemistry_dilatation_mf.Patch(lev,mfi);

        amrex::ParallelFor(bx, [=,this] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Model::Chemistry::SpeciesArray rhoY{};
            Set::Scalar gas_density = 0.0;
            for (int n = 0; n < ngas_species; ++n)
            {
                rhoY[n] = Util::Max(component_density(i,j,k,n), 0.0);
                gas_density += rhoY[n];
            }
            if (!(gas_density > density_floor) || !(T(i,j,k) > 0.0)) return;

            const Model::Gas::MoleFraction X =
                gas.MoleFractions(component_density, i, j, k);
            const Set::Scalar raw_gas_volume_fraction = Util::Clamp(
                gas_density * gas.R(X, i, j, k) * T(i,j,k) /
                    pressure_reference, 0.0, 1.0);
            Set::Scalar condensed_volume_fraction = 0.0;
            for (int n = ngas_species; n < nspecies; ++n)
                condensed_volume_fraction += Util::Max(
                    component_density(i,j,k,n), 0.0) *
                    condensed_inverse_reference_density[n];
            const Set::Scalar gas_accessibility =
                1.0 - Model::PhaseField::H(condensed_volume_fraction);
            const Set::Scalar chemistry_weight =
                raw_gas_volume_fraction > 0.0 ? Util::Min(
                    1.0, gas_accessibility /
                        raw_gas_volume_fraction) : 0.0;
            if (!(chemistry_weight > 0.0)) return;

            Set::Scalar temperature = T(i,j,k);
            auto [thermal_gas_volume_fraction, gas_heat_capacity, heat_capacity,
                conductivity, cp] = ComputeThermalState(
                    component_density, temperature, i, j, k);
            (void)thermal_gas_volume_fraction;
            (void)gas_heat_capacity;
            (void)conductivity;
            const Set::Scalar mixture_density = cp > 0.0 ?
                heat_capacity / cp : gas_density;

            auto result = chemistry.Advance(
                dt * chemistry_weight, pressure_reference, mixture_density,
                ngas_species,
                rhoY, temperature, &gas);
            if (!result.converged)
                Util::Abort(INFO, "Local chemistry integration failed at ",
                    i, ",", j, ",", k, " dt=", result.failed_dt,
                    " residual=", result.residual_norm);

            for (int n = 0; n < ngas_species; ++n)
                component_density(i,j,k,n) = rhoY[n];
            T(i,j,k) = temperature;
            integrated_dilatation(i,j,k) += result.volume_fraction_change;
        });
    }
}

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
std::tuple<Set::Scalar, Set::Scalar, Set::Scalar, Set::Scalar, Set::Scalar>
LowMach::ComputeThermalState(
    Set::Patch<const Set::Scalar> component_density,
    Set::Scalar temperature, int i, int j, int k) const
{
    Set::Scalar density = 0.0;
    Set::Scalar gas_density = 0.0;
    for (int n = 0; n < nspecies; ++n)
        density += Util::Max(component_density(i,j,k,n), 0.0);
    for (int n = 0; n < ngas_species; ++n)
        gas_density += Util::Max(component_density(i,j,k,n), 0.0);

    const Model::Gas::MoleFraction X =
        gas.MoleFractions(component_density, i, j, k);
    const Set::Scalar cp = gas.cp_mass(temperature, X, i, j, k);
    const Set::Scalar gas_conductivity =
        gas.thermal_conductivity(temperature, X, i, j, k);
    Set::Scalar gas_volume_fraction = 0.0;
    if (gas_density > density_floor && temperature > 0.0 &&
        pressure_reference > 0.0)
        gas_volume_fraction = Util::Clamp(
            gas_density * gas.R(X, i, j, k) * temperature /
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
    const Set::Scalar gas_constant = gas.R(X, i, j, k);
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

void
LowMach::ApplyImplicitDiffusion(Set::Scalar time, Set::Scalar dt)
{
    BL_PROFILE("Integrator::LowMach::ApplyImplicitDiffusion");
    if ((!implicit_momentum_diffusion && !implicit_thermal_diffusion &&
        !implicit_species_diffusion) || !(dt > 0.0)) return;

    const int nlev = finest_level + 1;
    const bool thermochemical_diffusion =
        implicit_thermal_diffusion || implicit_species_diffusion;
    diffusion.SetLayout(geom, refRatio(), temperature_mf, nlev, 1);
    diffusion.FillBoundary(temperature_mf, *temperature_bc, time, 1);
    diffusion.FillBoundary(
        component_density_mf, *component_density_bc, time, nspecies);
    if (implicit_momentum_diffusion)
        diffusion.FillBoundary(
            velocity_mf, *velocity_bc, time, AMREX_SPACEDIM);
    if (thermochemical_diffusion)
    {
        for (int lev = 0; lev < nlev; ++lev)
        {
            diffusion_dilatation_mf[lev]->setVal(0.0);

            for (amrex::MFIter mfi(*diffusion_dilatation_mf[lev],
                                    amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                Set::Patch<const Set::Scalar> component_density =
                    component_density_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> T = temperature_mf.Patch(lev,mfi);
                Set::Patch<Set::Scalar> volume_change =
                    diffusion_dilatation_mf.Patch(lev,mfi);
                const Set::Scalar p_reference = pressure_reference;
                const int ngas = ngas_species;

                amrex::ParallelFor(bx, [=,this] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar gas_density = 0.0;
                    for (int n = 0; n < ngas; ++n)
                        gas_density += component_density(i,j,k,n);
                    const Model::Gas::MoleFraction X =
                        gas.MoleFractions(component_density, i, j, k);
                    volume_change(i,j,k) = -gas_density *
                        gas.R(X, i, j, k) * T(i,j,k) / p_reference;
                });
            }
        }
    }

    if (implicit_species_diffusion)
    {
        diffusion.SetLayout(
            geom, refRatio(), component_density_mf, nlev, ngas_species);
        for (int lev = 0; lev < nlev; ++lev)
        {
            amrex::MultiFab& state = diffusion.State(lev, ngas_species);
            amrex::MultiFab& mass = diffusion.Mass(lev, ngas_species);
            amrex::MultiFab& mobility = diffusion.Mobility(lev, ngas_species);
            state.setVal(0.0);
            mass.setVal(0.0);
            mobility.setVal(0.0);

            for (amrex::MFIter mfi(state, amrex::TilingIfNotGPU());
                mfi.isValid(); ++mfi)
            {
                const amrex::Box& grown_box = mfi.growntilebox(1);
                const amrex::Box& valid_box = mfi.tilebox();
                Set::Patch<const Set::Scalar> component_density =
                    component_density_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> T = temperature_mf.Patch(lev,mfi);
                Set::Patch<Set::Scalar> Y = state.array(mfi);
                Set::Patch<Set::Scalar> a = mass.array(mfi);
                Set::Patch<Set::Scalar> b = mobility.array(mfi);
                const Set::Scalar p_reference = pressure_reference;
                const Set::Scalar rho_floor = density_floor;
                const int ngas = ngas_species;

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
                    valid_box, [=,this] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        Set::Scalar gas_density = 0.0;
                        for (int n = 0; n < ngas; ++n)
                            gas_density += component_density(i,j,k,n);
                        const Model::Gas::MoleFraction X =
                            gas.MoleFractions(component_density, i, j, k);
                        const Set::Scalar gas_volume_fraction = Util::Clamp(
                            gas_density * gas.R(X, i, j, k) * T(i,j,k) /
                                p_reference, 0.0, 1.0);
                        Set::Scalar condensed_volume_fraction = 0.0;
                        for (int n = ngas; n < nspecies; ++n)
                            condensed_volume_fraction += Util::Max(
                                component_density(i,j,k,n), 0.0) *
                                condensed_inverse_reference_density[n];
                        const Set::Scalar gas_accessibility =
                            1.0 - Model::PhaseField::H(
                                condensed_volume_fraction);
                        const Set::Scalar diffusion_weight =
                            gas_volume_fraction > 0.0 ? Util::Min(
                                1.0, gas_accessibility /
                                    gas_volume_fraction) : 0.0;
                        a(i,j,k) = Util::Max(gas_density, rho_floor);
                        b(i,j,k) = diffusion_weight * gas_density *
                            gas.diffusion_coefficient(
                                T(i,j,k), p_reference, X, i, j, k, 0);
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

    if (implicit_thermal_diffusion)
    {
        diffusion.SetLayout(geom, refRatio(), temperature_mf, nlev, 1);
        const auto inverse_reference_density =
            condensed_inverse_reference_density;
        const bool use_condensed_properties = condensed_thermal_transport;
        const bool use_thermal_bridge = thermal_bridge_peclet > 0.0;
        const int ngas = ngas_species;
        for (int lev = 0; lev < nlev; ++lev)
        {
            amrex::MultiFab& state = diffusion.State(lev, 1);
            amrex::MultiFab& mass = diffusion.Mass(lev, 1);
            amrex::MultiFab& mobility = diffusion.Mobility(lev, 1);
            amrex::MultiFab& tensor_mobility =
                diffusion.TensorMobility(lev, 1);
            amrex::MultiFab::Copy(
                state, *temperature_mf[lev], 0, 0, 1, state.nGrow());

            for (amrex::MFIter mfi(mass, amrex::TilingIfNotGPU());
                mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                Set::Patch<const Set::Scalar> component_density =
                    component_density_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> T = temperature_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> velocity =
                    velocity_mf.Patch(lev,mfi);
                Set::Patch<Set::Scalar> a = mass.array(mfi);
                Set::Patch<Set::Scalar> b = mobility.array(mfi);
                Set::Patch<Set::Scalar> B = tensor_mobility.array(mfi);
                const Set::Scalar rho_floor = density_floor;
                const int number_of_species = nspecies;
                const Set::Scalar bridge_width = thermal_bridge_width;
                const Set::Scalar inverse_bridge_peclet =
                    thermal_bridge_peclet > 0.0 ?
                    1.0 / thermal_bridge_peclet : 0.0;
                const auto DX = geom[lev].CellSizeArray();
                const amrex::Box domain = geom[lev].Domain();

                amrex::ParallelFor(bx, [=,this] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    auto [gas_volume_fraction, gas_heat_capacity,
                        heat_capacity, conductivity, cp] = ComputeThermalState(
                            component_density, T(i,j,k), i, j, k);
                    (void)gas_volume_fraction;
                    (void)gas_heat_capacity;
                    if (use_condensed_properties)
                    {
                        Set::Scalar condensed_volume_fraction = 0.0;
                        for (int n = ngas; n < number_of_species; ++n)
                        {
                            const Set::Scalar partial_density = Util::Max(
                                component_density(i,j,k,n), 0.0);
                            const Set::Scalar volume_fraction = partial_density *
                                inverse_reference_density[n];
                            condensed_volume_fraction += volume_fraction;
                        }
                        if (use_thermal_bridge)
                        {
                            const Set::Scalar eta = Model::PhaseField::H(
                                condensed_volume_fraction);
                            const Set::Scalar bridge_weight =
                                4.0 * eta * (1.0 - eta);
                            Set::Vector normal = Set::Vector::Zero();
                            const auto stencil =
                                Numeric::GetStencil(i,j,k,domain);
                            for (int n = ngas; n < number_of_species; ++n)
                                normal += inverse_reference_density[n] *
                                    Numeric::Gradient(component_density,
                                        i, j, k, n, DX.data(), stencil);
                            const Set::Scalar normal_norm = normal.norm();
                            if (normal_norm > 0.0) normal /= normal_norm;
                            const Set::Vector u(AMREX_D_DECL(
                                velocity(i,j,k,0), velocity(i,j,k,1),
                                velocity(i,j,k,2)));
                            const Set::Scalar bridge_conductivity =
                                bridge_weight * heat_capacity *
                                Util::Abs(u.dot(normal)) * bridge_width *
                                inverse_bridge_peclet;
                            for (int d = 0; d < AMREX_SPACEDIM; ++d)
                                for (int e = 0; e < AMREX_SPACEDIM; ++e)
                                    B(i,j,k,d * AMREX_SPACEDIM + e) =
                                        bridge_conductivity *
                                        normal(d) * normal(e);
                        }
                        b(i,j,k) = conductivity;
                        a(i,j,k) = Util::Max(
                            heat_capacity, rho_floor * cp);
                        return;
                    }
                    a(i,j,k) = Util::Max(heat_capacity, rho_floor * cp);
                    b(i,j,k) = conductivity;
                });
            }
        }

        diffusion.Solve(
            time, dt, temperature_bc->GetBCRec(), 1, use_thermal_bridge);
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

                amrex::ParallelFor(bx, [=,this] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar density = 0.0;
                    for (int n = 0; n < number_of_species; ++n)
                        density += component_density(i,j,k,n);
                    const Model::Gas::MoleFraction X =
                        gas.MoleFractions(component_density, i, j, k);
                    a(i,j,k) = Util::Max(density, rho_floor);
                    b(i,j,k) =
                        gas.dynamic_viscosity(T(i,j,k), X, i, j, k);
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

    if (thermochemical_diffusion)
    {
        for (int lev = 0; lev < nlev; ++lev)
        {
            for (amrex::MFIter mfi(*diffusion_dilatation_mf[lev],
                                    amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                Set::Patch<const Set::Scalar> component_density =
                    component_density_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> T = temperature_mf.Patch(lev,mfi);
                Set::Patch<Set::Scalar> volume_change =
                    diffusion_dilatation_mf.Patch(lev,mfi);
                const Set::Scalar p_reference = pressure_reference;
                const int ngas = ngas_species;

                amrex::ParallelFor(bx, [=,this] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar gas_density = 0.0;
                    for (int n = 0; n < ngas; ++n)
                        gas_density += component_density(i,j,k,n);
                    const Model::Gas::MoleFraction X =
                        gas.MoleFractions(component_density, i, j, k);
                    volume_change(i,j,k) += gas_density *
                        gas.R(X, i, j, k) * T(i,j,k) / p_reference;
                });
            }
        }
        diffusion.Synchronize(diffusion_dilatation_mf, 1);
    }
}

AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
std::tuple<Model::Chemistry::SpeciesArray, Set::Scalar /*temperature*/, Set::Scalar /*dilatation*/>
LowMach::ComputeThermochemicalSource(
    Set::Patch<const Set::Scalar> component_density,
    Set::Patch<const Set::Scalar> T,
    int i, int j, int k, const Set::Scalar* DX, Set::Scalar dt,
    bool include_reaction) const
{
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
        rhoY[n] = component_density(i,j,k,n);
        gas_density += rhoY[n];
        molar_density += rhoY[n] / gas.MW[n];
    }

    // Calculate and sanitize gas volume fractions
    const Model::Gas::MoleFraction X = gas.MoleFractions(component_density, i, j, k);
    Set::Scalar gas_volume_fraction = 0.0;
    if (gas_density > density_floor && T(i,j,k) > 0.0 && pressure_reference > 0.0)
        gas_volume_fraction = Util::Clamp(
            gas_density * gas.R(X, i, j, k) * T(i,j,k) / pressure_reference,
            0.0, 1.0);
    Set::Scalar condensed_volume_fraction = 0.0;
    for (int n = ngas_species; n < nspecies; ++n)
        condensed_volume_fraction += Util::Max(
            component_density(i,j,k,n), 0.0) *
            condensed_inverse_reference_density[n];
    const Set::Scalar gas_accessibility =
        1.0 - Model::PhaseField::H(condensed_volume_fraction);
    const Set::Scalar reacting_volume_fraction =
        Util::Min(gas_volume_fraction, gas_accessibility);
    const Set::Scalar chemistry_weight = gas_volume_fraction > 0.0 ?
        reacting_volume_fraction / gas_volume_fraction : 0.0;

    // Calculate relative density with respect to the gas volume fraction
    Model::Chemistry::SpeciesArray intrinsic_rhoY{};
    if (gas_volume_fraction > 0.0)
        for (int n = 0; n < ngas_species; ++n)
            intrinsic_rhoY[n] = rhoY[n] / gas_volume_fraction;

    Model::Chemistry::Source reaction{};
    if (include_reaction && reacting_volume_fraction > 0.0)
        reaction = chemistry.ComputeChemistrySources(
            pressure_reference, T(i,j,k), intrinsic_rhoY,
            dt * chemistry_weight, &gas);

    for (int n = 0; n < ngas_species; ++n)
        species[n] = reacting_volume_fraction * reaction.first[n];

    Set::Scalar enthalpy_diffusion = 0.0;
    if (ngas_species > 1 && !implicit_species_diffusion)
    {
        auto gas_density_at = [=,this] AMREX_GPU_DEVICE(int ii, int jj, int kk)
        {
            Set::Scalar value = 0.0;
            for (int n = 0; n < ngas_species; ++n)
                value += component_density(ii,jj,kk,n);
            return value;
        };
        auto mass_fraction_at = [=,this] AMREX_GPU_DEVICE(
            int ii, int jj, int kk, int n)
        {
            const Set::Scalar value = gas_density_at(ii,jj,kk);
            return value > density_floor ?
                component_density(ii,jj,kk,n) / value :
                (n == 0 ? 1.0 : 0.0);
        };
        auto transport_coefficient_at = [=,this] AMREX_GPU_DEVICE(
            int ii, int jj, int kk, int n)
        {
            const Set::Scalar local_gas_density =
                gas_density_at(ii,jj,kk);
            const Model::Gas::MoleFraction local_X =
                gas.MoleFractions(component_density, ii, jj, kk);
            const Set::Scalar local_gas_volume_fraction = Util::Clamp(
                local_gas_density * gas.R(local_X, ii, jj, kk) *
                    T(ii,jj,kk) / pressure_reference, 0.0, 1.0);
            Set::Scalar local_condensed_volume_fraction = 0.0;
            for (int m = ngas_species; m < nspecies; ++m)
                local_condensed_volume_fraction += Util::Max(
                    component_density(ii,jj,kk,m), 0.0) *
                    condensed_inverse_reference_density[m];
            const Set::Scalar local_accessibility =
                1.0 - Model::PhaseField::H(
                    local_condensed_volume_fraction);
            const Set::Scalar weight = local_gas_volume_fraction > 0.0 ?
                Util::Min(1.0, local_accessibility /
                    local_gas_volume_fraction) : 0.0;
            return weight * local_gas_density *
                gas.diffusion_coefficient(T(ii,jj,kk), pressure_reference,
                    local_X, ii, jj, kk, n);
        };
        const bool has_condensed_species = nspecies > ngas_species;
        auto face_coefficient = [=] AMREX_GPU_DEVICE(
            Set::Scalar a, Set::Scalar b)
        {
            if (!has_condensed_species) return 0.5 * (a + b);
            return a > 0.0 && b > 0.0 ? 2.0 * a * b / (a + b) : 0.0;
        };
        auto species_flux = [=,this] AMREX_GPU_DEVICE(
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
        auto conservative_species_flux = [=] AMREX_GPU_DEVICE(
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
                    gas.enthalpy_mol_species(T(i,j,k), n) / gas.MW[n];
                const Set::Scalar h_hi = 0.5 * (h +
                    gas.enthalpy_mol_species(T(i+di,j+dj,k+dk), n) /
                    gas.MW[n]);
                const Set::Scalar h_lo = 0.5 * (
                    gas.enthalpy_mol_species(T(i-di,j-dj,k-dk), n) /
                    gas.MW[n] + h);
                enthalpy_diffusion +=
                    ((h_hi - h) * flux_hi + (h - h_lo) * flux_lo) / DX[d];
            }
        }
    }

    auto [thermal_gas_volume_fraction, gas_heat_capacity, heat_capacity,
        conductivity, cp] = ComputeThermalState(
            component_density, T(i,j,k), i, j, k);
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
            molar_source += species[n] / gas.MW[n];
        dilatation += gas_volume_fraction * molar_source / molar_density;
    }
    return std::tie(species,temperature,dilatation);
}


void
LowMach::ProjectVelocity(Set::Scalar time, Set::Scalar dt)
{
    BL_PROFILE("Integrator::LowMach::ProjectVelocity");
    if (!projection_enabled || !(dt > 0.0)) return;

    const int nlev = finest_level + 1;
    const bool deformable_solid = deformable_solid_species >= 0;
    const bool rigid_solid = !rigid_solid_species.empty();
    const bool mixed_phase = deformable_solid || rigid_solid;
    const bool split_chemistry = chemistry.Split();
    const bool split_diffusion =
        implicit_thermal_diffusion || implicit_species_diffusion;
    if (!(pressure_reference == pressure_reference))
        pressure_reference = pressure_mf[0]->sum(0, false) /
                            static_cast<Set::Scalar>(geom[0].Domain().numPts());
    pressure_poisson.SetLayout(geom, refRatio(), velocity_mf, nlev);

    for (int lev = 0; lev < nlev; ++lev)
    {
        velocity_bc->define(geom[lev]);
        velocity_bc->FillBoundary(*velocity_mf[lev], 0, AMREX_SPACEDIM, time, 0);
        velocity_mf[lev]->FillBoundary(geom[lev].periodicity());
        UpdateComponentState(lev, *component_density_mf[lev]);

        if (rigid_solid)
        {
            // Eliminate the implicit Brinkman term locally; its mobility also
            // weights the pressure operator and correction below.
            const Set::Scalar inverse_relaxation_time = 1.0 / rigid_relaxation_time;
            const Set::Vector target_velocity = rigid_velocity;
            for (amrex::MFIter mfi(*velocity_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                Set::Patch<Set::Scalar> u = velocity_mf.Patch(lev,mfi);
                Set::Patch<const Set::Scalar> eta = rigid_eta_mf.Patch(lev,mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    const Set::Scalar weight =
                        Model::PhaseField::H(eta(i,j,k));
                    const Set::Scalar mobility = 1.0 / (1.0 + dt * weight * inverse_relaxation_time);
                    for (int d = 0; d < AMREX_SPACEDIM; ++d)
                        u(i,j,k,d) = mobility * u(i,j,k,d) +
                                    (1.0 - mobility) * target_velocity(d);
                });
            }
            velocity_bc->FillBoundary(*velocity_mf[lev], 0, AMREX_SPACEDIM, time, 0);
            velocity_mf[lev]->FillBoundary(geom[lev].periodicity());
        }
    }

    for (int lev = 0; lev < nlev; ++lev)
    {
        const Set::Scalar* DX = geom[lev].CellSize();
        amrex::MultiFab& u_mf = *velocity_mf[lev];
        amrex::MultiFab& beta_mf = pressure_poisson.Coefficient(lev);
        const Set::Scalar rho_floor = density_floor;
        const Set::Scalar inverse_relaxation_time = rigid_solid ?
            1.0 / rigid_relaxation_time : 0.0;

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
                Set::Patch<Set::Scalar> rhs = pressure_poisson.RHS(lev).array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    const Model::Mechanism::State state = {
                        component_density, rigid_eta, T(i,j,k), p_reference};
                    rhs(i,j,k) += mechanism.VolumeSource(state, i, j, k, DX);
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
            Set::Patch<const Set::Scalar> diffusion_dilatation;
            if (split_diffusion)
                diffusion_dilatation = diffusion_dilatation_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> eta = eta_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> rigid_eta = rigid_eta_mf.Patch(lev,mfi);
            Set::Patch<Set::Scalar> rhs = pressure_poisson.RHS(lev).array(mfi);
            const int ngas = ngas_species;
            const Set::Scalar p_reference = pressure_reference;
            const Set::Scalar inverse_dt = 1.0 / dt;

            amrex::ParallelFor(bx, [=,this] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                auto [species,temperature,dilatation] = ComputeThermochemicalSource(
                    component_density, T, i, j, k, DX, dt, !split_chemistry);

                rhs(i,j,k) += dilatation;
                if (split_chemistry)
                    rhs(i,j,k) += chemistry_dilatation(i,j,k) / dt;
                if (split_diffusion)
                    rhs(i,j,k) += diffusion_dilatation(i,j,k) / dt;

                if (mixed_phase)
                {
                    Set::Scalar gas_density = 0.0;
                    for (int n = 0; n < ngas; ++n)
                        gas_density += component_density(i,j,k,n);
                    const Model::Gas::MoleFraction X =
                        gas.MoleFractions(component_density, i, j, k);
                    Set::Scalar volume_fraction = gas_density *
                        gas.R(X, i, j, k) * T(i,j,k) / p_reference;
                    if (deformable_solid) volume_fraction += eta(i,j,k);
                    if (rigid_solid) volume_fraction += rigid_eta(i,j,k);

                    // Correct splitting and collocated-flux drift in the
                    // mixture volume constraint over one flow step.
                    rhs(i,j,k) += (volume_fraction - 1.0) * inverse_dt /
                        Util::Max(volume_fraction, 0.1);
                }
            });
        }

        for (amrex::MFIter mfi(pressure_poisson.RHS(lev), amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            Set::Patch<const Set::Scalar> rho = density_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> rigid_eta = rigid_eta_mf.Patch(lev,mfi);
            Set::Patch<Set::Scalar> beta = beta_mf.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                const Set::Scalar weight = rigid_solid ?
                    Model::PhaseField::H(rigid_eta(i,j,k)) : 0.0;
                const Set::Scalar mobility = 1.0 / (1.0 + dt * weight * inverse_relaxation_time);
                beta(i,j,k) = mobility / Util::Max(rho(i,j,k), rho_floor);
            });
        }
        pressure_poisson.PrepareRHS(lev, u_mf, dt);
    }

    pressure_poisson.Solve(time, pressure_bc->GetBCRec());

    for (int lev = 0; lev < nlev; ++lev)
        amrex::MultiFab::Copy(*pressure_correction_mf[lev],
                            pressure_poisson.Solution(lev), 0, 0, 1, 0);
    for (int lev = 0; lev < nlev; ++lev)
    {
        amrex::MultiFab& u_mf = *velocity_mf[lev];
        amrex::MultiFab& p_mf = *pressure_mf[lev];
        const Set::Scalar p_floor = pressure_floor;
        const Set::Scalar p_scale = pressure_scale;
        const Set::Scalar p_scale_inv = 1.0 / p_scale;
        const Set::Scalar p_reference = pressure_reference;
        const bool update_pressure = projection_update_pressure;

        pressure_poisson.ApplyCorrection(lev, u_mf, dt);

        for (amrex::MFIter mfi(u_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            Set::Patch<Set::Scalar> p = p_mf.array(mfi);
            Set::Patch<const Set::Scalar> phi = pressure_poisson.Solution(lev).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                if (update_pressure)
                    p(i,j,k) = Util::Max(p_reference + p_scale_inv * phi(i,j,k), p_floor);
            });
        }

        velocity_bc->define(geom[lev]);
        velocity_bc->FillBoundary(u_mf, 0, AMREX_SPACEDIM, time, 0);
        u_mf.FillBoundary(geom[lev].periodicity());
        pressure_bc->define(geom[lev]);
        pressure_bc->FillBoundary(p_mf, 0, 1, time, 0);
        p_mf.FillBoundary(geom[lev].periodicity());
    }

    for (int lev = nlev - 2; lev >= 0; --lev)
        amrex::average_down(*velocity_mf[lev + 1], *velocity_mf[lev],
                            geom[lev + 1], geom[lev], 0, AMREX_SPACEDIM, refRatio(lev));
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
    if (deformable_solid)
    {
        xi_ic->Initialize(lev, xi_mf, 0.0);
        xi_ic->Initialize(lev, xi_old_mf, 0.0);
    }
    pressure_ic->Initialize(lev, pressure_mf, 0.0);
    pressure_correction_mf[lev]->setVal(0.0);
    if (lev == 0 && !(pressure_reference == pressure_reference))
        pressure_reference = pressure_mf[0]->sum(0, false) /
                            static_cast<Set::Scalar>(geom[0].Domain().numPts());

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
LowMach::RHS(int lev, Set::Scalar /*time*/, Set::Scalar dt,
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
    if (deformable_solid && finite_solid_deviatoric_stress_divergence_sign != 0.0)
        UpdateSolidStress(lev, u_mf, *eta_mf[lev], *xi_mf);

    const Set::Scalar* DX = geom[lev].CellSize();
    amrex::Box domain = geom[lev].Domain();

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
            const Model::Gas::MoleFraction X =
                gas.MoleFractions(component_density, i, j, k);
            auto sten = Numeric::GetStencil(i, j, k, domain);
            Set::Scalar density = Util::Max(rho(i,j,k), density_floor);
            Set::Scalar mu = include_viscosity && !implicit_momentum_diffusion ?
                gas.dynamic_viscosity(T(i,j,k), X, i, j, k) : 0.0;

            Set::Vector div_sigma = Set::Vector::Zero();
            if (deformable_solid)
                div_sigma = Numeric::Divergence(solid_deviatoric_stress, i, j, k, DX, sten);

            Set::Vector adv_u = advect.Vector(
                u, u, i, j, k, 0, DX, {Numeric::Advect::Form::Advective}, sten);
            Set::Vector lap_u =
                Numeric::VectorLaplacian(u, i, j, k, 0, DX);

            // Pressure is applied once by the nonincremental projection after
            // this predictor; including the stored pressure here feeds the
            // projection solution back into the next time step.
            Set::Vector rhs_vec = adv_u
                                + g
                                + (mu / density) * lap_u
                                + (finite_solid_deviatoric_stress_divergence_sign / density) * div_sigma;

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
        const Set::Scalar p_reference = pressure_reference;
        for (amrex::MFIter mfi(component_density_rhs_mf, false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            Set::Patch<const Set::Scalar> T = T_mf.array(mfi);
            Set::Patch<const Set::Scalar> component_density = component_density_mf.array(mfi);
            Set::Patch<const Set::Scalar> rigid_eta = rigid_eta_mf.Patch(lev,mfi);
            Set::Patch<Set::Scalar> component_density_rhs = component_density_rhs_mf.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                const Model::Mechanism::State state = {
                    component_density, rigid_eta, T(i,j,k), p_reference};
                mechanism.Apply(component_density_rhs, state, i, j, k, DX);
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
        Set::Patch<const Set::Scalar> xi;
        if (deformable_solid)
            xi = xi_mf->array(mfi);
        Set::Patch<Set::Scalar> T_rhs = T_rhs_mf.array(mfi);
        Set::Patch<Set::Scalar> component_density_rhs = component_density_rhs_mf.array(mfi);
        Set::Patch<Set::Scalar> xi_rhs;
        if (deformable_solid)
            xi_rhs = xi_rhs_mf->array(mfi);
        const int ngas = ngas_species;
        const int solid = deformable_solid_species;
        const bool advect_T = advect_temperature;

        amrex::ParallelFor(bx, [=,this] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            auto sten = Numeric::GetStencil(i, j, k, domain);
            Set::Vector vel(AMREX_D_DECL(u(i,j,k,0),u(i,j,k,1),u(i,j,k,2)));

            auto [species,temperature,dilatation] =
                ComputeThermochemicalSource(component_density, T, i, j, k, DX, dt, !split_chemistry);
            (void)dilatation; // ignore unused

            T_rhs(i,j,k) = 0.0;
            if (advect_T)
            {
                auto [gas_volume_fraction, gas_heat_capacity, heat_capacity,
                    conductivity, cp] = ComputeThermalState(
                        component_density, T(i,j,k), i, j, k);
                (void)gas_volume_fraction;
                (void)conductivity;
                (void)cp;
                Set::Vector grad_T = Numeric::Gradient(T, i, j, k, 0, DX, sten);
                if (heat_capacity > 0.0)
                    T_rhs(i,j,k) -=
                        gas_heat_capacity / heat_capacity * vel.dot(grad_T);
            }
            T_rhs(i,j,k) += temperature;

            for (int n = 0; n < ngas; ++n)
            {
                const Set::Scalar mechanism_source = component_density_rhs(i,j,k,n);
                component_density_rhs(i,j,k,n) =
                    advect(component_density, u, i, j, k, n, DX,
                            {Numeric::Advect::Form::Conservative}, sten) + mechanism_source + species[n];
            }
            if (deformable_solid)
            {
                const Set::Scalar mechanism_source = component_density_rhs(i,j,k,solid);
                component_density_rhs(i,j,k,solid) =
                    advect(component_density, u, i, j, k, solid, DX,
                            {Numeric::Advect::Form::Conservative}, sten) + mechanism_source;
            }
            if (deformable_solid)
                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    xi_rhs(i,j,k,d) = advect(xi, u, i, j, k, d, DX,{Numeric::Advect::Form::Advective}, sten);
        });
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
LowMach::TimeStepBegin(Set::Scalar /*time*/, int /*iter*/)
{
    if (!dynamictimestep.on) return;

    Set::Scalar advmax = 0.0;
    Set::Scalar viscmax = 0.0;
    Set::Scalar elasticmax = 0.0;
    Set::Scalar phasefieldmax = 0.0;
    const bool deformable_solid = deformable_solid_species >= 0;
    const bool explicit_solid_deviatoric_stress = deformable_solid &&
        finite_solid_deviatoric_stress_divergence_sign != 0.0;
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        UpdateComponentState(lev, *component_density_mf[lev]);
        const Set::Scalar* DX = geom[lev].CellSize();
        Set::Scalar dxmin = std::min(DX[0], DX[1]);
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
                Set::Scalar speed = std::sqrt(u(i,j,k,0)*u(i,j,k,0) + u(i,j,k,1)*u(i,j,k,1));
                return {speed / dxmin};
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

            amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpMax,
                            amrex::ReduceOpMax> reduce_op;
            amrex::ReduceData<Set::Scalar, Set::Scalar, Set::Scalar> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;
            reduce_op.eval(bx, reduce_data, [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ReduceTuple
            {
                const Model::Gas::MoleFraction X = gas.MoleFractions(component_density, i, j, k);
                Set::Scalar nu = 0.0;
                if (viscous)
                {
                    Set::Scalar mu = gas.dynamic_viscosity(T(i,j,k), X, i, j, k);
                    nu = mu / Util::Max(rho(i,j,k), rho_floor);
                }
                if (conductive)
                {
                    const Set::Scalar cp = gas.cp_mass(T(i,j,k), X, i, j, k);
                    const Set::Scalar kappa =
                        gas.thermal_conductivity(T(i,j,k), X, i, j, k);
                    nu = Util::Max(nu, kappa /
                        (Util::Max(rho(i,j,k), rho_floor) * cp));
                }
                if (species_diffusive)
                    for (int n = 0; n < ngas; ++n)
                        nu = Util::Max(nu, gas.diffusion_coefficient(
                            T(i,j,k), p_reference, X, i, j, k, n));
                Set::Scalar elastic_rate = 0.0;
                if (elastic && Util::Clamp(eta(i,j,k), 0.0, 1.0) > eta_threshold)
                {
                    Set::Scalar density = Util::Max(rho(i,j,k), rho_floor);
                    Set::Scalar wave_speed = std::sqrt(Util::Max(kappa_solid + (4.0 / 3.0) * mu_solid, mu_solid) / density);
                    elastic_rate = wave_speed / dxmin;
                    nu = Util::Max(nu, (solid_viscosity + solid_interface_viscosity) / density);
                }
                return {nu / (dxmin * dxmin), elastic_rate, T(i,j,k)};
            });
            ReduceTuple hv = reduce_data.value();
            viscmax = std::max(viscmax, amrex::get<0>(hv));
            elasticmax = std::max(elasticmax, amrex::get<1>(hv));
            temperaturemax = std::max(temperaturemax, amrex::get<2>(hv));
        }
        for (const auto& mechanism : mechanisms)
            phasefieldmax = std::max(
                phasefieldmax, mechanism.StabilityRate(dxmin, temperaturemax));
    }
    amrex::ParallelDescriptor::ReduceRealMax(advmax);
    amrex::ParallelDescriptor::ReduceRealMax(viscmax);
    amrex::ParallelDescriptor::ReduceRealMax(elasticmax);
    amrex::ParallelDescriptor::ReduceRealMax(phasefieldmax);

    Set::Scalar adv_dt = cfl_v;
    if (advmax > 0.0) adv_dt = cfl / advmax;
    Set::Scalar visc_dt = viscmax > 0.0 ? 0.5 * cfl / viscmax : cfl_v;
    Set::Scalar elastic_dt = elasticmax > 0.0 ? cfl / elasticmax : cfl_v;
    Set::Scalar phasefield_dt = phasefieldmax > 0.0 ? cfl / phasefieldmax : cfl_v;
    DynamicTimestep_SyncTimeStep(0, std::min({adv_dt, visc_dt, elastic_dt, phasefield_dt}));
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
        const Set::Scalar* DX = geom[lev].CellSize();
        amrex::Box domain = geom[lev].Domain();

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
                auto sten = Numeric::GetStencil(i, j, k, domain);
                Set::Matrix grad_u = Numeric::Gradient(u, i, j, k, DX, sten);
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
                        << " divrms_interior " << divrms_interior << "\n";
}

void
LowMach::TimeStepComplete(Set::Scalar time, int iter)
{
    ApplyImplicitDiffusion(time + dt[0], dt[0]);
    ProjectVelocity(time + dt[0], dt[0]);
    PrintDiagnostics(time, iter);
}

void
LowMach::PreparePlotFile(Set::Scalar /*time*/, const amrex::Vector<int>& /*iter*/)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        UpdateComponentState(lev, *component_density_mf[lev]);
        if (deformable_solid_species >= 0)
            UpdateSolidStress(lev, *velocity_mf[lev], *eta_mf[lev], *xi_mf[lev], diagnostics_extended_fields);
        UpdateDerivedDiagnostics(lev, *velocity_mf[lev], *temperature_mf[lev]);
    }
}

void
LowMach::TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/)
{
    const Set::Scalar* DX = geom[lev].CellSize();
    Set::Scalar dr = std::sqrt(DX[0] * DX[0] + DX[1] * DX[1]);
    const bool deformable_solid = deformable_solid_species >= 0;
    const bool rigid_solid = !rigid_solid_species.empty();

    for (amrex::MFIter mfi(*temperature_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        Set::Patch<char> tag = tags.array(mfi);
        Set::Patch<const Set::Scalar> u = velocity_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> pressure = pressure_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> T = temperature_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> eta = eta_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> rigid_eta = rigid_eta_mf.Patch(lev,mfi);
        const Set::Scalar vcrit = velocity_refinement_criterion;
        const Set::Scalar pcrit = pressure_refinement_criterion;
        const Set::Scalar Tcrit = temperature_refinement_criterion;
        const Set::Scalar etacrit = eta_refinement_criterion;
        amrex::Box domain = geom[lev].Domain();

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            auto sten = Numeric::GetStencil(i, j, k, domain);
            Set::Matrix grad_u = Numeric::Gradient(u, i, j, k, DX, sten);
            Set::Vector grad_p = Numeric::Gradient(pressure, i, j, k, 0, DX, sten);
            Set::Vector grad_T = Numeric::Gradient(T, i, j, k, 0, DX, sten);
            bool refine_eta = false;
            if (deformable_solid)
            {
                Set::Vector grad_eta = Numeric::Gradient(eta, i, j, k, 0, DX, sten);
                const Set::Scalar eta_val = eta(i,j,k);
                refine_eta = grad_eta.lpNorm<2>() * dr * 2.0 > etacrit ||
                            (eta_val > etacrit && eta_val < 1.0 - etacrit);
            }
            if (rigid_solid)
            {
                Set::Vector grad_eta = Numeric::Gradient(rigid_eta, i, j, k, 0, DX, sten);
                const Set::Scalar eta_val = rigid_eta(i,j,k);
                refine_eta = refine_eta || grad_eta.lpNorm<2>() * dr * 2.0 > etacrit ||
                            (eta_val > etacrit && eta_val < 1.0 - etacrit);
            }
            if (grad_u.norm() * dr > vcrit ||
                grad_p.lpNorm<2>() * dr > pcrit ||
                grad_T.lpNorm<2>() * dr > Tcrit ||
                refine_eta)
                tag(i,j,k) = amrex::TagBox::SET;
        });
    }
}

}
