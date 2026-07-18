#include "LowMach.H"

#include "AMReX_MultiFabUtil.H"
#include "AMReX_TimeIntegrator.H"
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
    pp.query_default("include_conduction", value.include_conduction, true);
    pp.query_default("advect_temperature", value.advect_temperature, true);

    pp.queryclass<Model::Gas::Gas>("gas", value.gas);
    value.ngas_species = value.gas.nspecies;

    value.species.Parse(pp);
    value.nspecies = value.species.Size();
    value.reference_density.assign(value.nspecies, NAN);
    for (int n = 0; n < value.nspecies; ++n)
    {
        const Model::Species::Mechanics mechanics = value.species.GetMechanics(n);
        if (mechanics == Model::Species::Mechanics::Fluid)
        {
            if (n >= value.ngas_species)
                Util::Exception(INFO, "Fluid species ", value.species.Name(n),
                                " has no corresponding entry in gas.mw");
        }
        else if (mechanics == Model::Species::Mechanics::DeformableSolid)
        {
            if (value.deformable_solid_species >= 0)
                Util::Exception(INFO,
                                "LowMach currently supports one deformable solid species");
            value.deformable_solid_species = n;
        }
        else if (mechanics == Model::Species::Mechanics::RigidSolid)
            value.rigid_solid_species.push_back(n);
    }
    if (value.nspecies < value.ngas_species)
        Util::Exception(INFO, "species.names must contain one identifier for every gas species");
    for (int n = 0; n < value.ngas_species; ++n)
        if (value.species.GetMechanics(n) != Model::Species::Mechanics::Fluid)
            Util::Exception(INFO, "The first ", value.ngas_species,
                            " species must be fluid species described by the gas model");

    if (value.deformable_solid_species >= 0)
    {
        pp.queryclass("reference_map", value.reference_map_reconstruction);

        std::string solid_model_type;
        pp.query_required("solid.model.type", solid_model_type);
        if (solid_model_type != "finite.neohookean")
            Util::Exception(INFO, solid_model_type,
                            " is not a valid deformable solid model for LowMach");
        pp.query_default("solid.model.eta_threshold", value.finite_solid_eta_threshold, 0.5);
        pp.query_default("solid.model.J_floor", value.finite_solid_J_floor, 1.0e-6);
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
            pp.query_required(value.species.Name(n) + ".reference_density",
                              value.reference_density[n], Unit::Density());
            if (value.reference_density[n] <= 0.0)
                Util::Exception(INFO, value.species.Name(n),
                                ".reference_density must be positive");
        }
    }

    std::vector<std::string> mechanism_names;
    pp.queryarr_default("mechanisms.names", mechanism_names, {});
    value.mechanisms.resize(mechanism_names.size());
    for (int n = 0; n < static_cast<int>(mechanism_names.size()); ++n)
    {
        const std::string& id = mechanism_names[n];
        if (!Model::Species::Registry::IsIdentifier(id))
            Util::Exception(INFO, id, " is not a valid mechanism identifier");
        for (int m = 0; m < n; ++m)
            if (id == mechanism_names[m])
                Util::Exception(INFO, "Duplicate mechanism identifier ", id);
        value.mechanisms[n].Parse(pp, id, value.species, value.reference_density,
                                  value.gas.MW, value.gas.Rg);
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
    pp.select_default<IC::Constant,IC::Expression>("component_density.ic", value.component_density_ic, value.geom);
    pp.select_default<IC::Constant,IC::Expression>("pressure.ic", value.pressure_ic, value.geom);
    if (value.deformable_solid_species >= 0)
    {
        pp.select_default<BC::Constant::ZeroNeumann,BC::Constant,BC::Expression>("xi.bc", value.xi_bc, AMREX_SPACEDIM);
        pp.select_default<IC::Expression::X,IC::Constant,IC::Expression>("xi.ic", value.xi_ic, value.geom);
    }

    std::vector<std::string> species_suffix(value.nspecies);
    for (int n = 0; n < value.nspecies; ++n)
        species_suffix[n] = "_" + value.species.Name(n);
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
    if (value.deformable_solid_species >= 0)
        value.AddField<Set::Matrix,Set::HC::Cell>(value.solid_deviatoric_stress_mf, nullptr, 1, 1, "solid_deviatoric_stress", true, false);
    if (value.diagnostics_extended_fields)
    {
        value.AddField<Set::Scalar,Set::HC::Cell>(value.mass_fraction_mf, &value.bc_nothing, value.ngas_species, 1, "mass_fraction", true, false, gas_species_suffix);
        value.AddField<Set::Scalar,Set::HC::Cell>(value.mole_fraction_mf, &value.bc_nothing, value.ngas_species, 1, "mole_fraction", true, false, gas_species_suffix);
        value.AddField<Set::Scalar,Set::HC::Cell>(value.momentum_mf, &value.bc_nothing, AMREX_SPACEDIM, 1, "momentum", true, false, {"x","y"});
        value.AddField<Set::Scalar,Set::HC::Cell>(value.energy_mf, &value.bc_nothing, 1, 1, "energy", true, false);
        value.AddField<Set::Scalar,Set::HC::Cell>(value.vorticity_mf, &value.bc_nothing, 1, 1, "vorticity", true, false);
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
        amrex::Array4<Set::Matrix> F_field;
        if (write_diagnostics)
            F_field = deformation_gradient_mf[lev]->array(mfi);

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

void
LowMach::UpdateComponentState(int lev, const amrex::MultiFab& component_density_mf)
{
    density_mf[lev]->setVal(0.0);
    const bool deformable_solid = deformable_solid_species >= 0;
    const bool rigid_solid = !rigid_solid_species.empty();
    if (deformable_solid)
        eta_mf[lev]->setVal(0.0);
    if (rigid_solid)
        rigid_eta_mf[lev]->setVal(0.0);
    if (diagnostics_extended_fields)
    {
        mass_fraction_mf[lev]->setVal(0.0);
        mole_fraction_mf[lev]->setVal(0.0);
    }

    for (amrex::MFIter mfi(component_density_mf, true); mfi.isValid(); ++mfi)
    {
        amrex::Box bx = mfi.growntilebox(1);
        Set::Patch<const Set::Scalar> component_density = component_density_mf.array(mfi);
        Set::Patch<Set::Scalar> rho = density_mf.Patch(lev,mfi);
        amrex::Array4<Set::Scalar> eta;
        if (deformable_solid)
            eta = eta_mf[lev]->array(mfi);
        amrex::Array4<Set::Scalar> mass_fraction;
        amrex::Array4<Set::Scalar> mole_fraction;
        if (diagnostics_extended_fields)
        {
            mass_fraction = mass_fraction_mf[lev]->array(mfi);
            mole_fraction = mole_fraction_mf[lev]->array(mfi);
        }
        const int ngas = ngas_species;
        const int nsp = nspecies;
        const int solid = deformable_solid_species;
        const bool write_composition = diagnostics_extended_fields;
        const Set::Scalar solid_density = finite_solid_reference_density;
        const Set::Scalar composition_floor = small;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Scalar gas_partial_density = 0.0;
            Set::Scalar total_partial_density = 0.0;
            for (int n = 0; n < ngas; ++n)
                gas_partial_density += component_density(i,j,k,n);
            for (int n = 0; n < nsp; ++n)
                total_partial_density += component_density(i,j,k,n);

            // Partial densities are the conserved species state. Mechanical
            // volume fractions are reconstructed from them when needed.
            rho(i,j,k) = total_partial_density;
            if (deformable_solid)
                eta(i,j,k) = component_density(i,j,k,solid) / solid_density;

            if (write_composition)
            {
                const Model::Gas::MoleFraction X = gas.MoleFractions(component_density, i, j, k);
                for (int n = 0; n < ngas; ++n)
                {
                    mass_fraction(i,j,k,n) = gas_partial_density > composition_floor ?
                        component_density(i,j,k,n) / gas_partial_density : (n == 0 ? 1.0 : 0.0);
                    mole_fraction(i,j,k,n) = X(i,j,k,n);
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

void
LowMach::UpdateDerivedDiagnostics(int lev, const amrex::MultiFab& u_mf, const amrex::MultiFab& T_mf)
{
    if (!diagnostics_extended_fields) return;

    const Set::Scalar* DX = geom[lev].CellSize();
    amrex::Box domain = geom[lev].Domain();
    momentum_mf[lev]->setVal(0.0);
    energy_mf[lev]->setVal(0.0);
    vorticity_mf[lev]->setVal(0.0);

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

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            const Model::Gas::MoleFraction X = gas.MoleFractions(component_density, i, j, k);
            const Set::Scalar density = rho(i,j,k);
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
                M(i,j,k,d) = density * u(i,j,k,d);
            E(i,j,k) = gas.ComputeE(density, M(i,j,k,0), M(i,j,k,1), T(i,j,k), X, i, j, k);

            auto sten = Numeric::GetStencil(i, j, k, domain);
            Set::Matrix grad_u = Numeric::Gradient(u, i, j, k, DX, sten);
            omega(i,j,k) = grad_u(1,0) - grad_u(0,1);
        });
    }

    momentum_mf[lev]->FillBoundary(geom[lev].periodicity());
    energy_mf[lev]->FillBoundary(geom[lev].periodicity());
    vorticity_mf[lev]->FillBoundary(geom[lev].periodicity());
}


void
LowMach::ProjectVelocity(Set::Scalar time, Set::Scalar dt)
{
    BL_PROFILE("Integrator::LowMach::ProjectVelocity");
    if (!projection_enabled || !(dt > 0.0)) return;

    const int nlev = finest_level + 1;
    const bool rigid_solid = !rigid_solid_species.empty();
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
                    const Set::Scalar weight = Util::Clamp(eta(i,j,k), 0.0, 1.0);
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
        amrex::Box domain = geom[lev].Domain();
        amrex::MultiFab& u_mf = *velocity_mf[lev];
        amrex::MultiFab& T_mf = *temperature_mf[lev];
        amrex::MultiFab& component_density_mf = *this->component_density_mf[lev];
        amrex::MultiFab& rho_mf = *density_mf[lev];
        amrex::MultiFab& rhs_mf = pressure_poisson.RHS(lev);
        amrex::MultiFab& beta_mf = pressure_poisson.Coefficient(lev);
        const Set::Scalar inv_dt = 1.0 / dt;
        const Set::Scalar rho_floor = density_floor;
        const Set::Scalar inverse_relaxation_time = rigid_solid ?
            1.0 / rigid_relaxation_time : 0.0;

        rhs_mf.setVal(0.0);
        for (const auto& configured_mechanism : mechanisms)
        {
            const auto mechanism = configured_mechanism;
            const Set::Scalar p_reference = pressure_reference;
            for (amrex::MFIter mfi(rhs_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.tilebox();
                Set::Patch<const Set::Scalar> T = T_mf.array(mfi);
                Set::Patch<const Set::Scalar> component_density = component_density_mf.array(mfi);
                amrex::Array4<const Set::Scalar> rigid_eta;
                if (rigid_solid)
                    rigid_eta = rigid_eta_mf[lev]->array(mfi);
                Set::Patch<Set::Scalar> rhs = rhs_mf.array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    const Model::Mechanism::State state = {
                        component_density, rigid_eta, T(i,j,k), p_reference};
                    rhs(i,j,k) += mechanism.VolumeSource(state, i, j, k, DX);
                });
            }
        }

        for (amrex::MFIter mfi(rhs_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            Set::Patch<const Set::Scalar> u = u_mf.array(mfi);
            Set::Patch<const Set::Scalar> rho = rho_mf.array(mfi);
            amrex::Array4<const Set::Scalar> rigid_eta;
            if (rigid_solid)
                rigid_eta = rigid_eta_mf[lev]->array(mfi);
            Set::Patch<Set::Scalar> rhs = rhs_mf.array(mfi);
            Set::Patch<Set::Scalar> beta = beta_mf.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                auto sten = Numeric::GetStencil(i, j, k, domain);
                Set::Matrix grad_u = Numeric::Gradient(u, i, j, k, DX, sten);
                Set::Scalar div_u = 0.0;
                for (int d = 0; d < AMREX_SPACEDIM; ++d) div_u += grad_u(d,d);
                const Set::Scalar weight = rigid_solid ?
                    Util::Clamp(rigid_eta(i,j,k), 0.0, 1.0) : 0.0;
                const Set::Scalar mobility = 1.0 / (1.0 + dt * weight * inverse_relaxation_time);
                rhs(i,j,k) = (div_u - rhs(i,j,k)) * inv_dt;
                beta(i,j,k) = mobility / Util::Max(rho(i,j,k), rho_floor);
            });
        }
    }

    pressure_poisson.Solve(time, pressure_bc->GetBCRec());

    for (int lev = 0; lev < nlev; ++lev)
        amrex::MultiFab::Copy(*pressure_correction_mf[lev],
                              pressure_poisson.Solution(lev), 0, 0, 1, 0);
    for (int lev = 0; lev < nlev; ++lev)
    {
        amrex::MultiFab& u_mf = *velocity_mf[lev];
        amrex::MultiFab& p_mf = *pressure_mf[lev];
        amrex::MultiFab& beta_mf = pressure_poisson.Coefficient(lev);
        const Set::Scalar* DX = geom[lev].CellSize();
        amrex::Box domain = geom[lev].Domain();
        const Set::Scalar p_floor = pressure_floor;
        const Set::Scalar p_scale = pressure_scale;
        const Set::Scalar p_scale_inv = 1.0 / p_scale;
        const Set::Scalar p_reference = pressure_reference;
        const bool update_pressure = projection_update_pressure;

        for (amrex::MFIter mfi(u_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            Set::Patch<Set::Scalar> u = u_mf.array(mfi);
            Set::Patch<Set::Scalar> p = p_mf.array(mfi);
            Set::Patch<const Set::Scalar> beta = beta_mf.array(mfi);
            Set::Patch<const Set::Scalar> phi = pressure_poisson.Solution(lev).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                auto sten = Numeric::GetStencil(i, j, k, domain);
                Set::Vector grad_phi = Numeric::Gradient(phi, i, j, k, 0, DX, sten);
                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    u(i,j,k,d) -= dt * beta(i,j,k) * grad_phi(d);
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
    component_density_ic->Initialize(lev, component_density_mf, 0.0);
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
LowMach::RHS(int lev, Set::Scalar /*time*/,
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
    const bool rigid_solid = !rigid_solid_species.empty();
    if (deformable_solid && finite_solid_deviatoric_stress_divergence_sign != 0.0)
        UpdateSolidStress(lev, u_mf, *eta_mf[lev], *xi_mf);

    const Set::Scalar* DX = geom[lev].CellSize();
    amrex::Box domain = geom[lev].Domain();
    const Numeric::Advect::Options advective_options{Numeric::Advect::Form::Advective};

    u_rhs_mf.setVal(0.0, 0, AMREX_SPACEDIM, u_rhs_mf.nGrow());

    for (amrex::MFIter mfi(u_rhs_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        Set::Patch<const Set::Scalar> u = u_mf.array(mfi);
        Set::Patch<const Set::Scalar> T = T_mf.array(mfi);
        Set::Patch<const Set::Scalar> component_density = component_density_mf.array(mfi);
        Set::Patch<const Set::Scalar> rho = density_mf.Patch(lev,mfi);
        amrex::Array4<const Set::Matrix> solid_deviatoric_stress;
        if (deformable_solid)
            solid_deviatoric_stress = solid_deviatoric_stress_mf[lev]->array(mfi);
        Set::Patch<Set::Scalar> u_rhs = u_rhs_mf.array(mfi);
        const bool viscous = include_viscosity;
        const Set::Vector gravity = g;
        const Set::Scalar rho_floor = density_floor;
        const Set::Scalar deviatoric_stress_divergence_sign = finite_solid_deviatoric_stress_divergence_sign;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            const Model::Gas::MoleFraction X = gas.MoleFractions(component_density, i, j, k);
            auto sten = Numeric::GetStencil(i, j, k, domain);
            Set::Scalar density = Util::Max(rho(i,j,k), rho_floor);

            Set::Scalar mu = viscous ? gas.dynamic_viscosity(T(i,j,k), X, i, j, k) : 0.0;

            Set::Vector div_sigma = Set::Vector::Zero();
            if (deformable_solid)
                div_sigma = Numeric::Divergence(solid_deviatoric_stress, i, j, k, DX, sten);

            Set::Vector adv_u = advect.Vector(u, u, i, j, k, 0, DX, advective_options, sten);

            Set::Vector lap_u = Numeric::VectorLaplacian(u, i, j, k, 0, DX);

            // Pressure is applied once by the nonincremental projection after
            // this predictor; including the stored pressure here feeds the
            // projection solution back into the next time step.
            Set::Vector rhs_vec = adv_u
                                  + gravity
                                  + (mu / density) * lap_u
                                  + (deviatoric_stress_divergence_sign / density) * div_sigma;

            // Put the result into the multicomponent field
            for (int d = 0; d < AMREX_SPACEDIM; ++d)
                u_rhs(i,j,k,d) = rhs_vec(d);
        });
    }

    // Reuse the temperature RHS as temporary storage for the dilatation that
    // the phase-change mechanisms prescribe for the projection.
    T_rhs_mf.setVal(0.0);
    component_density_rhs_mf.setVal(0.0);
    for (const auto& configured_mechanism : mechanisms)
    {
        const auto mechanism = configured_mechanism;
        const Set::Scalar p_reference = pressure_reference;
        for (amrex::MFIter mfi(component_density_rhs_mf, false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            Set::Patch<const Set::Scalar> T = T_mf.array(mfi);
            Set::Patch<const Set::Scalar> component_density = component_density_mf.array(mfi);
            amrex::Array4<const Set::Scalar> rigid_eta;
            if (rigid_solid)
                rigid_eta = rigid_eta_mf[lev]->array(mfi);
            Set::Patch<Set::Scalar> T_rhs = T_rhs_mf.array(mfi);
            Set::Patch<Set::Scalar> component_density_rhs = component_density_rhs_mf.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                const Model::Mechanism::State state = {
                    component_density, rigid_eta, T(i,j,k), p_reference};
                mechanism.Apply(component_density_rhs, state, i, j, k, DX);
                T_rhs(i,j,k) += mechanism.VolumeSource(state, i, j, k, DX);
            });
        }
    }

    for (amrex::MFIter mfi(T_mf, false); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        Set::Patch<const Set::Scalar> u = u_mf.array(mfi);
        Set::Patch<const Set::Scalar> T = T_mf.array(mfi);
        Set::Patch<const Set::Scalar> component_density = component_density_mf.array(mfi);
        amrex::Array4<const Set::Scalar> xi;
        if (deformable_solid)
            xi = xi_mf->array(mfi);
        Set::Patch<const Set::Scalar> rho = density_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> T_rhs = T_rhs_mf.array(mfi);
        Set::Patch<Set::Scalar> component_density_rhs = component_density_rhs_mf.array(mfi);
        amrex::Array4<Set::Scalar> xi_rhs;
        if (deformable_solid)
            xi_rhs = xi_rhs_mf->array(mfi);
        const int ngas = ngas_species;
        const int solid = deformable_solid_species;
        const bool conductive = include_conduction;
        const bool advect_T = advect_temperature;
        const Set::Scalar rho_floor = density_floor;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            const Model::Gas::MoleFraction X = gas.MoleFractions(component_density, i, j, k);
            auto sten = Numeric::GetStencil(i, j, k, domain);
            Set::Vector vel = Set::Vector::Zero();
            vel(0) = u(i,j,k,0);
            vel(1) = u(i,j,k,1);
#if AMREX_SPACEDIM == 3
            vel(2) = u(i,j,k,2);
#endif
            Set::Scalar density = Util::Max(rho(i,j,k), rho_floor);
            const Set::Scalar volume_source = T_rhs(i,j,k);
            T_rhs(i,j,k) = 0.0;
            if (advect_T)
            {
                Set::Vector grad_T = Numeric::Gradient(T, i, j, k, 0, DX, sten);
                T_rhs(i,j,k) -= vel.dot(grad_T);
            }
            if (conductive)
            {
                Set::Scalar kappa = gas.thermal_conductivity(T(i,j,k), X, i, j, k);
                Set::Scalar cp = gas.cp_mass(T(i,j,k), X, i, j, k);
                if (kappa == kappa && cp == cp && cp > 0.0)
                    T_rhs(i,j,k) += kappa / (density * cp) * Numeric::Laplacian(T, i, j, k, 0, DX);
            }

            for (int n = 0; n < ngas; ++n)
            {
                const Set::Scalar mechanism_source = component_density_rhs(i,j,k,n);
                // The projection imposes div(u)=volume_source, so this is the
                // material form of partial-density conservation.
                component_density_rhs(i,j,k,n) =
                    advect(component_density, u, i, j, k, n, DX,
                           advective_options, sten) + mechanism_source -
                    component_density(i,j,k,n) * volume_source;
            }
            if (deformable_solid)
            {
                const Set::Scalar mechanism_source = component_density_rhs(i,j,k,solid);
                component_density_rhs(i,j,k,solid) =
                    advect(component_density, u, i, j, k, solid, DX,
                           advective_options, sten) + mechanism_source -
                    component_density(i,j,k,solid) * volume_source;
            }
            if (deformable_solid)
                for (int d = 0; d < AMREX_SPACEDIM; ++d)
                    xi_rhs(i,j,k,d) = advect(xi, u, i, j, k, d, DX, advective_options, sten);
        });
    }

    u_rhs_mf.FillBoundary(geom[lev].periodicity());
}

void
LowMach::Advance(int lev, Set::Scalar time, Set::Scalar dt)
{
    const bool deformable_solid = deformable_solid_species >= 0;
    std::swap(velocity_old_mf[lev], velocity_mf[lev]);
    std::swap(temperature_old_mf[lev], temperature_mf[lev]);
    std::swap(component_density_old_mf[lev], component_density_mf[lev]);
    if (deformable_solid)
        std::swap(xi_old_mf[lev], xi_mf[lev]);

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
        RHS(lev, rhs_time, rhs_mf[0], rhs_mf[1], rhs_mf[2],
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
        for (const auto& mechanism : mechanisms)
            phasefieldmax = std::max(phasefieldmax, mechanism.StabilityRate(dxmin));

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
            amrex::Array4<const Set::Scalar> eta;
            if (deformable_solid)
                eta = eta_mf[lev]->array(mfi);
            const Set::Scalar rho_floor = density_floor;
            const bool viscous = include_viscosity;
            const bool elastic = explicit_solid_deviatoric_stress;
            const Set::Scalar eta_threshold = Util::Clamp(finite_solid_eta_threshold, 0.0, 1.0);
            const Set::Scalar mu_solid = finite_solid_model.mu;
            const Set::Scalar kappa_solid = finite_solid_model.kappa;
            const Set::Scalar solid_viscosity = finite_solid_viscosity;
            const Set::Scalar solid_interface_viscosity = finite_solid_interface_viscosity;

            amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpMax> reduce_op;
            amrex::ReduceData<Set::Scalar, Set::Scalar> reduce_data(reduce_op);
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
                Set::Scalar elastic_rate = 0.0;
                if (elastic && Util::Clamp(eta(i,j,k), 0.0, 1.0) > eta_threshold)
                {
                    Set::Scalar density = Util::Max(rho(i,j,k), rho_floor);
                    Set::Scalar wave_speed = std::sqrt(Util::Max(kappa_solid + (4.0 / 3.0) * mu_solid, mu_solid) / density);
                    elastic_rate = wave_speed / dxmin;
                    nu = Util::Max(nu, (solid_viscosity + solid_interface_viscosity) / density);
                }
                return {nu / (dxmin * dxmin), elastic_rate};
            });
            ReduceTuple hv = reduce_data.value();
            viscmax = std::max(viscmax, amrex::get<0>(hv));
            elasticmax = std::max(elasticmax, amrex::get<1>(hv));
        }
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
        amrex::Array4<char> const& tag = tags.array(mfi);
        Set::Patch<const Set::Scalar> u = velocity_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> pressure = pressure_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> T = temperature_mf.Patch(lev,mfi);
        amrex::Array4<const Set::Scalar> eta;
        if (deformable_solid)
            eta = eta_mf[lev]->array(mfi);
        amrex::Array4<const Set::Scalar> rigid_eta;
        if (rigid_solid)
            rigid_eta = rigid_eta_mf[lev]->array(mfi);
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
