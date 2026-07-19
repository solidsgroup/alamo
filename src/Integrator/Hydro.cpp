
#include "Hydro.H"
#include "AMReX_MultiFab.H"
#include "IO/ParmParse.H"
#include "BC/Constant.H"
#include "BC/Expression.H"
#include "Numeric/Stencil.H"
#include "IC/Constant.H"
#include "IC/Laminate.H"
#include "IC/Expression.H"
#include "IC/BMP.H"
#include "IC/PNG.H"
#include "Solver/Local/Riemann/Roe.H"
#include "Solver/Local/Riemann/HLLE.H"
#include "Solver/Local/Riemann/HLLC.H"
#include "AMReX_TimeIntegrator.H"

#include "Model/Gas/Gas.H"
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
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar hydro_abs(Set::Scalar a)
{
    return a < 0.0 ? -a : a;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar hydro_min(Set::Scalar a, Set::Scalar b)
{
    return a < b ? a : b;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar hydro_max(Set::Scalar a, Set::Scalar b)
{
    return a > b ? a : b;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar hydro_minmod(Set::Scalar a, Set::Scalar b)
{
    if (a * b <= 0.0) return 0.0;
    Set::Scalar mag = hydro_min(hydro_abs(a), hydro_abs(b));
    return a < 0.0 ? -mag : mag;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar hydro_mc_slope(Set::Scalar lo, Set::Scalar center, Set::Scalar hi)
{
    Set::Scalar dl = center - lo;
    Set::Scalar dr = hi - center;
    return hydro_minmod(0.5 * (dl + dr), hydro_minmod(2.0 * dl, 2.0 * dr));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar hydro_clamp(Set::Scalar value, Set::Scalar lo, Set::Scalar hi)
{
    return hydro_min(hydro_max(value, hydro_min(lo, hi)), hydro_max(lo, hi));
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool hydro_is_finite(Set::Scalar value)
{
    constexpr Set::Scalar max_finite = 1.0e300;
    return value == value && value < max_finite && value > -max_finite;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar hydro_finite_or(Set::Scalar value, Set::Scalar fallback)
{
    return hydro_is_finite(value) ? value : fallback;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Scalar hydro_smootherstep(Set::Scalar value)
{
    Set::Scalar x = hydro_clamp(value, 0.0, 1.0);
    return x * x * x * (x * (x * 6.0 - 15.0) + 10.0);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Set::Matrix hydro_regularize_stress(Set::Matrix sigma, Set::Scalar cap)
{
    Set::Scalar norm2 = 0.0;
    for (int a = 0; a < AMREX_SPACEDIM; ++a)
        for (int b = 0; b < AMREX_SPACEDIM; ++b)
        {
            if (!hydro_is_finite(sigma(a,b))) return Set::Matrix::Zero();
            norm2 += sigma(a,b) * sigma(a,b);
        }

    if (cap > 0.0 && hydro_is_finite(cap) && norm2 > cap * cap)
    {
        Set::Scalar norm = std::sqrt(norm2);
        if (norm > 0.0 && hydro_is_finite(norm)) sigma *= cap / norm;
        else return Set::Matrix::Zero();
    }
    return sigma;
}
}

Hydro::Hydro(IO::ParmParse& pp) : Hydro()
{
    pp_queryclass(*this);
}

void
Hydro::Parse(Hydro& value, IO::ParmParse& pp)
{
    BL_PROFILE("Integrator::Hydro::Hydro()");
    {
        // pp.query_default("r_refinement_criterion",     value.r_refinement_criterion    , 0.01);
        // energy-based refinement
        // pp.query_default("e_refinement_criterion",     value.e_refinement_criterion    , 0.01);
        // momentum-based refinement
        // pp.query_default("m_refinement_criterion",     value.m_refinement_criterion    , 0.01);

        pp.forbid("scheme","use integration.type instead");

        // eta-based refinement
        pp.query_default("eta_refinement_criterion",   value.eta_refinement_criterion  , 0.01);
        // vorticity-based refinement
        pp.query_default("omega_refinement_criterion", value.omega_refinement_criterion, 0.01);
        // velocity gradient-based refinement
        pp.query_default("gradu_refinement_criterion", value.gradu_refinement_criterion, 0.01);
        // pressure-based refinement
        pp.query_default("p_refinement_criterion", value.p_refinement_criterion, 1e100);
        // density-based refinement
        pp.query_default("rho_refinement_criterion", value.rho_refinement_criterion, 1e100);

        pp_forbid("gamma", "replaced by gas->gamma(...)"); // gamma for gamma law
        pp_query_required("cfl", value.cfl); // cfl condition
        pp_query_default("cfl_v", value.cfl_v,1E100); // cfl condition
        pp_forbid("mu", "replaced with gas->dynamic_viscosity(...)"); // linear viscosity coefficient
        pp_forbid("Lfactor","replaced with mu");
        //pp_query_default("Lfactor", value.Lfactor,1.0); // (to be removed) test factor for viscous source
        pp_forbid("Pfactor","replaced with mu");
        //pp_query_default("Pfactor", value.Pfactor,1.0); // (to be removed) test factor for viscous source
        pp_forbid("pref", "deprecated - use absolute pressure"); // reference pressure for Roe solver

        pp_forbid("rho.bc","--> density.bc");
        pp_forbid("p.bc","--> pressure.bc");
        pp_forbid("v.bc", "--> velocity.bc");
        pp_forbid("pressure.bc","--> energy.bc");
        pp_forbid("velocity.bc","--> momentum.bc");

        // Boundary condition for density
        pp.select_default<BC::Constant,BC::Expression>("density.bc",value.density_bc,1);
        // Boundary condition for energy
        pp.select_default<BC::Constant,BC::Expression>("energy.bc",value.energy_bc,1);
        // Boundary condition for momentum
        pp.select_default<BC::Constant,BC::Expression>("momentum.bc",value.momentum_bc,2);

        if (!value.managed)
        {
            // Boundary condition for phase field order parameter
            pp.select_default<BC::Constant,BC::Expression>("pf.eta.bc",value.eta_bc,1);
        }

        pp_query_default("small",value.small,1E-8); // small regularization value
        pp_query_default("cutoff",value.cutoff,-1E100); // cutoff value
        pp_query_default("lagrange",value.lagrange,0.0); // lagrange no-penetration factor

        std::string eta_mode_str;
        pp.query_validate("eta.mode", eta_mode_str, {"static","evolving"});
        if (eta_mode_str == "static") value.eta_mode = EtaMode::Static;
        else if (eta_mode_str == "evolving") value.eta_mode = EtaMode::Evolving;
        if (value.managed && value.eta_mode == EtaMode::Evolving)
            Util::Exception(INFO,"Hydro eta.mode=evolving is only valid when Hydro owns eta; externally managed eta must use eta.mode=static.");

        pp_forbid("roefix","--> solver.roe.entropy_fix"); // Roe solver entropy fix

    }
    // Register FabFields:
    {
        int nghost = 1;
        int eta_nghost = value.eta_mode == EtaMode::Evolving ? 2 : nghost;

        if (!value.managed)
        {
            value.eta_mf = new Set::Field<Set::Scalar>();
            value.eta_old_mf = new Set::Field<Set::Scalar>();
            value.RegisterNewFab(*value.eta_mf,     value.eta_bc, 1, eta_nghost, "eta",     true, true);
            value.RegisterNewFab(*value.eta_old_mf, value.eta_bc, 1, eta_nghost, "eta_old", true, true);
        }
        value.RegisterNewFab(value.etadot_mf,  value.eta_bc, 1, nghost, "etadot",  true, false);

        value.RegisterNewFab(value.density_mf,     value.density_bc, 1, nghost, "density",     true , true);
        value.RegisterNewFab(value.density_old_mf, value.density_bc, 1, nghost, "density_old", false, true);

        value.RegisterNewFab(value.energy_mf,     value.energy_bc, 1, nghost, "energy",      true ,true);
        value.RegisterNewFab(value.energy_old_mf, value.energy_bc, 1, nghost, "energy_old" , false, true);

        value.RegisterNewFab(value.momentum_mf,     value.momentum_bc, 2, nghost, "momentum",     true ,true, {"x","y"});
        value.RegisterNewFab(value.momentum_old_mf, value.momentum_bc, 2, nghost, "momentum_old", false, true);

        value.RegisterNewFab(value.xi_mf,     &value.neumann_bc_D, 2, eta_nghost, "xi",     true, true, {"x","y"});
        value.RegisterNewFab(value.xi_old_mf, &value.neumann_bc_D, 2, eta_nghost, "xi_old", false, true);
        value.RegisterNewFab(value.F_mf,              &value.neumann_bc_DD, 4, nghost, "F",      true, false, {"_xx","_xy","_yx","_yy"});
        value.RegisterNewFab(value.P_mf,              &value.neumann_bc_DD, 4, nghost, "P",      true, false, {"_xx","_xy","_yx","_yy"});
        value.RegisterNewFab(value.elastic_stress_mf, &value.neumann_bc_DD, 4, nghost, "elastic_stress", false, false, {"_xx","_xy","_yx","_yy"});
        value.RegisterNewFab(value.fluid_stress_mf,   &value.neumann_bc_DD, 4, nghost, "fluid_stress",   false, false, {"_xx","_xy","_yx","_yy"});
        value.RegisterGeneralFab<Set::Matrix>(value.cauchy_stress_mf, 1, nghost, true, "stress", false);
 
        value.RegisterNewFab(value.pressure_mf,  &value.bc_nothing, 1, nghost, "pressure",  true, false);
        value.RegisterNewFab(value.temperature_mf,  &value.bc_nothing, 1, nghost, "temperature",  true, false);
        value.RegisterNewFab(value.velocity_mf,  &value.bc_nothing, 2, nghost, "velocity",  true, false,{"x","y"});
        value.RegisterNewFab(value.vorticity_mf, &value.bc_nothing, 1, nghost, "vorticity", true, false);

        if (value.UsesDiffuseBoundarySources())
        {
            value.RegisterNewFab(value.m0_mf, &value.bc_nothing, 1, 0, "m0", true, false);
            value.RegisterNewFab(value.u0_mf, &value.bc_nothing, 2, 0, "u0", true, false, {"x","y"});
            value.RegisterNewFab(value.q_mf,  &value.bc_nothing, 2, 0, "q",  true, false, {"x","y"});

            value.RegisterNewFab(value.solid.momentum_mf, &value.neumann_bc_D, 2, nghost, "solid.momentum", true, false, {"x","y"});
            value.RegisterNewFab(value.solid.density_mf,  &value.neumann_bc_1, 1, nghost, "solid.density",  true, false);
            value.RegisterNewFab(value.solid.energy_mf,   &value.neumann_bc_1, 1, nghost, "solid.energy",   true, false);

            value.RegisterNewFab(value.Source_mf, &value.bc_nothing, 4, 0, "Source", true, false);
        }

        value.RegisterNewFab(value.mass_fraction_mf,  &value.bc_nothing, 1, nghost, "mass_fraction",     true , true);
        value.RegisterNewFab(value.mole_fraction_mf,  &value.bc_nothing, 1, nghost, "mole_fraction",     true , true);
        value.RegisterNewFab(value.scratch_mf,  &value.bc_nothing, 1, nghost, "scratch",     false , false);
    }

    pp_forbid("Velocity.ic.type", "--> velocity.ic.type");
    pp_forbid("Pressure.ic", "--> pressure.ic");
    pp_forbid("SolidMomentum.ic", "--> solid.momentum.ic");
    pp_forbid("SolidDensity.ic.type", "--> solid.density.ic.type");
    pp_forbid("SolidEnergy.ic.type", "--> solid.energy.ic.type");
    pp_forbid("Density.ic.type", "--> density.ic.type");
    pp_forbid("rho_injected.ic.type","no longer using rho_injected use m0 instead");
    pp.forbid("mdot.ic.type", "replace mdot with u0");


    // ORDER PARAMETER

    if (!value.managed)
    {
        // eta initial condition
        pp.select_default<IC::Constant,IC::Laminate,IC::Expression,IC::BMP,IC::PNG>("eta.ic",value.eta_ic,value.geom);
    }

    // PRIMITIVE FIELD INITIAL CONDITIONS

    // velocity initial condition
    pp.select_default<IC::Constant,IC::Expression>("velocity.ic",value.velocity_ic,value.geom);
    // solid pressure initial condition
    pp.select_default<IC::Constant,IC::Expression>("pressure.ic",value.pressure_ic,value.geom);
    // density initial condition type
    pp.select_default<IC::Constant,IC::Expression>("density.ic",value.density_ic,value.geom);


    if (value.UsesDiffuseBoundarySources())
    {
        // SOLID FIELDS

        // solid momentum initial condition
        pp.select_default<IC::Constant,IC::Expression>("solid.momentum.ic",value.solid.momentum_ic,value.geom);
        // solid density initial condition
        pp.select_default<IC::Constant,IC::Expression>("solid.density.ic",value.solid.density_ic,value.geom);
        // solid energy initial condition
        pp.select_default<IC::Constant,IC::Expression>("solid.energy.ic",value.solid.energy_ic,value.geom);


        // DIFFUSE BOUNDARY SOURCES

        // diffuse boundary prescribed mass flux
        pp.select_default<IC::Constant,IC::Expression>("m0.ic",value.ic_m0,value.geom);
        // diffuse boundary prescribed velocity
        pp.select_default<IC::Constant,IC::Expression>("u0.ic",value.ic_u0,value.geom);
        // diffuse boundary prescribed heat flux
        pp.select_default<IC::Constant,IC::Expression>("q.ic",value.ic_q,value.geom);
    }

    std::string solid_model_type;
    pp.query_default("solid.model.type", solid_model_type, "none");
    pp.query_default("solid.model.eta_threshold", value.finite_solid_eta_threshold, 0.5);
    pp.query_default("solid.model.viscosity", value.finite_solid_viscosity, 0.0);
    pp.query_default("solid.model.bulk_viscosity", value.finite_solid_bulk_viscosity, 0.0);
    pp.query_default("solid.model.J_floor", value.finite_solid_J_floor, 1.0e-6);
    pp.query_default("solid.model.stress_cap", value.stress_cap, 1.0e100);
    pp.query_default("solid.model.stress_rhs_sign", value.stress_rhs_sign, 1.0);
    pp.query_default("solid.model.pressure_split", value.pressure_split_enabled, true);
    pp.query_default("solid.model.pressure_reference", value.pressure_reference, 0.0);
    pp.query_default("solid.model.reference_density", value.finite_solid_reference_density, 100.0);
    pp.query_default("solid.reference_map.relaxation", value.reference_map_relaxation, 0.0);
    pp.query_default("solid.reference_map.relaxation_eta_cutoff", value.reference_map_relaxation_eta_cutoff, 0.5);
    pp.query_default("solid.reference_map.extrapolation_cutoff", value.reference_map_extrapolation_cutoff, 0.5);
    pp.query_default("solid.reference_map.extrapolation_sweeps", value.reference_map_extrapolation_sweeps, 4);
    pp.query_default("solid.reference_map.smoothing_sweeps", value.reference_map_smoothing_sweeps, 2);
    pp.query_default("eta.ch.enabled", value.eta_ch_enabled, false);
    pp.query_default("eta.ch.mobility", value.eta_ch_mobility, 0.0);
    pp.query_default("eta.ch.kappa", value.eta_ch_kappa, 0.0);
    pp.query_default("eta.ch.barrier", value.eta_ch_barrier, 0.0);
    if (solid_model_type == "none")
    {
        value.finite_solid_enabled = false;
    }
    else if (solid_model_type == "finite.neohookean")
    {
        value.finite_solid_enabled = true;
        pp.queryclass<Model::Solid::Finite::NeoHookean>("solid.model.finite.neohookean", value.finite_solid_model);
    }
    else
    {
        Util::Exception(INFO, solid_model_type, " is not a valid Hydro solid.model.type");
    }

    // Riemann solver
    pp.select_default<  Solver::Local::Riemann::Roe,
                        Solver::Local::Riemann::HLLE,
                        Solver::Local::Riemann::HLLC>("solver",value.riemannsolver);

    // Gas model (Thermo, Transport, and EOS)
    pp.queryclass<Model::Gas::Gas>("gas", value.gas);
    value.nspecies = value.gas.nspecies;
    std::cout << value.nspecies << "\n";

    std::string prescribedflowmode_str;
    // 
    pp.query_validate("prescribedflowmode",prescribedflowmode_str,{"absolute","relative"});
    if (prescribedflowmode_str == "absolute") value.prescribedflowmode = PrescribedFlowMode::Absolute;
    else if (prescribedflowmode_str == "relative") value.prescribedflowmode = PrescribedFlowMode::Relative;

    // Gravitational acceleration vector
    pp.queryarr_default("g",value.g,Set::Vector::Zero());

    bool allow_unused;
    // Set this to true to allow unused inputs without error.
    // (Not recommended.)
    pp.query_default("allow_unused",allow_unused,false);
    if (!allow_unused && pp.AnyUnusedInputs(true, false))
    {
        Util::Warning(INFO,"The following inputs were specified but not used:");
        pp.AllUnusedInputs();
        Util::Exception(INFO,"Aborting. Specify 'allow_unused=True` to ignore this error.");
    }
}


void Hydro::Initialize(int lev)
{
    BL_PROFILE("Integrator::Hydro::Initialize");
 
    if (!managed)
    {
        eta_ic           ->Initialize(lev, *eta_mf,     0.0);
        eta_ic           ->Initialize(lev, *eta_old_mf, 0.0);
    }
    etadot_mf[lev]   ->setVal(0.0);

    //flux_mf[lev]   ->setVal(0.0);

    velocity_ic      ->Initialize(lev, velocity_mf, 0.0);
    pressure_ic      ->Initialize(lev, pressure_mf, 0.0);
    density_ic       ->Initialize(lev, density_mf,  0.0);

    density_ic       ->Initialize(lev, density_old_mf, 0.0);

    for (amrex::MFIter mfi(*xi_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox();
        Set::Patch<Set::Scalar> xi     = xi_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> xi_old = xi_old_mf.Patch(lev,mfi);
        amrex::Geometry const geom_lev = geom[lev];

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Vector pos = Set::Position(i, j, k, geom_lev, amrex::IndexType::TheCellType());
            xi(i,j,k,0) = pos(0);
            xi(i,j,k,1) = pos(1);
            xi_old(i,j,k,0) = pos(0);
            xi_old(i,j,k,1) = pos(1);
        });
    }

    if (UsesDiffuseBoundarySources())
    {
        solid.density_ic ->Initialize(lev, solid.density_mf,  0.0);
        solid.momentum_ic->Initialize(lev, solid.momentum_mf, 0.0);
        solid.energy_ic  ->Initialize(lev, solid.energy_mf,   0.0);

        ic_m0            ->Initialize(lev, m0_mf, 0.0);
        ic_u0            ->Initialize(lev, u0_mf, 0.0);
        ic_q             ->Initialize(lev, q_mf,  0.0);
        Source_mf[lev]->setVal(0.0);
    }

    if (managed)  { if (lev >= (int)mixed.size()) mixed.push_back(false);}
    else if (UsesDiffuseBoundarySources()) Mix(lev);
    else InitializeFluidState(lev);

    UpdateSolidKinematics(lev);
}


void Hydro::RebuildReferenceMapInFluid(int lev, amrex::MultiFab &eta_stage_mf, amrex::MultiFab &xi_stage_mf, Set::Scalar time)
{
    if (eta_mode != EtaMode::Evolving) return;

    const Set::Scalar cutoff_eta = hydro_clamp(reference_map_extrapolation_cutoff, 0.0, 1.0);
    const int extrap_sweeps = reference_map_extrapolation_sweeps < 0 ? 0 : reference_map_extrapolation_sweeps;
    const int smooth_sweeps = reference_map_smoothing_sweeps < 0 ? 0 : reference_map_smoothing_sweeps;
    amrex::Geometry const geom_lev = geom[lev];

    eta_bc->FillBoundary(eta_stage_mf,0,1,time,0);
    eta_stage_mf.FillBoundary(geom[lev].periodicity());
    neumann_bc_D.FillBoundary(xi_stage_mf,0,2,time,0);
    xi_stage_mf.FillBoundary(geom[lev].periodicity());

    for (amrex::MFIter mfi(xi_stage_mf, true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        amrex::Array4<const Set::Scalar> const& eta = eta_stage_mf.const_array(mfi);
        amrex::Array4<Set::Scalar> const& xi = xi_stage_mf.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            if (eta(i,j,k) >= cutoff_eta) return;
            Set::Vector pos = Set::Position(i, j, k, geom_lev, amrex::IndexType::TheCellType());
            for (int n = 0; n < AMREX_SPACEDIM; ++n) xi(i,j,k,n) = pos(n);
        });
    }
    neumann_bc_D.FillBoundary(xi_stage_mf,0,2,time,0);
    xi_stage_mf.FillBoundary(geom[lev].periodicity());

    for (int sweep = 0; sweep < extrap_sweeps; ++sweep)
    {
        for (amrex::MFIter mfi(xi_stage_mf, true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& eta = eta_stage_mf.const_array(mfi);
            amrex::Array4<Set::Scalar> const& xi = xi_stage_mf.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                if (eta(i,j,k) >= cutoff_eta) return;

                for (int n = 0; n < AMREX_SPACEDIM; ++n)
                {
                    Set::Scalar sum = 0.0;
                    Set::Scalar count = 0.0;

                    if (eta(i-1,j,k) >= cutoff_eta) { sum += xi(i-1,j,k,n); count += 1.0; }
                    if (eta(i+1,j,k) >= cutoff_eta) { sum += xi(i+1,j,k,n); count += 1.0; }
                    if (eta(i,j-1,k) >= cutoff_eta) { sum += xi(i,j-1,k,n); count += 1.0; }
                    if (eta(i,j+1,k) >= cutoff_eta) { sum += xi(i,j+1,k,n); count += 1.0; }
#if AMREX_SPACEDIM == 3
                    if (eta(i,j,k-1) >= cutoff_eta) { sum += xi(i,j,k-1,n); count += 1.0; }
                    if (eta(i,j,k+1) >= cutoff_eta) { sum += xi(i,j,k+1,n); count += 1.0; }
#endif
                    if (count > 0.0) xi(i,j,k,n) = sum / count;
                }
            });
        }
        neumann_bc_D.FillBoundary(xi_stage_mf,0,2,time,0);
        xi_stage_mf.FillBoundary(geom[lev].periodicity());
    }

    for (int sweep = 0; sweep < smooth_sweeps; ++sweep)
    {
        for (amrex::MFIter mfi(xi_stage_mf, true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const Set::Scalar> const& eta = eta_stage_mf.const_array(mfi);
            amrex::Array4<Set::Scalar> const& xi = xi_stage_mf.array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                if (eta(i,j,k) >= cutoff_eta) return;

                for (int n = 0; n < AMREX_SPACEDIM; ++n)
                {
                    Set::Scalar center = xi(i,j,k,n);
                    Set::Scalar avg = 0.25 * (xi(i-1,j,k,n) + xi(i+1,j,k,n) + xi(i,j-1,k,n) + xi(i,j+1,k,n));
                    xi(i,j,k,n) = 0.5 * center + 0.5 * avg;
                }
            });
        }
        neumann_bc_D.FillBoundary(xi_stage_mf,0,2,time,0);
        xi_stage_mf.FillBoundary(geom[lev].periodicity());
    }
}


void Hydro::ProjectReferenceMapDensity(int lev, amrex::MultiFab &rho_stage_mf, amrex::MultiFab &M_stage_mf, amrex::MultiFab &E_stage_mf, const amrex::MultiFab &eta_stage_mf, const amrex::MultiFab &xi_stage_mf)
{
    if (eta_mode != EtaMode::Evolving || !finite_solid_enabled) return;

    const Set::Scalar* DX = geom[lev].CellSize();
    amrex::Box domain = geom[lev].Domain();
    const Set::Scalar rho0_solid = finite_solid_reference_density;
    const Set::Scalar det_floor = hydro_max(small, finite_solid_J_floor);

    for (amrex::MFIter mfi(rho_stage_mf, true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        Set::Patch<Set::Scalar> rho = rho_stage_mf.array(mfi);
        Set::Patch<Set::Scalar> M   = M_stage_mf.array(mfi);
        Set::Patch<Set::Scalar> E   = E_stage_mf.array(mfi);
        amrex::Array4<const Set::Scalar> const& eta = eta_stage_mf.const_array(mfi);
        amrex::Array4<const Set::Scalar> const& xi = xi_stage_mf.const_array(mfi);
        Set::Patch<Set::Scalar> T = temperature_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> p = pressure_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> Y = mass_fraction_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> X = mole_fraction_mf.Patch(lev,mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Scalar eta_c = hydro_clamp(eta(i,j,k), 0.0, 1.0);
            Set::Scalar solid_weight = hydro_smootherstep(eta_c);
            if (solid_weight <= 0.0) return;

            auto sten = Numeric::GetStencil(i, j, k, domain);
            Set::Matrix grad_xi = Numeric::Gradient(xi, i, j, k, DX, sten);
            Set::Scalar det_grad_xi = grad_xi.determinant();
            if (!hydro_is_finite(det_grad_xi) || hydro_abs(det_grad_xi) <= det_floor) return;

            Set::Scalar rho_solid = hydro_max(rho0_solid * det_grad_xi, small);
            Set::Scalar rho_old = hydro_max(rho(i,j,k), small);
            Set::Scalar rho_new = hydro_max(solid_weight * rho_solid + (1.0 - solid_weight) * rho_old, small);

            Set::Scalar ux = M(i,j,k,0) / rho_old;
            Set::Scalar uy = M(i,j,k,1) / rho_old;
            Set::Scalar temp = gas.ComputeT(rho_old, M(i,j,k,0), M(i,j,k,1), E(i,j,k), T(i,j,k), X, i, j, k);
            if (!(temp > small) || !hydro_is_finite(temp)) temp = small;

            rho(i,j,k) = rho_new;
            M(i,j,k,0) = rho_new * ux;
            M(i,j,k,1) = rho_new * uy;
            gas.ComputeLocalFractions(rho, Y, X, i, j, k);
            E(i,j,k) = gas.ComputeE(rho_new, M(i,j,k,0), M(i,j,k,1), temp, X, i, j, k);
            T(i,j,k) = temp;
            p(i,j,k) = gas.ComputeP(rho_new, temp, X, i, j, k);
        });
    }
}

void Hydro::UpdateSolidKinematics(int lev)
{
    const Set::Scalar* DX = geom[lev].CellSize();
    amrex::Box domain = geom[lev].Domain();

    (*eta_mf)[lev]->FillBoundary(geom[lev].periodicity());
    xi_mf[lev]->FillBoundary(geom[lev].periodicity());

    for (amrex::MFIter mfi(*xi_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        amrex::Array4<const Set::Scalar> const& eta = (*(*eta_mf)[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& xi = (*xi_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& F_field = (*F_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& P_field = (*P_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& elastic_stress_field = (*elastic_stress_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& fluid_stress_field = (*fluid_stress_mf[lev]).array(mfi);
        Set::Patch<const Set::Scalar> velocity = velocity_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> pressure = pressure_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> T = temperature_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> molef = mole_fraction_mf.Patch(lev,mfi);
        const bool solid_enabled = finite_solid_enabled;
        const Model::Solid::Finite::NeoHookean solid_model = finite_solid_model;
        const Set::Scalar solid_viscosity = finite_solid_viscosity;
        const Set::Scalar solid_bulk_viscosity = finite_solid_bulk_viscosity;
        const Set::Scalar det_floor = hydro_max(small, finite_solid_J_floor);
        const Set::Scalar stress_limit = stress_cap;
        const bool use_pressure_split = pressure_split_enabled;
        const Set::Scalar p_ref = pressure_reference;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            auto sten = Numeric::GetStencil(i, j, k, domain);
            Set::Matrix F = Set::Matrix::Identity();
            Set::Matrix P = Set::Matrix::Zero();
            Set::Matrix sigma_elastic = Set::Matrix::Zero();
            Set::Matrix tau_fluid = Set::Matrix::Zero();

            Set::Scalar eta_weight = hydro_clamp(eta(i,j,k), 0.0, 1.0);
            Set::Scalar solid_weight = hydro_smootherstep(eta_weight);

            bool valid = false;
            if (solid_enabled && solid_weight > 0.0)
            {
                Set::Matrix grad_xi = Numeric::Gradient(xi, i, j, k, DX, sten);
                Set::Scalar det_grad_xi = grad_xi.determinant();
                valid = hydro_is_finite(det_grad_xi) && hydro_abs(det_grad_xi) > det_floor;
                if (valid)
                {
#if AMREX_SPACEDIM == 2
                    Set::Scalar inv_det = 1.0 / det_grad_xi;
                    F(0,0) =  grad_xi(1,1) * inv_det;
                    F(0,1) = -grad_xi(0,1) * inv_det;
                    F(1,0) = -grad_xi(1,0) * inv_det;
                    F(1,1) =  grad_xi(0,0) * inv_det;
#else
                    F = grad_xi.inverse();
#endif
                    valid = hydro_is_finite(F(0,0)) && hydro_is_finite(F(0,1)) &&
                            hydro_is_finite(F(1,0)) && hydro_is_finite(F(1,1));
                }
            }

            if (!valid)
            {
                F = Set::Matrix::Identity();
            }

            if (valid)
            {
                Set::Scalar J = F.determinant();
                if (hydro_is_finite(J) && hydro_abs(J) > det_floor)
                {
                    P = solid_model.DW(F);
                    sigma_elastic = (P * F.transpose()) / J;
                }
            }

            Set::Matrix grad_u = Numeric::Gradient(velocity, i, j, k, DX, sten);
            Set::Scalar div_u = grad_u.trace();
            if (solid_viscosity != 0.0 || solid_bulk_viscosity != 0.0)
            {
                sigma_elastic += solid_viscosity * (grad_u + grad_u.transpose()) +
                                 solid_bulk_viscosity * div_u * Set::Matrix::Identity();
            }

            sigma_elastic = hydro_regularize_stress(sigma_elastic, stress_limit);

            Set::Scalar mu = gas.dynamic_viscosity(T(i,j,k), molef, i, j, k);
            Set::Scalar lambda = 0.0;
            if (hydro_is_finite(mu))
            {
                tau_fluid = mu * (grad_u + grad_u.transpose()) + lambda * div_u * Set::Matrix::Identity();
            }
            if (use_pressure_split) tau_fluid -= hydro_finite_or(pressure(i,j,k) - p_ref, 0.0) * Set::Matrix::Identity();
            tau_fluid = hydro_regularize_stress(tau_fluid, stress_limit);

            F_field(i,j,k,0) = hydro_finite_or(F(0,0), 1.0);
            F_field(i,j,k,1) = hydro_finite_or(F(0,1), 0.0);
            F_field(i,j,k,2) = hydro_finite_or(F(1,0), 0.0);
            F_field(i,j,k,3) = hydro_finite_or(F(1,1), 1.0);
            P_field(i,j,k,0) = hydro_finite_or(P(0,0), 0.0);
            P_field(i,j,k,1) = hydro_finite_or(P(0,1), 0.0);
            P_field(i,j,k,2) = hydro_finite_or(P(1,0), 0.0);
            P_field(i,j,k,3) = hydro_finite_or(P(1,1), 0.0);
            elastic_stress_field(i,j,k,0) = hydro_finite_or(sigma_elastic(0,0), 0.0);
            elastic_stress_field(i,j,k,1) = hydro_finite_or(sigma_elastic(0,1), 0.0);
            elastic_stress_field(i,j,k,2) = hydro_finite_or(sigma_elastic(1,0), 0.0);
            elastic_stress_field(i,j,k,3) = hydro_finite_or(sigma_elastic(1,1), 0.0);
            fluid_stress_field(i,j,k,0) = hydro_finite_or(tau_fluid(0,0), 0.0);
            fluid_stress_field(i,j,k,1) = hydro_finite_or(tau_fluid(0,1), 0.0);
            fluid_stress_field(i,j,k,2) = hydro_finite_or(tau_fluid(1,0), 0.0);
            fluid_stress_field(i,j,k,3) = hydro_finite_or(tau_fluid(1,1), 0.0);
        });
    }

    F_mf[lev]->FillBoundary(geom[lev].periodicity());
    P_mf[lev]->FillBoundary(geom[lev].periodicity());
    elastic_stress_mf[lev]->FillBoundary(geom[lev].periodicity());
    fluid_stress_mf[lev]->FillBoundary(geom[lev].periodicity());

    for (amrex::MFIter mfi(*cauchy_stress_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        amrex::Array4<const Set::Scalar> const& eta = (*(*eta_mf)[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& xi = (*xi_mf[lev]).array(mfi);
        amrex::Array4<Set::Matrix> const& stress = (*cauchy_stress_mf[lev]).array(mfi);
        Set::Patch<const Set::Scalar> velocity = velocity_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> pressure = pressure_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> T = temperature_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> molef = mole_fraction_mf.Patch(lev,mfi);
        amrex::Box const cell_domain = domain;
        const bool solid_enabled = finite_solid_enabled;
        const Model::Solid::Finite::NeoHookean solid_model = finite_solid_model;
        const Set::Scalar solid_viscosity = finite_solid_viscosity;
        const Set::Scalar solid_bulk_viscosity = finite_solid_bulk_viscosity;
        const Set::Scalar det_floor = hydro_max(small, finite_solid_J_floor);
            const Set::Scalar stress_limit = stress_cap;
        const bool use_pressure_split = pressure_split_enabled;
        const Set::Scalar p_ref = pressure_reference;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            auto clamp_cell_i = [=] AMREX_GPU_DEVICE (int ii) -> int
            {
                return ii < cell_domain.smallEnd(0) ? cell_domain.smallEnd(0) : (ii > cell_domain.bigEnd(0) ? cell_domain.bigEnd(0) : ii);
            };
            auto clamp_cell_j = [=] AMREX_GPU_DEVICE (int jj) -> int
            {
                return jj < cell_domain.smallEnd(1) ? cell_domain.smallEnd(1) : (jj > cell_domain.bigEnd(1) ? cell_domain.bigEnd(1) : jj);
            };

            int il = clamp_cell_i(i - 1);
            int ir = clamp_cell_i(i);
            int jb = clamp_cell_j(j - 1);
            int jt = clamp_cell_j(j);

            Set::Matrix grad_xi = Set::Matrix::Zero();
            Set::Matrix grad_u = Set::Matrix::Zero();
            Set::Scalar eta_node = 0.25 * (eta(il,jb,k) + eta(ir,jb,k) + eta(il,jt,k) + eta(ir,jt,k));

            for (int n = 0; n < 2; ++n)
            {
                Set::Scalar xi_left  = 0.5 * (xi(il,jb,k,n) + xi(il,jt,k,n));
                Set::Scalar xi_right = 0.5 * (xi(ir,jb,k,n) + xi(ir,jt,k,n));
                Set::Scalar xi_bot   = 0.5 * (xi(il,jb,k,n) + xi(ir,jb,k,n));
                Set::Scalar xi_top   = 0.5 * (xi(il,jt,k,n) + xi(ir,jt,k,n));
                grad_xi(n,0) = (xi_right - xi_left) / DX[0];
                grad_xi(n,1) = (xi_top - xi_bot) / DX[1];

                Set::Scalar u_left  = 0.5 * (velocity(il,jb,k,n) + velocity(il,jt,k,n));
                Set::Scalar u_right = 0.5 * (velocity(ir,jb,k,n) + velocity(ir,jt,k,n));
                Set::Scalar u_bot   = 0.5 * (velocity(il,jb,k,n) + velocity(ir,jb,k,n));
                Set::Scalar u_top   = 0.5 * (velocity(il,jt,k,n) + velocity(ir,jt,k,n));
                grad_u(n,0) = (u_right - u_left) / DX[0];
                grad_u(n,1) = (u_top - u_bot) / DX[1];
            }

            Set::Matrix sigma_solid = Set::Matrix::Zero();
            Set::Scalar solid_weight = hydro_smootherstep(hydro_clamp(eta_node, 0.0, 1.0));
            if (solid_enabled && solid_weight > 0.0)
            {
                Set::Scalar det_grad_xi = grad_xi.determinant();
                if (hydro_is_finite(det_grad_xi) && hydro_abs(det_grad_xi) > det_floor)
                {
                    Set::Matrix F = Set::Matrix::Identity();
#if AMREX_SPACEDIM == 2
                    Set::Scalar inv_det = 1.0 / det_grad_xi;
                    F(0,0) =  grad_xi(1,1) * inv_det;
                    F(0,1) = -grad_xi(0,1) * inv_det;
                    F(1,0) = -grad_xi(1,0) * inv_det;
                    F(1,1) =  grad_xi(0,0) * inv_det;
#else
                    F = grad_xi.inverse();
#endif
                    Set::Scalar J = F.determinant();
                    if (hydro_is_finite(F(0,0)) && hydro_is_finite(F(0,1)) &&
                        hydro_is_finite(F(1,0)) && hydro_is_finite(F(1,1)) &&
                        hydro_is_finite(J) && hydro_abs(J) > det_floor)
                    {
                        Set::Matrix P = solid_model.DW(F);
                        sigma_solid = (P * F.transpose()) / J;
                    }
                }

                if (solid_viscosity != 0.0 || solid_bulk_viscosity != 0.0)
                {
                    Set::Scalar div_u = grad_u.trace();
                    sigma_solid += solid_viscosity * (grad_u + grad_u.transpose()) +
                                   solid_bulk_viscosity * div_u * Set::Matrix::Identity();
                }
            }

            Set::Matrix tau_fluid = Set::Matrix::Zero();
            Set::Scalar T_node = 0.25 * (T(il,jb,k) + T(ir,jb,k) + T(il,jt,k) + T(ir,jt,k));
            Set::Scalar p_node = 0.25 * (pressure(il,jb,k) + pressure(ir,jb,k) + pressure(il,jt,k) + pressure(ir,jt,k));
            Set::Scalar mu = gas.dynamic_viscosity(T_node, molef, il, jb, k);
            if (hydro_is_finite(mu))
            {
                Set::Scalar div_u = grad_u.trace();
                Set::Scalar lambda = 0.0;
                tau_fluid = mu * (grad_u + grad_u.transpose()) + lambda * div_u * Set::Matrix::Identity();
            }

            sigma_solid = hydro_regularize_stress(sigma_solid, stress_limit);
            if (use_pressure_split) tau_fluid -= hydro_finite_or(p_node - p_ref, 0.0) * Set::Matrix::Identity();
            tau_fluid = hydro_regularize_stress(tau_fluid, stress_limit);
            stress(i,j,k) = hydro_regularize_stress(solid_weight * sigma_solid + (1.0 - solid_weight) * tau_fluid, stress_limit);
        });
    }

    cauchy_stress_mf[lev]->FillBoundary(geom[lev].periodicity());
}

void Hydro::InitializeFluidState(int lev)
{
    for (amrex::MFIter mfi(*velocity_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox();

        Set::Patch<Set::Scalar> v       = velocity_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> p       = pressure_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> rho     = density_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> rho_old = density_old_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> M       = momentum_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> M_old   = momentum_old_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> E       = energy_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> E_old   = energy_old_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> Y       = mass_fraction_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> X       = mole_fraction_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> T       = temperature_mf.Patch(lev,mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            gas.ComputeLocalFractions(rho, Y, X, i, j, k);
            Set::Scalar density = gas.ComputeD(rho, i, j, k);
            T(i,j,k) = gas.ComputeT(p(i,j,k), density, X, i, j, k);

            M(i,j,k,0) = density * v(i,j,k,0);
            M(i,j,k,1) = density * v(i,j,k,1);
            E(i,j,k) = gas.ComputeE(density, M(i,j,k,0), M(i,j,k,1), T(i,j,k), X, i, j, k);

            rho_old(i,j,k) = rho(i,j,k);
            M_old(i,j,k,0) = M(i,j,k,0);
            M_old(i,j,k,1) = M(i,j,k,1);
            E_old(i,j,k) = E(i,j,k);
        });
    }
}

void Hydro::EnforceFluidStateBounds(int lev)
{
    EnforceFluidStateBounds(lev, *density_mf[lev], *momentum_mf[lev], *energy_mf[lev]);
}

void Hydro::EnforceFluidStateBounds(int lev, amrex::MultiFab &rho_mf, amrex::MultiFab &M_mf, amrex::MultiFab &E_mf)
{
    const Set::Scalar rho_floor = small;
    const Set::Scalar T_floor = small;

    for (amrex::MFIter mfi(rho_mf, true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        Set::Patch<Set::Scalar> rho = rho_mf.array(mfi);
        Set::Patch<Set::Scalar> M   = M_mf.array(mfi);
        Set::Patch<Set::Scalar> E   = E_mf.array(mfi);
        Set::Patch<Set::Scalar> v   = velocity_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> p   = pressure_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> T   = temperature_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> Y   = mass_fraction_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> X   = mole_fraction_mf.Patch(lev,mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            if (!(rho(i,j,k) > rho_floor)) rho(i,j,k) = rho_floor;

            gas.ComputeLocalFractions(rho, Y, X, i, j, k);
            Set::Scalar density = gas.ComputeD(rho, i, j, k);
            if (!(density > rho_floor))
            {
                rho(i,j,k) = rho_floor;
                density = rho_floor;
                gas.ComputeLocalFractions(rho, Y, X, i, j, k);
            }

            Set::Scalar temp = gas.ComputeT(density, M(i,j,k,0), M(i,j,k,1), E(i,j,k), T(i,j,k), X, i, j, k);
            if (!(temp > T_floor))
            {
                temp = T_floor;
                E(i,j,k) = gas.ComputeE(density, M(i,j,k,0), M(i,j,k,1), temp, X, i, j, k);
            }

            T(i,j,k) = temp;
            p(i,j,k) = gas.ComputeP(density, temp, X, i, j, k);
            v(i,j,k,0) = M(i,j,k,0) / density;
            v(i,j,k,1) = M(i,j,k,1) / density;
        });
    }
}

void Hydro::Mix(int lev)
{
    if (managed && mixed[lev]) return;

    for (amrex::MFIter mfi(*velocity_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox();

        Set::Patch<const Set::Scalar> eta_patch = eta_old_mf->Patch(lev,mfi);

        Set::Patch<Set::Scalar>       v         = velocity_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       p         = pressure_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       rho       = density_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       rho_old   = density_old_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       M         = momentum_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       M_old     = momentum_old_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       E         = energy_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       E_old     = energy_old_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> rho_solid = solid.density_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> M_solid   = solid.momentum_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> E_solid   = solid.energy_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       Y         = mass_fraction_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       X         = mole_fraction_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       T         = temperature_mf.Patch(lev,mfi);


        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {  
            Set::Scalar eta = invert ? 1.0-eta_patch(i,j,k)*eta_patch(i,j,k) : eta_patch(i,j,k);

            // Initially compute primitives (T,P,u) from given initial conditions
            // But from then on, compute them from mixed values to avoid zero T conditions
            // Except velocity - keep velocity from fluid values only
            gas.ComputeLocalFractions(rho, Y, X, i,j,k); // Get local mole/mass fractions from fluid densities
            Set::Scalar density = gas.ComputeD(rho, i, j, k); // If a gas mixture, this will compute the mixture density
            T(i,j,k) = gas.ComputeT(p(i,j,k), density, X, i, j, k);
            Set::Scalar E_fluid = gas.ComputeE(density, density*v(i,j,k,0), density*v(i,j,k,1), T(i,j,k), X, i, j, k);

            // Mix
            M(i, j, k, 0) = (rho(i, j, k)*v(i, j, k, 0))*eta +  M_solid(i, j, k, 0)*(1.0-eta);
            M(i, j, k, 1) = (rho(i, j, k)*v(i, j, k, 1))*eta +  M_solid(i, j, k, 1)*(1.0-eta);
            M_old(i, j, k, 0) = M(i, j, k, 0);
            M_old(i, j, k, 1) = M(i, j, k, 1);

            rho(i, j, k) = eta * rho(i, j, k) + (1.0 - eta) * rho_solid(i, j, k);
            rho_old(i, j, k) = rho(i, j, k);

            E(i, j, k) = E_fluid*eta + E_solid(i,j,k)*(1.0-eta);
            E_old(i, j, k) = E(i, j, k);
            //Util::Message(INFO,"Energy: ", E(i,j,k), " Pressure: ", p(i,j,k), " Temp: ", T(i,j,k), " Density: ",density, " R: ", gas.R(X,i,j,k), " MW: ", gas.GetMW(X,i,j,k), " Rg: ", Set::Constant::Rg);

            //gas.ComputeLocalFractions(rho, Y, X, i,j,k); // Get local mole/mass fractions from mixed densities
            //density = gas.ComputeD(rho, i, j, k);
            //T(i, j, k) = gas.ComputeT(density, M(i,j,k,0), M(i,j,k,1), E(i,j,k), T(i,j,k), X, i, j, k);
            //p(i, j, k) = gas.ComputeP(density, T(i,j,k), X, i, j, k);
            //v(i,j,k,0) = M(i,j,k,0)/density;
            //v(i,j,k,1) = M(i,j,k,1)/density;
        });
        //Util::Abort(INFO);
    }
    c_max = 0.0;
    vx_max = 0.0;
    vy_max = 0.0;
}

void Hydro::UpdateEta(int lev, Set::Scalar time)
{
    Util::Assert(INFO,TEST(!managed),"Should override this if Hydro is managed!");
    eta_ic->Initialize(lev, *eta_mf, time);
}

void Hydro::UpdateFluxes(int /*lev*/, Set::Scalar /*time*/, Set::Scalar /*dt*/)
{
    Util::Assert(INFO,TEST(!managed),"Should override this if Hydro is managed!");
}

void Hydro::TimeStepBegin(Set::Scalar, int /*iter*/)
{

}

void Hydro::TimeStepComplete(Set::Scalar, int lev)
{
    if (dynamictimestep.on)
        Integrator::DynamicTimestep_Update();
    return;

    const Set::Scalar* DX = geom[lev].CellSize();

    amrex::ParallelDescriptor::ReduceRealMax(c_max);
    amrex::ParallelDescriptor::ReduceRealMax(vx_max);
    amrex::ParallelDescriptor::ReduceRealMax(vy_max);

    Set::Scalar new_timestep = cfl / ((c_max + vx_max) / DX[0] + (c_max + vy_max) / DX[1]);

    Util::Assert(INFO, TEST(AMREX_SPACEDIM == 2));

    SetTimestep(new_timestep);
}

void Hydro::Advance(int lev, Set::Scalar time, Set::Scalar dt)
{

    if (!managed) std::swap(*eta_old_mf, *eta_mf);
    if (eta_mode == EtaMode::Evolving) std::swap(xi_old_mf[lev], xi_mf[lev]);
    std::swap(density_old_mf[lev],  density_mf[lev]);
    std::swap(momentum_old_mf[lev], momentum_mf[lev]);
    std::swap(energy_old_mf[lev],   energy_mf[lev]);
    
    //
    // UPDATE ETA AND CALCULATE ETADOT
    //

    if (eta_mode == EtaMode::Static)
    {
        if (!managed) UpdateEta(lev, time);
        if (managed) 
        {
            UpdateFluxes(lev,time,dt);
            Mix(lev);
        }
        for (amrex::MFIter mfi(*(velocity_mf)[lev], true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.growntilebox();
            amrex::Array4<const Set::Scalar> const& eta_new = (*(*eta_mf)[lev]).array(mfi);
            amrex::Array4<const Set::Scalar> const& eta = (*(*eta_old_mf)[lev]).array(mfi);
            amrex::Array4<Set::Scalar>       const& etadot = (*etadot_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {   

                etadot(i, j, k) = (eta_new(i, j, k) - eta(i, j, k)) / dt;
                if (invert) etadot(i,j,k) *= 1.0;

            });
        }
    }
    else
    {
        etadot_mf[lev]->setVal(0.0);
    }


    //
    // DO TIME INTEGRATION (driving the RHS function)
    //

    // Organize references to the "new" solution
    amrex::Vector<amrex::MultiFab> solution_new; 
    solution_new.emplace_back(*density_mf[lev].get(),amrex::MakeType::make_alias,0,1);
    solution_new.emplace_back(*momentum_mf[lev].get(),amrex::MakeType::make_alias,0,2);
    solution_new.emplace_back(*energy_mf[lev].get(),amrex::MakeType::make_alias,0,1);
    if (eta_mode == EtaMode::Evolving)
    {
        solution_new.emplace_back(*(*eta_mf)[lev].get(),amrex::MakeType::make_alias,0,1);
        solution_new.emplace_back(*xi_mf[lev].get(),amrex::MakeType::make_alias,0,2);
    }

    // Organize references to the "old" solution
    amrex::Vector<amrex::MultiFab> solution_old;
    solution_old.emplace_back(*density_old_mf[lev].get(),amrex::MakeType::make_alias,0,1);
    solution_old.emplace_back(*momentum_old_mf[lev].get(),amrex::MakeType::make_alias,0,2);
    solution_old.emplace_back(*energy_old_mf[lev].get(),amrex::MakeType::make_alias,0,1);
    if (eta_mode == EtaMode::Evolving)
    {
        solution_old.emplace_back(*(*eta_old_mf)[lev].get(),amrex::MakeType::make_alias,0,1);
        solution_old.emplace_back(*xi_old_mf[lev].get(),amrex::MakeType::make_alias,0,2);
    }

    // Create the time integrator
    amrex::TimeIntegrator timeintegrator(solution_new, time);

    // Set the time integrator RHS - in this case, just relay to our current RHS function
    timeintegrator.set_rhs([&](amrex::Vector<amrex::MultiFab> & rhs_mf, amrex::Vector<amrex::MultiFab> & solution_mf, const Set::Scalar time)
    {
        if (eta_mode == EtaMode::Evolving)
            RHS(lev, time,
                rhs_mf[0], rhs_mf[1], rhs_mf[2],
                solution_mf[0],solution_mf[1],solution_mf[2],
                &rhs_mf[3], &solution_mf[3],
                &rhs_mf[4], &solution_mf[4]);
        else
            RHS(lev, time,
                rhs_mf[0], rhs_mf[1], rhs_mf[2],
                solution_mf[0],solution_mf[1],solution_mf[2]);
    });

    // Take care of filling boundaries during stages
    timeintegrator.set_post_stage_action([&](amrex::Vector<amrex::MultiFab> & stage_mf, Set::Scalar time) 
    {
        if (!UsesDiffuseBoundarySources())
            EnforceFluidStateBounds(lev, stage_mf[0], stage_mf[1], stage_mf[2]);

        density_bc->FillBoundary(stage_mf[0],0,1,time,0);
        stage_mf[0].FillBoundary(true);
        momentum_bc->FillBoundary(stage_mf[1],0,2,time,0);
        stage_mf[1].FillBoundary(true);
        energy_bc->FillBoundary(stage_mf[2],0,1,time,0);
        stage_mf[2].FillBoundary(true);
        if (eta_mode == EtaMode::Evolving)
        {
            eta_bc->FillBoundary(stage_mf[3],0,1,time,0);
            stage_mf[3].FillBoundary(true);
            neumann_bc_D.FillBoundary(stage_mf[4],0,2,time,0);
            stage_mf[4].FillBoundary(true);
            RebuildReferenceMapInFluid(lev, stage_mf[3], stage_mf[4], time);
            ProjectReferenceMapDensity(lev, stage_mf[0], stage_mf[1], stage_mf[2], stage_mf[3], stage_mf[4]);
            density_bc->FillBoundary(stage_mf[0],0,1,time,0);
            stage_mf[0].FillBoundary(true);
            momentum_bc->FillBoundary(stage_mf[1],0,2,time,0);
            stage_mf[1].FillBoundary(true);
            energy_bc->FillBoundary(stage_mf[2],0,1,time,0);
            stage_mf[2].FillBoundary(true);
        }
    });
    
    // Do the update
    timeintegrator.advance(solution_old, solution_new, time, dt);
    if (!UsesDiffuseBoundarySources()) EnforceFluidStateBounds(lev);
    if (eta_mode == EtaMode::Evolving)
    {
        RebuildReferenceMapInFluid(lev, *(*eta_mf)[lev], *xi_mf[lev], time + dt);
        ProjectReferenceMapDensity(lev, *density_mf[lev], *momentum_mf[lev], *energy_mf[lev], *(*eta_mf)[lev], *xi_mf[lev]);
    }
    UpdateSolidKinematics(lev);


    //
    // APPLY CUTOFFS AND DO DYNAMIC TIMESTEP CALCULATION
    //

    Set::Scalar dt_max = std::numeric_limits<Set::Scalar>::max();
    if (eta_mode == EtaMode::Evolving)
    {
        for (amrex::MFIter mfi(*velocity_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            const Set::Scalar* DX = geom[lev].CellSize();

            Set::Patch<Set::Scalar> omega = vorticity_mf.Patch(lev,mfi);
            Set::Patch<Set::Scalar> u     = velocity_mf.Patch(lev,mfi);

            Set::Scalar *dt_max_handle = &dt_max;

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Matrix gradu = Numeric::Gradient(u, i, j, k, DX);
                omega(i, j, k) = gradu(1,0) - gradu(0,1);

                if (dynamictimestep.on)
                {
                    *dt_max_handle =                          std::fabs(cfl * DX[0] / (u(i,j,k,0) + small));
                    *dt_max_handle = std::min(*dt_max_handle, std::fabs(cfl * DX[1] / (u(i,j,k,1) + small)));
                }
            });
        }
    }
    else
    {
        for (amrex::MFIter mfi(*velocity_mf[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
            const Set::Scalar* DX = geom[lev].CellSize();
            
            Set::Patch<const Set::Scalar> eta_patch = eta_mf->Patch(lev,mfi);
            Set::Patch<const Set::Scalar> rho_solid = solid.density_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> M_solid   = solid.momentum_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> E_solid   = solid.energy_mf.Patch(lev,mfi);

            Set::Patch<Set::Scalar> rho_new       = density_mf.Patch(lev,mfi);
            Set::Patch<Set::Scalar> E_new         = energy_mf.Patch(lev,mfi);
            Set::Patch<Set::Scalar> M_new         = momentum_mf.Patch(lev,mfi);

            Set::Patch<Set::Scalar> omega         = vorticity_mf.Patch(lev,mfi);
            
            Set::Patch<Set::Scalar> u = velocity_mf.Patch(lev,mfi);
            Set::Patch<Set::Scalar> Source = Source_mf.Patch(lev,mfi);

            Set::Scalar *dt_max_handle = &dt_max;

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {   
                Set::Scalar eta = invert ? 1.0-eta_patch(i,j,k)*eta_patch(i,j,k) : eta_patch(i,j,k);

                if (eta < cutoff)
                {
                    rho_new(i,j,k,0) = rho_solid(i,j,k,0);
                    M_new(i,j,k,0)   = M_solid(i,j,k,0);
                    M_new(i,j,k,1)   = M_solid(i,j,k,1);
                    E_new(i,j,k,0)   = E_solid(i,j,k,0);
                }

                Set::Matrix gradu        = Numeric::Gradient(u, i, j, k, DX);
                omega(i, j, k) = eta * (gradu(1,0) - gradu(0,1));

                if (dynamictimestep.on)
                {
                    *dt_max_handle =                          std::fabs(cfl * DX[0] / (u(i,j,k,0)*eta + small));
                    *dt_max_handle = std::min(*dt_max_handle, std::fabs(cfl * DX[1] / (u(i,j,k,1)*eta + small)));
                    *dt_max_handle = std::min(*dt_max_handle, std::fabs(cfl_v * DX[0]*DX[0] / (Source(i,j,k,1)+small)));
                    *dt_max_handle = std::min(*dt_max_handle, std::fabs(cfl_v * DX[1]*DX[1] / (Source(i,j,k,2)+small)));
                }
            });
        }
    }


    if (dynamictimestep.on)
    {
        this->DynamicTimestep_SyncTimeStep(lev,dt_max);
    }

}//end Advance


void Hydro::RHS(int lev, Set::Scalar /*time*/, 
                amrex::MultiFab &rho_rhs_mf, 
                amrex::MultiFab &M_rhs_mf, 
                amrex::MultiFab &E_rhs_mf,
                const amrex::MultiFab &rho_mf,
                const amrex::MultiFab &M_mf,
                const amrex::MultiFab &E_mf,
                amrex::MultiFab *eta_rhs_mf,
                const amrex::MultiFab *eta_stage_mf,
                amrex::MultiFab *xi_rhs_mf,
                const amrex::MultiFab *xi_stage_mf)
{

    const bool evolve_eta = (eta_mode == EtaMode::Evolving);
    Util::Assert(INFO, TEST(!evolve_eta || (eta_rhs_mf != nullptr && eta_stage_mf != nullptr)));
    Util::Assert(INFO, TEST(!evolve_eta || (xi_rhs_mf != nullptr && xi_stage_mf != nullptr)));
    const bool use_diffuse_sources = UsesDiffuseBoundarySources();
    const amrex::MultiFab& eta_current_mf = evolve_eta ? *eta_stage_mf : *(*eta_old_mf)[lev];
    const amrex::MultiFab& eta_iter_mf = evolve_eta ? *eta_stage_mf : *(*eta_mf)[lev];

    if (evolve_eta)
    {
        for (amrex::MFIter mfi(*(velocity_mf)[lev], true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.growntilebox();

            Set::Patch<const Set::Scalar> rho = rho_mf.array(mfi);
            Set::Patch<const Set::Scalar> M   = M_mf.array(mfi);
            Set::Patch<const Set::Scalar> E   = E_mf.array(mfi);

            Set::Patch<Set::Scalar> v = velocity_mf.Patch(lev,mfi);
            Set::Patch<Set::Scalar> p = pressure_mf.Patch(lev,mfi);
            Set::Patch<Set::Scalar> T = temperature_mf.Patch(lev,mfi);
            Set::Patch<Set::Scalar> Y = mass_fraction_mf.Patch(lev,mfi);
            Set::Patch<Set::Scalar> X = mole_fraction_mf.Patch(lev,mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                gas.ComputeLocalFractions(rho, Y, X, i, j, k);
                Set::Scalar density = gas.ComputeD(rho, i, j, k);
                T(i,j,k) = gas.ComputeT(density, M(i,j,k,0), M(i,j,k,1), E(i,j,k), T(i,j,k), X, i, j, k);
                p(i,j,k) = gas.ComputeP(density, T(i,j,k), X, i, j, k);
                v(i,j,k,0) = M(i,j,k,0) / density;
                v(i,j,k,1) = M(i,j,k,1) / density;
            });
        }

        const Set::Scalar* DX = geom[lev].CellSize();
        amrex::Box domain = geom[lev].Domain();
        for (amrex::MFIter mfi(*eta_stage_mf, false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();

            Set::Patch<const Set::Scalar> rho = rho_mf.array(mfi);
            Set::Patch<const Set::Scalar> E   = E_mf.array(mfi);
            Set::Patch<const Set::Scalar> M   = M_mf.array(mfi);

            Set::Patch<Set::Scalar> rho_rhs = rho_rhs_mf.array(mfi);
            Set::Patch<Set::Scalar> M_rhs   = M_rhs_mf.array(mfi);
            Set::Patch<Set::Scalar> E_rhs   = E_rhs_mf.array(mfi);

            amrex::Array4<const Set::Scalar> const& eta_patch = eta_stage_mf->array(mfi);
            amrex::Array4<Set::Scalar> const& eta_rhs = eta_rhs_mf->array(mfi);
            amrex::Array4<const Set::Scalar> const& xi_patch = xi_stage_mf->array(mfi);
            amrex::Array4<Set::Scalar> const& xi_rhs = xi_rhs_mf->array(mfi);
            Set::Patch<const Set::Scalar> velocity = velocity_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> pstage   = pressure_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> Tstage   = temperature_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> molef    = mole_fraction_mf.Patch(lev,mfi);
            const bool solid_enabled = finite_solid_enabled;
            const Model::Solid::Finite::NeoHookean solid_model = finite_solid_model;
            const Set::Scalar solid_viscosity = finite_solid_viscosity;
            const Set::Scalar solid_bulk_viscosity = finite_solid_bulk_viscosity;
            const Set::Scalar det_floor = hydro_max(small, finite_solid_J_floor);
        const Set::Scalar stress_limit = stress_cap;
            const Set::Scalar rhs_stress_sign = stress_rhs_sign;
            const bool use_pressure_split = pressure_split_enabled;
            const Set::Scalar p_ref = pressure_reference;
            const Set::Scalar xi_relaxation = reference_map_relaxation;
            const Set::Scalar xi_relaxation_eta_cutoff = reference_map_relaxation_eta_cutoff;
            const Set::Scalar xi_evolve_cutoff = hydro_clamp(reference_map_extrapolation_cutoff, 0.0, 1.0);
            const bool ch_enabled = eta_ch_enabled && eta_ch_mobility > 0.0;
            const Set::Scalar ch_mobility = eta_ch_mobility;
            const Set::Scalar ch_kappa = eta_ch_kappa;
            const Set::Scalar ch_barrier = eta_ch_barrier;
            amrex::Geometry const geom_lev = geom[lev];

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar eta_c = eta_patch(i,j,k);
                Set::Scalar eta_xm = eta_patch(i-1,j,k);
                Set::Scalar eta_xp = eta_patch(i+1,j,k);
                Set::Scalar eta_ym = eta_patch(i,j-1,k);
                Set::Scalar eta_yp = eta_patch(i,j+1,k);

                Set::Scalar slope_xm = hydro_mc_slope(eta_patch(i-2,j,k), eta_xm, eta_c);
                Set::Scalar slope_x  = hydro_mc_slope(eta_xm, eta_c, eta_xp);
                Set::Scalar slope_xp = hydro_mc_slope(eta_c, eta_xp, eta_patch(i+2,j,k));
                Set::Scalar slope_ym = hydro_mc_slope(eta_patch(i,j-2,k), eta_ym, eta_c);
                Set::Scalar slope_y  = hydro_mc_slope(eta_ym, eta_c, eta_yp);
                Set::Scalar slope_yp = hydro_mc_slope(eta_c, eta_yp, eta_patch(i,j+2,k));

                Set::Scalar ux_xlo = 0.5 * (velocity(i-1,j,k,0) + velocity(i,j,k,0));
                Set::Scalar ux_xhi = 0.5 * (velocity(i,j,k,0) + velocity(i+1,j,k,0));
                Set::Scalar uy_ylo = 0.5 * (velocity(i,j-1,k,1) + velocity(i,j,k,1));
                Set::Scalar uy_yhi = 0.5 * (velocity(i,j,k,1) + velocity(i,j+1,k,1));

                Set::Scalar eta_xlo_l = hydro_clamp(eta_xm + 0.5 * slope_xm, eta_xm, eta_c);
                Set::Scalar eta_xlo_r = hydro_clamp(eta_c  - 0.5 * slope_x,  eta_xm, eta_c);
                Set::Scalar eta_xhi_l = hydro_clamp(eta_c  + 0.5 * slope_x,  eta_c, eta_xp);
                Set::Scalar eta_xhi_r = hydro_clamp(eta_xp - 0.5 * slope_xp, eta_c, eta_xp);

                Set::Scalar eta_ylo_l = hydro_clamp(eta_ym + 0.5 * slope_ym, eta_ym, eta_c);
                Set::Scalar eta_ylo_r = hydro_clamp(eta_c  - 0.5 * slope_y,  eta_ym, eta_c);
                Set::Scalar eta_yhi_l = hydro_clamp(eta_c  + 0.5 * slope_y,  eta_c, eta_yp);
                Set::Scalar eta_yhi_r = hydro_clamp(eta_yp - 0.5 * slope_yp, eta_c, eta_yp);

                Set::Scalar eta_flux_xlo = ux_xlo * (ux_xlo >= 0.0 ? eta_xlo_l : eta_xlo_r);
                Set::Scalar eta_flux_xhi = ux_xhi * (ux_xhi >= 0.0 ? eta_xhi_l : eta_xhi_r);
                Set::Scalar eta_flux_ylo = uy_ylo * (uy_ylo >= 0.0 ? eta_ylo_l : eta_ylo_r);
                Set::Scalar eta_flux_yhi = uy_yhi * (uy_yhi >= 0.0 ? eta_yhi_l : eta_yhi_r);

                Set::Scalar div_eta_u =
                    (eta_flux_xhi - eta_flux_xlo) / DX[0] +
                    (eta_flux_yhi - eta_flux_ylo) / DX[1];
                Set::Scalar div_u =
                    (ux_xhi - ux_xlo) / DX[0] +
                    (uy_yhi - uy_ylo) / DX[1];

                Set::Scalar eta_adv_rhs = -div_eta_u + eta_c * div_u;
                Set::Scalar eta_ch_rhs = 0.0;
                if (ch_enabled)
                {
                    auto eta_value = [=] AMREX_GPU_DEVICE (int ii, int jj) -> Set::Scalar
                    {
                        return hydro_clamp(eta_patch(ii,jj,k), 0.0, 1.0);
                    };
                    auto lap_eta = [=] AMREX_GPU_DEVICE (int ii, int jj) -> Set::Scalar
                    {
                        return (eta_value(ii+1,jj) - 2.0 * eta_value(ii,jj) + eta_value(ii-1,jj)) / (DX[0] * DX[0]) +
                               (eta_value(ii,jj+1) - 2.0 * eta_value(ii,jj) + eta_value(ii,jj-1)) / (DX[1] * DX[1]);
                    };
                    auto dfdeta = [=] AMREX_GPU_DEVICE (Set::Scalar a) -> Set::Scalar
                    {
                        return ch_barrier * 2.0 * a * (1.0 - a) * (1.0 - 2.0 * a);
                    };
                    Set::Scalar mu_c = dfdeta(eta_value(i,j)) - ch_kappa * lap_eta(i,j);
                    Set::Scalar mu_xp = dfdeta(eta_value(i+1,j)) - ch_kappa * lap_eta(i+1,j);
                    Set::Scalar mu_xm = dfdeta(eta_value(i-1,j)) - ch_kappa * lap_eta(i-1,j);
                    Set::Scalar mu_yp = dfdeta(eta_value(i,j+1)) - ch_kappa * lap_eta(i,j+1);
                    Set::Scalar mu_ym = dfdeta(eta_value(i,j-1)) - ch_kappa * lap_eta(i,j-1);
                    eta_ch_rhs = ch_mobility * ((mu_xp - 2.0 * mu_c + mu_xm) / (DX[0] * DX[0]) +
                                                (mu_yp - 2.0 * mu_c + mu_ym) / (DX[1] * DX[1]));
                }
                eta_rhs(i,j,k) = eta_adv_rhs + eta_ch_rhs;

                Set::Vector identity_xi = Set::Position(i, j, k, geom_lev, amrex::IndexType::TheCellType());
                Set::Scalar relax_cutoff = hydro_clamp(xi_relaxation_eta_cutoff, 1.0e-12, 1.0);
                Set::Scalar xi_relax_weight = (xi_relaxation > 0.0) ? (1.0 - hydro_smootherstep(eta_c / relax_cutoff)) : 0.0;

                for (int n = 0; n < 2; ++n)
                {
                    Set::Scalar xi_c = xi_patch(i,j,k,n);
                    Set::Scalar xi_xm = xi_patch(i-1,j,k,n);
                    Set::Scalar xi_xp = xi_patch(i+1,j,k,n);
                    Set::Scalar xi_ym = xi_patch(i,j-1,k,n);
                    Set::Scalar xi_yp = xi_patch(i,j+1,k,n);

                    Set::Scalar grad_xi_x = hydro_mc_slope(xi_xm, xi_c, xi_xp) / DX[0];
                    Set::Scalar grad_xi_y = hydro_mc_slope(xi_ym, xi_c, xi_yp) / DX[1];

                    if (eta_c >= xi_evolve_cutoff)
                    {
                        xi_rhs(i,j,k,n) = -velocity(i,j,k,0) * grad_xi_x - velocity(i,j,k,1) * grad_xi_y
                                          - xi_relaxation * xi_relax_weight * (xi_c - identity_xi(n));
                    }
                    else
                    {
                        xi_rhs(i,j,k,n) = 0.0;
                    }
                }

                #if AMREX_SPACEDIM == 2
                    Set::Vector u = Set::Vector(velocity(i, j, k, 0), velocity(i, j, k, 1));
                #endif

                #if AMREX_SPACEDIM == 3
                    Set::Vector u = Set::Vector(velocity(i, j, k, 0), velocity(i, j, k, 1), velocity(i, j, k, 2));
                #endif

                const int X = 0, Y = 1;
                Solver::Local::Riemann::State state_xlo(rho, M, E, i-1, j, k, X);
                Solver::Local::Riemann::State state_x  (rho, M, E, i  , j, k, X);
                Solver::Local::Riemann::State state_xhi(rho, M, E, i+1, j, k, X);

                Solver::Local::Riemann::State state_ylo(rho, M, E, i, j-1, k, Y);
                Solver::Local::Riemann::State state_y  (rho, M, E, i, j  , k, Y);
                Solver::Local::Riemann::State state_yhi(rho, M, E, i, j+1, k, Y);

                Solver::Local::Riemann::Flux flux_xlo, flux_ylo, flux_xhi, flux_yhi;

                try
                {
                    flux_xlo = riemannsolver->Solve(state_xlo, state_x, gas, molef, i, j, k, 0, small);
                    flux_ylo = riemannsolver->Solve(state_ylo, state_y, gas, molef, i, j, k, 2, small);
                    flux_xhi = riemannsolver->Solve(state_x, state_xhi, gas, molef, i, j, k, 1, small);
                    flux_yhi = riemannsolver->Solve(state_y, state_yhi, gas, molef, i, j, k, 3, small);
                }
                catch(...)
                {
                    Util::ParallelMessage(INFO,"lev=",lev);
                    Util::ParallelMessage(INFO,"i=",i,"j=",j);
                    Util::Abort(INFO);
                }

                rho_rhs(i,j,k) =
                    (flux_xlo.mass - flux_xhi.mass) / DX[0] +
                    (flux_ylo.mass - flux_yhi.mass) / DX[1];

                auto clamp_cell_i = [=] AMREX_GPU_DEVICE (int ii) -> int
                {
                    return ii < domain.smallEnd(0) ? domain.smallEnd(0) : (ii > domain.bigEnd(0) ? domain.bigEnd(0) : ii);
                };
                auto clamp_cell_j = [=] AMREX_GPU_DEVICE (int jj) -> int
                {
                    return jj < domain.smallEnd(1) ? domain.smallEnd(1) : (jj > domain.bigEnd(1) ? domain.bigEnd(1) : jj);
                };

                auto calc_nodal_stress = [=] AMREX_GPU_DEVICE (int ni, int nj) -> Set::Matrix
                {
                    int il = clamp_cell_i(ni - 1);
                    int ir = clamp_cell_i(ni);
                    int jb = clamp_cell_j(nj - 1);
                    int jt = clamp_cell_j(nj);

                    Set::Matrix grad_xi = Set::Matrix::Zero();
                    Set::Matrix grad_u = Set::Matrix::Zero();
                    Set::Scalar eta_node = 0.25 * (eta_patch(il,jb,k) + eta_patch(ir,jb,k) + eta_patch(il,jt,k) + eta_patch(ir,jt,k));

                    for (int n = 0; n < 2; ++n)
                    {
                        Set::Scalar xi_left  = 0.5 * (xi_patch(il,jb,k,n) + xi_patch(il,jt,k,n));
                        Set::Scalar xi_right = 0.5 * (xi_patch(ir,jb,k,n) + xi_patch(ir,jt,k,n));
                        Set::Scalar xi_bot   = 0.5 * (xi_patch(il,jb,k,n) + xi_patch(ir,jb,k,n));
                        Set::Scalar xi_top   = 0.5 * (xi_patch(il,jt,k,n) + xi_patch(ir,jt,k,n));
                        grad_xi(n,0) = (xi_right - xi_left) / DX[0];
                        grad_xi(n,1) = (xi_top - xi_bot) / DX[1];

                        Set::Scalar u_left  = 0.5 * (velocity(il,jb,k,n) + velocity(il,jt,k,n));
                        Set::Scalar u_right = 0.5 * (velocity(ir,jb,k,n) + velocity(ir,jt,k,n));
                        Set::Scalar u_bot   = 0.5 * (velocity(il,jb,k,n) + velocity(ir,jb,k,n));
                        Set::Scalar u_top   = 0.5 * (velocity(il,jt,k,n) + velocity(ir,jt,k,n));
                        grad_u(n,0) = (u_right - u_left) / DX[0];
                        grad_u(n,1) = (u_top - u_bot) / DX[1];
                    }

                    Set::Matrix sigma_solid = Set::Matrix::Zero();
                    Set::Scalar solid_weight = hydro_smootherstep(hydro_clamp(eta_node, 0.0, 1.0));
                    if (solid_enabled && solid_weight > 0.0)
                    {
                        Set::Scalar det_grad_xi = grad_xi.determinant();
                        if (hydro_is_finite(det_grad_xi) && hydro_abs(det_grad_xi) > det_floor)
                        {
                            Set::Matrix F = Set::Matrix::Identity();
#if AMREX_SPACEDIM == 2
                            Set::Scalar inv_det = 1.0 / det_grad_xi;
                            F(0,0) =  grad_xi(1,1) * inv_det;
                            F(0,1) = -grad_xi(0,1) * inv_det;
                            F(1,0) = -grad_xi(1,0) * inv_det;
                            F(1,1) =  grad_xi(0,0) * inv_det;
#else
                            F = grad_xi.inverse();
#endif
                            Set::Scalar J = F.determinant();
                            if (hydro_is_finite(F(0,0)) && hydro_is_finite(F(0,1)) &&
                                hydro_is_finite(F(1,0)) && hydro_is_finite(F(1,1)) &&
                                hydro_is_finite(J) && hydro_abs(J) > det_floor)
                            {
                                Set::Matrix P = solid_model.DW(F);
                                sigma_solid = (P * F.transpose()) / J;
                            }
                        }

                        if (solid_viscosity != 0.0 || solid_bulk_viscosity != 0.0)
                        {
                            Set::Scalar div_u_node = grad_u.trace();
                            sigma_solid += solid_viscosity * (grad_u + grad_u.transpose()) +
                                           solid_bulk_viscosity * div_u_node * Set::Matrix::Identity();
                        }
                    }

                    Set::Matrix tau_fluid = Set::Matrix::Zero();
                    Set::Scalar T_node = 0.25 * (Tstage(il,jb,k) + Tstage(ir,jb,k) + Tstage(il,jt,k) + Tstage(ir,jt,k));
                    Set::Scalar p_node = 0.25 * (pstage(il,jb,k) + pstage(ir,jb,k) + pstage(il,jt,k) + pstage(ir,jt,k));
                    Set::Scalar mu = gas.dynamic_viscosity(T_node, molef, i, j, k);
                    if (hydro_is_finite(mu))
                    {
                        Set::Scalar div_u_node = grad_u.trace();
                        Set::Scalar lambda = 0.0;
                        tau_fluid = mu * (grad_u + grad_u.transpose()) + lambda * div_u_node * Set::Matrix::Identity();
                    }

                    Set::Scalar p_mech = hydro_finite_or(p_node - p_ref, 0.0);
                    Set::Matrix fluid_actual = tau_fluid;
                    if (use_pressure_split) fluid_actual -= p_mech * Set::Matrix::Identity();
                    Set::Matrix actual_stress = solid_weight * sigma_solid + (1.0 - solid_weight) * fluid_actual;

                    sigma_solid = hydro_regularize_stress(sigma_solid, stress_limit);
                    tau_fluid = hydro_regularize_stress(tau_fluid, stress_limit);
                    actual_stress = hydro_regularize_stress(actual_stress, stress_limit);
                    if (use_pressure_split) actual_stress += p_mech * Set::Matrix::Identity();
                    return hydro_regularize_stress(actual_stress, stress_limit);
                };

                Set::Matrix sigma_ll = calc_nodal_stress(i,   j);
                Set::Matrix sigma_lr = calc_nodal_stress(i+1, j);
                Set::Matrix sigma_ul = calc_nodal_stress(i,   j+1);
                Set::Matrix sigma_ur = calc_nodal_stress(i+1, j+1);

                Set::Matrix sigma_xlo = 0.5 * (sigma_ll + sigma_ul);
                Set::Matrix sigma_xhi = 0.5 * (sigma_lr + sigma_ur);
                Set::Matrix sigma_ylo = 0.5 * (sigma_ll + sigma_lr);
                Set::Matrix sigma_yhi = 0.5 * (sigma_ul + sigma_ur);

                Set::Vector div_sigma = Set::Vector::Zero();
                div_sigma(0) = (sigma_xhi(0,0) - sigma_xlo(0,0)) / DX[0] +
                               (sigma_yhi(0,1) - sigma_ylo(0,1)) / DX[1];
                div_sigma(1) = (sigma_xhi(1,0) - sigma_xlo(1,0)) / DX[0] +
                               (sigma_yhi(1,1) - sigma_ylo(1,1)) / DX[1];

                M_rhs(i,j,k,0) =
                    (flux_xlo.momentum_normal  - flux_xhi.momentum_normal ) / DX[0] +
                    (flux_ylo.momentum_tangent - flux_yhi.momentum_tangent) / DX[1] +
                    rhs_stress_sign * div_sigma(0) +
                    g(0)*rho(i,j,k);

                M_rhs(i,j,k,1) =
                    (flux_xlo.momentum_tangent - flux_xhi.momentum_tangent) / DX[0] +
                    (flux_ylo.momentum_normal  - flux_yhi.momentum_normal ) / DX[1] +
                    rhs_stress_sign * div_sigma(1) +
                    g(1)*rho(i,j,k);

                E_rhs(i,j,k) =
                    (flux_xlo.energy - flux_xhi.energy) / DX[0] +
                    (flux_ylo.energy - flux_yhi.energy) / DX[1];
            });
        }

        return;
    }

    for (amrex::MFIter mfi(*(velocity_mf)[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox();
        amrex::Array4<const Set::Scalar> const& eta_patch = eta_current_mf.array(mfi);

        Set::Patch<const Set::Scalar> rho       = rho_mf.array(mfi);  // density
        Set::Patch<const Set::Scalar> M         = M_mf.array(mfi);    // momentum
        Set::Patch<const Set::Scalar> E         = E_mf.array(mfi);    // total energy (internal energy + kinetic energy) per unit volume (E/rho = e + 0.5*v^2)

        Set::Patch<const Set::Scalar> rho_solid = solid.density_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> M_solid   = solid.momentum_mf.Patch(lev,mfi);

        Set::Patch<Set::Scalar> scratch         = scratch_mf.Patch(lev,mfi);

        Set::Patch<Set::Scalar>       v         = velocity_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       p         = pressure_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       T         = temperature_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       Y         = mass_fraction_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       X         = mole_fraction_mf.Patch(lev,mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Scalar eta = invert ? 1.0-eta_patch(i,j,k)*eta_patch(i,j,k) : eta_patch(i,j,k);

            // Compute T and P primitives from mixed values
            Set::Scalar density = gas.ComputeD(rho, i, j, k);
            T(i,j,k) = gas.ComputeT(density, M(i,j,k,0), M(i,j,k,1), E(i,j,k), T(i,j,k), X, i, j, k);
            p(i,j,k) = gas.ComputeP(density, T(i,j,k), X, i, j, k);

            // Compute velocity from fluid values
            scratch(i,j,k) = (rho(i,j,k) - rho_solid(i,j,k)*(1.0 - eta))/(eta + small);
            gas.ComputeLocalFractions(scratch, Y, X, i, j, k);
            Set::Scalar density_fluid = gas.ComputeD(scratch, i, j, k);
            Set::Scalar Mx_fluid = (M(i,j,k,0) - M_solid(i,j,k,0)*(1.0 - eta))/(eta + small);
            Set::Scalar My_fluid = (M(i,j,k,1) - M_solid(i,j,k,1)*(1.0 - eta))/(eta + small);
            v(i,j,k,0) = Mx_fluid/density_fluid;
            v(i,j,k,1) = My_fluid/density_fluid;

            if (eta < small) 
            {
                v(i,j,k,0) *= eta;
                v(i,j,k,1) *= eta;

                #if AMREX_SPACEDIM == 3
                    v(i,j,k,2) *= eta;
                #endif
            }
        });
    }

    const Set::Scalar* DX = geom[lev].CellSize();
    amrex::Box domain = geom[lev].Domain();

    for (amrex::MFIter mfi(eta_iter_mf, false); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        
        // Inputs
        Set::Patch<const Set::Scalar> rho = rho_mf.array(mfi);
        Set::Patch<const Set::Scalar> E   = E_mf.array(mfi);
        Set::Patch<const Set::Scalar> M   = M_mf.array(mfi);

        // Outputs
        Set::Patch<Set::Scalar> rho_rhs = rho_rhs_mf.array(mfi);
        Set::Patch<Set::Scalar> M_rhs   = M_rhs_mf.array(mfi);
        Set::Patch<Set::Scalar> E_rhs   = E_rhs_mf.array(mfi);


        // Set::Patch<Set::Scalar>       rho_new = density_mf.Patch(lev,mfi);
        // Set::Patch<Set::Scalar>       E_new   = energy_mf.Patch(lev,mfi);
        // Set::Patch<Set::Scalar>       M_new   = momentum_mf.Patch(lev,mfi);

        Set::Patch<const Set::Scalar> rho_solid = solid.density_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> M_solid   = solid.momentum_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> E_solid   = solid.energy_mf.Patch(lev,mfi);

        Set::Patch<Set::Scalar>       omega     = vorticity_mf.Patch(lev,mfi);

        amrex::Array4<const Set::Scalar> const& eta_patch = eta_current_mf.array(mfi);
        Set::Patch<const Set::Scalar> etadot    = etadot_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> velocity  = velocity_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> T         = temperature_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> molef     = mole_fraction_mf.Patch(lev,mfi);

        Set::Patch<const Set::Scalar> m0        = m0_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> q         = q_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> _u0       = u0_mf.Patch(lev,mfi);

        amrex::Array4<Set::Scalar> const& Source = (*Source_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> eta_rhs;
        if (evolve_eta) eta_rhs = eta_rhs_mf->array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {   
            auto sten = Numeric::GetStencil(i, j, k, domain);

            Set::Scalar eta = invert ? 1.0-eta_patch(i,j,k)*eta_patch(i,j,k) : eta_patch(i,j,k);

            //Diffuse Sources
            Set::Vector grad_eta     = Numeric::Gradient(eta_patch, i, j, k, 0, DX);
            Set::Scalar grad_eta_mag = grad_eta.lpNorm<2>();
            Set::Matrix hess_eta     = Numeric::Hessian(eta_patch, i, j, k, 0, DX);
            if (invert) grad_eta *= -1.0;
            if (invert) hess_eta *= -1.0;
            
            #if AMREX_SPACEDIM == 2
                Set::Vector u            = Set::Vector(velocity(i, j, k, 0), velocity(i, j, k, 1)); // Velocity
                Set::Vector u0           = Set::Vector(_u0(i, j, k, 0), _u0(i, j, k, 1)); // Velocity
                Set::Vector q0           = Set::Vector(q(i,j,k,0), q(i,j,k,1));
            #endif

            #if AMREX_SPACEDIM == 3
                Set::Vector u            = Set::Vector(velocity(i, j, k, 0), velocity(i, j, k, 1), velocity(i, j, k, 2)); // Velocity
                Set::Vector u0           = Set::Vector(_u0(i, j, k, 0), _u0(i, j, k, 1), _u0(i, j, k, 2)); // Velocity
                Set::Vector q0           = Set::Vector(q(i,j,k,0), q(i,j,k,1), q(i,j,k,2));
            #endif

            if (evolve_eta) eta_rhs(i,j,k) = -u.dot(grad_eta);

            Set::Matrix gradM        = Numeric::Gradient(M, i, j, k, DX);
            Set::Vector gradrho      = Numeric::Gradient(rho,i,j,k,0,DX);
            Set::Matrix hess_rho     = Numeric::Hessian(rho,i,j,k,0,DX,sten);
            Set::Matrix gradu        = (gradM - u*gradrho.transpose()) / rho(i,j,k);

            if (prescribedflowmode == PrescribedFlowMode::Relative)
            {
                Set::Vector N = grad_eta / (grad_eta_mag + small);
                // Set::Vector T(N(1), -N(0));
                // u0 = N * u0(0) + T * u0(1);

                #if AMREX_SPACEDIM == 2
                    Set::Vector T(N(1), -N(0));
                    u0 = N * u0(0) + T * u0(1);
                #endif

                #if AMREX_SPACEDIM == 3
                    Set::Vector T;
                    T(0) = N(1);
                    T(1) = -N(0);
                    T(2) = 0;
                    u0 = N*u0(0) + T * u0(1);
                    // Might not be physcially accurate, need to find how to extend to 3 dimensions
                #endif
            }


            Set::Scalar mdot0 = m0(i,j,k)*grad_eta_mag;
            Set::Vector Pdot0 = Set::Vector::Zero(); // Linear momentum source term
            Set::Scalar qdot0 = q0.dot(grad_eta);

            Set::Scalar mu = gas.dynamic_viscosity(T(i,j,k), molef, i, j, k);

            // sten is necessary here because sometimes corner ghost
            // cells don't get filled
            Set::Matrix3 hess_M = Numeric::Hessian(M,i,j,k,DX);
            Set::Matrix3 hess_u = Set::Matrix3::Zero();
            for (int p = 0; p < 2; p++)
                for (int q = 0; q < 2; q++)
                    for (int r = 0; r < 2; r++)
                    {
                        hess_u(r,p,q) =
                            (hess_M(r,p,q) - gradu(r,q)*gradrho(p) - gradu(r,p)*gradrho(q) - u(r)*hess_rho(p,q))
                            / rho(i,j,k);
                    }

            Set::Vector Ldot0 = Set::Vector::Zero();
            Set::Vector div_tau = Set::Vector::Zero();
            Set::Scalar lambda = 0.0; //-2.0/3.0*mu_eff;
            for (int p = 0; p<2; p++)
                for (int q = 0; q<2; q++)
                    for (int r = 0; r<2; r++)
                        for (int s = 0; s<2; s++)
                        {
                            Ldot0(p) += 0.25 * (mu * ((p==r && q==s) + (p==s && q==r)) + lambda * (p==q && r==s)) * (u(r) - u0(r)) * hess_eta(q, s);
                            div_tau(p) += 0.5 * (mu * ((p==r && q==s) + (p==s && q==r)) + lambda * (p==q && r==s)) * (hess_u(r,q,s) + hess_u(s,q,r));

                        }

            if (use_diffuse_sources)
            {
                Source(i,j, k, 0) = mdot0;
                Source(i,j, k, 1) = Pdot0(0) - Ldot0(0);
                Source(i,j, k, 2) = Pdot0(1) - Ldot0(1);
                Source(i,j, k, 3) = qdot0;// - Ldot0(0)*v(i,j,k,0) - Ldot0(1)*v(i,j,k,1);

                // Lagrange terms to enforce no-penetration
                Source(i,j,k,1) -= lagrange*(u-u0).dot(grad_eta)*grad_eta(0);
                Source(i,j,k,2) -= lagrange*(u-u0).dot(grad_eta)*grad_eta(1);
            }
            else
            {
                Source(i,j,k,0) = 0.0;
                Source(i,j,k,1) = 0.0;
                Source(i,j,k,2) = 0.0;
                Source(i,j,k,3) = 0.0;
            }

            //Godunov flux
            //states of total fields
            const int X = 0, Y = 1;
            Solver::Local::Riemann::State state_xlo(rho, M, E, i-1, j, k, X);
            Solver::Local::Riemann::State state_x  (rho, M, E, i  , j, k, X); 
            Solver::Local::Riemann::State state_xhi(rho, M, E, i+1, j, k, X);

            Solver::Local::Riemann::State state_ylo(rho, M, E, i, j-1, k, Y);
            Solver::Local::Riemann::State state_y  (rho, M, E, i, j  , k, Y);
            Solver::Local::Riemann::State state_yhi(rho, M, E, i, j+1, k, Y);
            
            //states of solid fields
            Solver::Local::Riemann::State state_xlo_solid(rho_solid, M_solid, E_solid, i-1, j, k, X); 
            Solver::Local::Riemann::State state_x_solid  (rho_solid, M_solid, E_solid, i  , j, k, X); 
            Solver::Local::Riemann::State state_xhi_solid(rho_solid, M_solid, E_solid, i+1, j, k, X); 

            Solver::Local::Riemann::State state_ylo_solid(rho_solid, M_solid, E_solid, i, j-1, k, Y); 
            Solver::Local::Riemann::State state_y_solid  (rho_solid, M_solid, E_solid, i, j  , k, Y); 
            Solver::Local::Riemann::State state_yhi_solid(rho_solid, M_solid, E_solid, i, j+1, k, Y); 

            Solver::Local::Riemann::State state_xlo_fluid = invert ? 
                (state_xlo - (eta_patch(i-1,j,k))*state_xlo_solid) / (1.0 - eta_patch(i-1,j,k) + small) :
                (state_xlo - (1.0 - eta_patch(i-1,j,k))*state_xlo_solid) / (eta_patch(i-1,j,k) + small);
            Solver::Local::Riemann::State state_x_fluid   = invert ? 
                (state_x   - (eta_patch(i,j,k)  )*state_x_solid  )   / (1.0 - eta_patch(i,j,k)   + small): 
                (state_x   - (1.0 - eta_patch(i,j,k)  )*state_x_solid  ) / (eta_patch(i,j,k)   + small);
            Solver::Local::Riemann::State state_xhi_fluid = invert ? 
                (state_xhi - (eta_patch(i+1,j,k))*state_xhi_solid) / (1.0 - eta_patch(i+1,j,k) + small) : 
                (state_xhi - (1.0 - eta_patch(i+1,j,k))*state_xhi_solid) / (eta_patch(i+1,j,k) + small);
            Solver::Local::Riemann::State state_ylo_fluid = invert ? 
                (state_ylo - (eta_patch(i,j-1,k))*state_ylo_solid) / (1.0 - eta_patch(i,j-1,k) + small): 
                (state_ylo - (1.0 - eta_patch(i,j-1,k))*state_ylo_solid) / (eta_patch(i,j-1,k) + small);
            Solver::Local::Riemann::State state_y_fluid =   invert ? 
                (state_y   - (eta_patch(i,j,k)  )*state_y_solid  )  / (1.0 - eta_patch(i,j,k)   + small): 
                (state_y   - (1.0 - eta_patch(i,j,k)  )*state_y_solid  ) / (eta_patch(i,j,k)   + small);
            Solver::Local::Riemann::State state_yhi_fluid = invert ? 
                (state_yhi - (eta_patch(i,j+1,k))*state_yhi_solid) / (1.0 - eta_patch(i,j+1,k) + small): 
                (state_yhi - (1.0 - eta_patch(i,j+1,k))*state_yhi_solid) / (eta_patch(i,j+1,k) + small);

            Solver::Local::Riemann::Flux flux_xlo, flux_ylo, flux_xhi, flux_yhi;

            try
            {
                //lo interface fluxes
                flux_xlo = riemannsolver->Solve(state_xlo_fluid, state_x_fluid, gas, molef, i, j, k, 0, small) * eta;
                flux_ylo = riemannsolver->Solve(state_ylo_fluid, state_y_fluid, gas, molef, i, j, k, 2, small) * eta;

                //hi interface fluxes
                flux_xhi = riemannsolver->Solve(state_x_fluid, state_xhi_fluid, gas, molef, i, j, k, 1, small) * eta;
                flux_yhi = riemannsolver->Solve(state_y_fluid, state_yhi_fluid, gas, molef, i, j, k, 3, small) * eta;
            }
            catch(...)
            {
                Util::ParallelMessage(INFO,"lev=",lev);
                Util::ParallelMessage(INFO,"i=",i,"j=",j);
                Util::Abort(INFO);
            }
                

            Set::Scalar drhof_dt = 
                (flux_xlo.mass - flux_xhi.mass) / DX[0] +
                (flux_ylo.mass - flux_yhi.mass) / DX[1] +
                (use_diffuse_sources ? Source(i, j, k, 0) : 0.0);

            rho_rhs(i,j,k) = 
                // rho_new(i, j, k) = rho(i, j, k) + 
                //(
                    drhof_dt +
                    // todo add drhos_dt term if want time-evolving rhos
                    (use_diffuse_sources ? etadot(i,j,k) * (rho(i,j,k) - rho_solid(i,j,k)) / (eta + small) : 0.0)
                // ) * dt;
                ;


                
            Set::Scalar dMxf_dt =
                (flux_xlo.momentum_normal  - flux_xhi.momentum_normal ) / DX[0] +
                (flux_ylo.momentum_tangent - flux_yhi.momentum_tangent) / DX[1] +
                div_tau(0) * eta +
                g(0)*rho(i,j,k) +
                (use_diffuse_sources ? Source(i, j, k, 1) : 0.0);

            M_rhs(i,j,k,0) = 
                //M_new(i, j, k, 0) = M(i, j, k, 0) +
                // ( 
                    dMxf_dt + 
                    // todo add dMs_dt term if want time-evolving Ms
                    (use_diffuse_sources ? etadot(i,j,k)*(M(i,j,k,0) - M_solid(i,j,k,0)) / (eta + small) : 0.0)
                // ) * dt;
                ;

            Set::Scalar dMyf_dt =
                (flux_xlo.momentum_tangent - flux_xhi.momentum_tangent) / DX[0] +
                (flux_ylo.momentum_normal  - flux_yhi.momentum_normal ) / DX[1] +
                div_tau(1) * eta + 
                g(1)*rho(i,j,k) +
                (use_diffuse_sources ? Source(i, j, k, 2) : 0.0);

            M_rhs(i,j,k,1) = 
                //M_new(i, j, k, 1) = M(i, j, k, 1) +
                //( 
                    dMyf_dt +
                    // todo add dMs_dt term if want time-evolving Ms
                    (use_diffuse_sources ? etadot(i,j,k)*(M(i,j,k,1) - M_solid(i,j,k,1)) / (eta+small) : 0.0)
                // )*dt;
                ;

            Set::Scalar dEf_dt =
                (flux_xlo.energy - flux_xhi.energy) / DX[0] +
                (flux_ylo.energy - flux_yhi.energy) / DX[1] +
                (use_diffuse_sources ? Source(i, j, k, 3) : 0.0);

            E_rhs(i,j,k) = 
            // E_new(i, j, k) = E(i, j, k) + 
            //     ( 
                    dEf_dt +
                    // todo add dEs_dt term if want time-evolving Es
                    (use_diffuse_sources ? etadot(i,j,k)*(E(i,j,k) - E_solid(i,j,k)) / (eta+small) : 0.0)
                // ) * dt;
                ;
            
#ifdef AMREX_DEBUG
            if ((rho_rhs(i,j,k) != rho_rhs(i,j,k)) ||
                (M_rhs(i,j,k,0) != M_rhs(i,j,k,0)) ||
                (M_rhs(i,j,k,1) != M_rhs(i,j,k,1)) ||
                (E_rhs(i,j,k) != E_rhs(i,j,k)))
            {
                Util::ParallelMessage(INFO,"rho_rhs=",rho_rhs(i,j,k));
                Util::ParallelMessage(INFO,"Mx_rhs=",M_rhs(i,j,k,0));
                Util::ParallelMessage(INFO,"Mx_rhs=",M_rhs(i,j,k,1));
                Util::ParallelMessage(INFO,"E_rhs=",E_rhs(i,j,k));

                Util::ParallelMessage(INFO,"lev=",lev);
                Util::ParallelMessage(INFO,"i=",i," j=",j);
                Util::ParallelMessage(INFO,"drhof_dt ",drhof_dt); // dies
                Util::ParallelMessage(INFO,"flux_xlo.mass ",flux_xlo.mass);
                Util::ParallelMessage(INFO,"flux_xhi.mass ",flux_xhi.mass); // dies, depends on state_xx, state_xhi, state_x_solid, state_xhi_solid, eta, small
                Util::ParallelMessage(INFO,"flux_ylo.mass ",flux_ylo.mass);
                Util::ParallelMessage(INFO,"flux_xhi.mass ",flux_yhi.mass);
                Util::ParallelMessage(INFO,"eta ",eta);
                Util::ParallelMessage(INFO,"etadot ",etadot(i,j,k));
                Util::ParallelMessage(INFO,"Source ",Source(i,j,k,0));
                Util::ParallelMessage(INFO,"state_x ",state_x); // <<<<
                Util::ParallelMessage(INFO,"state_y ",state_y);
                Util::ParallelMessage(INFO,"state_x_solid ",state_x_solid); // <<<<
                Util::ParallelMessage(INFO,"state_y_solid ",state_y_solid);
                Util::ParallelMessage(INFO,"state_xhi ",state_xhi); // <<<<
                Util::ParallelMessage(INFO,"state_yhi ",state_yhi);
                Util::ParallelMessage(INFO,"state_xhi_solid ",state_xhi_solid);
                Util::ParallelMessage(INFO,"state_yhi_solids ",state_yhi_solid);
                Util::ParallelMessage(INFO,"state_xlo ",state_xlo);
                Util::ParallelMessage(INFO,"state_ylo ",state_ylo);
                Util::ParallelMessage(INFO,"state_xlo_solid ",state_xlo_solid);
                Util::ParallelMessage(INFO,"state_ylo_solid ",state_ylo_solid);

                Util::ParallelMessage(INFO,"Mx_solid ",M_solid(i,j,k,0));
                Util::ParallelMessage(INFO,"My_solid ",M_solid(i,j,k,1));
                Util::ParallelMessage(INFO,"small ",small);
                Util::ParallelMessage(INFO,"Mx ",M(i,j,k,0));
                Util::ParallelMessage(INFO,"My ",M(i,j,k,1));
                Util::ParallelMessage(INFO,"dMx/dt ",dMxf_dt);
                Util::ParallelMessage(INFO,"dMy/dt ",dMyf_dt);


                Util::Message(INFO,flux_xlo.momentum_tangent);
                Util::Message(INFO,flux_xhi.momentum_tangent);
                Util::Message(INFO,DX[0]);
                Util::Message(INFO,flux_ylo.momentum_normal);
                Util::Message(INFO,flux_yhi.momentum_normal);
                Util::Message(INFO,DX[1]);
                Util::Message(INFO,div_tau);
                Util::Message(INFO,Source(i, j, k, 2));
                
                Util::Message(INFO,hess_eta);
                Util::Message(INFO,velocity(i,j,k,0));
                Util::Message(INFO,velocity(i,j,k,1));

                Util::Exception(INFO);
            }
#endif



            // todo - may need to move this for higher order schemes...
            omega(i, j, k) = eta * (gradu(1,0) - gradu(0,1));
        });
    }
}

void Hydro::Regrid(int lev, Set::Scalar /* time */)
{
    BL_PROFILE("Integrator::Hydro::Regrid");
    if (UsesDiffuseBoundarySources()) Source_mf[lev]->setVal(0.0);
    if (lev < finest_level) return;

    Util::Message(INFO, "Regridding on level", lev);
}//end regrid

//void Hydro::TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, Set::Scalar time, int ngrow)
void Hydro::TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar, int)
{
    BL_PROFILE("Integrator::Flame::TagCellsForRefinement");

    const Set::Scalar* DX = geom[lev].CellSize();
    Set::Scalar dr = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));

    // Eta criterion for refinement
    for (amrex::MFIter mfi(*(*eta_mf)[lev], true); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<char> const& tags = a_tags.array(mfi);
        amrex::Array4<const Set::Scalar> const& eta = (*(*eta_mf)[lev]).array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Set::Vector grad_eta = Numeric::Gradient(eta, i, j, k, 0, DX);
            if (grad_eta.lpNorm<2>() * dr * 2 > eta_refinement_criterion) tags(i, j, k) = amrex::TagBox::SET;
        });
    }

    // Vorticity criterion for refinement
    for (amrex::MFIter mfi(*vorticity_mf[lev], true); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<char> const& tags = a_tags.array(mfi);
        amrex::Array4<const Set::Scalar> const& omega = (*vorticity_mf[lev]).array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            auto sten = Numeric::GetStencil(i, j, k, bx);
            Set::Vector grad_omega = Numeric::Gradient(omega, i, j, k, 0, DX, sten);
            if (grad_omega.lpNorm<2>() * dr * 2 > omega_refinement_criterion) tags(i, j, k) = amrex::TagBox::SET;
        });
    }

    // Gradu criterion for refinement
    for (amrex::MFIter mfi(*velocity_mf[lev], true); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<char> const& tags = a_tags.array(mfi);
        amrex::Array4<const Set::Scalar> const& v = (*velocity_mf[lev]).array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            auto sten = Numeric::GetStencil(i, j, k, bx);
            Set::Matrix grad_u = Numeric::Gradient(v, i, j, k, DX, sten);
            if (grad_u.lpNorm<2>() * dr * 2 > gradu_refinement_criterion) tags(i, j, k) = amrex::TagBox::SET;
        });
    }

    // Pressure criterion for refinement
    for (amrex::MFIter mfi(*pressure_mf[lev], true); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<char> const& tags = a_tags.array(mfi);
        amrex::Array4<const Set::Scalar> const& p = (*pressure_mf[lev]).array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            auto sten = Numeric::GetStencil(i, j, k, bx);
            Set::Vector grad_p = Numeric::Gradient(p, i, j, k, 0, DX, sten);
            if (grad_p.lpNorm<2>() * dr * 2 > p_refinement_criterion) tags(i, j, k) = amrex::TagBox::SET;
        });
    }

    // Density criterion for refinement
    for (amrex::MFIter mfi(*density_mf[lev], true); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<char> const& tags = a_tags.array(mfi);
        amrex::Array4<const Set::Scalar> const& rho = (*density_mf[lev]).array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            auto sten = Numeric::GetStencil(i, j, k, bx);
            Set::Vector grad_rho = Numeric::Gradient(rho, i, j, k, 0, DX, sten);
            if (grad_rho.lpNorm<2>() * dr * 2 > rho_refinement_criterion) tags(i, j, k) = amrex::TagBox::SET;
        });
    }

}

}
