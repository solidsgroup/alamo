
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
#include "AMReX_FabArrayUtility.H"
#include <string>

#include "Model/Gas/Gas.H"
#include "Model/Gas/Thermo/Thermo.H"
#include "Model/Gas/Thermo/CpConstant.H"
#include "Model/Gas/Transport/Transport.H"
#include "Model/Gas/Transport/Mixture_Averaged.H"
#include "Model/Gas/EOS/EOS.H"
#include "Model/Gas/EOS/CPG.H"

#include "Model/Chemistry/Chemistry.H"
#include "Model/Chemistry/Frozen.H"
#include "Model/Chemistry/Equilibrium.H"

namespace Integrator
{

namespace
{
    // Check that all partial densities are finite and greater than or equal to 0.0.
    // If returns true: negative densities are clipped and remaining densities are scaled to
    // conserve the original density
    // If returns false: no changes are made to the partial densities, and run quietly continues
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool TryProjectSpeciesDensities(
        std::array<Set::Scalar, NSPECIES>& rhoY,
        Set::Scalar& raw_density,
        Set::Scalar& positive_density,
        Set::Scalar& scale,
        int& nonfinite_species)
    {
        raw_density = 0.0;
        positive_density = 0.0;
        scale = 0.0;
        nonfinite_species = -1;
        for (int n = 0; n < NSPECIES; ++n)
        {
            if (!std::isfinite(rhoY[n]))
            {
                nonfinite_species = n;
                return false;
            }
            raw_density += rhoY[n];
            if (rhoY[n] > 0.0) positive_density += rhoY[n];
        }

        if (!std::isfinite(raw_density) || raw_density <= 0.0 || positive_density <= 0.0)
        {
            return false;
        }

        scale = raw_density / positive_density;
        if (!std::isfinite(scale) || scale <= 0.0)
        {
            return false;
        }

        for (int n = 0; n < NSPECIES; ++n)
        {
            rhoY[n] = (rhoY[n] > 0.0) ? rhoY[n] * scale : 0.0;
        }
        return true;
    }

    // Helper function to test if densities are scaled or not
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool TryProjectSpeciesDensities(std::array<Set::Scalar, NSPECIES>& rhoY)
    {
        Set::Scalar raw_density = 0.0;
        Set::Scalar positive_density = 0.0;
        Set::Scalar scale = 0.0;
        int nonfinite_species = -1;
        return TryProjectSpeciesDensities(
            rhoY, raw_density, positive_density, scale, nonfinite_species);
    }

    // Noisy density projection that exits loudly if it fails
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void ProjectSpeciesDensities(std::array<Set::Scalar, NSPECIES>& rhoY, int i, int j, int k)
    {
        Set::Scalar raw_density = 0.0;
        Set::Scalar positive_density = 0.0;
        Set::Scalar scale = 0.0;
        int nonfinite_species = -1;
        if (TryProjectSpeciesDensities(rhoY, raw_density, positive_density, scale, nonfinite_species)) return;

        if (nonfinite_species >= 0)
        {
            Util::Abort(INFO, "Non-finite species density before projection at (",
                i, ",", j, ",", k, "), species ", nonfinite_species, ": ", rhoY[nonfinite_species]);
        }
        if (std::isfinite(raw_density) && raw_density > 0.0 && positive_density > 0.0)
        {
            Util::Abort(INFO, "Invalid species projection scale at (",
                i, ",", j, ",", k, "): scale=", scale,
                ", raw density=", raw_density, ", positive density=", positive_density);
        }
        Util::Abort(INFO, "Cannot project invalid species state at (",
            i, ",", j, ",", k, "): raw density=", raw_density,
            ", positive density=", positive_density);
    }

    // Stencil to apply central differencing even at bounds of validbox
    // This ensures ghost cells are correctly utilized to enforce BCs
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    std::array<Numeric::StencilType, AMREX_SPACEDIM> BCGhostStencil()
    {
        return { AMREX_D_DECL(  Numeric::StencilType::Central,
                                Numeric::StencilType::Central,
                                Numeric::StencilType::Central) };
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int OutsideDirection(
        int i, int j, int k,
        int xlo, int xhi,
        int ylo, int yhi,
        int zlo, int zhi)
    {
        int dir = -1;
        if (i < xlo || i > xhi) dir = 0;
#if AMREX_SPACEDIM >= 2
        if (j < ylo || j > yhi)
        {
            if (dir >= 0) return -2;
            dir = 1;
        }
#endif
#if AMREX_SPACEDIM == 3
        if (k < zlo || k > zhi)
        {
            if (dir >= 0) return -2;
            dir = 2;
        }
#else
        (void)k; (void)zlo; (void)zhi;
#endif
        return dir;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int OutsideDistance(
        int i, int j, int k,
        int xlo, int xhi,
        int ylo, int yhi,
        int zlo, int zhi,
        int dir)
    {
        if (dir == 0) return i < xlo ? xlo - i : i - xhi;
#if AMREX_SPACEDIM >= 2
        if (dir == 1) return j < ylo ? ylo - j : j - yhi;
#endif
#if AMREX_SPACEDIM == 3
        if (dir == 2) return k < zlo ? zlo - k : k - zhi;
#else
        (void)k; (void)zlo; (void)zhi;
#endif
        return 0;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool IsValidAxisGhostOrCornerWithin(
        int i, int j, int k,
        int xlo, int xhi,
        int ylo, int yhi,
        int zlo, int zhi,
        int axis_depth,
        int corner_depth)
    {
        int outside = 0;
        int max_distance = 0;
        if (i < xlo || i > xhi)
        {
            ++outside;
            const int distance = i < xlo ? xlo - i : i - xhi;
            if (distance > max_distance) max_distance = distance;
        }
#if AMREX_SPACEDIM >= 2
        if (j < ylo || j > yhi)
        {
            ++outside;
            const int distance = j < ylo ? ylo - j : j - yhi;
            if (distance > max_distance) max_distance = distance;
        }
#endif
#if AMREX_SPACEDIM == 3
        if (k < zlo || k > zhi)
        {
            ++outside;
            const int distance = k < zlo ? zlo - k : k - zhi;
            if (distance > max_distance) max_distance = distance;
        }
#else
        (void)k; (void)zlo; (void)zhi;
#endif
        if (outside == 0) return true;
        if (outside == 1) return max_distance <= axis_depth;
        return max_distance <= corner_depth;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool IsValidOrAxisGhostWithin(
        int i, int j, int k,
        int xlo, int xhi,
        int ylo, int yhi,
        int zlo, int zhi,
        int max_depth)
    {
        const int dir = OutsideDirection(i, j, k, xlo, xhi, ylo, yhi, zlo, zhi);
        if (dir == -2) return false;
        if (dir < 0) return true;
        return OutsideDistance(i, j, k, xlo, xhi, ylo, yhi, zlo, zhi, dir) <= max_depth;
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
        // Gas model (Thermo, Transport, and EOS)
        pp.queryclass<Model::Gas::Gas>("gas", value.gas);
        Util::Message(INFO, "Model::Gas::EOS       = ", value.gas.eos.model_name());
        Util::Message(INFO, "Model::Gas::Thermo    = ", value.gas.thermo.model_name());
        Util::Message(INFO, "Model::Gas::Transport = ", value.gas.transport.model_name());
        Util::Message(INFO, "NSPECIES = ", NSPECIES);

        // Riemann solver
        pp.select_default<  Solver::Local::Riemann::Roe,
                            Solver::Local::Riemann::HLLE,
                            Solver::Local::Riemann::HLLC>("solver",value.riemannsolver);

        // Chemistry Model
        pp.select<  Model::Chemistry::Frozen,
                    Model::Chemistry::Equilibrium,
                    Model::Chemistry::FiniteRate,
                    Model::Chemistry::Rocfire>("chemistry",value.chemistry);

        // Advection Scheme
        pp.select<Numeric::Advect::MUSCL,
                  Numeric::Advect::Upwind,
                  Numeric::Advect::Centered,
                  Numeric::Advect::QUICK,
                  Numeric::Advect::WENO5>("advection",value.advect);

        std::string flux_scheme_str;
        pp.query_validate("flux_scheme", flux_scheme_str, {"riemann","advect"});
        if (flux_scheme_str == "riemann") value.flux_scheme = FluxScheme::Riemann;
        else if (flux_scheme_str == "advect") value.flux_scheme = FluxScheme::Advect;

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
        // density-based refinement
        pp.query_default("temp_refinement_criterion", value.temp_refinement_criterion, 1e100);

        pp_forbid("gamma", "replaced by gas->gamma(...)"); // gamma for gamma law
        pp_query_required("cfl", value.cfl); // cfl condition
        pp_query_default("cfl_v", value.cfl_v,1E100); // cfl condition
        pp_forbid("mu", "replaced with gas->dynamic_viscosity(...)"); // linear viscosity coefficient

        // Boundary condition for density
        pp.select_default<BC::Constant,BC::Expression>("density.bc",value.density_bc, NSPECIES);
        // Boundary condition for energy
        pp.select_default<BC::Constant,BC::Expression>("energy.bc",value.energy_bc,1);
        // Boundary condition for momentum
        pp.select_default<BC::Constant,BC::Expression>("momentum.bc",value.momentum_bc,AMREX_SPACEDIM);

        if (!value.managed)
        {
            // Boundary condition for phase field order parameter
            pp.select_default<BC::Constant,BC::Expression>("pf.eta.bc",value.eta_bc,1);
        }

        pp_query_default("small",value.small,1E-8); // small regularization value
        pp_query_default("cutoff",value.cutoff,-1E100); // cutoff value
        pp_query_default("lagrange",value.lagrange,0.0); // lagrange no-penetration factor
        pp_query_default("details",value.details,false); // save detailed data (viscosity, heat conductivity, etc.)
        // Physical solid/gas interfacial contact conductance [W/m^2/K] - see
        // AdvanceSolidEnergy. Only used when managed (Flame-driven); required
        // rather than defaulted so switching to this eps-independent
        // formulation can't silently go unnoticed by existing inputs.
        if (value.managed) pp_query_required("solid.h_interface", value.solid.h_interface, Unit::HeatTransferCoefficient());
    }
    // Register FabFields:
    {
        int nghost = value.flux_scheme == FluxScheme::Advect ? value.advect.NGhost() : 1;

        if (!value.managed)
        {
            value.eta_mf = new Set::Field<Set::Scalar>();
            value.eta_old_mf = new Set::Field<Set::Scalar>();
            value.RegisterNewFab(*value.eta_mf,     value.eta_bc, 1, nghost, "eta",     true, true);
            value.RegisterNewFab(*value.eta_old_mf, value.eta_bc, 1, nghost, "eta_old", true, true);
        }
        value.RegisterNewFab(value.etadot_mf,  value.eta_bc, 1, nghost, "etadot",  true, false);

        value.RegisterNewFab(value.density_mf,     value.density_bc, NSPECIES, nghost, "density",     true , true);
        value.RegisterNewFab(value.density_old_mf, value.density_bc, NSPECIES, nghost, "density_old", false, true);

        value.RegisterNewFab(value.energy_mf,     value.energy_bc, 1, nghost, "energy",      true ,true);
        value.RegisterNewFab(value.energy_old_mf, value.energy_bc, 1, nghost, "energy_old" , false, true);

        std::vector<std::string> vector_suffix = {AMREX_D_DECL("x","y","z")};
        value.RegisterNewFab(value.momentum_mf,     value.momentum_bc, AMREX_SPACEDIM, nghost, "momentum",     true ,true, vector_suffix);
        value.RegisterNewFab(value.momentum_old_mf, value.momentum_bc, AMREX_SPACEDIM, nghost, "momentum_old", false, true);

        value.RegisterNewFab(value.pressure_mf,  &value.bc_nothing, 1, nghost, "pressure",  true, false);
        value.RegisterNewFab(value.temperature_mf,  &value.bc_nothing, 1, nghost, "temperature",  true, false);
        value.RegisterNewFab(value.velocity_mf,  &value.bc_nothing, AMREX_SPACEDIM, nghost, "velocity",  true, false, vector_suffix);
        #if AMREX_SPACEDIM == 2
        value.RegisterNewFab(value.vorticity_mf, &value.bc_nothing, 1, nghost, "vorticity", true, false);
        #elif AMREX_SPACEDIM == 3
        value.RegisterNewFab(value.vorticity_mf, &value.bc_nothing, 3, nghost, "vorticity", true, false);
        #endif

        value.RegisterNewFab(value.m0_mf,           &value.bc_nothing, NSPECIES, 0, "m0",  true, false);
        value.RegisterNewFab(value.u0_mf,           &value.bc_nothing, AMREX_SPACEDIM, 0, "u0",  true, false, vector_suffix);
        value.RegisterNewFab(value.q_mf,            &value.bc_nothing, AMREX_SPACEDIM, 0, "q",   true, false, vector_suffix);

        value.neumann_bc_D = new BC::Constant(BC::Constant::ZeroNeumann(AMREX_SPACEDIM));
        value.neumann_bc_1 = new BC::Constant(BC::Constant::ZeroNeumann(1));
        value.neumann_bc_N = new BC::Constant(BC::Constant::ZeroNeumann(NSPECIES));

        value.RegisterNewFab(value.solid.density_mf,  value.neumann_bc_N, NSPECIES, nghost, "solid.density", true, false);
        value.RegisterNewFab(value.solid.momentum_mf, value.neumann_bc_D, AMREX_SPACEDIM, nghost, "solid.momentum", true, false, vector_suffix);
        value.RegisterNewFab(value.solid.energy_mf,   value.neumann_bc_1, 1, nghost, "solid.energy",   true, false);
        value.RegisterNewFab(value.solid.rho_phys_mf, value.neumann_bc_1, 1, nghost, "solid.rho_phys", true, false);
        value.RegisterNewFab(value.solid.cp_mf,       value.neumann_bc_1, 1, nghost, "solid.cp",       true, false);
        value.RegisterNewFab(value.solid.k_mf,        value.neumann_bc_1, 1, nghost, "solid.k",        true, false);
        value.RegisterNewFab(value.solid.laser_mf,    value.neumann_bc_1, 1, nghost, "solid.laser",    true, false);

        value.RegisterNewFab(value.Source_mf, &value.bc_nothing, NSPECIES+AMREX_SPACEDIM+1, 0, "Source", true, false);

        value.RegisterNewFab(value.mass_fraction_mf,    &value.bc_nothing, NSPECIES, nghost, "mass_fraction",     true , true);
        value.RegisterNewFab(value.mole_fraction_mf,    &value.bc_nothing, NSPECIES, nghost, "mole_fraction",     true , true);
        value.RegisterNewFab(value.scratch_mf,          &value.bc_nothing, NSPECIES, nghost, "scratch",           false, false);

        if ( value.details )
        {
            value.RegisterNewFab(value.viscosity_mf, &value.bc_nothing, 1, nghost, "viscosity", true, true);
            value.RegisterNewFab(value.thermal_conductivity_coeff_mf, &value.bc_nothing, 1, nghost, "thermal_conductivity_coeff", true, true);
            value.RegisterNewFab(value.diffusion_coeff_mf, &value.bc_nothing, NSPECIES, nghost, "diffusion_coeff", true, true);
            value.RegisterNewFab(value.wdot_mf, &value.bc_nothing, NSPECIES, nghost, "wdot", true, true);
            value.RegisterNewFab(value.qdot_mf, &value.bc_nothing, 1, nghost, "qdot", true, true);
        }
    }

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

    solid.density_ic ->Initialize(lev, solid.density_mf,  0.0);
    solid.momentum_ic->Initialize(lev, solid.momentum_mf, 0.0);
    solid.energy_ic  ->Initialize(lev, solid.energy_mf,   0.0);

    if (!managed)
    {
        eta_bc->define(geom[lev]);
        eta_bc->FillBoundary(*(*eta_mf)[lev], 0, 1, 0.0, 0);
        eta_bc->FillBoundary(*(*eta_old_mf)[lev], 0, 1, 0.0, 0);
    }
    neumann_bc_N->define(geom[lev]);
    neumann_bc_D->define(geom[lev]);
    neumann_bc_1->define(geom[lev]);
    neumann_bc_N->FillBoundary(*solid.density_mf[lev], 0, NSPECIES, 0.0, 0);
    neumann_bc_D->FillBoundary(*solid.momentum_mf[lev], 0, AMREX_SPACEDIM, 0.0, 0);
    neumann_bc_1->FillBoundary(*solid.energy_mf[lev], 0, 1, 0.0, 0);

    ic_m0            ->Initialize(lev, m0_mf, 0.0);
    ic_u0            ->Initialize(lev, u0_mf, 0.0);
    ic_q             ->Initialize(lev, q_mf,  0.0);

    Source_mf[lev]   ->setVal(0.0);

    // momentum_mf/energy_mf (the fluid conserved fields) are never set from an IC
    // directly - they are only ever populated by Mix(), which derives them from
    // density_ic/velocity_ic/pressure_ic. Standalone (!managed) Hydro always ran
    // Mix() here so momentum/energy were valid from t=0. The managed (Flame) path
    // used to skip this and only mix on the first Advance() call, one full
    // timestep late, leaving momentum_mf/energy_mf at zero at t=0.
    if (managed) { if (lev >= (int)mixed.size()) mixed.push_back(false); }
    Mix(lev);
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
            Set::Scalar eta = invert ? 1.0-eta_patch(i,j,k): eta_patch(i,j,k);

            // Initially compute primitives (T,P,u) from given initial conditions
            // But from then on, compute them from mixed values to avoid zero T conditions
            // Except velocity - keep velocity from fluid values only
            Set::Scalar density = gas.ComputeLocalFractions(rho, Y, X, i,j,k); // Get local mole/mass fractions from fluid densities

            T(i,j,k) = gas.ComputeT_from_primitives(p(i,j,k), density, X, i, j, k);
            #if AMREX_SPACEDIM == 2
            Set::Scalar E_fluid = gas.ComputeE(density, density*v(i,j,k,0), density*v(i,j,k,1), T(i,j,k), X, i, j, k);
            #elif AMREX_SPACEDIM == 3
            Set::Scalar E_fluid = gas.ComputeE(density, density*v(i,j,k,0), density*v(i,j,k,1), density*v(i,j,k,2), T(i,j,k), X, i, j, k);
            #endif

            // Mix
            M(i, j, k, 0) = density*v(i, j, k, 0)*eta +  M_solid(i, j, k, 0)*(1.0-eta);
            M(i, j, k, 1) = density*v(i, j, k, 1)*eta +  M_solid(i, j, k, 1)*(1.0-eta);
            M_old(i, j, k, 0) = M(i, j, k, 0);
            M_old(i, j, k, 1) = M(i, j, k, 1);

            #if AMREX_SPACEDIM == 3
            M(i, j, k, 2) = density*v(i, j, k, 2)*eta +  M_solid(i, j, k, 2)*(1.0-eta);
            M_old(i, j, k, 2) = M(i, j, k, 2);
            #endif


            for (int n=0; n<NSPECIES; ++n)
            {
                rho(i, j, k, n) = eta * rho(i, j, k, n) + (1.0 - eta) * rho_solid(i, j, k, n);
                rho_old(i, j, k, n) = rho(i, j, k, n);
            }

            E(i, j, k) = E_fluid*eta + E_solid(i,j,k)*(1.0-eta);
            E_old(i, j, k) = E(i, j, k);
            //Util::Message(INFO,"Energy: ", E(i,j,k), " Pressure: ", p(i,j,k), " Temp: ", T(i,j,k), " Density: ",density, " R: ", gas.R(X,i,j,k), " MW: ", gas.GetMW(X,i,j,k), " Rg: ", Set::Constant::Rg);

        });
        //Util::Abort(INFO);
    }
    if (managed) { if (lev < (int)mixed.size()) mixed[lev] = true; }
    c_max = 0.0;
    vx_max = 0.0;
    vy_max = 0.0;
    vz_max = 0.0;
}

void Hydro::UpdateEta(int lev, Set::Scalar time)
{
    Util::Assert(INFO,TEST(!managed),"Should override this if Hydro is managed!");
    eta_ic->Initialize(lev, *eta_mf, time);
    eta_bc->define(geom[lev]);
    eta_bc->FillBoundary(*(*eta_mf)[lev], 0, 1, time, 0);
}

void Hydro::UpdateFluxes(int /*lev*/, Set::Scalar /*time*/, Set::Scalar /*dt*/)
{
    Util::Assert(INFO,TEST(!managed),"Should override this if Hydro is managed!");
}

// When cutoff is configured, keep eta <= cutoff cells on the solid state and
// reconstruct the retained mixed cells from the fluid/solid blend.
void Hydro::ApplyCutoffToConserved(int lev, amrex::MultiFab& rho_mf, amrex::MultiFab& M_mf, amrex::MultiFab& E_mf, bool include_ghost, bool use_old_eta)
{
    if (!(cutoff >= 0.0 && cutoff < 1.0)) return;

    Set::Field<Set::Scalar>* eta_field = use_old_eta ? eta_old_mf : eta_mf;
    for (amrex::MFIter mfi(rho_mf, true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = include_ghost ? mfi.growntilebox() : mfi.tilebox();

        Set::Patch<const Set::Scalar> eta_patch = eta_field->Patch(lev,mfi);
        Set::Patch<const Set::Scalar> rho_solid = solid.density_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> M_solid   = solid.momentum_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> E_solid   = solid.energy_mf.Patch(lev,mfi);

        amrex::Array4<Set::Scalar> const& rho = rho_mf.array(mfi);
        amrex::Array4<Set::Scalar> const& M   = M_mf.array(mfi);
        amrex::Array4<Set::Scalar> const& E   = E_mf.array(mfi);

        Set::Scalar cutoff_local = cutoff;
        Set::Scalar small_local = small;
        bool invert_local = invert;
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Scalar eta = invert_local ? 1.0 - eta_patch(i,j,k): eta_patch(i,j,k);
            if (eta < 0.0) eta = 0.0;
            if (eta > 1.0) eta = 1.0;

            auto set_solid_state = [&]() AMREX_GPU_DEVICE
            {
                for (int n = 0; n < NSPECIES; ++n)
                {
                    rho(i,j,k,n) = rho_solid(i,j,k,n);
                }
                M(i,j,k,0) = M_solid(i,j,k,0);
                M(i,j,k,1) = M_solid(i,j,k,1);
                #if AMREX_SPACEDIM == 3
                M(i,j,k,2) = M_solid(i,j,k,2);
                #endif
                E(i,j,k,0) = E_solid(i,j,k,0);
            };

            if (!std::isfinite(eta) || eta <= cutoff_local)
            {
                set_solid_state();
                return;
            }

            std::array<Set::Scalar, NSPECIES> rhoY_fluid;
            for (int n = 0; n < NSPECIES; ++n)
            {
                rhoY_fluid[n] =
                    (rho(i,j,k,n) - (1.0 - eta) * rho_solid(i,j,k,n)) / (eta + small_local);
            }

            if (!TryProjectSpeciesDensities(rhoY_fluid))
            {
                set_solid_state();
                return;
            }

            Set::Scalar density_fluid = 0.0;
            for (int n = 0; n < NSPECIES; ++n) density_fluid += rhoY_fluid[n];

            const Set::Scalar Mx_fluid =
                (M(i,j,k,0) - (1.0 - eta) * M_solid(i,j,k,0)) / (eta + small_local);
            const Set::Scalar My_fluid =
                (M(i,j,k,1) - (1.0 - eta) * M_solid(i,j,k,1)) / (eta + small_local);
            #if AMREX_SPACEDIM == 3
            const Set::Scalar Mz_fluid =
                (M(i,j,k,2) - (1.0 - eta) * M_solid(i,j,k,2)) / (eta + small_local);
            #endif
            const Set::Scalar E_fluid =
                (E(i,j,k,0) - (1.0 - eta) * E_solid(i,j,k,0)) / (eta + small_local);
            #if AMREX_SPACEDIM == 2
            const Set::Scalar kinetic_fluid =
                0.5 * (Mx_fluid * Mx_fluid + My_fluid * My_fluid) / (density_fluid + small_local);
            #elif AMREX_SPACEDIM == 3
            const Set::Scalar kinetic_fluid =
                0.5 * (Mx_fluid * Mx_fluid + My_fluid * My_fluid + Mz_fluid * Mz_fluid) / (density_fluid + small_local);
            #endif

            if (!std::isfinite(Mx_fluid) || !std::isfinite(My_fluid) ||
                #if AMREX_SPACEDIM == 3
                !std::isfinite(Mz_fluid) ||
                #endif
                !std::isfinite(E_fluid) || !std::isfinite(kinetic_fluid) ||
                E_fluid <= kinetic_fluid)
            {
                set_solid_state();
                return;
            }

            for (int n = 0; n < NSPECIES; ++n)
            {
                rho(i,j,k,n) = eta * rhoY_fluid[n] + (1.0 - eta) * rho_solid(i,j,k,n);
            }
        });
    }
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
    amrex::ParallelDescriptor::ReduceRealMax(vz_max);

    Set::Scalar new_timestep = cfl / (
        AMREX_D_TERM( (c_max + vx_max) / DX[0],
                    + (c_max + vy_max) / DX[1],
                    + (c_max + vz_max) / DX[2]));

    SetTimestep(new_timestep);
}

// Advances solid.energy_mf via its own domain-wide conduction. This is the sole
// update path for solid.energy_mf after its initial IC seed - Flame no longer
// rewrites it, so there is no round trip through Hydro's own T. Since
// temperature_mf here still holds the previous step's value (this runs before
// this step's RK/RHS stages recompute it), the scheme is a standard
// one-step-lagged explicit update, not an algebraic loop.
//
// Two physically distinct mechanisms, kept separate on purpose:
//
//  1. Bulk conduction WITHIN the solid: div(phi_s * alpha_solid * grad(T_solid)),
//     diffusing the solid-only caloric temperature (T_solid = E_solid/(rho_phys*cp)
//     + T_ref), never the shared/blended temperature_mf. Differentiating the
//     blended field here was the bug found in practice: at the diffuse
//     interface, grad(temperature_mf) inherits however steep the (separately
//     unstable) reacting gas temperature happens to be, and even a tiny
//     physical alpha multiplying an unbounded gradient produced a growing,
//     runaway solid-energy drift that fed back into the fluid-side temperature
//     reconstruction. T_solid is smooth and non-reactive, so its own gradient
//     stays physically bounded regardless of what the gas is doing.
//
//  2. Interfacial exchange with the gas: a bounded Robin-type flux using the
//     actual (un-blended) T_gas - T_solid difference,
//     q_interface = h_interface * |grad(phi_s)| * (T_gas - T_solid), where
//     h_interface (solid.h_interface) is a physical contact conductance
//     [W/m^2/K], independent of the phase field's diffuse width eps. This is
//     what lets heat genuinely conduct in from the gas side. |grad(phi_s)|
//     (first power, not squared) is the standard phase-field surface delta
//     function: its integral across the interface is exactly 1 for any eps
//     (fundamental theorem of calculus, since phi_s runs from 1 to 0), so
//     as eps->0 this term converges to a genuine sharp-interface Robin
//     condition h_interface*(T_gas(0)-T_solid(0)) for a fixed h_interface.
//     (An earlier version used |grad(phi_s)|^2 with k_solid, i.e. an
//     eps-dependent effective h=k_solid/(eps*sqrt(pi)) - that term's areal
//     integral diverges as eps->0, so refining eps silently strengthened
//     the coupling every time rather than converging to a fixed answer;
//     confirmed empirically via a mesh/eps refinement study that showed the
//     error growing, not shrinking, with refinement.)
void Hydro::AdvanceSolidEnergy(int lev, Set::Scalar /*time*/, Set::Scalar dt)
{
    const Set::Scalar* DX = geom[lev].CellSize();
    const amrex::BoxArray &ba = energy_mf[lev]->boxArray();
    const amrex::DistributionMapping &dm = energy_mf[lev]->DistributionMap();
    amrex::MultiFab alpha_solid_mf(ba, dm, 1, 1);
    amrex::MultiFab T_solid_mf(ba, dm, 1, 1);

    // T_solid's energy datum must match the active gas EOS - see Hydro::RHS for
    // why (same reasoning: E_solid and the gas energy must share a zero point).
    const Set::Scalar T_ref_energy =
        (std::string(gas.eos.model_name()) == "tpg")
            ? Model::Gas::EOS::sensible_reference_temperature : 0.0;

    for (amrex::MFIter mfi(*(*eta_mf)[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(1);
        Set::Patch<const Set::Scalar> rho_phys = solid.rho_phys_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> cp_solid = solid.cp_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> k_solid  = solid.k_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> E_solid  = solid.energy_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       alpha    = alpha_solid_mf.array(mfi);
        Set::Patch<Set::Scalar>       T_solid  = T_solid_mf.array(mfi);
        Set::Scalar small_local = small;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            alpha(i,j,k) = k_solid(i,j,k) / (rho_phys(i,j,k)*cp_solid(i,j,k) + small_local);
            T_solid(i,j,k) = T_ref_energy + E_solid(i,j,k) / (rho_phys(i,j,k)*cp_solid(i,j,k) + small_local);
        });
    }
    alpha_solid_mf.FillBoundary(geom[lev].periodicity());
    T_solid_mf.FillBoundary(geom[lev].periodicity());

    bool invert_local = invert;
    Set::Scalar eta_cutoff_local = (cutoff >= 0.0 && cutoff < 1.0) ? cutoff : small;
    Set::Scalar h_interface_local = solid.h_interface;

    for (amrex::MFIter mfi(*(*eta_mf)[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();

        Set::Patch<const Set::Scalar> eta_patch = eta_mf->Patch(lev,mfi);
        Set::Patch<const Set::Scalar> T         = temperature_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> rho_phys  = solid.rho_phys_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> cp_solid  = solid.cp_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> alpha     = alpha_solid_mf.array(mfi);
        Set::Patch<const Set::Scalar> T_solid_p = T_solid_mf.array(mfi);
        Set::Patch<const Set::Scalar> laser     = solid.laser_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       E_solid   = solid.energy_mf.Patch(lev,mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            // Always use the central stencil here, even at the box's own
            // edges: alpha_solid_mf/T_solid_mf ghost cells were just filled
            // above (periodic wrap or same-level neighbor, whichever
            // applies), so they're always valid - unlike Numeric::GetStencil's
            // one-sided Lo/Hi fallback (meant for true, unfilled physical
            // boundaries), which would otherwise silently downgrade to a
            // lower-order derivative at the domain edge. For a periodic
            // domain that one-sided fallback is actively wrong (ghost data
            // is valid there too) and was confirmed to seed a growing
            // error concentrated at the domain seam in testing.
            Set::Vector grad_raw    = Numeric::Gradient(eta_patch, i, j, k, 0, DX);
            Set::Vector grad_Tsolid = Numeric::Gradient(T_solid_p, i, j, k, 0, DX);
            Set::Scalar lap_Tsolid  = Numeric::Laplacian(T_solid_p, i, j, k, 0, DX);
            Set::Vector grad_alpha  = Numeric::Gradient(alpha, i, j, k, 0, DX);

            // eta_patch is Flame's own phase field (1=solid,0=gas under invert);
            // when invert is set, eta_patch already *is* the solid fraction, so
            // its gradient needs no sign flip.
            Set::Vector grad_phis = invert_local ? grad_raw : (-1.0*grad_raw);
            Set::Scalar eta_gas   = invert_local ? 1.0 - eta_patch(i,j,k) : eta_patch(i,j,k);
            Set::Scalar phi_s     = 1.0 - eta_gas;

            // (1) Bulk solid conduction - product-rule expansion of
            // div(phi_s * alpha * grad(T_solid)); T_solid only, never the
            // reactive blended field.
            Set::Scalar dTsolid_dt =
                grad_phis.dot(grad_Tsolid * alpha(i,j,k)) +
                grad_alpha.dot(phi_s * grad_Tsolid) +
                phi_s * alpha(i,j,k) * lap_Tsolid;

            Set::Scalar dEsolid_dt = rho_phys(i,j,k) * cp_solid(i,j,k) * dTsolid_dt;

            // (2) Interfacial exchange with the gas - bounded Robin-type flux.
            // T_gas is recovered by algebraically un-mixing the previous step's
            // blended temperature_mf (T = eta_gas*T_gas + phi_s*T_solid), the
            // same eta_recon floor RHS uses to keep the division bounded as
            // eta_gas -> 0.
            const Set::Scalar eta_recon = eta_gas > eta_cutoff_local ? eta_gas : eta_cutoff_local;
            Set::Scalar T_gas_local = (T(i,j,k) - (1.0 - eta_recon) * T_solid_p(i,j,k)) / eta_recon;
            Set::Scalar q_interface = h_interface_local * grad_phis.lpNorm<2>() * (T_gas_local - T_solid_p(i,j,k));
            dEsolid_dt += q_interface;

            // Laser flux (W/m^2) is localized to the regressing interface via
            // |grad(solid fraction)| (1/m), giving a volumetric source (W/m^3) -
            // same convention as the retired Flame loop's alpha*heatflux*grad_eta_mag.
            dEsolid_dt += laser(i,j,k) * grad_phis.lpNorm<2>();

            E_solid(i,j,k) += dt * dEsolid_dt;
        });
    }
}

void Hydro::Advance(int lev, Set::Scalar time, Set::Scalar dt)
{
    if (!managed) std::swap((*eta_old_mf)[lev], (*eta_mf)[lev]);
    std::swap(density_old_mf[lev],  density_mf[lev]);
    std::swap(momentum_old_mf[lev], momentum_mf[lev]);
    std::swap(energy_old_mf[lev],   energy_mf[lev]);

    //
    // UPDATE ETA AND CALCULATE ETADOT
    //

    if (!managed) UpdateEta(lev, time);
    if (managed)
    {
        UpdateFluxes(lev,time,dt);
        Mix(lev);
    }

    if (!managed)
    {
        eta_bc->define(geom[lev]);
        eta_bc->FillBoundary(*(*eta_mf)[lev], 0, 1, time, 0);
        eta_bc->FillBoundary(*(*eta_old_mf)[lev], 0, 1, time, 0);
    }
    neumann_bc_N->define(geom[lev]);
    neumann_bc_D->define(geom[lev]);
    neumann_bc_1->define(geom[lev]);
    neumann_bc_N->FillBoundary(*solid.density_mf[lev], 0, NSPECIES, time, 0);
    neumann_bc_D->FillBoundary(*solid.momentum_mf[lev], 0, AMREX_SPACEDIM, time, 0);
    neumann_bc_1->FillBoundary(*solid.energy_mf[lev], 0, 1, time, 0);
    neumann_bc_1->FillBoundary(*solid.rho_phys_mf[lev], 0, 1, time, 0);
    neumann_bc_1->FillBoundary(*solid.cp_mf[lev], 0, 1, time, 0);
    neumann_bc_1->FillBoundary(*solid.k_mf[lev], 0, 1, time, 0);
    neumann_bc_1->FillBoundary(*solid.laser_mf[lev], 0, 1, time, 0);

    // Hydro is the sole owner of solid.energy_mf's time evolution: this reads the
    // still-previous-step temperature_mf (not yet touched this Advance call) so
    // there is no same-step round trip through the freshly-computed T. Must run
    // before the RK stages below, which will recompute temperature_mf from the
    // conserved E this field feeds into via Mix()/the reconstruction branch.
    if (managed) AdvanceSolidEnergy(lev, time, dt);

    for (amrex::MFIter mfi(*(velocity_mf)[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox();
        amrex::Array4<const Set::Scalar> const& eta_new = (*(*eta_mf)[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& eta = (*(*eta_old_mf)[lev]).array(mfi);
        amrex::Array4<Set::Scalar>       const& etadot = (*etadot_mf[lev]).array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {

            etadot(i, j, k) = (eta_new(i, j, k) - eta(i, j, k)) / dt;
            if (invert) etadot(i,j,k) *= -1.0;

        });
    }


    //
    // DO TIME INTEGRATION (driving the RHS function)
    //

    // Organize references to the "new" solution
    amrex::Vector<amrex::MultiFab> solution_new;
    solution_new.emplace_back(*density_mf[lev].get(),amrex::MakeType::make_alias,0,NSPECIES);
    solution_new.emplace_back(*momentum_mf[lev].get(),amrex::MakeType::make_alias,0,AMREX_SPACEDIM);
    solution_new.emplace_back(*energy_mf[lev].get(),amrex::MakeType::make_alias,0,1);

    // Organize references to the "old" solution
    amrex::Vector<amrex::MultiFab> solution_old;
    solution_old.emplace_back(*density_old_mf[lev].get(),amrex::MakeType::make_alias,0,NSPECIES);
    solution_old.emplace_back(*momentum_old_mf[lev].get(),amrex::MakeType::make_alias,0,AMREX_SPACEDIM);
    solution_old.emplace_back(*energy_old_mf[lev].get(),amrex::MakeType::make_alias,0,1);

    // Create the time integrator
    amrex::TimeIntegrator timeintegrator(solution_new, time);

    // Set the time integrator RHS - in this case, just relay to our current RHS function
    timeintegrator.set_rhs([&](amrex::Vector<amrex::MultiFab> & rhs_mf, amrex::Vector<amrex::MultiFab> & solution_mf, const Set::Scalar time)
    {
        RHS(lev, time, dt,
            rhs_mf[0], rhs_mf[1], rhs_mf[2],
            solution_mf[0],solution_mf[1],solution_mf[2]);
    });

    auto fill_conserved_boundaries = [&](   amrex::MultiFab& rho, amrex::MultiFab& M,
                                            amrex::MultiFab& E, Set::Scalar fill_time,
                                            bool use_old_eta)
    {
        ApplyCutoffToConserved(lev, rho, M, E, false, use_old_eta);
        density_bc->define(geom[lev]);
        momentum_bc->define(geom[lev]);
        energy_bc->define(geom[lev]);
        density_bc->FillBoundary(rho, 0, NSPECIES, fill_time, 0);
        momentum_bc->FillBoundary(M, 0, AMREX_SPACEDIM, fill_time, 0);
        energy_bc->FillBoundary(E, 0, 1, fill_time, 0);
        ApplyCutoffToConserved(lev, rho, M, E, true, use_old_eta);
    };

    // Integrator::TimeStep fills AMR coarse/fine ghost cells for evolving
    // fields before this call. The transient RK stage aliases still need their
    // same-level and physical ghost cells refreshed after every stage update.
    timeintegrator.set_post_stage_action([&](amrex::Vector<amrex::MultiFab> & stage_mf, Set::Scalar time)
    {
        fill_conserved_boundaries(stage_mf[0], stage_mf[1], stage_mf[2], time, true);
    });

    // Do the update
    timeintegrator.advance(solution_old, solution_new, time, dt);

    fill_conserved_boundaries(  *density_mf[lev], *momentum_mf[lev],
                                *energy_mf[lev], time + dt, false);

    //
    // APPLY CUTOFFS AND DO DYNAMIC TIMESTEP CALCULATION
    //

    Set::Scalar dt_max = std::numeric_limits<Set::Scalar>::max();
    const amrex::Box domain = geom[lev].Domain();
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
            Set::Scalar eta = invert ? 1.0-eta_patch(i,j,k): eta_patch(i,j,k);

            if (cutoff >= 0.0 && cutoff < 1.0 && eta < cutoff)
            {
                for (int n=0; n<NSPECIES; ++n)
                {
                    rho_new(i,j,k,n) = rho_solid(i,j,k,n);
                }
                M_new(i,j,k,0)   = M_solid(i,j,k,0);
                M_new(i,j,k,1)   = M_solid(i,j,k,1);
                #if AMREX_SPACEDIM == 3
                M_new(i,j,k,2)   = M_solid(i,j,k,2);
                #endif
                E_new(i,j,k,0)   = E_solid(i,j,k,0);
            }

            auto sten = Numeric::GetStencil(i, j, k, domain);
            Set::Matrix gradu        = Numeric::Gradient(u, i, j, k, DX, sten);
            #if AMREX_SPACEDIM == 2
            omega(i, j, k) = eta * (gradu(1,0) - gradu(0,1));
            #elif AMREX_SPACEDIM == 3
            omega(i, j, k, 0) = eta * (gradu(2,1) - gradu(1,2));
            omega(i, j, k, 1) = eta * (gradu(0,2) - gradu(2,0));
            omega(i, j, k, 2) = eta * (gradu(1,0) - gradu(0,1));
            #endif

            if (dynamictimestep.on)
            {
                *dt_max_handle =                          std::fabs(cfl * DX[0] / (u(i,j,k,0)*eta + small));
                *dt_max_handle = std::min(*dt_max_handle, std::fabs(cfl * DX[1] / (u(i,j,k,1)*eta + small)));
                *dt_max_handle = std::min(*dt_max_handle, std::fabs(cfl_v * DX[0]*DX[0] / (Source(i,j,k,NSPECIES)+small)));
                *dt_max_handle = std::min(*dt_max_handle, std::fabs(cfl_v * DX[1]*DX[1] / (Source(i,j,k,NSPECIES+1)+small)));
            }
        });
    }


    if (dynamictimestep.on)
    {
        this->DynamicTimestep_SyncTimeStep(lev,dt_max);
    }

}//end Advance

// Plotfiles show correct conserved quantities, but derived quantities are one time step behind
// This function computes the derived quantities at all AMR levels for correct time matching.
void Hydro::PreparePlotFileData()
{
    BL_PROFILE("Integrator::Hydro::PreparePlotFileData");
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        RefreshDerivedPlotFields(lev);
    }
}

// Compute the derived quantities from conservatives for a given AMR level
void Hydro::RefreshDerivedPlotFields(int lev)
{
    BL_PROFILE("Integrator::Hydro::RefreshDerivedPlotFields");

    // See Hydro::RHS for why this datum must match the active gas EOS.
    const Set::Scalar T_ref_energy =
        (std::string(gas.eos.model_name()) == "tpg")
            ? Model::Gas::EOS::sensible_reference_temperature : 0.0;

    const Set::Scalar* DX = geom[lev].CellSize();
    const amrex::Box domain = geom[lev].Domain();
    for (amrex::MFIter mfi(*velocity_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();

        amrex::Array4<const Set::Scalar> const& eta_patch = (*(*eta_mf)[lev]).array(mfi);

        Set::Patch<const Set::Scalar> rho       = density_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> M         = momentum_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> E         = energy_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> rho_solid = solid.density_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> M_solid   = solid.momentum_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> E_solid   = solid.energy_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> rho_phys  = solid.rho_phys_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> cp_solid  = solid.cp_mf.Patch(lev,mfi);

        Set::Patch<Set::Scalar> scratch         = scratch_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> v               = velocity_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> p               = pressure_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> T               = temperature_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> Y               = mass_fraction_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> X               = mole_fraction_mf.Patch(lev,mfi);

        amrex::Array4<Set::Scalar> mu_arr;
        amrex::Array4<Set::Scalar> k_arr;
        amrex::Array4<Set::Scalar> D_arr;
        amrex::Array4<Set::Scalar> wdot_arr;
        amrex::Array4<Set::Scalar> qdot_arr;
        if (details)
        {
            mu_arr    = viscosity_mf.Patch(lev,mfi);
            k_arr     = thermal_conductivity_coeff_mf.Patch(lev,mfi);
            D_arr     = diffusion_coeff_mf.Patch(lev,mfi);
            wdot_arr  = wdot_mf.Patch(lev,mfi);
            qdot_arr  = qdot_mf.Patch(lev,mfi);
        }

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Scalar eta = invert ? 1.0-eta_patch(i,j,k): eta_patch(i,j,k);
            const Set::Scalar eta_cutoff = (cutoff >= 0.0 && cutoff < 1.0) ? cutoff : small;

            std::array<Set::Scalar, NSPECIES> rhoY_fluid;
            Set::Scalar Mx_fluid = 0.0;
            Set::Scalar My_fluid = 0.0;
            #if AMREX_SPACEDIM == 3
            Set::Scalar Mz_fluid = 0.0;
            #endif
            if (eta <= eta_cutoff)
            {
                for (int n=0; n<NSPECIES; ++n) rhoY_fluid[n] = rho_solid(i,j,k,n);
                Mx_fluid = M_solid(i,j,k,0);
                My_fluid = M_solid(i,j,k,1);
                #if AMREX_SPACEDIM == 3
                Mz_fluid = M_solid(i,j,k,2);
                #endif
            }
            else
            {
                for (int n=0; n<NSPECIES; ++n)
                {
                    rhoY_fluid[n] = (rho(i,j,k,n) - rho_solid(i,j,k,n)*(1.0 - eta))/(eta + small);
                }
                Mx_fluid = (M(i,j,k,0) - M_solid(i,j,k,0)*(1.0 - eta))/(eta + small);
                My_fluid = (M(i,j,k,1) - M_solid(i,j,k,1)*(1.0 - eta))/(eta + small);
                #if AMREX_SPACEDIM == 3
                Mz_fluid = (M(i,j,k,2) - M_solid(i,j,k,2)*(1.0 - eta))/(eta + small);
                #endif
            }
            ProjectSpeciesDensities(rhoY_fluid, i, j, k);
            for (int n=0; n<NSPECIES; ++n) scratch(i,j,k,n) = rhoY_fluid[n];

            Set::Scalar density_fluid = gas.ComputeLocalFractions(scratch, Y, X, i, j, k);

            // Solid-only caloric temperature - always well-defined, since
            // E_solid is never touched by the eta<=cutoff conserved-state
            // forcing below.
            Set::Scalar T_solid_caloric =
                T_ref_energy + E_solid(i,j,k) / (rho_phys(i,j,k) * cp_solid(i,j,k) + small);

            // Gas-side reconstruction is only meaningful where the conserved
            // state actually still holds a mixed fluid/solid blend. For
            // eta<=eta_cutoff, ApplyCutoffToConserved has already overwritten
            // E(i,j,k) (and rho, M) to equal the pure-solid state exactly -
            // the real fluid energy that would be needed to reconstruct a
            // "gas" temperature there no longer exists. Previously this branch
            // ran unconditionally with an eta_recon=max(eta,eta_cutoff) floor,
            // intended to keep the division bounded - but once E==E_solid
            // exactly, that reconstruction algebraically collapses to
            // Ef_forT==E_solid (built on the *physical* solid density/cp
            // scale) divided into density_fluid (the *tame* Riemann-blend
            // density, ~195x smaller) inside gas.ComputeT: a physically
            // nonsensical, unbounded specific-energy mismatch that produced
            // spurious temperatures in the tens-to-hundreds-of-thousands of
            // Kelvin range - visible even after being weighted by the small
            // true eta in the blend below, and confirmed to grow worse (not
            // better) under mesh/eta refinement. Skip it entirely below
            // cutoff and just report the solid caloric value directly.
            Set::Scalar T_gas_inversion = T_solid_caloric;
            if (eta > eta_cutoff)
            {
                Set::Scalar Mx_forT = (M(i,j,k,0) - M_solid(i,j,k,0)*(1.0 - eta))/(eta + small);
                Set::Scalar My_forT = (M(i,j,k,1) - M_solid(i,j,k,1)*(1.0 - eta))/(eta + small);
                #if AMREX_SPACEDIM == 3
                Set::Scalar Mz_forT = (M(i,j,k,2) - M_solid(i,j,k,2)*(1.0 - eta))/(eta + small);
                #endif
                Set::Scalar Ef_forT = (E(i,j,k) - E_solid(i,j,k)*(1.0 - eta))/(eta + small);

                #if AMREX_SPACEDIM == 2
                T_gas_inversion = gas.ComputeT(density_fluid, Mx_forT, My_forT, Ef_forT, T(i,j,k), X, i, j, k);
                #elif AMREX_SPACEDIM == 3
                T_gas_inversion = gas.ComputeT(density_fluid, Mx_forT, My_forT, Mz_forT, Ef_forT, T(i,j,k), X, i, j, k);
                #endif
            }
            T(i,j,k) = eta*T_gas_inversion + (1.0-eta)*T_solid_caloric;

            // Pressure keeps its previous eta<=eta_cutoff hard cutoff (separate
            // design choice - see Hydro::RHS) rather than following T's smooth
            // blend.
            if (eta > eta_cutoff)
            {
                p(i,j,k) = gas.ComputeP(density_fluid, T_gas_inversion, X, i, j, k);
            }
            // else: leave p(i,j,k) as the already-prescribed local gas pressure
            // (see Hydro::RHS).
            v(i,j,k,0) = Mx_fluid / density_fluid;
            v(i,j,k,1) = My_fluid / density_fluid;
            #if AMREX_SPACEDIM == 3
            v(i,j,k,2) = Mz_fluid / density_fluid;
            #endif

            if (eta < small)
            {
                v(i,j,k,0) *= eta;
                v(i,j,k,1) *= eta;
                #if AMREX_SPACEDIM == 3
                v(i,j,k,2) *= eta;
                #endif
            }

            if (details)
            {
                D_arr(i,j,k,0) = 0.0;
                gas.diffusion_coeffs(D_arr, T(i,j,k), p(i,j,k), X, i, j, k);
                k_arr(i,j,k) = gas.thermal_conductivity(T(i,j,k), X, i, j, k);
                mu_arr(i,j,k) = gas.dynamic_viscosity(T(i,j,k), X, i, j, k);

                std::array<Set::Scalar, NSPECIES> rhoY;
                for (int n=0; n<NSPECIES; ++n) rhoY[n] = rho(i,j,k,n);

                std::array<Set::Scalar, NSPECIES> wdot;
                Set::Scalar qdot = 0.0;
                std::tie(wdot, qdot) = chemistry.ComputeChemistrySources(p(i,j,k), T(i,j,k), rhoY, dt[lev], &gas);
                qdot_arr(i,j,k) = qdot;
                for (int n=0; n<NSPECIES; ++n)
                {
                    wdot_arr(i,j,k,n) = wdot[n];
                }
            }
        });
    }

    for (amrex::MFIter mfi(*vorticity_mf[lev], false); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        Set::Patch<const Set::Scalar> eta_patch = eta_mf->Patch(lev,mfi);
        Set::Patch<const Set::Scalar> u         = velocity_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar> omega          = vorticity_mf.Patch(lev,mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Scalar eta = invert ? 1.0-eta_patch(i,j,k): eta_patch(i,j,k);
            auto sten = Numeric::GetStencil(i, j, k, domain);
            Set::Matrix gradu = Numeric::Gradient(u, i, j, k, DX, sten);
            #if AMREX_SPACEDIM == 2
            omega(i, j, k) = eta * (gradu(1,0) - gradu(0,1));
            #elif AMREX_SPACEDIM == 3
            omega(i, j, k, 0) = eta * (gradu(2,1) - gradu(1,2));
            omega(i, j, k, 1) = eta * (gradu(0,2) - gradu(2,0));
            omega(i, j, k, 2) = eta * (gradu(1,0) - gradu(0,1));
            #endif
        });
    }

    vorticity_mf[lev]->FillBoundary(true);
}


void Hydro::RHS(int lev, Set::Scalar time, Set::Scalar dt,
                amrex::MultiFab &rho_rhs_mf,
                amrex::MultiFab &M_rhs_mf,
                amrex::MultiFab &E_rhs_mf,
                amrex::MultiFab &rho_mf,
                amrex::MultiFab &M_mf,
                amrex::MultiFab &E_mf)
{
    ApplyCutoffToConserved(lev, rho_mf, M_mf, E_mf, false, true);

    // RHS evaluation samples ghost cells directly for boundary fluxes and
    // centered gradients. Refresh physical and same-level ghosts here as a guard
    // for the initial RK evaluation as well as subsequent stage evaluations.
    density_bc->define(geom[lev]);
    momentum_bc->define(geom[lev]);
    energy_bc->define(geom[lev]);
    density_bc->FillBoundary(rho_mf, 0, NSPECIES, time, 0);
    momentum_bc->FillBoundary(M_mf, 0, AMREX_SPACEDIM, time, 0);
    energy_bc->FillBoundary(E_mf, 0, 1, time, 0);
    if (!managed)
    {
        eta_bc->define(geom[lev]);
        eta_bc->FillBoundary(*(*eta_mf)[lev], 0, 1, time, 0);
        eta_bc->FillBoundary(*(*eta_old_mf)[lev], 0, 1, time, 0);
    }
    neumann_bc_N->define(geom[lev]);
    neumann_bc_D->define(geom[lev]);
    neumann_bc_1->define(geom[lev]);
    neumann_bc_N->FillBoundary(*solid.density_mf[lev], 0, NSPECIES, time, 0);
    neumann_bc_D->FillBoundary(*solid.momentum_mf[lev], 0, AMREX_SPACEDIM, time, 0);
    neumann_bc_1->FillBoundary(*solid.energy_mf[lev], 0, 1, time, 0);
    neumann_bc_1->FillBoundary(*solid.rho_phys_mf[lev], 0, 1, time, 0);
    neumann_bc_1->FillBoundary(*solid.cp_mf[lev], 0, 1, time, 0);
    neumann_bc_1->FillBoundary(*solid.k_mf[lev], 0, 1, time, 0);
    neumann_bc_1->FillBoundary(*solid.laser_mf[lev], 0, 1, time, 0);
    ApplyCutoffToConserved(lev, rho_mf, M_mf, E_mf, true, true);

    // The solid caloric energy datum must match whichever gas EOS is active (see
    // Flame::UpdateFluxes for why); look it up once here rather than per cell.
    const Set::Scalar T_ref_energy =
        (std::string(gas.eos.model_name()) == "tpg")
            ? Model::Gas::EOS::sensible_reference_temperature : 0.0;

    int nghost = flux_scheme == FluxScheme::Advect ? advect.NGhost() : 1;
    int primitive_nghost = nghost < 2 ? nghost : 2;
    int diffusion_nghost = 1;
    const amrex::BoxArray &ba = energy_mf[lev]->boxArray();
    const amrex::DistributionMapping &dm = energy_mf[lev]->DistributionMap();
    amrex::MultiFab rho_fluid_mf(ba,dm,NSPECIES,nghost);    // fluid species densities
    amrex::MultiFab M_fluid_mf(ba,dm,AMREX_SPACEDIM,nghost); // fluid momentum
    amrex::MultiFab E_fluid_mf(ba,dm,1,nghost);              // fluid energy
    amrex::MultiFab rho_sum_mf(ba,dm,1,primitive_nghost);             // sum_k[rhoY_k]
    amrex::MultiFab mixed_k_mf(ba,dm,1,primitive_nghost);             // mixture averaged thermal conductivity coefficient
    amrex::MultiFab mixed_kT_mf(ba,dm,AMREX_SPACEDIM,diffusion_nghost);  // mixture averaged thermal conductivity * temperature gradient
    amrex::MultiFab mixed_mu_mf(ba,dm,1,primitive_nghost);            // mixture averaged dynamic viscosity
    amrex::MultiFab mixed_H_mf(ba,dm,1,primitive_nghost);             // Perfect gas mixture enthalpy, H=cp_mix*T
    amrex::MultiFab DKM_mf(ba,dm,NSPECIES,primitive_nghost);      // Diffusion coefficent for species k into mixture
    amrex::MultiFab rhoHDYx_mf(ba,dm,NSPECIES,diffusion_nghost);  // species enthalpy diffusion, rho*H*D*dY/dx
    amrex::MultiFab rhoHDYy_mf(ba,dm,NSPECIES,diffusion_nghost);  // species enthalpy diffusion, rho*H*D*dY/dy
    amrex::MultiFab rhoDYx_mf(ba,dm,NSPECIES,diffusion_nghost);   // Fickian diffusion, rho*D*dY/dx
    amrex::MultiFab rhoDYy_mf(ba,dm,NSPECIES,diffusion_nghost);   // Fickian diffusion, rho*D*dY/dy
    #if AMREX_SPACEDIM == 3
    amrex::MultiFab rhoHDYz_mf(ba,dm,NSPECIES,diffusion_nghost);  // species enthalpy diffusion, rho*H*D*dY/dz
    amrex::MultiFab rhoDYz_mf(ba,dm,NSPECIES,diffusion_nghost);   // Fickian diffusion, rho*D*dY/dz
    #endif

    // Values only to be written if details=true
    amrex::Array4<Set::Scalar> mu_arr;
    amrex::Array4<Set::Scalar> k_arr;
    amrex::Array4<Set::Scalar> D_arr;
    amrex::Array4<Set::Scalar> wdot_arr;
    amrex::Array4<Set::Scalar> qdot_arr;

    const Set::Scalar* DX = geom[lev].CellSize();

    Set::Scalar diffuse_source_norm = 1.0;
    {
        // Normalize against the same cell-centered eta mask used by cutoff,
        // so retained cells carry the full diffuse source integral.
        const Set::Scalar cutoff_local = cutoff;
        const bool invert_local = invert;

        Set::Scalar reference_source_weight = amrex::ReduceSum(
            *(*eta_old_mf)[lev], amrex::IntVect(0),
            [=] AMREX_GPU_HOST_DEVICE (const amrex::Box& bx,
                                       const amrex::Array4<const Set::Scalar>& eta_arr) -> Set::Scalar
            {
                Set::Scalar sum = 0.0;
                const auto lo = amrex::lbound(bx);
                const auto hi = amrex::ubound(bx);
                auto sten = BCGhostStencil();
                for (int k = lo.z; k <= hi.z; ++k)
                for (int j = lo.y; j <= hi.y; ++j)
                for (int i = lo.x; i <= hi.x; ++i)
                {
                    Set::Vector grad_eta = Numeric::Gradient(eta_arr, i, j, k, 0, DX, sten);
                    if (invert_local) grad_eta *= -1.0;
                    sum += grad_eta.lpNorm<2>();
                }
                return sum;
            });

        Set::Scalar active_source_weight = amrex::ReduceSum(
            *(*eta_old_mf)[lev], amrex::IntVect(0),
            [=] AMREX_GPU_HOST_DEVICE (const amrex::Box& bx,
                                       const amrex::Array4<const Set::Scalar>& eta_arr) -> Set::Scalar
            {
                Set::Scalar sum = 0.0;
                const auto lo = amrex::lbound(bx);
                const auto hi = amrex::ubound(bx);
                auto sten = BCGhostStencil();
                for (int k = lo.z; k <= hi.z; ++k)
                for (int j = lo.y; j <= hi.y; ++j)
                for (int i = lo.x; i <= hi.x; ++i)
                {
                    Set::Scalar eta = invert_local ?
                        1.0 - eta_arr(i,j,k,0):
                        eta_arr(i,j,k,0);
                    if (cutoff_local >= 0.0 && cutoff_local < 1.0 && eta <= cutoff_local) continue;
                    Set::Vector grad_eta = Numeric::Gradient(eta_arr, i, j, k, 0, DX, sten);
                    if (invert_local) grad_eta *= -1.0;
                    sum += grad_eta.lpNorm<2>();
                }
                return sum;
            });

        amrex::ParallelDescriptor::ReduceRealSum(reference_source_weight);
        amrex::ParallelDescriptor::ReduceRealSum(active_source_weight);
        if (active_source_weight > small && reference_source_weight > 0.0)
        {
            diffuse_source_norm = reference_source_weight / active_source_weight;
        }
    }

    for (amrex::MFIter mfi(*(velocity_mf)[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox();
        const amrex::Box valid_bx = mfi.validbox();
        const auto valid_lo = amrex::lbound(valid_bx);
        const auto valid_hi = amrex::ubound(valid_bx);
        const int first_pass_ghost_depth = primitive_nghost;
        amrex::Array4<const Set::Scalar> const& eta_patch = (*(*eta_old_mf)[lev]).array(mfi);

        Set::Patch<const Set::Scalar> rho       = rho_mf.array(mfi);  // density
        Set::Patch<const Set::Scalar> M         = M_mf.array(mfi);    // momentum
        Set::Patch<const Set::Scalar> E         = E_mf.array(mfi);    // total energy (internal energy + kinetic energy) per unit volume (E/rho = e + 0.5*v^2)

        Set::Patch<const Set::Scalar> rho_solid = solid.density_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> M_solid   = solid.momentum_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> E_solid   = solid.energy_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> rho_phys  = solid.rho_phys_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> cp_solid  = solid.cp_mf.Patch(lev,mfi);

        Set::Patch<Set::Scalar> scratch         = scratch_mf.Patch(lev,mfi);

        Set::Patch<Set::Scalar>       v         = velocity_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       p         = pressure_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       T         = temperature_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       Y         = mass_fraction_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       X         = mole_fraction_mf.Patch(lev,mfi);

        Set::Patch<Set::Scalar>       rho_fluid = rho_fluid_mf.array(mfi);
        Set::Patch<Set::Scalar>       M_fluid   = M_fluid_mf.array(mfi);
        Set::Patch<Set::Scalar>       E_fluid   = E_fluid_mf.array(mfi);
        Set::Patch<Set::Scalar>       DKM       = DKM_mf.array(mfi);
        Set::Patch<Set::Scalar>       rho_sum   = rho_sum_mf.array(mfi);
        Set::Patch<Set::Scalar>       mixed_k   = mixed_k_mf.array(mfi);
        Set::Patch<Set::Scalar>       mixed_mu  = mixed_mu_mf.array(mfi);
        Set::Patch<Set::Scalar>       mixed_H   = mixed_H_mf.array(mfi);
        Set::Patch<Set::Scalar>       mixed_kT  = mixed_kT_mf.array(mfi);
        Set::Patch<Set::Scalar>       rhoHDYx   = rhoHDYx_mf.array(mfi);
        Set::Patch<Set::Scalar>       rhoHDYy   = rhoHDYy_mf.array(mfi);
        Set::Patch<Set::Scalar>       rhoDYx    = rhoDYx_mf.array(mfi);
        Set::Patch<Set::Scalar>       rhoDYy    = rhoDYy_mf.array(mfi);
        #if AMREX_SPACEDIM == 3
        Set::Patch<Set::Scalar>       rhoHDYz   = rhoHDYz_mf.array(mfi);
        Set::Patch<Set::Scalar>       rhoDYz    = rhoDYz_mf.array(mfi);
        #endif

        if (details)
        {
            mu_arr    = viscosity_mf.Patch(lev,mfi);
            k_arr     = thermal_conductivity_coeff_mf.Patch(lev,mfi);
            D_arr     = diffusion_coeff_mf.Patch(lev,mfi);
            wdot_arr  = wdot_mf.Patch(lev,mfi);
            qdot_arr  = qdot_mf.Patch(lev,mfi);
        }

        // First ParallelFor loop to get initial values needed for gradients
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            if (!IsValidAxisGhostOrCornerWithin(
                    i, j, k,
                    valid_lo.x, valid_hi.x,
                    valid_lo.y, valid_hi.y,
                    valid_lo.z, valid_hi.z,
                    first_pass_ghost_depth,
                    1))
            {
                return;
            }

            Set::Scalar eta = invert ? 1.0-eta_patch(i,j,k) : eta_patch(i,j,k);
            const Set::Scalar eta_cutoff = (cutoff >= 0.0 && cutoff < 1.0) ? cutoff : small;

            // Recover the gas state once before computing primitive or
            // transport quantities; the flux branches reuse these fields.
            std::array<Set::Scalar, NSPECIES> rhoY_fluid;
            Set::Scalar Mx_fluid = 0.0;
            Set::Scalar My_fluid = 0.0;
            #if AMREX_SPACEDIM == 3
            Set::Scalar Mz_fluid = 0.0;
            #endif
            Set::Scalar Ef_fluid = 0.0;
            if (eta <= eta_cutoff)
            {
                for (int n=0; n<NSPECIES; ++n) rhoY_fluid[n] = rho_solid(i,j,k,n);
                Mx_fluid = M_solid(i,j,k,0);
                My_fluid = M_solid(i,j,k,1);
                #if AMREX_SPACEDIM == 3
                Mz_fluid = M_solid(i,j,k,2);
                #endif
                Ef_fluid = E_solid(i,j,k);
            }
            else
            {
                for (int n=0; n<NSPECIES; ++n)
                {
                    rhoY_fluid[n] = (rho(i,j,k,n) - rho_solid(i,j,k,n)*(1.0 - eta))/(eta + small);
                }
                Mx_fluid = (M(i,j,k,0) - M_solid(i,j,k,0)*(1.0 - eta))/(eta + small);
                My_fluid = (M(i,j,k,1) - M_solid(i,j,k,1)*(1.0 - eta))/(eta + small);
                #if AMREX_SPACEDIM == 3
                Mz_fluid = (M(i,j,k,2) - M_solid(i,j,k,2)*(1.0 - eta))/(eta + small);
                #endif
                Ef_fluid = (E(i,j,k) - E_solid(i,j,k)*(1.0 - eta))/(eta + small);
            }

            ProjectSpeciesDensities(rhoY_fluid, i, j, k);
            for (int n=0; n<NSPECIES; ++n)
            {
                scratch(i,j,k,n) = rhoY_fluid[n];
                rho_fluid(i,j,k,n) = rhoY_fluid[n];
            }
            M_fluid(i,j,k,0) = Mx_fluid;
            M_fluid(i,j,k,1) = My_fluid;
            #if AMREX_SPACEDIM == 3
            M_fluid(i,j,k,2) = Mz_fluid;
            #endif
            E_fluid(i,j,k) = Ef_fluid;

            Set::Scalar density_fluid = gas.ComputeLocalFractions(scratch, Y, X, i, j, k);

            // Solid-only caloric temperature - always well-defined, since
            // E_solid is never touched by the eta<=cutoff conserved-state
            // forcing in ApplyCutoffToConserved.
            Set::Scalar T_solid_caloric =
                T_ref_energy + E_solid(i,j,k) / (rho_phys(i,j,k) * cp_solid(i,j,k) + small);

            // Gas-side reconstruction is only meaningful where the conserved
            // state still holds a genuine mixed fluid/solid blend. For
            // eta<=eta_cutoff, ApplyCutoffToConserved (called at the top of
            // RHS) has already overwritten E(i,j,k)/rho/M to the pure-solid
            // state exactly - the real fluid energy needed to reconstruct a
            // "gas" temperature there no longer exists. Previously this ran
            // unconditionally with an eta_recon=max(eta,eta_cutoff) floor to
            // keep the division bounded - but once E==E_solid exactly, that
            // reconstruction algebraically collapses to Ef_forT==E_solid
            // (built on the *physical* solid density/cp scale) divided by
            // density_fluid (the *tame* Riemann-blend density, ~195x
            // smaller) inside gas.ComputeT: a physically nonsensical,
            // unbounded specific-energy mismatch that produced spurious
            // temperatures in the tens-to-hundreds-of-thousands of Kelvin
            // range - visible even after being weighted by the small true
            // eta in the blend below, and confirmed (via a mesh/eta
            // refinement study) to grow worse, not better, under refinement.
            // Skip it entirely below cutoff and just report the solid
            // caloric value directly.
            Set::Scalar T_gas_inversion = T_solid_caloric;
            if (eta > eta_cutoff)
            {
                Set::Scalar Mx_forT = (M(i,j,k,0) - M_solid(i,j,k,0)*(1.0 - eta))/(eta + small);
                Set::Scalar My_forT = (M(i,j,k,1) - M_solid(i,j,k,1)*(1.0 - eta))/(eta + small);
                #if AMREX_SPACEDIM == 3
                Set::Scalar Mz_forT = (M(i,j,k,2) - M_solid(i,j,k,2)*(1.0 - eta))/(eta + small);
                #endif
                Set::Scalar Ef_forT = (E(i,j,k) - E_solid(i,j,k)*(1.0 - eta))/(eta + small);

                #if AMREX_SPACEDIM == 2
                T_gas_inversion = gas.ComputeT(density_fluid, Mx_forT, My_forT, Ef_forT, T(i,j,k), X, i, j, k);
                #elif AMREX_SPACEDIM == 3
                T_gas_inversion = gas.ComputeT(density_fluid, Mx_forT, My_forT, Mz_forT, Ef_forT, T(i,j,k), X, i, j, k);
                #endif
            }
            T(i,j,k) = eta*T_gas_inversion + (1.0-eta)*T_solid_caloric;

            // Pressure keeps its previous eta<=eta_cutoff hard cutoff (a separate,
            // deliberate design choice - see below) rather than following T's new
            // smooth blend: deriving it from density_fluid (the tame solid.density)
            // times a gas-side T would reintroduce the original t=0 pressure-spike
            // bug. Only update p(i,j,k) above cutoff, using the unblended
            // T_gas_inversion (not the blended T) paired with the pure
            // reconstructed gas density - the blended T is reserved for
            // conduction/kinetics only.
            if (eta > eta_cutoff)
            {
                p(i,j,k) = gas.ComputeP(density_fluid, T_gas_inversion, X, i, j, k);
            }
            // else: leave p(i,j,k) as already prescribed.

            v(i,j,k,0) = Mx_fluid/density_fluid;
            v(i,j,k,1) = My_fluid/density_fluid;
            #if AMREX_SPACEDIM == 3
            v(i,j,k,2) = Mz_fluid/density_fluid;
            #endif

            if (eta < small)
            {
                v(i,j,k,0) *= eta;
                v(i,j,k,1) *= eta;

            #if AMREX_SPACEDIM == 3
                v(i,j,k,2) *= eta;
            #endif
            }

            rho_sum(i,j,k) = density_fluid;
            gas.diffusion_coeffs(DKM, T(i,j,k), p(i,j,k), X, i, j, k);
            mixed_k(i,j,k) = gas.thermal_conductivity(T(i,j,k), X, i, j, k);
            mixed_mu(i,j,k) = gas.dynamic_viscosity(T(i,j,k), X, i, j, k);
            mixed_H(i,j,k) = gas.enthalpy_mass(T(i,j,k), X, i, j, k);

            if ( details )
            {
                mu_arr(i,j,k) = mixed_mu(i,j,k);
                k_arr(i,j,k) = mixed_k(i,j,k);
                for (int n=0; n<NSPECIES; ++n) D_arr(i,j,k,n) = DKM(i,j,k,n);
            }
        });

        // Second ParallelFor loop to get first gradients
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            if (!IsValidOrAxisGhostWithin(
                    i, j, k,
                    valid_lo.x, valid_hi.x,
                    valid_lo.y, valid_hi.y,
                    valid_lo.z, valid_hi.z,
                    1))
            {
                return;
            }

            const int outside_dir = OutsideDirection(
                i, j, k,
                valid_lo.x, valid_hi.x,
                valid_lo.y, valid_hi.y,
                valid_lo.z, valid_hi.z);
            auto sten = Numeric::GetStencil(i,j,k,bx);

            auto compute_x = [&]() AMREX_GPU_DEVICE
            {
                for (int n=0; n<NSPECIES; ++n)
                {
                    Set::Scalar grad_Yx = Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(Y,i,j,k,n,DX,sten);
                    rhoHDYx(i,j,k,n) = rho_sum(i,j,k)*mixed_H(i,j,k)*DKM(i,j,k,n)*grad_Yx;
                    rhoDYx(i,j,k,n)  = rho_sum(i,j,k)*DKM(i,j,k,n)*grad_Yx;
                }
                Set::Scalar grad_Tx = Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(T,i,j,k,0,DX,sten);
                mixed_kT(i,j,k,0) = mixed_k(i,j,k)*grad_Tx;
            };

#if AMREX_SPACEDIM >= 2
            auto compute_y = [&]() AMREX_GPU_DEVICE
            {
                for (int n=0; n<NSPECIES; ++n)
                {
                    Set::Scalar grad_Yy = Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(Y,i,j,k,n,DX,sten);
                    rhoHDYy(i,j,k,n) = rho_sum(i,j,k)*mixed_H(i,j,k)*DKM(i,j,k,n)*grad_Yy;
                    rhoDYy(i,j,k,n)  = rho_sum(i,j,k)*DKM(i,j,k,n)*grad_Yy;
                }
                Set::Scalar grad_Ty = Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(T,i,j,k,0,DX,sten);
                mixed_kT(i,j,k,1) = mixed_k(i,j,k)*grad_Ty;
            };
#endif

#if AMREX_SPACEDIM == 3
            auto compute_z = [&]() AMREX_GPU_DEVICE
            {
                for (int n=0; n<NSPECIES; ++n)
                {
                    Set::Scalar grad_Yz = Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(Y,i,j,k,n,DX,sten);
                    rhoHDYz(i,j,k,n) = rho_sum(i,j,k)*mixed_H(i,j,k)*DKM(i,j,k,n)*grad_Yz;
                    rhoDYz(i,j,k,n)  = rho_sum(i,j,k)*DKM(i,j,k,n)*grad_Yz;
                }
                Set::Scalar grad_Tz = Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(T,i,j,k,0,DX,sten);
                mixed_kT(i,j,k,2) = mixed_k(i,j,k)*grad_Tz;
            };
#endif

            if (outside_dir < 0 || outside_dir == 0) compute_x();
#if AMREX_SPACEDIM >= 2
            if (outside_dir < 0 || outside_dir == 1) compute_y();
#endif
#if AMREX_SPACEDIM == 3
            if (outside_dir < 0 || outside_dir == 2) compute_z();
#endif
        });
    }

    density_bc->FillBoundary(rho_fluid_mf, 0, NSPECIES, time, 0);
    momentum_bc->FillBoundary(M_fluid_mf, 0, AMREX_SPACEDIM, time, 0);
    energy_bc->FillBoundary(E_fluid_mf, 0, 1, time, 0);
    rho_sum_mf.FillBoundary(true);
    mixed_kT_mf.FillBoundary(true);
    rhoHDYx_mf.FillBoundary(true);
    rhoHDYy_mf.FillBoundary(true);
    rhoDYx_mf.FillBoundary(true);
    rhoDYy_mf.FillBoundary(true);
    #if AMREX_SPACEDIM == 3
    rhoHDYz_mf.FillBoundary(true);
    rhoDYz_mf.FillBoundary(true);
    #endif

    const auto advect_op = advect;
    const Numeric::Advect::Options conservative_options{Numeric::Advect::Form::Conservative};

    for (amrex::MFIter mfi(*(*eta_mf)[lev], false); mfi.isValid(); ++mfi)
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
        // Physical (not the tame Riemann-blend) solid density, needed below to
        // recover the injected mass's specific internal energy from E_solid.
        Set::Patch<const Set::Scalar> rho_phys  = solid.rho_phys_mf.Patch(lev,mfi);

        Set::Patch<Set::Scalar>       omega     = vorticity_mf.Patch(lev,mfi);

        Set::Patch<const Set::Scalar> eta_patch = eta_old_mf->Patch(lev,mfi);
        Set::Patch<const Set::Scalar> etadot    = etadot_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> velocity  = velocity_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> T         = temperature_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> molef     = mole_fraction_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> pressure  = pressure_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> temp      = temperature_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> scratch   = scratch_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> rho_fluid = rho_fluid_mf.array(mfi);
        Set::Patch<const Set::Scalar> M_fluid   = M_fluid_mf.array(mfi);
        Set::Patch<const Set::Scalar> E_fluid   = E_fluid_mf.array(mfi);

        Set::Patch<const Set::Scalar> m0        = m0_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> q         = q_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> _u0       = u0_mf.Patch(lev,mfi);

        amrex::Array4<Set::Scalar> const& Source = (*Source_mf[lev]).array(mfi);

        Set::Patch<Set::Scalar>       rho_sum   = rho_sum_mf.array(mfi);
        Set::Patch<Set::Scalar>       mixed_kT  = mixed_kT_mf.array(mfi);
        Set::Patch<Set::Scalar>       rhoHDYx   = rhoHDYx_mf.array(mfi);
        Set::Patch<Set::Scalar>       rhoHDYy   = rhoHDYy_mf.array(mfi);
        Set::Patch<Set::Scalar>       rhoDYx    = rhoDYx_mf.array(mfi);
        Set::Patch<Set::Scalar>       rhoDYy    = rhoDYy_mf.array(mfi);
        #if AMREX_SPACEDIM == 3
        Set::Patch<Set::Scalar>       rhoHDYz   = rhoHDYz_mf.array(mfi);
        Set::Patch<Set::Scalar>       rhoDYz    = rhoDYz_mf.array(mfi);
        #endif

        if (details)
        {
            wdot_arr = wdot_mf.Patch(lev,mfi);
            qdot_arr = qdot_mf.Patch(lev,mfi);
        }

        // Third and final ParallelFor loop to get 2nd gradients and compute fluxes
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            // Conservative ghost cells and derived ghost cells were filled for
            // this RHS stage, so boundary valid cells can use centered stencils
            // and Riemann states that sample the specified physical BCs.
            auto sten = BCGhostStencil();

            Set::Scalar eta = invert ? 1.0-eta_patch(i,j,k) : eta_patch(i,j,k);

            //Diffuse Sources
            Set::Vector grad_eta     = Numeric::Gradient(eta_patch, i, j, k, 0, DX, sten);
            Set::Scalar grad_eta_mag = grad_eta.lpNorm<2>();
            Set::Matrix hess_eta     = Numeric::Hessian(eta_patch, i, j, k, 0, DX, sten);
            if (invert) grad_eta *= -1.0;
            if (invert) hess_eta *= -1.0;
            const bool cutoff_enabled = cutoff >= 0.0 && cutoff < 1.0;
            const Set::Scalar eta_cutoff = cutoff_enabled ? cutoff : small;
            const Set::Scalar source_cutoff = cutoff_enabled ? cutoff : 0.0;
            Set::Scalar source_delta = eta > source_cutoff ? diffuse_source_norm * grad_eta_mag : 0.0;
            Set::Scalar source_scale = source_delta / (grad_eta_mag + small);

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
            // The cell-centered eta flux scaling omits the advective
            // F dot grad(eta) term. Add it for transported conserved quantities.
            Set::Scalar eta_transport_weight = cutoff_enabled ? 1.0 - cutoff : 1.0;
            Set::Scalar eta_transport_rate = source_delta > 0.0 ? eta_transport_weight * u.dot(grad_eta) : 0.0;

            Set::Matrix gradM        = Numeric::Gradient(M, i, j, k, DX, sten);
            Set::Vector gradrho      = Numeric::Gradient(rho_sum,i,j,k,0,DX, sten);
            Set::Matrix hess_rho     = Numeric::Hessian(rho_sum,i,j,k,0,DX,sten);
            Set::Matrix gradu        = (gradM - u*gradrho.transpose()) / rho_sum(i,j,k);

            Set::Scalar div_mixed_kT =
                Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(mixed_kT,i,j,k,0,DX,sten)
#if AMREX_SPACEDIM >= 2
                + Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(mixed_kT,i,j,k,1,DX,sten)
#endif
#if AMREX_SPACEDIM == 3
                + Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(mixed_kT,i,j,k,2,DX,sten)
#endif
                ;
            // Species diffusion divergences are computed direction-by-direction
            // below so transverse ghost-line values are not read.

            if (prescribedflowmode == PrescribedFlowMode::Relative)
            {
                Set::Vector N = grad_eta / (grad_eta_mag + small);
                // Set::Vector T(N(1), -N(0));
                // u0 = N * u0(0) + T * u0(1);
                // Normal is u0(0), tangential is u0(1)

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


            std::vector<Set::Scalar> mdot0(NSPECIES);
            Set::Scalar mdot0_total = 0.0;
            for (int n=0; n<NSPECIES; ++n )
            {
                mdot0[n] = m0(i,j,k,n)*source_delta;
                mdot0_total += mdot0[n];
            }

            Set::Vector Pdot0 = mdot0_total * u0;
            Set::Scalar qdot0 = q0.dot(grad_eta) * source_scale;
            // E_solid = rho_phys*cp*(T-T_ref) is built from the *physical* solid
            // density (see Flame::UpdateFluxes / Hydro::AdvanceSolidEnergy), so its
            // specific internal energy is E_solid/rho_phys - not E_solid summed over
            // species and divided by the numerically-tame Riemann-blend density
            // (solid.density_mf, hydro.rho_ap/htpb) used only for the fluid/solid
            // conserved-state blend. Dividing by the tame density here over-injected
            // internal energy by ~rho_phys/rho_tame (two orders of magnitude). This
            // matches the same E_solid/(rho_phys*cp) convention used by the C1
            // temperature blend above.
            if (rho_phys(i,j,k) > small)
            {
                qdot0 += mdot0_total * E_solid(i,j,k) / rho_phys(i,j,k);
            }
            // The momentum source Pdot0 = mdot0_total*u0 injects mass moving at u0,
            // which raises the fluid's kinetic energy (M_rhs picks up Pdot0 via
            // Source(momentum_source_comp)). Without a matching kinetic-energy flux
            // here, the energy budget at the interface is unbalanced by an amount
            // that depends on the local flow velocity u (0.5*|u|^2 - u.u0), which
            // flips sign across the diffuse interface and was found to drive a
            // growing, sign-flipping internal-energy (temperature) oscillation
            // there. Adding the kinetic energy of the injected stream itself here
            // makes the net induced internal-energy source
            // mdot0_total*(e_int + 0.5*|u-u0|^2) - a bounded, non-negative mixing
            // dissipation instead of a sign-indefinite drain.
            qdot0 += 0.5 * mdot0_total * u0.dot(u0);

            Set::Scalar mu = gas.dynamic_viscosity(T(i,j,k), molef, i, j, k);

            Set::Matrix3 hess_M = Numeric::Hessian(M,i,j,k,DX,sten);
            Set::Matrix3 hess_u = Set::Matrix3::Zero();
            for (int p = 0; p < AMREX_SPACEDIM; p++)
                for (int q = 0; q < AMREX_SPACEDIM; q++)
                    for (int r = 0; r < AMREX_SPACEDIM; r++)
                    {
                        hess_u(r,p,q) =
                            (hess_M(r,p,q) - gradu(r,q)*gradrho(p) - gradu(r,p)*gradrho(q) - u(r)*hess_rho(p,q))
                            / rho_sum(i,j,k);
                    }

            Set::Vector Ldot0 = Set::Vector::Zero();
            Set::Vector div_tau = Set::Vector::Zero();
            Set::Scalar lambda = 0.0; //-2.0/3.0*mu_eff;
            for (int p = 0; p<AMREX_SPACEDIM; p++)
                for (int q = 0; q<AMREX_SPACEDIM; q++)
                    for (int r = 0; r<AMREX_SPACEDIM; r++)
                        for (int s = 0; s<AMREX_SPACEDIM; s++)
                        {
                            Ldot0(p) += 0.25 * (mu * ((p==r && q==s) + (p==s && q==r)) + lambda * (p==q && r==s)) * (u(r) - u0(r)) * hess_eta(q, s);
                            div_tau(p) += 0.5 * (mu * ((p==r && q==s) + (p==s && q==r)) + lambda * (p==q && r==s)) * (hess_u(r,q,s) + hess_u(s,q,r));

                        }

            Set::Scalar E_fluid_transport = E_fluid(i,j,k);

            // Transport flux
            Solver::Local::Riemann::FluxDivergence transport_flux;
            if (flux_scheme == FluxScheme::Riemann)
            {
                const int X = 0, Y = 1;
#if AMREX_SPACEDIM == 3
                const int Z = 2;
#endif

                Solver::Local::Riemann::State state_xlo_fluid(rho_fluid, M_fluid, E_fluid, i-1, j, k, X);
                Solver::Local::Riemann::State state_x_fluid  (rho_fluid, M_fluid, E_fluid, i  , j, k, X);
                Solver::Local::Riemann::State state_xhi_fluid(rho_fluid, M_fluid, E_fluid, i+1, j, k, X);

                Solver::Local::Riemann::State state_ylo_fluid(rho_fluid, M_fluid, E_fluid, i, j-1, k, Y);
                Solver::Local::Riemann::State state_y_fluid  (rho_fluid, M_fluid, E_fluid, i, j  , k, Y);
                Solver::Local::Riemann::State state_yhi_fluid(rho_fluid, M_fluid, E_fluid, i, j+1, k, Y);
#if AMREX_SPACEDIM == 3
                Solver::Local::Riemann::State state_zlo_fluid(rho_fluid, M_fluid, E_fluid, i, j, k-1, Z);
                Solver::Local::Riemann::State state_z_fluid  (rho_fluid, M_fluid, E_fluid, i, j, k  , Z);
                Solver::Local::Riemann::State state_zhi_fluid(rho_fluid, M_fluid, E_fluid, i, j, k+1, Z);
#endif

                try
                {
                    transport_flux = riemannsolver->ComputeFluxDivergence(
                        state_xlo_fluid, state_x_fluid, state_xhi_fluid,
                        state_ylo_fluid, state_y_fluid, state_yhi_fluid,
#if AMREX_SPACEDIM == 3
                        state_zlo_fluid, state_z_fluid, state_zhi_fluid,
#endif
                        gas, molef, i, j, k, eta, DX, small);
                }
                catch(...)
                {
                    Util::ParallelMessage(INFO,"lev=",lev);
                    Util::ParallelMessage(INFO,"i=",i,"j=",j,"k=",k);
                    Util::Abort(INFO);
                }
            }
            else if (flux_scheme == FluxScheme::Advect)
            {
                for (int n=0; n<NSPECIES; ++n)
                {
                    transport_flux.mass[n] = eta * advect_op.Scalar(rho_fluid, velocity, i, j, k, n, DX, conservative_options, sten);
                }

                Set::Vector grad_pressure = Numeric::Gradient(pressure, i, j, k, 0, DX, sten);
                Set::Scalar flux_pressure = advect_op.Scalar(pressure, velocity, i, j, k, 0, DX, conservative_options, sten);
                transport_flux.momentum = eta * (
                    advect_op.Vector(M_fluid, velocity, i, j, k, 0, DX, conservative_options, sten) -
                    grad_pressure);
                transport_flux.energy = eta * (
                    advect_op.Scalar(E_fluid, velocity, i, j, k, 0, DX, conservative_options, sten) +
                    flux_pressure);
            }
            else
            {
                Util::Abort(INFO, "Unknown Hydro flux scheme");
            }

            const int momentum_source_comp = NSPECIES;
            const int energy_source_comp = NSPECIES + AMREX_SPACEDIM;
            std::array<Set::Scalar, NSPECIES> drhof_dt_hydro;
            for (int n=0; n<NSPECIES; ++n)
            {
                Source(i,j, k, n) = mdot0[n];
                drhof_dt_hydro[n] =
                    transport_flux.mass[n] +
                    Source(i, j, k, n) -
                    scratch(i,j,k,n) * eta_transport_rate;
                if (NSPECIES > 1)
                {
                    // species diffusion term, d/dx_i(rho*DKM*Y,i)
                    Set::Scalar div_rhoDY =
                        Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(rhoDYx,i,j,k,n,DX,sten)
#if AMREX_SPACEDIM >= 2
                        + Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(rhoDYy,i,j,k,n,DX,sten)
#endif
#if AMREX_SPACEDIM == 3
                        + Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(rhoDYz,i,j,k,n,DX,sten)
#endif
                        ;
                    drhof_dt_hydro[n] += eta * div_rhoDY;
                }
            }

            Set::Vector interface_momentum_source = -Ldot0;

            Source(i,j, k, momentum_source_comp  ) = Pdot0(0) + interface_momentum_source(0);
            Source(i,j, k, momentum_source_comp+1) = Pdot0(1) + interface_momentum_source(1);
#if AMREX_SPACEDIM == 3
            Source(i,j, k, momentum_source_comp+2) = Pdot0(2) + interface_momentum_source(2);
#endif
            Source(i,j, k, energy_source_comp) = qdot0 + interface_momentum_source.dot(u);

            // Lagrange terms to enforce no-penetration
            Set::Vector lagrange_momentum_source =
                -lagrange*(u-u0).dot(grad_eta)*grad_eta;
            Source(i,j,k,momentum_source_comp  ) += lagrange_momentum_source(0);
            Source(i,j,k,momentum_source_comp+1) += lagrange_momentum_source(1);
#if AMREX_SPACEDIM == 3
            Source(i,j,k,momentum_source_comp+2) += lagrange_momentum_source(2);
#endif
            Source(i,j,k,energy_source_comp) += lagrange_momentum_source.dot(u);

            std::array<Set::Scalar, NSPECIES> rhoY_intermediate;
            for (int n=0; n<NSPECIES; ++n)
            {
                rhoY_intermediate[n] = scratch(i,j,k,n);
            }
            if (eta > eta_cutoff)
            {
                for (int n=0; n<NSPECIES; ++n)
                {
                    rhoY_intermediate[n] += drhof_dt_hydro[n] * dt / (eta + small);
                }
                if (!TryProjectSpeciesDensities(rhoY_intermediate))
                {
                    for (int n=0; n<NSPECIES; ++n)
                    {
                        rhoY_intermediate[n] = scratch(i,j,k,n);
                    }
                    ProjectSpeciesDensities(rhoY_intermediate, i, j, k);
                }
            }

            std::array<Set::Scalar, NSPECIES> wdot;
            Set::Scalar qdot = 0.0;
            std::tie(wdot, qdot) = chemistry.ComputeChemistrySources(pressure(i,j,k), temp(i,j,k), rhoY_intermediate, dt, &gas);

            if (details)
            {
                qdot_arr(i,j,k) = qdot;
                for (int n=0; n<NSPECIES; ++n) wdot_arr(i,j,k,n) = wdot[n];
            }

            for (int n=0; n<NSPECIES; ++n)
            {
                rho_rhs(i,j,k,n) =
                        drhof_dt_hydro[n] + eta * wdot[n] +
                        etadot(i,j,k) * (rho(i,j,k,n) - rho_solid(i,j,k,n)) / (eta + small)
                    ;
            }

            Set::Scalar dMxf_dt =
                transport_flux.momentum(0) +
                div_tau(0) * eta +
                g(0)*rho_sum(i,j,k) +
                Source(i, j, k, momentum_source_comp) -
                rho_sum(i,j,k) * u(0) * eta_transport_rate;

            M_rhs(i,j,k,0) =
                //M_new(i, j, k, 0) = M(i, j, k, 0) +
                // (
                    dMxf_dt +
                    // todo add dMs_dt term if want time-evolving Ms
                    etadot(i,j,k)*(M(i,j,k,0) - M_solid(i,j,k,0)) / (eta + small)
                // ) * dt;
                ;

            Set::Scalar dMyf_dt =
                transport_flux.momentum(1) +
                div_tau(1) * eta +
                g(1)*rho_sum(i,j,k) +
                Source(i, j, k, momentum_source_comp+1) -
                rho_sum(i,j,k) * u(1) * eta_transport_rate;

            M_rhs(i,j,k,1) =
                //M_new(i, j, k, 1) = M(i, j, k, 1) +
                //(
                    dMyf_dt +
                    // todo add dMs_dt term if want time-evolving Ms
                    etadot(i,j,k)*(M(i,j,k,1) - M_solid(i,j,k,1)) / (eta+small)
                // )*dt;
                ;

#if AMREX_SPACEDIM == 3
            Set::Scalar dMzf_dt =
                transport_flux.momentum(2) +
                div_tau(2) * eta +
                g(2)*rho_sum(i,j,k) +
                Source(i, j, k, momentum_source_comp+2) -
                rho_sum(i,j,k) * u(2) * eta_transport_rate;

            M_rhs(i,j,k,2) =
                    dMzf_dt +
                    etadot(i,j,k)*(M(i,j,k,2) - M_solid(i,j,k,2)) / (eta+small)
                ;
#endif

            Set::Scalar dEf_dt =
                transport_flux.energy +
                eta * (div_tau.dot(u) + div_mixed_kT) +
                rho_sum(i,j,k)*g.dot(u) +
                Source(i, j, k, energy_source_comp) +
                eta * qdot -
                E_fluid_transport * eta_transport_rate;

            if (NSPECIES > 1)
            {
                for (int n=0; n<NSPECIES; ++n)
                {
                    // Species energy diffusion term: d/dx_i(rho*H*DKM*Y,i)
                    Set::Scalar div_rhoHDY =
                        Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(rhoHDYx,i,j,k,n,DX,sten)
#if AMREX_SPACEDIM >= 2
                        + Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(rhoHDYy,i,j,k,n,DX,sten)
#endif
#if AMREX_SPACEDIM == 3
                        + Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(rhoHDYz,i,j,k,n,DX,sten)
#endif
                        ;
                    dEf_dt += eta * div_rhoHDY;
                }
            }

            E_rhs(i,j,k) =
            // E_new(i, j, k) = E(i, j, k) +
            //     (
                    dEf_dt +
                    // todo add dEs_dt term if want time-evolving Es
                    etadot(i,j,k)*(E(i,j,k) - E_solid(i,j,k)) / (eta+small)
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
                //Util::ParallelMessage(INFO,"drhof_dt ",drhof_dt); // dies
                Util::ParallelMessage(INFO,"eta ",eta);
                Util::ParallelMessage(INFO,"etadot ",etadot(i,j,k));
                Util::ParallelMessage(INFO,"Source ",Source(i,j,k,0));
                Util::ParallelMessage(INFO,"transport mass[0] ",transport_flux.mass[0]);
                Util::ParallelMessage(INFO,"transport Mx ",transport_flux.momentum(0));
                Util::ParallelMessage(INFO,"transport My ",transport_flux.momentum(1));
                Util::ParallelMessage(INFO,"transport E ",transport_flux.energy);

                Util::ParallelMessage(INFO,"Mx_solid ",M_solid(i,j,k,0));
                Util::ParallelMessage(INFO,"My_solid ",M_solid(i,j,k,1));
                Util::ParallelMessage(INFO,"small ",small);
                Util::ParallelMessage(INFO,"Mx ",M(i,j,k,0));
                Util::ParallelMessage(INFO,"My ",M(i,j,k,1));
                Util::ParallelMessage(INFO,"dMx/dt ",dMxf_dt);
                Util::ParallelMessage(INFO,"dMy/dt ",dMyf_dt);


                Util::Message(INFO,transport_flux.momentum(0));
                Util::Message(INFO,transport_flux.momentum(1));
                Util::Message(INFO,DX[0]);
                Util::Message(INFO,transport_flux.energy);
                Util::Message(INFO,transport_flux.mass[0]);
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
            #if AMREX_SPACEDIM == 2
            omega(i, j, k) = eta * (gradu(1,0) - gradu(0,1));
            #elif AMREX_SPACEDIM == 3
            omega(i, j, k, 0) = eta * (gradu(2,1) - gradu(1,2));
            omega(i, j, k, 1) = eta * (gradu(0,2) - gradu(2,0));
            omega(i, j, k, 2) = eta * (gradu(1,0) - gradu(0,1));
            #endif
        });
    }
}

void Hydro::Regrid(int lev, Set::Scalar /* time */)
{
    BL_PROFILE("Integrator::Hydro::Regrid");
    Source_mf[lev]->setVal(0.0);
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
            #if AMREX_SPACEDIM == 2
            Set::Vector grad_omega = Numeric::Gradient(omega, i, j, k, 0, DX, sten);
            #elif AMREX_SPACEDIM == 3
            Set::Matrix grad_omega = Numeric::Gradient(omega, i, j, k, DX, sten);
            #endif
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

    // Temperature criterion for refinement
    for (amrex::MFIter mfi(*temperature_mf[lev], true); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<char> const& tags = a_tags.array(mfi);
        amrex::Array4<const Set::Scalar> const& temp = (*temperature_mf[lev]).array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            auto sten = Numeric::GetStencil(i, j, k, bx);
            Set::Vector grad_temp = Numeric::Gradient(temp, i, j, k, 0, DX, sten);
            if (grad_temp.lpNorm<2>() * dr * 2 > temp_refinement_criterion) tags(i, j, k) = amrex::TagBox::SET;
        });
    }
}

}
