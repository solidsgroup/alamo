
#include "Hydro.H"
#include "IO/ParmParse.H"
#include "BC/Constant.H"
#include "BC/Expression.H"
#include "Numeric/Stencil.H"
#include "IC/Constant.H"
#include "IC/Laminate.H"
#include "IC/Expression.H"
#include "IC/BMP.H"
#include "IC/PNG.H"
#include <typeinfo>
//#include "Solver/Local/Riemann/Roe.H"
#include "Solver/Local/Riemann/HLLE.H"
#include "Solver/Local/Riemann/HLLC.H"
#if AMREX_SPACEDIM == 2

namespace Integrator
{

double collision_integral( double T )
{
    return 1.06036/pow(T, 0.15610) + 0.19300/exp(0.47635*T) + 1.03587/exp(1.52996*T) + 1.76474/exp(3.89411*T);
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

        std::string scheme_str;
        // time integration scheme to use
        pp.query_validate("scheme",scheme_str, {"forwardeuler","ssprk3","rk4"});
        if (scheme_str == "forwardeuler") value.scheme = IntegrationScheme::ForwardEuler;
        else if (scheme_str == "ssprk3") value.scheme = IntegrationScheme::SSPRK3;
        else if (scheme_str == "rk4") value.scheme = IntegrationScheme::RK4;

        if (pp.contains("restart")) value.restart_found = true;
        std::cout << "restart " << value.restart_found << std::endl;

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

        pp_query_required("gamma", value.gamma); // gamma for gamma law
        pp_query_required("cfl", value.cfl); // cfl condition
        pp_query_default("cfl_v", value.cfl_v,1E100); // cfl condition
        pp_query_required("mu", value.mu); // linear viscosity coefficient
        pp_forbid("Lfactor","replaced with mu");
        //pp_query_default("Lfactor", value.Lfactor,1.0); // (to be removed) test factor for viscous source
        pp_forbid("Pfactor","replaced with mu");
        //pp_query_default("Pfactor", value.Pfactor,1.0); // (to be removed) test factor for viscous source
        pp_query_default("pref", value.pref,1.0); // reference pressure for Roe solver
        pp_query_default("tref", value.tref,1.0); // reference pressure for Roe solver

        pp_forbid("rho.bc","--> density.bc");
        pp_forbid("p.bc","--> pressure.bc");
        pp_forbid("v.bc", "--> velocity.bc");
        pp_forbid("pressure.bc","--> energy.bc");
        pp_forbid("velocity.bc","--> momentum.bc");

        // Species inputs
        pp_queryarr_default("species", value.species, {"N2", "O2"});
        pp_queryarr_default("species_mw", value.species_mw, {28.0134, 31.9988}); // g/mol or kg/kmol
        value.nspecies = value.species.size();

        pp_queryarr_default("species_k", value.species_k, {2.623368E-2, 2.560608E-2}); // W/m-K, thermal conductivity
        pp_queryarr_default("species_cp", value.species_cp, {1040.0, 920.0}); // J/kg-K, specific heat by mass
        pp_queryarr_default("species_mu", value.species_mu, {175.4E-7, 203.1E-7}); //kg/m-s, dynamic viscosity
        pp_queryarr_default("species_LJdiameter", value.species_LJdiameter, {3.667, 3.433}); // Angstroms, Lennard-Jones potential collision diameter
        pp_queryarr_default("species_LJwelldepth", value.species_LJwelldepth, {99.8, 113.0}); // K, Lennard-Jones potential e/k (k: Boltzmann constant)

        // Boundary condition for density
        pp.select_default<BC::Constant,BC::Expression>("density.bc",value.density_bc,value.nspecies);
        // Boundary condition for energy
        pp.select_default<BC::Constant,BC::Expression>("energy.bc",value.energy_bc,1);
        // Boundary condition for momentum
        pp.select_default<BC::Constant,BC::Expression>("momentum.bc",value.momentum_bc,2);

        if (!value.managed)
        {
            // Boundary condition for phase field order parameter
            pp.select_default<BC::Constant,BC::Expression>("pf.eta.bc",value.eta_bc,1);
        }

        // // Boundary condition for tracer field
        // pp.select_default<BC::Constant,BC::Expression>("tracer.bc",value.tracer_bc,1);

        pp_query_default("small",value.small,1E-8); // small regularization value
        pp_query_default("cutoff",value.cutoff,-1E100); // cutoff value
        pp_query_default("lagrange",value.lagrange,0.0); // lagrange no-penetration factor
        pp_query_default("lagrange_m0",value.lagrange_m0,1.0); // lagrange no-penetration factor

        pp_forbid("roefix","--> solver.roe.entropy_fix"); // Roe solver entropy fix

    }
    // Register FabFields:
    {
        int nghost = 1;

        if (!value.managed)
        {
            value.eta_mf = new Set::Field<Set::Scalar>();
            value.eta_old_mf = new Set::Field<Set::Scalar>();
            value.RegisterNewFab(*value.eta_mf,     value.eta_bc, 1, nghost, "eta",     true, true);
            value.RegisterNewFab(*value.eta_old_mf, value.eta_bc, 1, nghost, "eta_old", true, true);
        }
        value.RegisterNewFab(value.etadot_mf,  value.eta_bc, 1, nghost, "etadot",  true, false);

        // If nspecies > 1, then replace continuity equation with N multicomponent continuity equations where rho[i] := rho*Y[i]
        value.RegisterNewFab(value.density_mf,     value.density_bc, value.nspecies, nghost, "density",     true , true);
        value.RegisterNewFab(value.density_old_mf, value.density_bc, value.nspecies, nghost, "density_old", false, true);

        value.RegisterNewFab(value.energy_mf,     value.energy_bc, 1, nghost, "energy",      true ,true);
        value.RegisterNewFab(value.energy_old_mf, value.energy_bc, 1, nghost, "energy_old" , false, true);

        value.RegisterNewFab(value.momentum_mf,     value.momentum_bc, 2, nghost, "momentum",     true ,true, {"x","y"});
        value.RegisterNewFab(value.momentum_old_mf, value.momentum_bc, 2, nghost, "momentum_old", false, true);
 
        // value.RegisterNewFab(value.tracer_mf,     value.tracer_bc, 1, nghost, "tracer",     true ,true);
        // value.RegisterNewFab(value.tracer_old_mf, value.tracer_bc, 1, nghost, "tracer_old", false);

        value.RegisterNewFab(value.pressure_mf,     &value.bc_nothing, 1, nghost, "pressure",  true, false);
        value.RegisterNewFab(value.temperature_mf,  &value.bc_nothing, 1, nghost, "temperature",  true, false);
        value.RegisterNewFab(value.velocity_mf,     &value.bc_nothing, 2, nghost, "velocity",  true, false,{"x","y"});
        value.RegisterNewFab(value.vorticity_mf,    &value.bc_nothing, 1, nghost, "vorticity", true, false);

        // If nspecies > 1, then m0 also needs nspecies source terms
        value.RegisterNewFab(value.m0_mf,           &value.bc_nothing, value.nspecies, 0, "m0",  true, false);
        value.RegisterNewFab(value.u0_mf,           &value.bc_nothing, 2, 0, "u0",  true, false, {"x","y"});
        value.RegisterNewFab(value.q_mf,            &value.bc_nothing, 2, 0, "q",   true, false, {"x","y"});

        // Solid density terms for nspecies (needed in diffuse region for smooth convergence)
        value.neumann_bc_N = new BC::Constant(value.nspecies);
      	*value.neumann_bc_N = BC::Constant::ZeroNeumann(value.nspecies);
        value.RegisterNewFab(value.solid.density_mf,  value.neumann_bc_N, value.nspecies, nghost, "solid.density", true, false);
        value.RegisterNewFab(value.solid.momentum_mf, &value.neumann_bc_D, 2, nghost, "solid.momentum", true, false, {"x","y"});
        value.RegisterNewFab(value.solid.energy_mf,   &value.neumann_bc_1, 1, nghost, "solid.energy",   true, false);

        // nspecies source terms for continuity, plus 2 for 2D momentum, plus 1 for energy
        value.RegisterNewFab(value.Source_mf, &value.bc_nothing, value.nspecies+3, 0, "Source", true, false);
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

    // density initial condition type
    // pp.select_default<IC::Constant,IC::Expression>("tracer.ic",value.tracer_ic,value.geom);


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

    // Riemann solver
    pp.select_default</*Solver::Local::Riemann::Roe,*/
                      Solver::Local::Riemann::HLLE,
                      Solver::Local::Riemann::HLLC>("solver",value.riemannsolver);


    std::string prescribedflowmode_str;
    // 
    pp.query_validate("prescribedflowmode",prescribedflowmode_str,{"absolute","relative"});
    if (prescribedflowmode_str == "absolute") value.prescribedflowmode = PrescribedFlowMode::Absolute;
    else if (prescribedflowmode_str == "relative") value.prescribedflowmode = PrescribedFlowMode::Relative;



    Util::Message(INFO);
    pp.queryarr_default("g",value.g,Set::Vector::Zero());
    Util::Message(INFO);


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
    density_ic       ->Initialize(lev, density_mf, 0.0);

    density_ic       ->Initialize(lev, density_old_mf, 0.0);

    // tracer_ic       ->Initialize(lev, tracer_old_mf, 0.0);
    // tracer_ic       ->Initialize(lev, tracer_mf, 0.0);

    solid.density_ic ->Initialize(lev, solid.density_mf, 0.0);
    solid.momentum_ic->Initialize(lev, solid.momentum_mf, 0.0);
    solid.energy_ic  ->Initialize(lev, solid.energy_mf, 0.0);

    ic_m0            ->Initialize(lev, m0_mf, 0.0);
    ic_u0            ->Initialize(lev, u0_mf, 0.0);
    ic_q             ->Initialize(lev, q_mf,  0.0);

    Source_mf[lev]   ->setVal(0.0);


    if (lev >= (int)mixed.size()) mixed.push_back(false);
}

void Hydro::Mix(int lev)
{
    if (restart_found) {
        if (lev >= (int)mixed.size()) mixed.push_back(true);
        else mixed[lev] = true;
    }
    if (mixed[lev])  return;

    Util::Message(INFO,"MIXING");
    for (amrex::MFIter mfi(*velocity_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox();

        Set::Patch<const Set::Scalar> eta_patch = eta_old_mf->Patch(lev,mfi);

        Set::Patch<const Set::Scalar> v         = velocity_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> p         = pressure_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       rho       = density_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       rho_old   = density_old_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       M         = momentum_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       M_old     = momentum_old_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       E         = energy_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       E_old     = energy_old_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> rho_solid = solid.density_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> M_solid   = solid.momentum_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> E_solid   = solid.energy_mf.Patch(lev,mfi);


        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {  
            Set::Scalar eta = invert ? 1.0-eta_patch(i,j,k) : eta_patch(i,j,k);
            double rho_sum = 0.0;

            rho_sum = 0.0;
            for (int n=0; n<nspecies; ++n) {
                rho(i, j, k, n) = eta * rho(i, j, k, n) + (1.0 - eta) * rho_solid(i, j, k, n);
                rho_old(i, j, k, n) = rho(i, j, k, n);
                rho_sum += rho(i,j,k,n);
            }
            double cp_mix = 0.0;
            double cv_mix = 0.0;
            for (int n=0; n<nspecies; ++n) {
                cp_mix += rho(i,j,k,n)/rho_sum * species_cp[n];
                cv_mix += rho(i,j,k,n)/rho_sum * (species_cp[n] - 8314.45/species_mw[n]);
            }
            gamma = cp_mix/cv_mix;

            M(i, j, k, 0) = (rho_sum*v(i, j, k, 0))*eta  +  M_solid(i, j, k, 0)*(1.0-eta);
            M(i, j, k, 1) = (rho_sum*v(i, j, k, 1))*eta  +  M_solid(i, j, k, 1)*(1.0-eta);
            M_old(i, j, k, 0) = M(i, j, k, 0);
            M_old(i, j, k, 1) = M(i, j, k, 1);

            E(i, j, k) =
                (0.5 * (v(i, j, k, 0) * v(i, j, k, 0) + v(i, j, k, 1) * v(i, j, k, 1)) * rho_sum + p(i, j, k) / (gamma - 1.0)) * eta 
                + 
                E_solid(i, j, k) * (1.0 - eta);
            E_old(i, j, k) = E(i, j, k);

        });
    }
    c_max = 0.0;
    vx_max = 0.0;
    vy_max = 0.0;
    mixed[lev] = true;
}

void Hydro::UpdateEta(int lev, Set::Scalar time)
{
    Util::Assert(INFO,TEST(!managed),"Should override this if Hydro is managed!");
    eta_ic->Initialize(lev, *eta_mf, time);
}

void Hydro::UpdateFluxes(int lev, Set::Scalar time)
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

    if (!managed)
    {
        std::swap(*eta_old_mf, *eta_mf);
        UpdateEta(lev, time);
    }
    UpdateFluxes(lev,time);

    Mix(lev);

    std::swap(density_old_mf[lev],  density_mf[lev]);
    std::swap(momentum_old_mf[lev], momentum_mf[lev]);
    std::swap(energy_old_mf[lev],   energy_mf[lev]);
    // std::swap(tracer_old_mf[lev],   tracer_mf[lev]);
    Set::Scalar dt_max = std::numeric_limits<Set::Scalar>::max();

    //ic_u0->Initialize(lev, u0_mf,    time);

    for (amrex::MFIter mfi(*(velocity_mf)[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox();
        amrex::Array4<const Set::Scalar> const& eta_new = (*(*eta_mf)[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& eta = (*(*eta_old_mf)[lev]).array(mfi);
        amrex::Array4<Set::Scalar>       const& etadot = (*etadot_mf[lev]).array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            // etadot(i, j, k) = (eta_new(i, j, k) - eta(i, j, k)) / dt;
            // if (invert) etadot(i,j,k) *= 1.0;
            etadot(i,j,k) = 0.0;
        });
    }
    if (!lev) Util::Warning(INFO,"zeroing out etadot for the moment");


    if (scheme == IntegrationScheme::ForwardEuler) // forward euler
    {

        RHS(lev, time, 
            *density_mf[lev],     *momentum_mf[lev],     *energy_mf[lev],   
            *density_old_mf[lev], *momentum_old_mf[lev], *energy_old_mf[lev]);

        for (amrex::MFIter mfi(*(*eta_mf)[lev], false); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox();
        
            Set::Patch<const Set::Scalar> rho_rhs  = density_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> E_rhs    = energy_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> M_rhs    = momentum_mf.Patch(lev,mfi);

            Set::Patch<const Set::Scalar> rho_old  = density_old_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> E_old    = energy_old_mf.Patch(lev,mfi);
            Set::Patch<const Set::Scalar> M_old    = momentum_old_mf.Patch(lev,mfi);

            Set::Patch<Set::Scalar> rho_new       = density_mf.Patch(lev,mfi);
            Set::Patch<Set::Scalar> E_new         = energy_mf.Patch(lev,mfi);
            Set::Patch<Set::Scalar> M_new         = momentum_mf.Patch(lev,mfi);
        
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {   
                for (int n=0; n<nspecies; ++n) {
                    rho_new(i,j,k,n) = rho_old(i,j,k,n) + dt * rho_rhs(i,j,k,n);
                }
                M_new(i,j,k,0) = M_old(i,j,k,0)     + dt * M_rhs(i,j,k,0);
                M_new(i,j,k,1) = M_old(i,j,k,1)     + dt * M_rhs(i,j,k,1);
                E_new(i,j,k)     = E_old(i,j,k)     + dt * E_rhs(i,j,k);
            });
        }
    }


//    else if (scheme == IntegrationScheme::SSPRK3)
//    {
//        // Butcher Tableau
//        //     |
//        //  1  |  1    
//        // 1/2 | 1/4  1/4
//        // ---------------------
//        //     | 1/6  1/6  2/3
//        
//        Set::Scalar 
//            /* */  
//            /* */  c2 = 1.0 ,    a21 = 1.0,  
//            /* */  c3 = 0.5,     a31 = 0.25, a32 = 0.25, 
//            /*     ---------------------------------------------    */
//            /* */                b1 = 1./6,  b2 = 1./6., b3 = 2./3.;
//
//
//        const amrex::BoxArray &ba = density_mf[lev]->boxArray();
//        const amrex::DistributionMapping &dm = density_mf[lev]->DistributionMap();
//        const int ng = density_mf[lev]->nGrow();
//
//        // handles to old solution
//        const amrex::MultiFab &density_old = *density_old_mf[lev];
//        const amrex::MultiFab &momentum_old = *momentum_old_mf[lev];
//        const amrex::MultiFab &energy_old = *energy_old_mf[lev];
//
//        
//        // temporary storage
//        amrex::MultiFab density_k1(ba,dm,1,0), momentum_k1(ba,dm,2,0), energy_k1(ba,dm,1,0);
//        amrex::MultiFab density_k2(ba,dm,1,0), momentum_k2(ba,dm,2,0), energy_k2(ba,dm,1,0);
//        amrex::MultiFab density_k3(ba,dm,1,0), momentum_k3(ba,dm,2,0), energy_k3(ba,dm,1,0);
//            
//        // buffer to hold combs of k1
//        amrex::MultiFab density_temp(ba,dm,1,ng), momentum_temp(ba,dm,2,ng), energy_temp(ba,dm,1,ng);
//
//        // fill the ghost cells from the _old fields, which were updated from the coarse patch.
//        density_temp.ParallelCopyToGhost(*density_old_mf[lev],0,0,1,amrex::IntVect(1),amrex::IntVect(1));
//        momentum_temp.ParallelCopyToGhost(*momentum_old_mf[lev],0,0,2,amrex::IntVect(1),amrex::IntVect(1));
//        energy_temp.ParallelCopyToGhost(*energy_old_mf[lev],0,0,1,amrex::IntVect(1),amrex::IntVect(1));
//
//        // handles to new solution
//        amrex::MultiFab &density_new = *density_mf[lev];
//        amrex::MultiFab &momentum_new = *momentum_mf[lev];
//        amrex::MultiFab &energy_new = *energy_mf[lev];
//
//
//        //
//        // Calculate K1
//        //
//        // k1 = RHS(t, yold)
//
//        RHS(lev,time,
//            density_k1,momentum_k1,energy_k1,
//            *density_old_mf[lev], *momentum_old_mf[lev], *energy_old_mf[lev]);
//
//
//        //
//        // Calculate K2
//        // 
//        // ytemp = yold + dt*( a21*k1 )
//        amrex::MultiFab::LinComb(density_temp,  1.0, density_old,  0, dt*a21, density_k1,  0, 0, 1, 0);
//        amrex::MultiFab::LinComb(momentum_temp, 1.0, momentum_old, 0, dt*a21, momentum_k1, 0, 0, 2, 0);
//        amrex::MultiFab::LinComb(energy_temp,   1.0, energy_old,   0, dt*a21, energy_k1,   0, 0, 1, 0);
//        // fill boundary
//        density_bc ->FillBoundary(density_temp,  0, 1, time, 0); density_temp.FillBoundary(true);
//        momentum_bc->FillBoundary(momentum_temp, 0, 2, time, 0); momentum_temp.FillBoundary(true);
//        energy_bc  ->FillBoundary(energy_temp,   0, 1, time, 0); energy_temp.FillBoundary(true);
//        // k2 = RHS(t + c2*dt, ytemp)
//        RHS(lev,time + c2*dt,
//            density_k2, momentum_k2, energy_k2,
//            density_temp, momentum_temp, energy_temp);
//
//        //
//        // Calculate K3
//        //
//        // ytemp = yold + dt*( a31*k1 + a32*k2 )
//        //
//        // 1. ytemp = yold + dt*a31*k1
//        amrex::MultiFab::LinComb(density_temp,  1.0, density_old,  0, dt*a31, density_k1,  0, 0, 1, 0);
//        amrex::MultiFab::LinComb(momentum_temp, 1.0, momentum_old, 0, dt*a31, momentum_k1, 0, 0, 2, 0);
//        amrex::MultiFab::LinComb(energy_temp,   1.0, energy_old,   0, dt*a31, energy_k1,   0, 0, 1, 0);
//        // 2. ytemp += dt*a32*k2
//        amrex::MultiFab::Saxpy(density_temp,  dt*a32, density_k2,  0, 0, 1, 0);
//        amrex::MultiFab::Saxpy(momentum_temp, dt*a32, momentum_k2, 0, 0, 2, 0);
//        amrex::MultiFab::Saxpy(energy_temp,   dt*a32, energy_k2,   0, 0, 1, 0);
//        // 3. fill boundary
//        density_bc ->FillBoundary(density_temp,  0, 1, time+c2*dt, 0); density_temp.FillBoundary(true);
//        momentum_bc->FillBoundary(momentum_temp, 0, 2, time+c2*dt, 0); momentum_temp.FillBoundary(true);
//        energy_bc  ->FillBoundary(energy_temp,   0, 1, time+c2*dt, 0); energy_temp.FillBoundary(true);
//        // 4. k3 = RHS(t + c3*dt, ytemp)
//        RHS(lev,time + c3*dt,
//            density_k3, momentum_k3, energy_k3,
//            density_temp, momentum_temp, energy_temp);
//        
//        //
//        // Assemble to get ynew
//        //
//        // ynew = yold + dt*(b1*k1 + b2*k2 + b3*k3)
//        //
//        // 1. ynew = yold + dt*b1*k1
//        amrex::MultiFab::LinComb(density_new,  1.0, density_old,  0, dt*b1, density_k1,  0, 0, 1, 0);
//        amrex::MultiFab::LinComb(momentum_new, 1.0, momentum_old, 0, dt*b1, momentum_k1, 0, 0, 2, 0);
//        amrex::MultiFab::LinComb(energy_new,   1.0, energy_old,   0, dt*b1, energy_k1,   0, 0, 1, 0);
//        // 2. ynew += dt*b2*k2
//        amrex::MultiFab::Saxpy(density_new,  dt*b2, density_k2,  0, 0, 1, 0);
//        amrex::MultiFab::Saxpy(momentum_new, dt*b2, momentum_k2, 0, 0, 2, 0);
//        amrex::MultiFab::Saxpy(energy_new,   dt*b2, energy_k2,   0, 0, 1, 0);
//        // 2. ynew += dt*b3*k3
//        amrex::MultiFab::Saxpy(density_new,  dt*b3, density_k3,  0, 0, 1, 0);
//        amrex::MultiFab::Saxpy(momentum_new, dt*b3, momentum_k3, 0, 0, 2, 0);
//        amrex::MultiFab::Saxpy(energy_new,   dt*b3, energy_k3,   0, 0, 1, 0);
//
//    }
//
//
//
//    else if (scheme == IntegrationScheme::RK4)
//    {
//        //
//        // RK4 time integration scheme:
//        //
//        
//        const amrex::BoxArray &ba = density_mf[lev]->boxArray();
//        const amrex::DistributionMapping &dm = density_mf[lev]->DistributionMap();
//        const int ng = density_mf[lev]->nGrow();
//
//        // handles to old solution
//        const amrex::MultiFab &density_old = *density_old_mf[lev];
//        const amrex::MultiFab &momentum_old = *momentum_old_mf[lev];
//        const amrex::MultiFab &energy_old = *energy_old_mf[lev];
//
//        // runge kutta stages
//        amrex::MultiFab density_k1(ba,dm,1,0), momentum_k1(ba,dm,2,0), energy_k1(ba,dm,1,0);
//        amrex::MultiFab density_k2(ba,dm,1,0), momentum_k2(ba,dm,2,0), energy_k2(ba,dm,1,0);
//        amrex::MultiFab density_k3(ba,dm,1,0), momentum_k3(ba,dm,2,0), energy_k3(ba,dm,1,0);
//        amrex::MultiFab density_k4(ba,dm,1,0), momentum_k4(ba,dm,2,0), energy_k4(ba,dm,1,0);
//        
//        // temporary storage
//        amrex::MultiFab density_st(ba,dm,1,ng), momentum_st(ba,dm,2,ng), energy_st(ba,dm,1,ng);
//            
//        // fill the ghost cells from the _old fields, which were updated from the coarse patch.
//        density_st.ParallelCopyToGhost(*density_old_mf[lev],0,0,1,amrex::IntVect(1),amrex::IntVect(1));
//        momentum_st.ParallelCopyToGhost(*momentum_old_mf[lev],0,0,2,amrex::IntVect(1),amrex::IntVect(1));
//        energy_st.ParallelCopyToGhost(*energy_old_mf[lev],0,0,1,amrex::IntVect(1),amrex::IntVect(1));
//
//
//        // handles to new solution
//        amrex::MultiFab &density_new = *density_mf[lev];
//        amrex::MultiFab &momentum_new = *momentum_mf[lev];
//        amrex::MultiFab &energy_new = *energy_mf[lev];
//
//
//        //
//        // K1
//        // 
//
//        RHS(lev,time,
//            density_k1,momentum_k1,energy_k1,
//            density_old, momentum_old, energy_old);
//
//        //
//        // K2
//        //
//
//        // [state] = [old] + (dt/2)[k1]
//        amrex::MultiFab::LinComb(density_st,  1.0, density_old,  0, dt/2.0, density_k1,  0, 0, 1, 0);
//        amrex::MultiFab::LinComb(momentum_st, 1.0, momentum_old, 0, dt/2.0, momentum_k1, 0, 0, 2, 0);
//        amrex::MultiFab::LinComb(energy_st,   1.0, energy_old,   0, dt/2.0, energy_k1,   0, 0, 1, 0);
//
//        density_bc->FillBoundary(density_st, 0, 1, time, 0);   density_st.FillBoundary(true);
//        momentum_bc->FillBoundary(momentum_st, 0, 2, time, 0); momentum_st.FillBoundary(true);
//        energy_bc->FillBoundary(energy_st,0,1,time,0);         energy_st.FillBoundary(true);
//        
//
//        RHS(lev,time,
//            density_k2, momentum_k2, energy_k2,
//            density_st, momentum_st, energy_st);
//
//        //
//        // K3
//        //
//
//        // [state] = [old] + (dt/2)[k2]
//        amrex::MultiFab::LinComb(density_st,  1.0, density_old,  0, dt/2.0, density_k2,  0, 0, 1, 0);
//        amrex::MultiFab::LinComb(momentum_st, 1.0, momentum_old, 0, dt/2.0, momentum_k2, 0, 0, 2, 0);
//        amrex::MultiFab::LinComb(energy_st,   1.0, energy_old,   0, dt/2.0, energy_k2,   0, 0, 1, 0);
//
//        density_bc->FillBoundary(density_st, 0, 1, time, 0);   density_st.FillBoundary(true);
//        momentum_bc->FillBoundary(momentum_st, 0, 2, time, 0); momentum_st.FillBoundary(true);
//        energy_bc->FillBoundary(energy_st,0,1,time,0);         energy_st.FillBoundary(true);
//
//        RHS(lev,time,
//            density_k3, momentum_k3, energy_k3,
//            density_st, momentum_st, energy_st);
//        
//        //
//        // K4
//        //
//        
//        // [state] = [old] + (dt/2)[k3]
//        amrex::MultiFab::LinComb(density_st,  1.0, density_old,  0, dt, density_k3,  0, 0, 1, 0);
//        amrex::MultiFab::LinComb(momentum_st, 1.0, momentum_old, 0, dt, momentum_k3, 0, 0, 2, 0);
//        amrex::MultiFab::LinComb(energy_st,   1.0, energy_old,   0, dt, energy_k3,   0, 0, 1, 0);
//
//        density_bc-> FillBoundary(density_st,  0, 1, time, 0); density_st.FillBoundary(true);
//        momentum_bc->FillBoundary(momentum_st, 0, 2, time, 0); momentum_st.FillBoundary(true);
//        energy_bc->  FillBoundary(energy_st,   0, 1, time, 0); energy_st.FillBoundary(true);
//
//
//        RHS(lev,time,
//            density_k4, momentum_k4, energy_k4,
//            density_st, momentum_st, energy_st);
//
//        
//        // [new] = [old] + (dt/6)k1
//        amrex::MultiFab::LinComb(density_new,  1.0, density_old,  0, (dt/6.0), density_k1,  0, 0, 1, 0);
//        amrex::MultiFab::LinComb(momentum_new, 1.0, momentum_old, 0, (dt/6.0), momentum_k1, 0, 0, 2, 0);
//        amrex::MultiFab::LinComb(energy_new,   1.0, energy_old,   0, (dt/6.0), energy_k1,   0, 0, 1, 0);
//
//        // [new] += (2 dt/6)k2
//        amrex::MultiFab::Saxpy(density_new,  (dt/3.0), density_k2,  0, 0, 1, 0);
//        amrex::MultiFab::Saxpy(momentum_new, (dt/3.0), momentum_k2, 0, 0, 2, 0);
//        amrex::MultiFab::Saxpy(energy_new,   (dt/3.0), energy_k2,   0, 0, 1, 0);
//
//        // [new] += (2 dt/6)k3
//        amrex::MultiFab::Saxpy(density_new,  (dt/3.0), density_k3,  0, 0, 1, 0);
//        amrex::MultiFab::Saxpy(momentum_new, (dt/3.0), momentum_k3, 0, 0, 2, 0);
//        amrex::MultiFab::Saxpy(energy_new,   (dt/3.0), energy_k3,   0, 0, 1, 0);
//                                 
//        // [new] += (dt/6)k4
//        amrex::MultiFab::Saxpy(density_new,  (dt/6.0), density_k4,  0, 0, 1, 0);
//        amrex::MultiFab::Saxpy(momentum_new, (dt/6.0), momentum_k4, 0, 0, 2, 0);
//        amrex::MultiFab::Saxpy(energy_new,   (dt/6.0), energy_k4,   0, 0, 1, 0);
//
//    }


    for (amrex::MFIter mfi(*(*eta_mf)[lev], false); mfi.isValid(); ++mfi)
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

        // Set::Patch<Set::Scalar> tracer_new = tracer_mf.Patch(lev,mfi);
        // Set::Patch<const Set::Scalar> tracer   = tracer_old_mf.Patch(lev,mfi);

        Set::Scalar *dt_max_handle = &dt_max;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {   
            Set::Scalar eta = invert ? 1.0-eta_patch(i,j,k) : eta_patch(i,j,k);

            if (eta < cutoff)
            {
                for (int n=0; n<nspecies; ++n) {
                    rho_new(i,j,k,n) = rho_solid(i,j,k,n);
                }
                M_new(i,j,k,0)   = M_solid(i,j,k,0);
                M_new(i,j,k,1)   = M_solid(i,j,k,1);
                E_new(i,j,k,0)   = E_solid(i,j,k,0);
            }

            Set::Vector vel(u(i,j,k,0),u(i,j,k,1));
            Set::Matrix gradu        = Numeric::Gradient(u, i, j, k, DX);
            omega(i, j, k) = eta * (gradu(1,0) - gradu(0,1));

            // Set::Vector gradtracer        = Numeric::Gradient(tracer, i, j, k,0, DX);
            // tracer_new(i,j,k) = tracer(i,j,k) - vel.dot(gradtracer)*dt;

            if (dynamictimestep.on)
            {
                *dt_max_handle =                          std::fabs(cfl * DX[0] / (u(i,j,k,0)*eta + small));
                *dt_max_handle = std::min(*dt_max_handle, std::fabs(cfl * DX[1] / (u(i,j,k,1)*eta + small)));
                *dt_max_handle = std::min(*dt_max_handle, std::fabs(cfl_v * DX[0]*DX[0] / (Source(i,j,k,nspecies+1)+small)));
                *dt_max_handle = std::min(*dt_max_handle, std::fabs(cfl_v * DX[1]*DX[1] / (Source(i,j,k,nspecies+2)+small)));
            }
        });
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
                const amrex::MultiFab &E_mf)
{

    for (amrex::MFIter mfi(*(velocity_mf)[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox();
        amrex::Array4<const Set::Scalar> const& eta_patch = (*(*eta_old_mf)[lev]).array(mfi);

        Set::Patch<const Set::Scalar> rho       = rho_mf.array(mfi);  // density
        Set::Patch<const Set::Scalar> M         = M_mf.array(mfi);    // momentum
        Set::Patch<const Set::Scalar> E         = E_mf.array(mfi);    // total energy (internal energy + kinetic energy) per unit volume (E/rho = e + 0.5*v^2)

        Set::Patch<const Set::Scalar> rho_solid = solid.density_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> M_solid   = solid.momentum_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> E_solid   = solid.energy_mf.Patch(lev,mfi);

        Set::Patch<Set::Scalar>       v         = velocity_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       p         = pressure_mf.Patch(lev,mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Scalar eta = invert ? 1.0-eta_patch(i,j,k) : eta_patch(i,j,k);

            Set::Scalar etarho_fluid  = 0.0;
            double rho_sum = 0.0;
            double rho_solid_sum = 0.0;
            for (int n=0; n<nspecies; ++n) {
                rho_sum += rho(i,j,k,n);
                rho_solid_sum += rho_solid(i,j,k,n);
            }
            double cp_mix = 0.0;
            double cv_mix = 0.0;
            for (int n=0; n<nspecies; ++n) {
                cp_mix += rho(i,j,k,n)/rho_sum * species_cp[n];
                cv_mix += rho(i,j,k,n)/rho_sum * (species_cp[n] - 8314.45/species_mw[n]);
            }
            gamma = cp_mix/cv_mix;

            etarho_fluid += rho_sum - (1.-eta) * rho_solid_sum;
            Set::Scalar etaE_fluid    = E(i,j,k)   - (1.-eta) * E_solid(i,j,k);

            Set::Vector etaM_fluid( M(i,j,k,0) - (1.-eta) * M_solid(i,j,k,0),
                                    M(i,j,k,1) - (1.-eta) * M_solid(i,j,k,1) );

            //THESE ARE FLUID VELOCITY AND PRESSURE

            v(i,j,k,0) = etaM_fluid(0) / (etarho_fluid + small);
            v(i,j,k,1) = etaM_fluid(1) / (etarho_fluid + small);
            p(i,j,k) = (etaE_fluid - 0.5 * (etaM_fluid(0)*etaM_fluid(0) + etaM_fluid(1)*etaM_fluid(1)) / (etarho_fluid + small)) * ((gamma - 1.0) / (eta + small))+pref;

            if (eta < small) 
            {
                v(i,j,k,0) *= eta;
                v(i,j,k,1) *= eta;
            }
        });
    }

    const Set::Scalar* DX = geom[lev].CellSize();
    amrex::Box domain = geom[lev].Domain();

    int nghost = 1;
    const amrex::BoxArray &ba = energy_mf[lev]->boxArray();
    const amrex::DistributionMapping &dm = energy_mf[lev]->DistributionMap();
    amrex::MultiFab density_sum_mf(ba,dm,1,nghost); // sum_k[rhoY_k]
    amrex::MultiFab mixed_k_mf(ba,dm,1,nghost);     // mixture averaged thermal conductivity coefficient
    amrex::MultiFab mixed_kT_mf(ba,dm,2,nghost);     // mixture averaged thermal conductivity
    amrex::MultiFab mixed_mu_mf(ba,dm,1,nghost);    // mixture averaged dynamic viscosity
    amrex::MultiFab mixed_H_mf(ba,dm,1,nghost);     // Perfect gas mixture enthalpy, H=cp_mix*T
    amrex::MultiFab DKM_mf(ba,dm,nspecies,nghost);  // Diffusion coefficent for species k into mixture
    amrex::MultiFab MF_mf(ba,dm,nspecies,nghost);   // mass fraction
    amrex::MultiFab rhoHDYx_mf(ba,dm,nspecies,nghost);   // species enthalpy diffusion, rho*H*D*dY/dx
    amrex::MultiFab rhoHDYy_mf(ba,dm,nspecies,nghost);   // species enthalpy diffusion, rho*H*D*dY/dy
    amrex::MultiFab rhoDYx_mf(ba,dm,nspecies,nghost);   // Fickian diffusion, rho*D*dY/dx
    amrex::MultiFab rhoDYy_mf(ba,dm,nspecies,nghost);   // Fickian diffusion, rho*D*dY/dy
    amrex::MultiFab omega(ba,dm, nspecies,nghost);  // mass generation/destruction (reaction rate, kg/m^3/s)

    for (amrex::MFIter mfi(*(*eta_mf)[lev], false); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        const amrex::Box& bx_ghost = mfi.growntilebox();

        // Inputs
        Set::Patch<const Set::Scalar> rho  = rho_mf.array(mfi);
        Set::Patch<const Set::Scalar> E    = E_mf.array(mfi);
        Set::Patch<const Set::Scalar> M    = M_mf.array(mfi);

        // Outputs
        Set::Patch<Set::Scalar> rho_rhs  = rho_rhs_mf.array(mfi);
        Set::Patch<Set::Scalar> M_rhs    = M_rhs_mf.array(mfi);
        Set::Patch<Set::Scalar> E_rhs    = E_rhs_mf.array(mfi);


        // Set::Patch<Set::Scalar>       rho_new = density_mf.Patch(lev,mfi);
        // Set::Patch<Set::Scalar>       E_new   = energy_mf.Patch(lev,mfi);
        // Set::Patch<Set::Scalar>       M_new   = momentum_mf.Patch(lev,mfi);

        Set::Patch<const Set::Scalar> rho_solid = solid.density_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> M_solid   = solid.momentum_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> E_solid   = solid.energy_mf.Patch(lev,mfi);

        Set::Patch<Set::Scalar>       omega     = vorticity_mf.Patch(lev,mfi);
        //Set::Patch<Set::Scalar>       pressure  = pressure_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>    temperature  = temperature_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>        rho_sum  = density_sum_mf.array(mfi);
        Set::Patch<Set::Scalar>        mixed_k  = mixed_k_mf.array(mfi);
        Set::Patch<Set::Scalar>       mixed_kT  = mixed_kT_mf.array(mfi);
        Set::Patch<Set::Scalar>       mixed_mu  = mixed_mu_mf.array(mfi);
        Set::Patch<Set::Scalar>        mixed_H  = mixed_H_mf.array(mfi);
        Set::Patch<Set::Scalar>            DKM  = DKM_mf.array(mfi);
        Set::Patch<Set::Scalar>             MF  = MF_mf.array(mfi);
        Set::Patch<Set::Scalar>        rhoHDYx  = rhoHDYx_mf.array(mfi);
        Set::Patch<Set::Scalar>        rhoHDYy  = rhoHDYy_mf.array(mfi);
        Set::Patch<Set::Scalar>         rhoDYx  = rhoDYx_mf.array(mfi);
        Set::Patch<Set::Scalar>         rhoDYy  = rhoDYy_mf.array(mfi);

        Set::Patch<const Set::Scalar> eta_patch = eta_old_mf->Patch(lev,mfi);
        Set::Patch<const Set::Scalar> etadot    = etadot_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> velocity  = velocity_mf.Patch(lev,mfi);

        Set::Patch<const Set::Scalar> m0        = m0_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> q         = q_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> _u0       = u0_mf.Patch(lev,mfi);

        amrex::Array4<Set::Scalar> const& Source = (*Source_mf[lev]).array(mfi);

        // Need to compute temperature field in order to get gradients for next ParallelFor loop
        amrex::ParallelFor(bx_ghost, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {   
            double rhoijk_sum = 0.0;
            for (int n=0; n<nspecies; ++n) {
                rhoijk_sum += rho(i,j,k,n);
            }
            rho_sum(i,j,k) = rhoijk_sum;
            std::vector<double> species_massf(nspecies);
            std::vector<double> species_molef(nspecies);
            for (int n=0; n<nspecies; ++n) {
                species_massf[n] = rho(i,j,k,n)/rhoijk_sum;
            }
            
            // Normalize given values so sum = 1; Get mole fractions
            double mass_total = std::reduce(species_massf.begin(), species_massf.end());
            double moles = 0.0;
            for (int n=0; n<nspecies; ++n) {
                species_massf[n] /= mass_total;
                moles += species_massf[n] / species_mw[n];
            }
            for (int n=0; n<nspecies; ++n) {
                species_molef[n] = species_massf[n] / species_mw[n] / moles;
            }
        
            double mixed_cp = 0.0;
            double mixed_mw = 0.0;
            for (int n=0; n<nspecies; ++n)
            {
                mixed_cp += species_massf[n] * species_cp[n];
                mixed_mw += species_molef[n] * species_mw[n];
            }
            double mixed_cv = mixed_cp - 8314.45/mixed_mw;
            gamma = mixed_cp/mixed_cv;
            Set::Vector u            = Set::Vector(velocity(i, j, k, 0), velocity(i, j, k, 1));
            double internal_energy = (E(i,j,k) - 0.5 * rhoijk_sum * (pow(u(0), 2.0) + pow(u(1), 2.0))) / rhoijk_sum;
            double pressure = (gamma - 1.0) * rho_sum(i,j,k) * internal_energy + pref;
            temperature(i,j,k) = internal_energy / mixed_cv + tref;
            mixed_H(i,j,k) = mixed_cp*temperature(i,j,k);


            // Multispecies effects
            // Mixture viscosity and heat conduction coefficient
            // Multicomponent extension of Chapman-Enskog
            // Transport Phenomena, Revised 2nd Edition / Bird, Stewart, Lightfoot // Ch. 1.4 & 9.3
            //
            // Binary diffusion constants
            // Chapman-Enskog
            // Transport PHenomena, Revised 2nd Edition / Bird, Stewart, Lightfoot // Ch. 17.3
            //
            // Upper triangular matrix for species row and species column since D_AB = D_BA, i.e.
            //
            //                species A | species B | species C
            //              ------------------------------------
            //    species A |    D_AA   |   D_AB    |   D_AC   |
            //              ------------------------------------
            //    species B |           |   D_BB    |   D_BC   |
            //              ------------------------------------
            //    species C |           |           |   D_CC   |
            //              ------------------------------------

            double sigmaAB = 0.0;
            double epsAB = 0.0;
            double nondimT = 0.0;
            double omegaAB = 0.0;
            double DAB = 0.0;
            mixed_mu(i,j,k) = 0.0;
            mixed_k(i,j,k) = 0.0;
            double phi = 0.0;
            for (int a=0; a<nspecies; ++a)
            {
                phi = 0.0;
                DKM(i,j,k,a) = 0.0;
                for (int b=0; b<nspecies; ++b)
                {
                    phi += species_molef[b] * 1.0/sqrt(8.0) *
                           pow(1.0 + species_mw[a]/species_mw[b], -0.5) *
                           pow(1.0 + sqrt(species_mu[a]/species_mu[b]) *
                           pow(species_mw[b]/species_mw[a], 0.25), 2.0);
                    if ( b != a ) 
                    {
                        sigmaAB = 0.5*(species_LJdiameter[a] + species_LJdiameter[b]);
                        epsAB = sqrt(species_LJwelldepth[a] * species_LJwelldepth[b]);
                        nondimT = temperature(i,j,k)/epsAB;
                        omegaAB = collision_integral(nondimT);
                        DAB = 0.0018583*sqrt(pow(temperature(i,j,k), 3.0) * (1.0/species_mw[a] + 1.0/species_mw[b])) / 
                                    (pressure/101325.0*pow(sigmaAB,2.0)*omegaAB);
                        DAB /= 10000.0; // convert from cm^2/s to m^2/s
                        DKM(i,j,k,a) += species_molef[b]/DAB;
                    }
                }
                DKM(i,j,k,a) = (1.0 - species_molef[a])/DKM(i,j,k,a);
                MF(i,j,k,a) = species_massf[a];
                if ( DKM(i,j,k,a) != DKM(i,j,k,a) ) {DKM(i,j,k,a) = 0.0;} // Set to zero if nan or inf (pure species locally)

                mixed_mu(i,j,k) += species_molef[a] * species_mu[a] / phi;
                mixed_k(i,j,k)  += species_molef[a] * species_k[a] / phi;
            }
            // Calculate reaction rates using modified arrhenius (3-body optional)
            // Per Cantera documentation and source code, cantera.org
            //
            // Arrhenius: aA + bB <-> cC + dD
            // R = kf * ( [A]^a * [B]^b - 1/Kc * [C]^c * [D]^d )
            // kf = a*T^b*exp(-E/R0/T)
            // Kc = Kp * (Pref/Ru/T)^(dnu), Kp = exp(-G0/Ru/T), G0 = Sum(nu_j * (H0_j - T*S0_j))
            //
            // Third-Body: A + B + M <-> C + D + M
            // R = kf * [M] * ([A]*[B] - 1/Kc * [C]*[D])
            // [M] = Sum_j( e_j*[C_j] )
            //
            // Get mass net generation rate from reaction rate
            // omega_k = MW_k * Sum_i( (nu''_k - nu'_k) * R_i )
            for (int a=0; a<nspecies; ++a)
            {
                omega(i,j,k,a) = 0.0;
                for (int b=0; b<nreactions; ++b)
                {
                    double nu_diff = rxn[b].products[a] - rxn[b].reactants[a];
                    double kf = rxn[b].factor * pow( temperature(i,j,k), rxn[b].exponent ) * exp( rxn[b].energy/Ru/temperature(i,j,k) );
                    double C_react = 1.0;
                    double C_prod = 1.0;
                    double third_body = 0.0;
                    double Kc = 1.0;
                    for (int c=0; c<nspecies; ++c)
                    {
                        if (rxn[b].type == "arrhenius")
                        // Elementary reaction with modified arrhenius rate
                        {
                            C_react *= pow(species_molef[c]*pressure/Ru/temperature(i,j,k), rxn[b].reactants[c]);
                            if (rxn[b].reverse) {
                                C_prod *= pow(species_molef[c]*pressure/Ru/temperature(i,j,k), rxn[b].products[c]);
                            }
                            third_body = 1.0; // third_body of unity removes impact of third_body which is nonexistent in this case
                        }
                        else if (rxn[b].type == "third-body")
                        {
                            // Need to figure out how to determine if species[c] is in list of rxn.reactants
                            if (species[c] in rxn[b].reactants) {
                                C_react *= species_molef[c]*pressure/Ru/temperature(i,j,k);
                            }
                            // Need to figure out how to determine if species[c] is in list of rxn.products
                            if ((rxn[b].reverse) and (species[c] in rxn[b].products)) {
                                C_prod *= species_molef[c]*pressure/Ru/temperature(i,j,k);
                            }
                            third_body += rxn[b].efficiency[c] * species_molef[c]*pressure/Ru/temperature(i,j,k);
                        }
                        else
                        {
                            Util::ParallelMessage(INFO,"Reaction type is invalid or unspecified.");
                            Util::Abort(INFO);
                        }
                    }
                    if (rxn[b].reverse) {
                        // TODO: calculate Kc for reverse rate
                        Util::ParallelMessage(INFO,"Reverse reactions not yet fully implemented.");
                        Util::Abort(INFO);
                    }
                    omega(i,j,k,a) += nu_diff * kf * third_body * (C_react - rxn[b].reverse*C_prod/Kc);
                }
                omega(i,j,k,a) *= species_mw(a); // Convert molar concentration to mass concentration (i.e. density)
            }
        });
        amrex::ParallelFor(bx_ghost, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {   
            for (int n=0; n<nspecies; ++n) 
            {
                Set::Vector grad_MF = Numeric::Gradient(MF,i,j,k,n,DX);
                if (i == bx_ghost.smallEnd(0)) {
                    grad_MF(0) = (MF(i+1,j,k,n) - MF(i,j,k,n))/DX[0];
                } else if (i == bx_ghost.bigEnd(0)) {
                    grad_MF(0) = (MF(i,j,k,n) - MF(i-1,j,k,n))/DX[0];
                }
                if (j == bx_ghost.smallEnd(1)) {
                    grad_MF(1) = (MF(i,j+1,k,n) - MF(i,j,k,n))/DX[1];
                } else if (j == bx_ghost.bigEnd(1)) {
                    grad_MF(1) = (MF(i,j,k,n) - MF(i,j-1,k,n))/DX[1];
                }
                rhoHDYx(i,j,k,n)    = rho_sum(i,j,k)*mixed_H(i,j,k)*DKM(i,j,k,n)*grad_MF(0);
                rhoHDYy(i,j,k,n)    = rho_sum(i,j,k)*mixed_H(i,j,k)*DKM(i,j,k,n)*grad_MF(1);
                rhoDYx(i,j,k,n)     = rho_sum(i,j,k)*DKM(i,j,k,n)*grad_MF(0);
                rhoDYy(i,j,k,n)     = rho_sum(i,j,k)*DKM(i,j,k,n)*grad_MF(1);
            }
            Set::Vector gradT       = Numeric::Gradient(temperature,i,j,k,0,DX);
            if (i == bx_ghost.smallEnd(0)) {
                gradT(0) = (temperature(i+1,j,k) - temperature(i,j,k))/DX[0];
            } else if (i == bx_ghost.bigEnd(0)) {
                gradT(0) = (temperature(i,j,k) - temperature(i-1,j,k))/DX[0];
            }
            if (j == bx_ghost.smallEnd(1)) {
                gradT(1) = (temperature(i,j+1,k) - temperature(i,j,k))/DX[1];
            } else if (j == bx_ghost.bigEnd(1)) {
                gradT(1) = (temperature(i,j,k) - temperature(i,j-1,k))/DX[1];
            }
            mixed_kT(i,j,k,0)       = mixed_k(i,j,k)*gradT(0);
            mixed_kT(i,j,k,1)       = mixed_k(i,j,k)*gradT(1);
            //Util::ParallelMessage(INFO,"\nijk: ", i, " ", j, " ", k, " DX: ", DX[0]);
            //Util::ParallelMessage(INFO,"bx_ghost: ", bx_ghost, " ", bx_ghost.smallEnd(0));
            //Util::ParallelMessage(INFO,"mixed_kT: ", mixed_kT(i,j,k,0), " ", mixed_kT(i,j,k,1), " mixed_k: ", mixed_k(i,j,k), " gradT: ", gradT(0), " ", gradT(1));
            //Util::ParallelMessage(INFO,"temperature [i-1,i,i+1,j-1,j+1]: ", temperature(i-1,j,k), ",", temperature(i,j,k), ",", temperature(i+1,j,k), ",", temperature(i,j-1,k), ",", temperature(i,j+1,k));
            //for (int n=0; n<nspecies; ++n) {
            //    Util::ParallelMessage(INFO,"rhoHDYx[", n, "]: ", rhoHDYx(i,j,k,n));
            //    Util::ParallelMessage(INFO,"rhoDYx[", n, "]: ", rhoDYx(i,j,k,n));
            //    Util::ParallelMessage(INFO,"rhoHDYy[", n, "]: ", rhoHDYy(i,j,k,n));
            //    Util::ParallelMessage(INFO,"rhoDYy[", n, "]: ", rhoDYy(i,j,k,n));
            //}
        });

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {   
            auto sten = Numeric::GetStencil(i, j, k, domain);

            Set::Scalar eta = invert ? 1.0-eta_patch(i,j,k) : eta_patch(i,j,k);

            //Diffuse Sources
            Set::Vector grad_eta     = Numeric::Gradient(eta_patch, i, j, k, 0, DX);
            Set::Scalar grad_eta_mag = grad_eta.lpNorm<2>();
            Set::Matrix hess_eta     = Numeric::Hessian(eta_patch, i, j, k, 0, DX);
            if (invert) grad_eta *= -1.0;
            if (invert) hess_eta *= -1.0;

            Set::Vector u            = Set::Vector(velocity(i, j, k, 0), velocity(i, j, k, 1));
            Set::Vector u0           = Set::Vector(_u0(i, j, k, 0), _u0(i, j, k, 1));

            Set::Matrix gradM        = Numeric::Gradient(M, i, j, k, DX);
            Set::Vector gradrho      = Numeric::Gradient(rho_sum,i,j,k,0,DX);
            Set::Matrix hess_rho     = Numeric::Hessian(rho_sum,i,j,k,0,DX);
            Set::Matrix gradu        = (gradM - u*gradrho.transpose()) / rho_sum(i,j,k);

            Set::Vector grad_mixed_kTx  = Numeric::Gradient(mixed_kT,i,j,k,0,DX);
            Set::Vector grad_mixed_kTy  = Numeric::Gradient(mixed_kT,i,j,k,1,DX);

            // These need to be recalculated for each species
            //Set::Vector grad_rhoHDYx    = Numeric::Gradient(rhoHDYx,i,j,k,0,DX);
            //Set::Vector grad_rhoHDYy    = Numeric::Gradient(rhoHDYy,i,j,k,0,DX);
            //Set::Vector grad_rhoDYx     = Numeric::Gradient(rhoDYx,i,j,k,0,DX);
            //Set::Vector grad_rhoDYy     = Numeric::Gradient(rhoDYy,i,j,k,0,DX);
            // End

            Set::Vector q0           = Set::Vector(q(i,j,k,0),q(i,j,k,1));

            std::vector<double> species_massf(nspecies);
            std::vector<double> species_molef(nspecies);
            for (int n=0; n<nspecies; ++n) {
                species_massf[n] = rho(i,j,k,n)/rho_sum(i,j,k);
            }
            // Get mole fractions
            double moles = 0.0;
            for (int n=0; n<nspecies; ++n) {
                moles += species_massf[n] / species_mw[n];
            }
            for (int n=0; n<nspecies; ++n) {
                species_molef[n] = species_massf[n] / species_mw[n] / moles;
            }
        
            double mixed_cp = 0.0;
            double mixed_mw = 0.0;
            for (int a=0; a<nspecies; ++a)
            {
                mixed_cp += species_massf[a] * species_cp[a];
                mixed_mw += species_molef[a] * species_mw[a];
            }
            double mixed_cv = mixed_cp - 8314.45/mixed_mw;
            gamma = mixed_cp/mixed_cv;
            double internal_energy = (E(i,j,k) - 0.5 * rho_sum(i,j,k) * (pow(u(0), 2.0) + pow(u(1), 2.0))) / rho_sum(i,j,k);
            double pressure = (gamma - 1.0) * rho_sum(i,j,k) * internal_energy + pref;
            temperature(i,j,k) = internal_energy / mixed_cv + tref;


            if (prescribedflowmode == PrescribedFlowMode::Relative)
            {
                Set::Vector N = grad_eta / (grad_eta_mag + small);
                Set::Vector T(N(1), -N(0));
                u0 = N * u0(0) + T * u0(1);
            }


            std::vector<double> mdot0 = std::vector<double>(nspecies);
            for (int n=0; n<nspecies; ++n) {
                mdot0[n] = -lagrange_m0 * (rho(i,j,k,n)-m0(i,j,k,n)) * grad_eta_mag;
            }
            Set::Vector Pdot0 = Set::Vector::Zero(); 
            Set::Scalar qdot0 = q0.dot(grad_eta);

            Set::Matrix3 hess_M = Numeric::Hessian(M,i,j,k,DX);
            Set::Matrix3 hess_u = Set::Matrix3::Zero();
            for (int p = 0; p < 2; p++)
                for (int q = 0; q < 2; q++)
                    for (int r = 0; r < 2; r++)
                    {
                        hess_u(r,p,q) =
                            (hess_M(r,p,q) - gradu(r,q)*gradrho(p) - gradu(r,p)*gradrho(q) - u(r)*hess_rho(p,q))
                            / rho_sum(i,j,k);
                    }



            Set::Vector Ldot0 = Set::Vector::Zero();
            Set::Vector div_tau = Set::Vector::Zero();
            for (int p = 0; p<2; p++)
                for (int q = 0; q<2; q++)
                    for (int r = 0; r<2; r++)
                        for (int s = 0; s<2; s++)
                        {
                            Set::Scalar Mpqrs = 0.0;
                            if (p==r && q==s) Mpqrs += 0.5 * mixed_mu(i,j,k);

                            Ldot0(p) += 0.5*Mpqrs * (u(r) - u0(r)) * hess_eta(q, s);
                            div_tau(p) += 2.0*Mpqrs * hess_u(r,s,q);
                        }
           
            for (int n=0; n<nspecies; ++n) { 
                Source(i,j, k, n) = mdot0[n];
            }
            Source(i,j, k, nspecies+0) = Pdot0(0) - Ldot0(0);
            Source(i,j, k, nspecies+1) = Pdot0(1) - Ldot0(1);
            Source(i,j, k, nspecies+2) = qdot0;// - Ldot0(0)*v(i,j,k,0) - Ldot0(1)*v(i,j,k,1);

            // Lagrange terms to enforce no-penetration
            Source(i,j,k,nspecies+0) -= lagrange*(u-u0).dot(grad_eta)*grad_eta(0);
            Source(i,j,k,nspecies+1) -= lagrange*(u-u0).dot(grad_eta)*grad_eta(1);

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
                flux_xlo = riemannsolver->Solve(state_xlo_fluid, state_x_fluid, pref, small, species_mw, species_cp) * eta;
                flux_ylo = riemannsolver->Solve(state_ylo_fluid, state_y_fluid, pref, small, species_mw, species_cp) * eta;

                //hi interface fluxes
                flux_xhi = riemannsolver->Solve(state_x_fluid, state_xhi_fluid, pref, small, species_mw, species_cp) * eta;
                flux_yhi = riemannsolver->Solve(state_y_fluid, state_yhi_fluid, pref, small, species_mw, species_cp) * eta;
            }
            catch(...)
            {
                Util::ParallelMessage(INFO,"lev=",lev);
                Util::ParallelMessage(INFO,"i=",i,"j=",j);
                Util::Abort(INFO);
            }

            for (int n=0; n<nspecies; ++n) {
                Set::Scalar drhof_dt = 
                    (flux_xlo.mass[n] - flux_xhi.mass[n]) / DX[0] +
                    (flux_ylo.mass[n] - flux_yhi.mass[n]) / DX[1] +
                    Source(i, j, k, n);
                if (nspecies > 1)
                {
                    // species diffusion term, d/dx_i(rho*DKM*Y,i)
                    Set::Vector grad_rhoDYx     = Numeric::Gradient(rhoDYx,i,j,k,n,DX);
                    Set::Vector grad_rhoDYy     = Numeric::Gradient(rhoDYy,i,j,k,n,DX);
                    drhof_dt += eta * (grad_rhoDYx[0] + grad_rhoDYy[1]);
                }

                rho_rhs(i,j,k,n) = 
                    // rho_new(i, j, k) = rho(i, j, k) + 
                    //(
                        drhof_dt +
                        // todo add drhos_dt term if want time-evolving rhos
                        etadot(i,j,k) * (rho(i,j,k,n) - rho_solid(i,j,k,n)) / (eta + small)
                    // ) * dt;
                    ;
            }
                    
            Set::Scalar dMxf_dt =
                (flux_xlo.momentum_normal  - flux_xhi.momentum_normal ) / DX[0] +
                (flux_ylo.momentum_tangent - flux_yhi.momentum_tangent) / DX[1] +
                eta * (div_tau(0) + g(0)*rho_sum(i,j,k)) +
                Source(i, j, k, nspecies+0);

            M_rhs(i,j,k,0) = 
                //M_new(i, j, k, 0) = M(i, j, k, 0) +
                // ( 
                    dMxf_dt + 
                    // todo add dMs_dt term if want time-evolving Ms
                    etadot(i,j,k)*(M(i,j,k,0) - M_solid(i,j,k,0)) / (eta + small)
                // ) * dt;
                ;

            Set::Scalar dMyf_dt =
                (flux_xlo.momentum_tangent - flux_xhi.momentum_tangent) / DX[0] +
                (flux_ylo.momentum_normal  - flux_yhi.momentum_normal ) / DX[1] +
                eta * (div_tau(1) + g(1)*rho_sum(i,j,k)) +
                Source(i, j, k, nspecies+1);
                
            M_rhs(i,j,k,1) = 
                //M_new(i, j, k, 1) = M(i, j, k, 1) +
                //( 
                    dMyf_dt +
                    // todo add dMs_dt term if want time-evolving Ms
                    etadot(i,j,k)*(M(i,j,k,1) - M_solid(i,j,k,1)) / (eta+small)
                // )*dt;
                ;

            Set::Scalar dEf_dt =
                (flux_xlo.energy - flux_xhi.energy) / DX[0] +
                (flux_ylo.energy - flux_yhi.energy) / DX[1] +
                eta * (div_tau.dot(u) + (grad_mixed_kTx[0] + grad_mixed_kTy[1]) + rho_sum(i,j,k)*g.dot(u)) +
                Source(i, j, k, nspecies+2);

            if (nspecies > 1)
            {
                for (int n=0; n<nspecies; ++n) {
                    // Species energy diffusion term: d/dx_i(rho*H*DKM*Y,i)
                    Set::Vector grad_rhoHDYx     = Numeric::Gradient(rhoHDYx,i,j,k,n,DX);
                    Set::Vector grad_rhoHDYy     = Numeric::Gradient(rhoHDYy,i,j,k,n,DX);
                    dEf_dt += eta * (grad_rhoHDYx[0] + grad_rhoHDYy[1]);
                              
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
            if ((rho_rhs(i,j,k,0) != rho_rhs(i,j,k,0)) ||
                (M_rhs(i,j,k,0) != M_rhs(i,j,k,0)) ||
                (M_rhs(i,j,k,1) != M_rhs(i,j,k,1)) ||
                (E_rhs(i,j,k) != E_rhs(i,j,k)))
            {
                Util::ParallelMessage(INFO,"i,j,k=",i, ",", j, ",", k);
                Util::ParallelMessage(INFO,"rho_rhs=",rho_rhs(i,j,k,0));
                Util::ParallelMessage(INFO,"Mx_rhs=",M_rhs(i,j,k,0));
                Util::ParallelMessage(INFO,"Mx_rhs=",M_rhs(i,j,k,1));
                Util::ParallelMessage(INFO,"E_rhs=",E_rhs(i,j,k));

                Util::ParallelMessage(INFO,"lev=",lev);
                Util::ParallelMessage(INFO,"i=",i," j=",j);
                //Util::ParallelMessage(INFO,"drhof_dt ",drhof_dt); // dies
                Util::ParallelMessage(INFO,"flux_xlo.mass ",flux_xlo.mass[0]);
                Util::ParallelMessage(INFO,"flux_xhi.mass ",flux_xhi.mass[0]); // dies, depends on state_xx, state_xhi, state_x_solid, state_xhi_solid, gamma, eta, pref, small
                Util::ParallelMessage(INFO,"flux_ylo.mass ",flux_ylo.mass[0]);
                Util::ParallelMessage(INFO,"flux_xhi.mass ",flux_yhi.mass[0]);
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
                Util::Message(INFO,Source(i, j, k, nspecies+2));
                
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


#endif
