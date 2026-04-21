
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

#include "Model/Chemistry/Chemistry.H"
#include "Model/Chemistry/Frozen.H"
#include "Model/Chemistry/Equilibrium.H"

namespace Integrator
{

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
        std::cout << value.gas.thermo->model_name() << "\n";
        std::cout << value.gas.transport->model_name() << "\n";
        std::cout << value.gas.eos->model_name() << "\n";
        std::cout << NSPECIES << "\n";

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
        pp.select_default<BC::Constant,BC::Expression>("density.bc",value.density_bc, NSPECIES);
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
        pp_query_default("details",value.details,false); // save detailed data (viscosity, heat conductivity, etc.)

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

        value.RegisterNewFab(value.density_mf,     value.density_bc, NSPECIES, nghost, "density",     true , true);
        value.RegisterNewFab(value.density_old_mf, value.density_bc, NSPECIES, nghost, "density_old", false, true);

        value.RegisterNewFab(value.energy_mf,     value.energy_bc, 1, nghost, "energy",      true ,true);
        value.RegisterNewFab(value.energy_old_mf, value.energy_bc, 1, nghost, "energy_old" , false, true);

        value.RegisterNewFab(value.momentum_mf,     value.momentum_bc, 2, nghost, "momentum",     true ,true, {"x","y"});
        value.RegisterNewFab(value.momentum_old_mf, value.momentum_bc, 2, nghost, "momentum_old", false, true);
 
        value.RegisterNewFab(value.pressure_mf,  &value.bc_nothing, 1, nghost, "pressure",  true, false);
        value.RegisterNewFab(value.temperature_mf,  &value.bc_nothing, 1, nghost, "temperature",  true, false);
        value.RegisterNewFab(value.velocity_mf,  &value.bc_nothing, 2, nghost, "velocity",  true, false,{"x","y"});
        value.RegisterNewFab(value.vorticity_mf, &value.bc_nothing, 1, nghost, "vorticity", true, false);

        value.RegisterNewFab(value.m0_mf,           &value.bc_nothing, NSPECIES, 0, "m0",  true, false);
        value.RegisterNewFab(value.u0_mf,           &value.bc_nothing, 2, 0, "u0",  true, false, {"x","y"});
        value.RegisterNewFab(value.q_mf,            &value.bc_nothing, 2, 0, "q",   true, false, {"x","y"});

        pp.select_default<BC::Constant,BC::Expression>("bc_D",value.neumann_bc_D, AMREX_SPACEDIM);
        pp.select_default<BC::Constant,BC::Expression>("bc_1",value.neumann_bc_1, 1);
        pp.select_default<BC::Constant,BC::Expression>("bc_N",value.neumann_bc_N, NSPECIES);

        value.RegisterNewFab(value.solid.density_mf,  value.neumann_bc_N, NSPECIES, nghost, "solid.density", true, false);
        value.RegisterNewFab(value.solid.momentum_mf, value.neumann_bc_D, 2, nghost, "solid.momentum", true, false, {"x","y"});
        value.RegisterNewFab(value.solid.energy_mf,   value.neumann_bc_1, 1, nghost, "solid.energy",   true, false);

        value.RegisterNewFab(value.Source_mf, &value.bc_nothing, NSPECIES+3, 0, "Source", true, false);

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
    pp.select_default<  Solver::Local::Riemann::Roe,
                        Solver::Local::Riemann::HLLE,
                        Solver::Local::Riemann::HLLC>("solver",value.riemannsolver);

    // Chemistry Model
    pp.select<  Model::Chemistry::Frozen,
                Model::Chemistry::Equilibrium,
                Model::Chemistry::FiniteRateKinetics,
                Model::Chemistry::Rocfire>("chemistry",value.chemistry);

    std::string prescribedflowmode_str;
    // 
    pp.query_validate("prescribedflowmode",prescribedflowmode_str,{"absolute","relative"});
    if (prescribedflowmode_str == "absolute") value.prescribedflowmode = PrescribedFlowMode::Absolute;
    else if (prescribedflowmode_str == "relative") value.prescribedflowmode = PrescribedFlowMode::Relative;

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

    ic_m0            ->Initialize(lev, m0_mf, 0.0);
    ic_u0            ->Initialize(lev, u0_mf, 0.0);
    ic_q             ->Initialize(lev, q_mf,  0.0);

    Source_mf[lev]   ->setVal(0.0);

    if (managed)  { if (lev >= (int)mixed.size()) mixed.push_back(false);}
    else  Mix(lev);
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
            T(i,j,k) = gas.ComputeT_from_primitives(p(i,j,k), density, X, i, j, k);
            Set::Scalar E_fluid = gas.ComputeE(density, density*v(i,j,k,0), density*v(i,j,k,1), T(i,j,k), X, i, j, k);

            // Mix
            M(i, j, k, 0) = density*v(i, j, k, 0)*eta +  M_solid(i, j, k, 0)*(1.0-eta);
            M(i, j, k, 1) = density*v(i, j, k, 1)*eta +  M_solid(i, j, k, 1)*(1.0-eta);
            M_old(i, j, k, 0) = M(i, j, k, 0);
            M_old(i, j, k, 1) = M(i, j, k, 1);

            for (int n=0; n<NSPECIES; ++n)
            {
                rho(i, j, k, n) = eta * rho(i, j, k, n) + (1.0 - eta) * rho_solid(i, j, k, n);
                rho_old(i, j, k, n) = rho(i, j, k, n);
            }

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


    //
    // DO TIME INTEGRATION (driving the RHS function)
    //

    // Organize references to the "new" solution
    amrex::Vector<amrex::MultiFab> solution_new; 
    solution_new.emplace_back(*density_mf[lev].get(),amrex::MakeType::make_alias,0,NSPECIES);
    solution_new.emplace_back(*momentum_mf[lev].get(),amrex::MakeType::make_alias,0,2);
    solution_new.emplace_back(*energy_mf[lev].get(),amrex::MakeType::make_alias,0,1);

    // Organize references to the "old" solution
    amrex::Vector<amrex::MultiFab> solution_old;
    solution_old.emplace_back(*density_old_mf[lev].get(),amrex::MakeType::make_alias,0,NSPECIES);
    solution_old.emplace_back(*momentum_old_mf[lev].get(),amrex::MakeType::make_alias,0,2);
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

    // Take care of filling boundaries during stages
    timeintegrator.set_post_stage_action([&](amrex::Vector<amrex::MultiFab> & stage_mf, Set::Scalar time) 
    {
        density_bc->FillBoundary(stage_mf[0],0,NSPECIES,time,0);   
        stage_mf[0].FillBoundary(true);
        momentum_bc->FillBoundary(stage_mf[1],0,2,time,0);  
        stage_mf[1].FillBoundary(true);
        energy_bc->FillBoundary(stage_mf[2],0,1,time,0);    
        stage_mf[2].FillBoundary(true);
    });
    
    // Do the update
    timeintegrator.advance(solution_old, solution_new, time, dt);


    //
    // APPLY CUTOFFS AND DO DYNAMIC TIMESTEP CALCULATION
    //

    Set::Scalar dt_max = std::numeric_limits<Set::Scalar>::max();
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
                for (int n=0; n<NSPECIES; ++n)
                {
                    rho_new(i,j,k,n) = rho_solid(i,j,k,n);
                }
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


void Hydro::RHS(int lev, Set::Scalar /*time*/, Set::Scalar dt,
                amrex::MultiFab &rho_rhs_mf, 
                amrex::MultiFab &M_rhs_mf, 
                amrex::MultiFab &E_rhs_mf,
                const amrex::MultiFab &rho_mf,
                const amrex::MultiFab &M_mf,
                const amrex::MultiFab &E_mf)
{
    int nghost = 1;
    const amrex::BoxArray &ba = energy_mf[lev]->boxArray();
    const amrex::DistributionMapping &dm = energy_mf[lev]->DistributionMap();
    amrex::MultiFab rho_sum_mf(ba,dm,1,nghost);             // sum_k[rhoY_k]
    amrex::MultiFab mixed_k_mf(ba,dm,1,nghost);             // mixture averaged thermal conductivity coefficient
    amrex::MultiFab mixed_kT_mf(ba,dm,2,nghost);            // mixture averaged thermal conductivity
    amrex::MultiFab mixed_mu_mf(ba,dm,1,nghost);            // mixture averaged dynamic viscosity
    amrex::MultiFab mixed_H_mf(ba,dm,1,nghost);             // Perfect gas mixture enthalpy, H=cp_mix*T
    amrex::MultiFab DKM_mf(ba,dm,NSPECIES,nghost);      // Diffusion coefficent for species k into mixture
    amrex::MultiFab rhoHDYx_mf(ba,dm,NSPECIES,nghost);  // species enthalpy diffusion, rho*H*D*dY/dx
    amrex::MultiFab rhoHDYy_mf(ba,dm,NSPECIES,nghost);  // species enthalpy diffusion, rho*H*D*dY/dy
    amrex::MultiFab rhoDYx_mf(ba,dm,NSPECIES,nghost);   // Fickian diffusion, rho*D*dY/dx
    amrex::MultiFab rhoDYy_mf(ba,dm,NSPECIES,nghost);   // Fickian diffusion, rho*D*dY/dy

    // Values only to be written if details=true
    amrex::Array4<Set::Scalar> mu_arr;
    amrex::Array4<Set::Scalar> k_arr;
    amrex::Array4<Set::Scalar> D_arr;
    amrex::Array4<Set::Scalar> wdot_arr;
    amrex::Array4<Set::Scalar> qdot_arr;

    const Set::Scalar* DX = geom[lev].CellSize();
    amrex::Box domain = geom[lev].Domain();

    for (amrex::MFIter mfi(*(velocity_mf)[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox();
        amrex::Array4<const Set::Scalar> const& eta_patch = (*(*eta_old_mf)[lev]).array(mfi);

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
            Set::Scalar eta = invert ? 1.0-eta_patch(i,j,k)*eta_patch(i,j,k) : eta_patch(i,j,k);

            // Compute T and P primitives from mixed values
            gas.ComputeLocalFractions(rho, Y, X, i, j, k);
            Set::Scalar density = gas.ComputeD(rho, i, j, k);
            T(i,j,k) = gas.ComputeT(density, M(i,j,k,0), M(i,j,k,1), E(i,j,k), T(i,j,k), X, i, j, k);
            p(i,j,k) = gas.ComputeP(density, T(i,j,k), X, i, j, k);

            // Compute velocity from fluid values
            for (int n=0; n<NSPECIES; ++n)
            {
                scratch(i,j,k,n) = (rho(i,j,k,n) - rho_solid(i,j,k,n)*(1.0 - eta))/(eta + small);
            }
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

            rho_sum(i,j,k) = density;
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
            for (int n=0; n<NSPECIES; ++n)
            {
                Set::Vector grad_Y;
                if (i == bx.smallEnd(0)) {
                    grad_Y(0) = (Y(i+1,j,k,n) - Y(i,j,k,n))/DX[0];
                } else if (i == bx.bigEnd(0)) {
                    grad_Y(0) = (Y(i,j,k,n) - Y(i-1,j,k,n))/DX[0];
                } else {
                    grad_Y(0) = (Y(i+1,j,k,n) - Y(i-1,j,k,n))/(2.0*DX[0]);
                }
                if (j == bx.smallEnd(1)) {
                    grad_Y(1) = (Y(i,j+1,k,n) - Y(i,j,k,n))/DX[1];
                } else if (j == bx.bigEnd(1)) {
                    grad_Y(1) = (Y(i,j,k,n) - Y(i,j-1,k,n))/DX[1];
                } else {
                    grad_Y(1) = (Y(i,j+1,k,n) - Y(i,j-1,k,n))/(2.0*DX[1]);
                }
                rhoHDYx(i,j,k,n)    = rho_sum(i,j,k)*mixed_H(i,j,k)*DKM(i,j,k,n)*grad_Y(0);
                rhoHDYy(i,j,k,n)    = rho_sum(i,j,k)*mixed_H(i,j,k)*DKM(i,j,k,n)*grad_Y(1);
                rhoDYx(i,j,k,n)     = rho_sum(i,j,k)*DKM(i,j,k,n)*grad_Y(0);
                rhoDYy(i,j,k,n)     = rho_sum(i,j,k)*DKM(i,j,k,n)*grad_Y(1);
            }
            Set::Vector gradT;
            if (i == bx.smallEnd(0)) {
                gradT(0) = (T(i+1,j,k) - T(i,j,k))/DX[0];
            } else if (i == bx.bigEnd(0)) {
                gradT(0) = (T(i,j,k) - T(i-1,j,k))/DX[0];
            } else {
                gradT(0) = (T(i+1,j,k) - T(i-1,j,k))/(2.0*DX[0]);
            }
            if (j == bx.smallEnd(1)) {
                gradT(1) = (T(i,j+1,k) - T(i,j,k))/DX[1];
            } else if (j == bx.bigEnd(1)) {
                gradT(1) = (T(i,j,k) - T(i,j-1,k))/DX[1];
            } else {
                gradT(1) = (T(i,j+1,k) - T(i,j-1,k))/(2.0*DX[1]);
            }
            mixed_kT(i,j,k,0)       = mixed_k(i,j,k)*gradT(0);
            mixed_kT(i,j,k,1)       = mixed_k(i,j,k)*gradT(1);
        });
    }

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

        Set::Patch<Set::Scalar>       omega     = vorticity_mf.Patch(lev,mfi);

        Set::Patch<const Set::Scalar> eta_patch = eta_old_mf->Patch(lev,mfi);
        Set::Patch<const Set::Scalar> etadot    = etadot_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> velocity  = velocity_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> T         = temperature_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> molef     = mole_fraction_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> pressure  = pressure_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> temp      = temperature_mf.Patch(lev,mfi);

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

        // Third and final ParallelFor loop to get 2nd gradients and compute fluxes
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

            Set::Matrix gradM        = Numeric::Gradient(M, i, j, k, DX);
            Set::Vector gradrho      = Numeric::Gradient(rho_sum,i,j,k,0,DX);
            Set::Matrix hess_rho     = Numeric::Hessian(rho_sum,i,j,k,0,DX,sten);
            Set::Matrix gradu        = (gradM - u*gradrho.transpose()) / rho_sum(i,j,k);

            Set::Vector grad_mixed_kTx  = Numeric::Gradient(mixed_kT,i,j,k,0,DX);
            Set::Vector grad_mixed_kTy  = Numeric::Gradient(mixed_kT,i,j,k,1,DX);
            // Gradients of rhoHDY and rhoDY are computed for individual species in for loops later

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


            std::vector<double> mdot0(NSPECIES);
            for (int n=0; n<NSPECIES; ++n ) { mdot0[n] = m0(i,j,k,n)*grad_eta_mag; }
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
                            / rho_sum(i,j,k);
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


            std::array<double, NSPECIES> rhoY;
            for (int n=0; n<NSPECIES; ++n)
            {
                Source(i,j, k, n) = mdot0[n];
                rhoY[n] = rho(i,j,k,n);
            }
            Source(i,j, k, NSPECIES  ) = Pdot0(0) - Ldot0(0);
            Source(i,j, k, NSPECIES+1) = Pdot0(1) - Ldot0(1);
            Source(i,j, k, NSPECIES+2) = qdot0;// - Ldot0(0)*v(i,j,k,0) - Ldot0(1)*v(i,j,k,1);

            std::array<double, NSPECIES> wdot;
            Set::Scalar qdot = 0.0;
            std::tie(wdot, qdot) = chemistry.compute(pressure(i,j,k), temp(i,j,k), rhoY, dt, &gas);

            if (details)
            {
                qdot_arr(i,j,k) = qdot;
                for (int n=0; n<NSPECIES; ++n) wdot_arr(i,j,k,n) = wdot[n];
            }

            // Lagrange terms to enforce no-penetration
            Source(i,j,k,NSPECIES  ) -= lagrange*(u-u0).dot(grad_eta)*grad_eta(0);
            Source(i,j,k,NSPECIES+1) -= lagrange*(u-u0).dot(grad_eta)*grad_eta(1);

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
                
            Set::Scalar drhof_dt = 0.0;
            for (int n=0; n<NSPECIES; ++n)
            {
                drhof_dt = 
                    (flux_xlo.mass[n] - flux_xhi.mass[n]) / DX[0] +
                    (flux_ylo.mass[n] - flux_yhi.mass[n]) / DX[1] +
                    Source(i, j, k, n);
                if (NSPECIES > 1)
                {
                    // species diffusion term, d/dx_i(rho*DKM*Y,i)
                    Set::Vector grad_rhoDYx     = Numeric::Gradient(rhoDYx,i,j,k,n,DX);
                    Set::Vector grad_rhoDYy     = Numeric::Gradient(rhoDYy,i,j,k,n,DX);
                    drhof_dt += eta * (grad_rhoDYx[0] + grad_rhoDYy[1]);
                }
                drhof_dt += eta * wdot[n];

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
                div_tau(0) * eta +
                g(0)*rho_sum(i,j,k) +
                Source(i, j, k, NSPECIES);

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
                div_tau(1) * eta + 
                g(1)*rho_sum(i,j,k) +
                Source(i, j, k, NSPECIES+1);

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
                eta * (div_tau.dot(u) + (grad_mixed_kTx[0] + grad_mixed_kTy[1])) +
                rho_sum(i,j,k)*g.dot(u) +
                Source(i, j, k, NSPECIES+2) +
                eta * qdot;

            if (NSPECIES > 1)
            {
                for (int n=0; n<NSPECIES; ++n)
                {
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
                //Util::ParallelMessage(INFO,"flux_xlo.mass ",flux_xlo.mass);
                //Util::ParallelMessage(INFO,"flux_xhi.mass ",flux_xhi.mass); // dies, depends on state_xx, state_xhi, state_x_solid, state_xhi_solid, eta, small
                //Util::ParallelMessage(INFO,"flux_ylo.mass ",flux_ylo.mass);
                //Util::ParallelMessage(INFO,"flux_xhi.mass ",flux_yhi.mass);
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
