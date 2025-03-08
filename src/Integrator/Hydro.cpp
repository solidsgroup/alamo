
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
#include "Solver/Local/Riemann/Roe.H"

#if AMREX_SPACEDIM == 2

namespace Integrator
{

Hydro::Hydro(IO::ParmParse& pp) : Hydro()
{
    pp.queryclass(*this);
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

        pp_forbid("rho.bc","--> density.bc");
        pp_forbid("p.bc","--> pressure.bc");
        pp_forbid("v.bc","--> velocity.bc");
        value.density_bc = new BC::Expression(1, pp, "density.bc");
        pp_forbid("pressure.bc","--> energy.bc");
        value.energy_bc = new BC::Constant(1, pp, "energy.bc");
        pp_forbid("velocity.bc","--> momentum.bc");
        value.momentum_bc = new BC::Expression(2, pp, "momentum.bc");
        value.eta_bc = new BC::Constant(1, pp, "pf.eta.bc");

        pp_query_default("small",value.small,1E-8); // small regularization value
        pp_query_default("cutoff",value.cutoff,-1E100); // cutoff value
        pp_query_default("lagrange",value.lagrange,0.0); // lagrange no-penetration factor

        pp_forbid("roefix","--> solver.roe.entropy_fix"); // Roe solver entropy fix

    }
    // Register FabFields:
    {
        int nghost = 2;

        value.RegisterNewFab(value.eta_mf,     value.eta_bc, 1, nghost, "eta",     true );
        value.RegisterNewFab(value.eta_old_mf, value.eta_bc, 1, nghost, "eta_old", true);
        value.RegisterNewFab(value.etadot_mf,  value.eta_bc, 1, nghost, "etadot",  true );

        value.RegisterNewFab(value.density_mf,     value.density_bc, 1, nghost, "density",     true );
        value.RegisterNewFab(value.density_old_mf, value.density_bc, 1, nghost, "density_old", false);

        value.RegisterNewFab(value.energy_mf,     value.energy_bc, 1, nghost, "energy",      true );
        value.RegisterNewFab(value.energy_old_mf, value.energy_bc, 1, nghost, "energy_old" , false);

        value.RegisterNewFab(value.momentum_mf,     value.momentum_bc, 2, nghost, "momentum",     true , {"x","y"});
        value.RegisterNewFab(value.momentum_old_mf, value.momentum_bc, 2, nghost, "momentum_old", false);
 
        value.RegisterNewFab(value.pressure_mf,  &value.bc_nothing, 1, nghost, "pressure",  true);
        value.RegisterNewFab(value.velocity_mf,  &value.bc_nothing, 2, nghost, "velocity",  true, {"x","y"});
        value.RegisterNewFab(value.vorticity_mf, &value.bc_nothing, 1, nghost, "vorticity", true);

        value.RegisterNewFab(value.m0_mf,           &value.bc_nothing, 1, 0, "m0",  true);
        value.RegisterNewFab(value.u0_mf,           &value.bc_nothing, 2, 0, "u0",  true, {"x","y"});
        value.RegisterNewFab(value.q_mf,            &value.bc_nothing, 2, 0, "q",   true, {"x","y"});

        value.RegisterNewFab(value.solid.momentum_mf, &value.neumann_bc_D, 2, nghost, "solid.momentum", true, {"x","y"});
        value.RegisterNewFab(value.solid.density_mf,  &value.neumann_bc_1,  1, nghost, "solid.density", true);
        value.RegisterNewFab(value.solid.energy_mf,   &value.neumann_bc_1, 1, nghost, "solid.energy",   true);

        value.RegisterNewFab(value.Source_mf, &value.bc_nothing, 4, 0, "Source", true);
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

    // eta initial condition
    pp.select_default<IC::Constant,IC::Laminate,IC::Expression,IC::BMP,IC::PNG>("eta.ic",value.eta_ic,value.geom);

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
    pp.select_default<Solver::Local::Riemann::Roe>("solver",value.roesolver);
}


void Hydro::Initialize(int lev)
{
    BL_PROFILE("Integrator::Hydro::Initialize");
 
    eta_ic           ->Initialize(lev, eta_mf,     0.0);
    eta_ic           ->Initialize(lev, eta_old_mf, 0.0);
    etadot_mf[lev]   ->setVal(0.0);

    //flux_mf[lev]   ->setVal(0.0);

    velocity_ic      ->Initialize(lev, velocity_mf, 0.0);
    pressure_ic      ->Initialize(lev, pressure_mf, 0.0);
    density_ic       ->Initialize(lev, density_mf, 0.0);

    density_ic       ->Initialize(lev, density_old_mf, 0.0);


    solid.density_ic ->Initialize(lev, solid.density_mf, 0.0);
    solid.momentum_ic->Initialize(lev, solid.momentum_mf, 0.0);
    solid.energy_ic  ->Initialize(lev, solid.energy_mf, 0.0);

    ic_m0            ->Initialize(lev, m0_mf, 0.0);
    ic_u0            ->Initialize(lev, u0_mf,    0.0);
    ic_q             ->Initialize(lev, q_mf,            0.0);

    Source_mf[lev]   ->setVal(0.0);

    Mix(lev);
}

void Hydro::Mix(int lev)
{
    for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox();

        Set::Patch<const Set::Scalar> eta       = eta_mf.Patch(lev,mfi);

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
            rho(i, j, k) = eta(i, j, k) * rho(i, j, k) + (1.0 - eta(i, j, k)) * rho_solid(i, j, k);
            rho_old(i, j, k) = rho(i, j, k);

            M(i, j, k, 0) = (rho(i, j, k)*v(i, j, k, 0))*eta(i, j, k)  +  M_solid(i, j, k, 0)*(1.0-eta(i, j, k));
            M(i, j, k, 1) = (rho(i, j, k)*v(i, j, k, 1))*eta(i, j, k)  +  M_solid(i, j, k, 1)*(1.0-eta(i, j, k));
            M_old(i, j, k, 0) = M(i, j, k, 0);
            M_old(i, j, k, 1) = M(i, j, k, 1);

            E(i, j, k) =
                (0.5 * (v(i, j, k, 0) * v(i, j, k, 0) + v(i, j, k, 1) * v(i, j, k, 1)) * rho(i, j, k) + p(i, j, k) / (gamma - 1.0)) * eta(i, j, k) 
                + 
                E_solid(i, j, k) * (1.0 - eta(i, j, k));
            E_old(i, j, k) = E(i, j, k);
        });
    }
    c_max = 0.0;
    vx_max = 0.0;
    vy_max = 0.0;
}

void Hydro::UpdateEta(int lev, Set::Scalar time)
{
    eta_ic->Initialize(lev, eta_mf, time);
}

void Hydro::TimeStepBegin(Set::Scalar, int /*iter*/)
{

}

void Hydro::TimeStepComplete(Set::Scalar, int lev)
{
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

    std::swap(eta_old_mf, eta_mf);
    std::swap(density_old_mf[lev],  density_mf[lev]);
    std::swap(momentum_old_mf[lev], momentum_mf[lev]);
    std::swap(energy_old_mf[lev],   energy_mf[lev]);
    Set::Scalar dt_max = std::numeric_limits<Set::Scalar>::max();
    
    UpdateEta(lev, time);

    for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox();
        amrex::Array4<const Set::Scalar> const& eta_new = (*eta_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& eta = (*eta_old_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar>       const& etadot = (*etadot_mf[lev]).array(mfi);

        Set::Patch<const Set::Scalar> rho       = density_old_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> M         = momentum_old_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> E         = energy_old_mf.Patch(lev,mfi); // total energy (internal energy + kinetic energy) per unit volume
                                                                                // E/rho = e + 0.5*v^2

        Set::Patch<const Set::Scalar> rho_solid = solid.density_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> M_solid   = solid.momentum_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> E_solid   = solid.energy_mf.Patch(lev,mfi);

        Set::Patch<Set::Scalar>       v         = velocity_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       p         = pressure_mf.Patch(lev,mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            etadot(i, j, k) = (eta_new(i, j, k) - eta(i, j, k)) / dt;

            Set::Scalar etarho_fluid  = rho(i,j,k) - (1.-eta(i,j,k)) * rho_solid(i,j,k);
            Set::Scalar etaE_fluid    = E(i,j,k)   - (1.-eta(i,j,k)) * E_solid(i,j,k);

            Set::Vector etaM_fluid( M(i,j,k,0) - (1.-eta(i,j,k)) * M_solid(i,j,k,0),
                                    M(i,j,k,1) - (1.-eta(i,j,k)) * M_solid(i,j,k,1) );

            //THESE ARE FLUID VELOCITY AND PRESSURE

            v(i,j,k,0) = etaM_fluid(0) / (etarho_fluid + small);
            v(i,j,k,1) = etaM_fluid(1) / (etarho_fluid + small);

            p(i,j,k)   = (etaE_fluid / (eta(i, j, k) + small) - 0.5 * (etaM_fluid(0)*etaM_fluid(0) + etaM_fluid(1)*etaM_fluid(1)) / (etarho_fluid + small)) * ((gamma - 1.0) / (eta(i, j, k) + small))-pref;
        });
    }

    const Set::Scalar* DX = geom[lev].CellSize();
    amrex::Box domain = geom[lev].Domain();

    for (amrex::MFIter mfi(*eta_mf[lev], false); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        
        Set::Patch<const Set::Scalar> rho = density_old_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> E   = energy_old_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> M   = momentum_old_mf.Patch(lev,mfi);

        Set::Patch<Set::Scalar>       rho_new = density_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       E_new   = energy_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       M_new   = momentum_mf.Patch(lev,mfi);

        Set::Patch<const Set::Scalar> rho_solid = solid.density_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> M_solid   = solid.momentum_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> E_solid   = solid.energy_mf.Patch(lev,mfi);

        Set::Patch<Set::Scalar>       omega     = vorticity_mf.Patch(lev,mfi);

        Set::Patch<const Set::Scalar> eta       = eta_old_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> etadot    = etadot_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> velocity  = velocity_mf.Patch(lev,mfi);

        Set::Patch<const Set::Scalar> m0        = m0_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> q         = q_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> _u0       = u0_mf.Patch(lev,mfi);

        amrex::Array4<Set::Scalar> const& Source = (*Source_mf[lev]).array(mfi);

        Set::Scalar *dt_max_handle = &dt_max;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {   
            auto sten = Numeric::GetStencil(i, j, k, domain);

            //Diffuse Sources
            Set::Vector grad_eta     = Numeric::Gradient(eta, i, j, k, 0, DX);
            Set::Scalar grad_eta_mag = grad_eta.lpNorm<2>();
            Set::Matrix hess_eta     = Numeric::Hessian(eta, i, j, k, 0, DX);

            Set::Vector u            = Set::Vector(velocity(i, j, k, 0), velocity(i, j, k, 1));
            Set::Vector u0           = Set::Vector(_u0(i, j, k, 0), _u0(i, j, k, 1));

            Set::Matrix gradM   = Numeric::Gradient(M, i, j, k, DX);
            Set::Vector gradrho = Numeric::Gradient(rho,i,j,k,0,DX);
            Set::Matrix hess_rho = Numeric::Hessian(rho,i,j,k,0,DX,sten);
            Set::Matrix gradu   = (gradM - u*gradrho.transpose()) / rho(i,j,k);

            Set::Vector q0           = Set::Vector(q(i,j,k,0),q(i,j,k,1));


            Set::Scalar mdot0 = -m0(i,j,k) * grad_eta_mag;
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
                            / rho(i,j,k);
                    }

            Set::Vector Ldot0 = Set::Vector::Zero();
            Set::Vector div_tau = Set::Vector::Zero();
            for (int p = 0; p<2; p++)
                for (int q = 0; q<2; q++)
                    for (int r = 0; r<2; r++)
                        for (int s = 0; s<2; s++)
                        {
                            Set::Scalar Mpqrs = 0.0;
                            if (p==r && q==s) Mpqrs += 0.5 * mu;

                            Ldot0(p) += 0.5*Mpqrs * (u(r) - u0(r)) * hess_eta(q, s);
                            div_tau(p) += 2.0*Mpqrs * hess_u(r,s,q);
                        }
            
            Source(i,j, k, 0) = mdot0;
            Source(i,j, k, 1) = Pdot0(0) - Ldot0(0);
            Source(i,j, k, 2) = Pdot0(1) - Ldot0(1);
            Source(i,j, k, 3) = qdot0;// - Ldot0(0)*v(i,j,k,0) - Ldot0(1)*v(i,j,k,1);

            // Lagrange terms to enforce no-penetration
            Source(i,j,k,1) -= lagrange*u.dot(grad_eta)*grad_eta(0);
            Source(i,j,k,2) -= lagrange*u.dot(grad_eta)*grad_eta(1);

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
            
            Solver::Local::Riemann::State state_xlo_fluid = (state_xlo - (1.0 - eta(i-1,j,k))*state_xlo_solid) / (eta(i-1,j,k) + small);
            Solver::Local::Riemann::State state_x_fluid   = (state_x   - (1.0 - eta(i,j,k)  )*state_x_solid  ) / (eta(i,j,k)   + small);
            Solver::Local::Riemann::State state_xhi_fluid = (state_xhi - (1.0 - eta(i+1,j,k))*state_xhi_solid) / (eta(i+1,j,k) + small);

            Solver::Local::Riemann::State state_ylo_fluid = (state_ylo - (1.0 - eta(i,j-1,k))*state_ylo_solid) / (eta(i,j-1,k) + small);
            Solver::Local::Riemann::State state_y_fluid =   (state_y   - (1.0 - eta(i,j,k)  )*state_y_solid  ) / (eta(i,j,k)   + small);
            Solver::Local::Riemann::State state_yhi_fluid = (state_yhi - (1.0 - eta(i,j+1,k))*state_yhi_solid) / (eta(i,j+1,k) + small);

            Solver::Local::Riemann::Flux flux_xlo, flux_ylo, flux_xhi, flux_yhi;

            try
            {
                //lo interface fluxes
                flux_xlo = roesolver->Solve(state_xlo_fluid, state_x_fluid, gamma, pref, small) * eta(i,j,k);
                flux_ylo = roesolver->Solve(state_ylo_fluid, state_y_fluid, gamma, pref, small) * eta(i,j,k);

                //hi interface fluxes
                flux_xhi = roesolver->Solve(state_x_fluid, state_xhi_fluid, gamma, pref, small) * eta(i,j,k);
                flux_yhi = roesolver->Solve(state_y_fluid, state_yhi_fluid, gamma, pref, small) * eta(i,j,k);
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
                Source(i, j, k, 0);

            rho_new(i, j, k) = rho(i, j, k) + 
                (
                    drhof_dt +
                    // todo add drhos_dt term if want time-evolving rhos
                    etadot(i,j,k) * (rho(i,j,k) - rho_solid(i,j,k)) / (eta(i,j,k) + small)
                    ) * dt;

            if (rho_new(i,j,k) != rho_new(i,j,k))
            {
                Util::ParallelMessage(INFO,"lev=",lev);
                Util::ParallelMessage(INFO,"i=",i,"j=",j);
                Util::ParallelMessage(INFO,"drhof_dt",drhof_dt); // dies
                Util::ParallelMessage(INFO,"flux_xlo.mass",flux_xlo.mass);
                Util::ParallelMessage(INFO,"flux_xhi.mass",flux_xhi.mass); // dies, depends on state_xx, state_xhi, state_x_solid, state_xhi_solid, gamma, eta, pref, small
                Util::ParallelMessage(INFO,"flux_ylo.mass",flux_ylo.mass);
                Util::ParallelMessage(INFO,"flux_xhi.mass",flux_yhi.mass);
                Util::ParallelMessage(INFO,"eta",eta(i,j,k));
                Util::ParallelMessage(INFO,"Source",Source(i,j,k,0));
                Util::ParallelMessage(INFO,"state_x",state_x); // <<<<
                Util::ParallelMessage(INFO,"state_y",state_y);
                Util::ParallelMessage(INFO,"state_x_solid",state_x_solid); // <<<<
                Util::ParallelMessage(INFO,"state_y_solid",state_y_solid);
                Util::ParallelMessage(INFO,"state_xhi",state_xhi); // <<<<
                Util::ParallelMessage(INFO,"state_yhi",state_yhi);
                Util::ParallelMessage(INFO,"state_xhi_solid",state_xhi_solid);
                Util::ParallelMessage(INFO,"state_yhi_solids",state_yhi_solid);
                Util::ParallelMessage(INFO,"state_xlo",state_xlo);
                Util::ParallelMessage(INFO,"state_ylo",state_ylo);
                Util::ParallelMessage(INFO,"state_xlo_solid",state_xlo_solid);
                Util::ParallelMessage(INFO,"state_ylo_solid",state_ylo_solid);
                Util::Exception(INFO);
            }

                
            Set::Scalar dMxf_dt =
                (flux_xlo.momentum_normal  - flux_xhi.momentum_normal ) / DX[0] +
                (flux_ylo.momentum_tangent - flux_yhi.momentum_tangent) / DX[1] +
                div_tau(0) * eta(i,j,k) +
                //(mu * (lap_ux * eta(i, j, k))) +
                Source(i, j, k, 1);

            M_new(i, j, k, 0) = M(i, j, k, 0) +
                ( 
                    dMxf_dt + 
                    // todo add dMs_dt term if want time-evolving Ms
                    etadot(i,j,k)*(M(i,j,k,0) - M_solid(i,j,k,0)) / (eta(i,j,k) + small)
                    ) * dt;

            Set::Scalar dMyf_dt =
                (flux_xlo.momentum_tangent - flux_xhi.momentum_tangent) / DX[0] +
                (flux_ylo.momentum_normal  - flux_yhi.momentum_normal ) / DX[1] +
                div_tau(1) * eta(i,j,k) + 
                //(mu * (lap_uy * eta(i, j, k))) +
                Source(i, j, k, 2);
                
            M_new(i, j, k, 1) = M(i, j, k, 1) +
                ( 
                    dMyf_dt +
                    // todo add dMs_dt term if want time-evolving Ms
                    etadot(i,j,k)*(M(i,j,k,1) - M_solid(i,j,k,1)) / (eta(i,j,k)+small)
                    )*dt;

            Set::Scalar dEf_dt =
                (flux_xlo.energy - flux_xhi.energy) / DX[0] +
                (flux_ylo.energy - flux_yhi.energy) / DX[1] +
                Source(i, j, k, 3);
                
            E_new(i, j, k) = (E(i, j, k) + pref/(gamma - 1.0))*pow(rho_new(i,j,k)/rho(i,j,k),gamma-1.0) - pref/(gamma - 1.0) +
                ( 
                    dEf_dt +
                    // todo add dEs_dt term if want time-evolving Es
                    etadot(i,j,k)*(E(i,j,k) - E_solid(i,j,k)) / (eta(i,j,k)+small)
                    ) * dt;


            if (eta(i,j,k) < cutoff)
            {
                rho_new(i,j,k,0) = rho_solid(i,j,k,0);
                M_new(i,j,k,0)   = M_solid(i,j,k,0);
                M_new(i,j,k,1)   = M_solid(i,j,k,1);
                E_new(i,j,k,0)   = E_solid(i,j,k,0);
            }


            //Set::Vector grad_ux = Numeric::Gradient(v, i, j, k, 0, DX);
            //Set::Vector grad_uy = Numeric::Gradient(v, i, j, k, 1, DX);

            *dt_max_handle =                          std::fabs(cfl * DX[0] / (u(0)*eta(i,j,k) + small));
            *dt_max_handle = std::min(*dt_max_handle, std::fabs(cfl * DX[1] / (u(1)*eta(i,j,k) + small)));
            *dt_max_handle = std::min(*dt_max_handle, std::fabs(cfl_v * DX[0]*DX[0] / (Source(i,j,k,1)+small)));
            *dt_max_handle = std::min(*dt_max_handle, std::fabs(cfl_v * DX[1]*DX[1] / (Source(i,j,k,2)+small)));

            // Compute vorticity
            omega(i, j, k) = eta(i, j, k) * (gradu(1,0) - gradu(0,1));

        });

    }
    this->DynamicTimestep_SyncTimeStep(lev,dt_max);

}//end Advance

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
    for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<char> const& tags = a_tags.array(mfi);
        amrex::Array4<const Set::Scalar> const& eta = (*eta_mf[lev]).array(mfi);

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
