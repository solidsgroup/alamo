
#include "Hydro.H"
#include "IO/ParmParse.H"
#include "BC/Constant.H"
#include "BC/Nothing.H"
#include "BC/Expression.H"
#include "Numeric/Stencil.H"
#include "Numeric/Function.H"
#include "IC/Laminate.H"
#include "IC/PSRead.H"
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
        pp.query_default("gradu_refinement_criterion", value.gradu_refinement_criterion, 0.01);
        // pressure-based refinement
        pp.query_default("p_refinement_criterion", value.p_refinement_criterion, 1e100);
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

        value.RegisterNewFab(value.momentum_mf,     value.momentum_bc, 2, nghost, "momentum",     true );
        value.RegisterNewFab(value.momentum_old_mf, value.momentum_bc, 2, nghost, "momentum_old", false);
 
        value.RegisterNewFab(value.pressure_mf,  &value.bc_nothing, 1, nghost, "pressure",  true);
        value.RegisterNewFab(value.velocity_mf,  &value.bc_nothing, 2, nghost, "velocity",  true);
        value.RegisterNewFab(value.vorticity_mf, &value.bc_nothing, 1, nghost, "vorticity", true);

        value.RegisterNewFab(value.rho_injected_mf, &value.bc_nothing, 1, 0, "rho_injected", true);
        value.RegisterNewFab(value.mdot_mf,         &value.bc_nothing, 2, 0, "mdot",         true);
        value.RegisterNewFab(value.q_mf,            &value.bc_nothing, 2, 0, "q",            true);
        value.RegisterNewFab(value.flux_mf,         &value.bc_nothing, 1, 0, "flux",            true);

        value.RegisterNewFab(value.solid.momentum_mf, &value.neumann_bc_D, 2, nghost, "solid.momentum", true);
        value.RegisterNewFab(value.solid.density_mf,  &value.neumann_bc_1,  1, nghost, "solid.density", true);
        value.RegisterNewFab(value.solid.energy_mf,   &value.neumann_bc_1, 1, nghost, "solid.energy",   true);

        value.RegisterNewFab(value.Source_mf, &value.bc_nothing, 4, 0, "Source", true);
    }
    {
        std::string type = "constant";
        // eta initial condition
        pp_query_validate("eta.ic.type", type, {"constant","laminate","expression","bmp","png"});
        if (type == "constant") value.eta_ic = new IC::Constant(value.geom, pp, "eta.ic.constant");
        else if (type == "laminate") value.eta_ic = new IC::Laminate(value.geom, pp, "eta.ic.laminate");
        else if (type == "expression") value.eta_ic = new IC::Expression(value.geom, pp, "eta.ic.expression");
        else if (type == "bmp") value.eta_ic = new IC::BMP(value.geom, pp, "eta.ic.bmp");
        else if (type == "png") value.eta_ic = new IC::PNG(value.geom, pp, "eta.ic.png");
        else Util::Abort(INFO, "Invalid eta.ic: ", type);
    }
    {
        std::string type = "constant";
        pp_forbid("Velocity.ic.type", "--> velocity.ic.type");
        // velocity initial condition
        pp_query_validate("velocity.ic.type", type,{"constant","expression"});
        if (type == "constant") value.velocity_ic = new IC::Constant(value.geom, pp, "velocity.ic.constant");
        else if (type == "expression") value.velocity_ic = new IC::Expression(value.geom, pp, "velocity.ic.expression");
    }
    {
        std::string type = "constant";
        pp_forbid("Pressure.ic", "--> pressure.ic");
        // solid pressure IC type
        pp_query_validate("pressure.ic.type", type, {"constant","expression"});
        if (type == "constant") value.pressure_ic = new IC::Constant(value.geom, pp, "pressure.ic.constant");
        else if (type == "expression") value.pressure_ic = new IC::Expression(value.geom, pp, "pressure.ic.expression");
        else Util::Abort(INFO, "Invalid Pressure.ic: ", type);
    }
    {
        std::string type = "constant";
        pp_forbid("SolidMomentum.ic", "--> solid.momentum.ic");
        // solid momentum IC type
        pp_query_validate("solid.momentum.ic.type", type, {"constant","expression"});
        if (type == "constant") value.solid.momentum_ic = new IC::Constant(value.geom, pp, "solid.momentum.ic.constant");
        else if (type == "expression") value.solid.momentum_ic = new IC::Expression(value.geom, pp, "solid.momentum.ic.expression");
        else Util::Abort(INFO, "Invalid solid.momentum.ic: ", type);
    }
    {
        std::string type = "constant";
        pp_forbid("SolidDensity.ic.type", "--> solid.density.ic.type");
        // solid density IC type
        pp_query_validate("solid.density.ic.type", type, {"constant","expression"});
        if (type == "constant") value.solid.density_ic = new IC::Constant(value.geom, pp, "solid.density.ic.constant");
        else if (type == "expression") value.solid.density_ic = new IC::Expression(value.geom, pp, "solid.density.ic.expression");
        else Util::Abort(INFO, "Invalid solid.density.ic: ", type);
    }
    {
        std::string type = "constant";
        pp_forbid("SolidEnergy.ic.type", "--> solid.energy.ic.type");
        // solid energy IC type
        pp_query_validate("solid.energy.ic.type", type, {"constant","expression"});
        if (type == "constant") value.solid.energy_ic = new IC::Constant(value.geom, pp, "solid.energy.ic.constant");
        else if (type == "expression") value.solid.energy_ic = new IC::Expression(value.geom, pp, "solid.energy.ic.expression");
        else Util::Abort(INFO, "Invalid solid.energy.ic: ", type);
    }
    {
        std::string type = "constant";
        pp_forbid("Density.ic.type", "--> density.ic.type");
        // density initial condition type
        pp_query_validate("density.ic.type", type, {"constant","expression"});
        if (type == "constant") value.density_ic = new IC::Constant(value.geom, pp, "density.ic.constant");
        else if (type == "expression") value.density_ic = new IC::Expression(value.geom, pp, "density.ic.expression");
        else Util::Abort(INFO, "Invalid density.ic: ", type);
    }
    {
        std::string type = "constant";
        // injected density initial condition type
        pp_query_validate("rho_injected.ic.type", type, {"constant","expression"});
        if (type == "constant") value.ic_rho_injected = new IC::Constant(value.geom, pp, "rho_injected.ic.constant");
        else if (type == "expression") value.ic_rho_injected = new IC::Expression(value.geom, pp, "rho_injected.ic.expression");
        else Util::Abort(INFO, "Invalid rho_injected.ic: ", type);
    }
    {
        std::string type = "constant";
        // mass flux initial condition
        pp.query("mdot.ic.type", type);
        if (type == "constant") value.ic_mdot = new IC::Constant(value.geom, pp, "mdot.ic.constant");
        else if (type == "expression") value.ic_mdot = new IC::Expression(value.geom, pp, "mdot.ic.expression");
        else Util::Abort(INFO, "Invalid mdot.ic: ", type);
    }
    {
        std::string type = "constant";
        // injected heat initial condition
        pp.query("q.ic.type", type);
        if (type == "constant") value.ic_q = new IC::Constant(value.geom, pp, "q.ic.constant");
        else if (type == "expression") value.ic_q = new IC::Expression(value.geom, pp, "q.ic.expression");
        else Util::Abort(INFO, "Invalid q.ic: ", type);
    }
}


void Hydro::Initialize(int lev)
{
    BL_PROFILE("Integrator::Hydro::Initialize");
 
    eta_ic           ->Initialize(lev, eta_mf,     0.0);
    eta_ic           ->Initialize(lev, eta_old_mf, 0.0);
    etadot_mf[lev]   ->setVal(0.0);

    flux_mf[lev]   ->setVal(0.0);

    velocity_ic      ->Initialize(lev, velocity_mf, 0.0);
    pressure_ic      ->Initialize(lev, pressure_mf, 0.0);
    density_ic       ->Initialize(lev, density_mf, 0.0);

    density_ic       ->Initialize(lev, density_old_mf, 0.0);


    solid.density_ic ->Initialize(lev, solid.density_mf, 0.0);
    solid.momentum_ic->Initialize(lev, solid.momentum_mf, 0.0);
    solid.energy_ic  ->Initialize(lev, solid.energy_mf, 0.0);

    ic_rho_injected  ->Initialize(lev, rho_injected_mf, 0.0);
    ic_mdot          ->Initialize(lev, mdot_mf,    0.0);
    ic_q             ->Initialize(lev, q_mf,            0.0);

    Source_mf[lev]   ->setVal(0.0);

    Mix(lev);
}

void Hydro::Mix(int lev)
{
    Util::Message(INFO, eta_mf[lev]->nComp());

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
        Set::Patch<const Set::Scalar> E         = energy_old_mf.Patch(lev,mfi);

        Set::Patch<const Set::Scalar> rho_solid = solid.density_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> M_solid   = solid.momentum_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> E_solid   = solid.energy_mf.Patch(lev,mfi);

        Set::Patch<Set::Scalar>       v         = velocity_mf.Patch(lev,mfi);
        Set::Patch<Set::Scalar>       p         = pressure_mf.Patch(lev,mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            etadot(i, j, k) = (eta_new(i, j, k) - eta(i, j, k)) / dt;

            Set::Scalar etarho_fluid  = (rho(i,j,k) - (1.-eta(i,j,k)) * rho_solid(i,j,k));
            Set::Scalar etaE_fluid    = E(i,j,k)   - (1.-eta(i,j,k)) * E_solid(i,j,k);

            Set::Vector etaM_fluid( M(i,j,k,0) - (1.-eta(i,j,k)) * M_solid(i,j,k,0),
                                    M(i,j,k,1) - (1.-eta(i,j,k)) * M_solid(i,j,k,1) );

            //THESE ARE FLUID VELOCITY AND PRESSURE

            v(i,j,k,0) = etaM_fluid(0) / (etarho_fluid + small);
            v(i,j,k,1) = etaM_fluid(1) / (etarho_fluid + small);

            p(i,j,k)   = (etaE_fluid - 0.5 * (etaM_fluid(0)*etaM_fluid(0) + etaM_fluid(1)*etaM_fluid(1)) / (etarho_fluid + small)) * ((gamma - 1.0) * (eta(i, j, k) + small));
        });
    }

    const Set::Scalar* DX = geom[lev].CellSize();
    //amrex::Box domain = geom[lev].Domain();

    for (amrex::MFIter mfi(*eta_mf[lev], false); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        amrex::Array4<const Set::Scalar> const& rho = (*density_old_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& E   = (*energy_old_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& M   = (*momentum_old_mf[lev]).array(mfi);

        amrex::Array4<Set::Scalar> const& rho_new = (*density_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& E_new   = (*energy_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& M_new   = (*momentum_mf[lev]).array(mfi);

        amrex::Array4<const Set::Scalar> const& rho_solid = (*solid.density_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& M_solid   = (*solid.momentum_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& E_solid   = (*solid.energy_mf[lev]).array(mfi);

        amrex::Array4<Set::Scalar> const& omega   = (*vorticity_mf[lev]).array(mfi);

        amrex::Array4<Set::Scalar> const& flux   = (*flux_mf[lev]).array(mfi);

        amrex::Array4<const Set::Scalar> const& eta    = (*eta_old_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& etadot = (*etadot_mf[lev]).array(mfi);

        amrex::Array4<const Set::Scalar> const& v          = (*velocity_mf[lev]).array(mfi);

        amrex::Array4<const Set::Scalar> const& rho_injected = (*rho_injected_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& q            = (*q_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& mdot         = (*mdot_mf[lev]).array(mfi);


        amrex::Array4<Set::Scalar> const& Source = (*Source_mf[lev]).array(mfi);

        Set::Scalar *dt_max_handle = &dt_max;

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {   
            //Diffuse Sources
            Set::Vector grad_eta     = Numeric::Gradient(eta, i, j, k, 0, DX);
            Set::Scalar grad_eta_mag = grad_eta.lpNorm<2>();
            Set::Matrix hess_eta     = Numeric::Hessian(eta, i, j, k, 0, DX);
            Set::Scalar lap_eta      = Numeric::Laplacian(eta, i, j, k, 0, DX);

            Set::Scalar boundary_energy = 0.5 * (mdot(i, j, k, 0) * mdot(i, j, k, 0) + mdot(i, j, k, 1) * mdot(i, j, k, 1)) / rho_injected(i, j, k);
            Set::Vector u_injected      = Set::Vector(mdot(i, j, k, 0), mdot(i, j, k, 1)) / rho_injected(i, j, k);
            //Set::Vector u_interface     = -etadot(i, j, k)*grad_eta/(grad_eta_mag * grad_eta_mag + small);
            Set::Vector u_interface     = Set::Vector::Zero();

            Set::Matrix I            = Set::Matrix::Identity();
            Set::Vector u_applied    = u_injected + u_interface /*+ Set::Vector(0.1, 0.0)*/;
            Set::Vector u            = Set::Vector(v(i, j, k, 0), v(i, j, k, 1));
            //Set::Matrix gradu        = Numeric::Gradient(v, i, j, k, DX);


            Set::Matrix gradM   = Numeric::Gradient(M, i, j, k, DX);
            Set::Vector gradrho = Numeric::Gradient(rho,i,j,k,0,DX);
            Set::Matrix gradu   = (gradM - u*gradrho.transpose()) / rho(i,j,k);
            Set::Scalar divu    = gradu.trace();

            //Set::Scalar divu         = 0.0;//Numeric::Divergence(v, i, j, k, 2, DX);
            Set::Matrix T            = 0.5 * mu * (gradu + gradu.transpose());    
            Set::Vector q0           = Set::Vector(q(i,j,k,0),q(i,j,k,1));

            //Boundary flux
            //states of solid/fluid boundary cells, eta = 1.0
            Solver::Local::Riemann::Roe::State solid_state(rho_solid(i,j,k), M_solid(i,j,k,0), M_solid(i,j,k,1), E_solid(i,j,k), eta(i,j,k));
            //Solver::Local::Riemann::Roe::State boundary_state(rho_injected(i, j, k), mdot(i, j, k, 0), mdot(i, j, k, 1), boundary_energy, 1.0);
            Solver::Local::Riemann::Roe::State boundary_state(0.0,0.0,0.0,0.0, 1.0);
            Solver::Local::Riemann::Roe::State empty_state(small, 0.0, 0.0, 0.0, 1.0);
            Solver::Local::Riemann::Roe::State current_state(rho(i,j,k), M(i,j,k,0), M(i,j,k,1), E(i,j,k), eta(i,j,k));

            //Solver::Local::Riemann::Roe::Flux prescribed_boundary_flux = Solver::Local::Riemann::Roe::Solve(boundary_state, boundary_state, empty_state, empty_state, gamma, 1.0, 0.0, small);
            //Solver::Local::Riemann::Roe::Flux prescribed_boundary_flux = Solver::Local::Riemann::Roe::Solve(current_state, boundary_state, solid_state, solid_state, gamma, eta(i,j,k), pref, small);

            //Active Source Terms
            // Set::Scalar mdot0 = -1000.*grad_eta.dot(u); //-prescribed_boundary_flux.mass * grad_eta_mag;
            // Set::Vector Pdot0 = Set::Vector::Zero();; //-Set::Vector(prescribed_boundary_flux.momentum_normal * grad_eta(0), prescribed_boundary_flux.momentum_tangent * grad_eta(1)) - T*grad_eta;
            // Set::Scalar qdot0 = 0.0; //-prescribed_boundary_flux.energy * grad_eta_mag - (T*u).dot(grad_eta) + q0.dot(grad_eta); 
            //Set::Vector Ldot0 = rho(i, j, k) * (u*u.transpose() - u_applied*u_applied.transpose())*grad_eta;

            Set::Scalar mdot0 = 0.0; //-prescribed_boundary_flux.mass * grad_eta_mag;
            //Set::Vector Pdot0 = rho(i,j,k)*(u_applied - u)*grad_eta_mag; //-Set::Vector(prescribed_boundary_flux.momentum_normal * grad_eta(0), prescribed_boundary_flux.momentum_tangent * grad_eta(1)) - T*grad_eta;
            Set::Vector Pdot0 = Set::Vector::Zero(); //- Pfactor*(u.dot(grad_eta)) * grad_eta/(grad_eta_mag + small); //-Set::Vector(prescribed_boundary_flux.momentum_normal * grad_eta(0), prescribed_boundary_flux.momentum_tangent * grad_eta(1)) - T*grad_eta;
            Set::Scalar qdot0 = 0.0; //-prescribed_boundary_flux.energy * grad_eta_mag - (T*u).dot(grad_eta) + q0.dot(grad_eta); 



            Set::Vector Ldot0 = Set::Vector::Zero();

            for (int p = 0; p<2; p++)
            for (int q = 0; q<2; q++)
            for (int r = 0; r<2; r++)
            for (int s = 0; s<2; s++)
            {
                Set::Scalar Mpqrs = 0.0;
                if (p==r && q==s) Mpqrs += 0.5 * mu;
                if (p==s && q==r) Mpqrs += 0.5 * mu;
                if (p==q && r==s) Mpqrs += 0.5 * mu;
                Ldot0(p) += 0.5 * Mpqrs * (v(i, j, k, r) - u_applied(r)) * hess_eta(q, s);
            }
            
            Source(i,j, k, 0) = mdot0;
            Source(i,j, k, 1) = Pdot0(0) - Ldot0(0);
            Source(i,j, k, 2) = Pdot0(1) - Ldot0(1);
            Source(i,j, k, 3) = qdot0;// - Ldot0(0)*v(i,j,k,0) - Ldot0(1)*v(i,j,k,1);

            //Viscous Terms
            Set::Scalar lapMx  = Numeric::Laplacian(M, i, j, k, 0, DX);
            Set::Scalar lapMy  = Numeric::Laplacian(M, i, j, k, 1, DX);
            Set::Scalar laprho = Numeric::Laplacian(rho, i, j, k, 0, DX);

            Set::Scalar lap_ux = (lapMx - laprho*u(0) - 2.0*gradrho(0)*gradu(0,0)) / rho(i,j,k);
            Set::Scalar lap_uy = (lapMy - laprho*u(1) - 2.0*gradrho(1)*gradu(1,1)) / rho(i,j,k);
            
            // Set::Scalar lap_ux = Numeric::Laplacian(v, i, j, k, 0, DX);
            // Set::Scalar lap_uy = Numeric::Laplacian(v, i, j, k, 1, DX);
            //Set::Scalar div_u  = Numeric::Divergence(v, i, j, k, 2, DX); // currently causes error!
            // Set::Vector grad_ux = Numeric::Gradient(v, i, j, k, 0, DX);
            // Set::Vector grad_uy = Numeric::Gradient(v, i, j, k, 1, DX);
            //Set::Vector grad_eta = Numeric::Gradient(eta, i, j, k, 0, DX);

            //Godunov flux
            //states of total fields
            Solver::Local::Riemann::Roe::State   state_x(rho(i, j, k), M(i, j, k, 0), M(i, j, k, 1), E(i, j, k), eta(i, j, k));
            Solver::Local::Riemann::Roe::State   state_y(rho(i, j, k), M(i, j, k, 1), M(i, j, k, 0), E(i, j, k), eta(i, j, k));

            Solver::Local::Riemann::Roe::State lo_statex(rho(i - 1, j, k), M(i - 1, j, k, 0), M(i - 1, j, k, 1), E(i - 1, j, k), eta(i - 1, j, k));
            Solver::Local::Riemann::Roe::State hi_statex(rho(i + 1, j, k), M(i + 1, j, k, 0), M(i + 1, j, k, 1), E(i + 1, j, k), eta(i + 1, j, k));

            Solver::Local::Riemann::Roe::State lo_statey(rho(i, j - 1, k), M(i, j - 1, k, 1), M(i, j - 1, k, 0), E(i, j - 1, k), eta(i, j - 1, k));
            Solver::Local::Riemann::Roe::State hi_statey(rho(i, j + 1, k), M(i, j + 1, k, 1), M(i, j + 1, k, 0), E(i, j + 1, k), eta(i, j + 1, k));
            
            //states of solid fields
            Solver::Local::Riemann::Roe::State    statex_solid(rho_solid(i, j, k), M_solid(i, j, k, 0), M_solid(i, j, k, 1), E_solid(i, j, k), eta(i, j, k));
            Solver::Local::Riemann::Roe::State    statey_solid(rho_solid(i, j, k), M_solid(i, j, k, 1), M_solid(i, j, k, 0), E_solid(i, j, k), eta(i, j, k));

            Solver::Local::Riemann::Roe::State lo_statex_solid(rho_solid(i - 1, j, k), M_solid(i - 1, j, k, 0), M_solid(i - 1, j, k, 1), E_solid(i - 1, j, k), eta(i - 1, j, k));
            Solver::Local::Riemann::Roe::State hi_statex_solid(rho_solid(i + 1, j, k), M_solid(i + 1, j, k, 0), M_solid(i + 1, j, k, 1), E_solid(i + 1, j, k), eta(i + 1, j, k));

            Solver::Local::Riemann::Roe::State lo_statey_solid(rho_solid(i, j - 1, k), M_solid(i, j - 1, k, 1), M_solid(i, j - 1, k, 0), E_solid(i, j - 1, k), eta(i, j - 1, k));
            Solver::Local::Riemann::Roe::State hi_statey_solid(rho_solid(i, j + 1, k), M_solid(i, j + 1, k, 1), M_solid(i, j + 1, k, 0), E_solid(i, j + 1, k), eta(i, j + 1, k));

            Solver::Local::Riemann::Roe::Flux flux_xlo, flux_ylo, flux_xhi, flux_yhi;

            try
            {
                //lo interface fluxes
                flux_xlo = Solver::Local::Riemann::Roe::Solve(lo_statex, state_x, lo_statex_solid, statex_solid, gamma, eta(i, j, k), pref, small);
                flux_ylo = Solver::Local::Riemann::Roe::Solve(lo_statey, state_y, lo_statey_solid, statey_solid, gamma, eta(i, j, k), pref, small);

                //hi interface fluxes
                flux_xhi = Solver::Local::Riemann::Roe::Solve(state_x, hi_statex, statex_solid, hi_statex_solid, gamma, eta(i, j, k), pref, small);
                flux_yhi = Solver::Local::Riemann::Roe::Solve(state_y, hi_statey, statey_solid, hi_statey_solid, gamma, eta(i, j, k), pref, small);
            }
            catch(...)
            {
                Util::ParallelMessage(INFO,"lev=",lev);
                Util::ParallelMessage(INFO,"i=",i,"j=",j);
                Util::Abort(INFO);
                
            }


            flux(i,j,k) = flux_xhi.mass - flux_xlo.mass;

            Set::Scalar smallmod = small; //(1.0 - eta(i,j,k))*0.001;


            Set::Scalar drhof_dt = 
                (flux_xlo.mass - flux_xhi.mass) / DX[0] +
                (flux_ylo.mass - flux_yhi.mass) / DX[1] +
                Source(i, j, k, 0);

            rho_new(i, j, k) = rho(i, j, k) + 
                (
                    drhof_dt +
                    // todo add drhos_dt term
                    etadot(i,j,k) * (rho(i,j,k) - rho_solid(i,j,k)) / (eta(i,j,k) + smallmod)
                ) * dt;

            if (rho_new(i,j,k) != rho_new(i,j,k))
            {
                Util::ParallelMessage(INFO,"lev=",lev);
                Util::ParallelMessage(INFO,"i=",i,"j=",j);
                Util::ParallelMessage(INFO,"drhof_dt",drhof_dt); // dies
                Util::ParallelMessage(INFO,"flux_xlo.mass",flux_xlo.mass);
                Util::ParallelMessage(INFO,"flux_xhi.mass",flux_xhi.mass); // dies, depends on state_xx, hi_statex, statex_solid, hi_statex_solid, gamma, eta, pref, small
                Util::ParallelMessage(INFO,"flux_ylo.mass",flux_ylo.mass);
                Util::ParallelMessage(INFO,"flux_xhi.mass",flux_yhi.mass);
                Util::ParallelMessage(INFO,"eta",eta(i,j,k));
                Util::ParallelMessage(INFO,"Source",Source(i,j,k,0));
                Util::ParallelMessage(INFO,"state_x",state_x); // <<<<
                Util::ParallelMessage(INFO,"state_y",state_y);
                Util::ParallelMessage(INFO,"statex_solid",statex_solid); // <<<<
                Util::ParallelMessage(INFO,"statey_solid",statey_solid);
                Util::ParallelMessage(INFO,"hi_statex",hi_statex); // <<<<
                Util::ParallelMessage(INFO,"hi_statey",hi_statey);
                Util::ParallelMessage(INFO,"hi_statex_solid",hi_statex_solid);
                Util::ParallelMessage(INFO,"hi_statey_solids",hi_statey_solid);
                Util::ParallelMessage(INFO,"lo_statex",lo_statex);
                Util::ParallelMessage(INFO,"lo_statey",lo_statey);
                Util::ParallelMessage(INFO,"lo_statex_solid",lo_statex_solid);
                Util::ParallelMessage(INFO,"lo_statey_solid",lo_statey_solid);
                Util::Exception(INFO);
            }

                
            Set::Scalar dMxf_dt =
                (flux_xlo.momentum_normal  - flux_xhi.momentum_normal ) / DX[0] +
                (flux_ylo.momentum_tangent - flux_yhi.momentum_tangent) / DX[1] +
                (mu * (lap_ux * eta(i, j, k))) +
                Source(i, j, k, 1);

            M_new(i, j, k, 0) = M(i, j, k, 0) +
                ( 
                    dMxf_dt + 
                    // todo add dMs_dt term
                    etadot(i,j,k)*(M(i,j,k,0) - M_solid(i,j,k,0)) / (eta(i,j,k) + smallmod)
                ) * dt;
                

            Set::Scalar dMyf_dt =
                (flux_xlo.momentum_tangent - flux_xhi.momentum_tangent) / DX[0] +
                (flux_ylo.momentum_normal  - flux_yhi.momentum_normal ) / DX[1] +
                (mu * (lap_uy * eta(i, j, k))) +
                Source(i, j, k, 2);
                
            M_new(i, j, k, 1) = M(i, j, k, 1) +
                ( 
                    dMyf_dt +
                    // todo add dMs_dt term
                    etadot(i,j,k)*(M(i,j,k,1) - M_solid(i,j,k,1)) / (eta(i,j,k)+smallmod)
                )*dt;

            Set::Scalar dEf_dt =
                (flux_xlo.energy - flux_xhi.energy) / DX[0] +
                (flux_ylo.energy - flux_yhi.energy) / DX[1] +
                Source(i, j, k, 3);
                
            E_new(i, j, k) = E(i, j, k) + 
                ( 
                    dEf_dt +
                    // todo add dEs_dt term
                    etadot(i,j,k)*(E(i,j,k) - E_solid(i,j,k)) / (eta(i,j,k)+smallmod)
                ) * dt;

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

}//end TagCells

// void Hydro::Integrate(int amrlev, Set::Scalar /*time*/, int /*step*/, const amrex::MFIter &mfi, const amrex::Box &box)
// {
//   BL_PROFILE("Hydro::Integrate");
//   const Set::Scalar *DX = geom[amrlev].CellSize();
//   Set::Scalar dv = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);
//   amrex::Array4<amrex::Real> const &eta = (*eta_mf[amrlev]).array(mfi);
//   amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k){
//   volume += eta(i, j, k, 0) * dv;
//   Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX);
//   Set::Scalar normgrad = grad.lpNorm<2>();
//   Set::Scalar da = normgrad * dv;
//   area += da;
//   });


//  }//end Integrate

//  void Hydro::UpdateModel(int /*a_step*/)
//  {
//    for (int lev = 0; lev <= finest_level; ++lev)
//    {
//      eta_mf[lev] -> FillBoundary();
//      density_mf[lev] -> FillBoundary();
//      energy_mf[lev] -> FillBoundary();
//      Momentum[lev] -> FillBoundary();
//      
//      for (MFIter mfi(*model_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
//      {
//  amrex::Box bx = mfi.nodaltilebox();
//  amrex::Array4<model_type> const &model = model_mf[lev]->array(mfi);
//
//  amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k){
//  // TODO
//
//  });
//
//
//     } // end For2

//     Util::RealFillBoundary(*model_mf[lev], geom[lev]);
//     amrex::MultiFab::Copy(*psi_mf[lev], *eta_mf[lev], 0, 0, 1, psi_mf[lev]-> nGrow());

//    } //end For1
//}//end update


}//end code


#endif
