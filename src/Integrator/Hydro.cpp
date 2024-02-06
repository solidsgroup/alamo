#include "Hydro.H"
#include "IO/ParmParse.H"
#include "BC/Constant.H"
#include "Numeric/Stencil.H"
#include "IC/Laminate.H"
#include "IC/PSRead.H"
#include "IC/Expression.H"
#include "Solver/Local/Riemann/Roe.H"

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
        pp.query("r_refinement_criterion", value.r_refinement_criterion);
        pp.query("e_refinement_criterion", value.e_refinement_criterion);
        pp.query("m_refinement_criterion", value.m_refinement_criterion);
        pp.query("eta_refinement_criterion", value.eta_refinement_criterion);

        pp.query("gamma", value.gamma);
        pp.query("cfl", value.cfl);
        pp.query("mu", value.mu);

        pp.query("rho_solid", value.rho_solid);
        pp.query("rho_fluid", value.rho_fluid);

	pp.query("E_solid", value.E_solid);

        pp.query("rho_interface", value.rho_interface);
        pp.query("v_x_interface", value.v_x_interface);
        pp.query("v_y_interface", value.v_y_interface);
        pp.query("delta_p_interface", value.delta_p_interface);

        value.bc_eta = new BC::Constant(1, pp, "pf.eta.bc");
        value.bc_rho = new BC::Constant(1, pp, "rho.bc");
        value.bc_p = new BC::Constant(1, pp, "p.bc");
        value.bc_v = new BC::Constant(2, pp, "v.bc");
	
	//value.bc_M = new BC::Constant(2, pp, "M.bc");
        //value.bc_E = new BC::Constant(1, pp, "E.bc");
	
	//value.bc_omega = new BC::Constant(2, pp, "omega.bc");
    }
    // Register FabFields:
    {
        int nghost = 2;
        value.RegisterNewFab(value.eta_mf, value.bc_eta, 1, nghost, "eta", true);
        value.RegisterNewFab(value.eta_old_mf, value.bc_eta, 1, nghost, "eta_old", false);

        value.RegisterNewFab(value.etadot_mf, value.bc_eta, 1, nghost, "etadot", true);

        value.RegisterNewFab(value.etaDensity_mf, value.bc_rho, 1, nghost, "etaDensity", true);
        value.RegisterNewFab(value.etaDensity_old_mf, value.bc_rho, 1, nghost, "etarho_old", false);
	
        value.RegisterNewFab(value.DensityMix_mf, value.bc_rho, 1, nghost, "DensityMix", true);

        value.RegisterNewFab(value.etaEnergy_mf, value.bc_p, 1, nghost, "etaEnergy", true);
        value.RegisterNewFab(value.etaEnergy_old_mf, value.bc_p, 1, nghost, "etaE_old", false);
	
        value.RegisterNewFab(value.EnergyMix_mf, value.bc_p, 1, nghost, "EnergyMix", true);

        value.RegisterNewFab(value.etaMomentum_mf, value.bc_v, 2, nghost, "etaMomentum", true);
        value.RegisterNewFab(value.etaMomentum_old_mf, value.bc_v, 2, nghost, "etaM_old", false);
	
        value.RegisterNewFab(value.MomentumMix_mf, value.bc_v, 2, nghost, "MomentumMix", true);

        value.RegisterNewFab(value.Velocity_mf, value.bc_v, 2, nghost, "Velocity", true);

        value.RegisterNewFab(value.Vorticity_mf, value.bc_eta, 1, nghost, "Vorticity", true);

        value.RegisterNewFab(value.Pressure_mf, value.bc_p, 1, nghost, "Pressure", true);
    }
    {
        std::string type = "constant";
        pp.query("eta.ic.type", type);
        if (type == "constant") value.ic_eta = new IC::Constant(value.geom, pp, "eta.ic.constant");
        else if (type == "laminate") value.ic_eta = new IC::Laminate(value.geom, pp, "eta.ic.laminate");
        else if (type == "expression") value.ic_eta = new IC::Expression(value.geom, pp, "eta.ic.expression");
        else Util::Abort(INFO, "Invalid eta.ic: ", type);
    }
    {
        std::string type = "constant";
        pp.query("etadot.ic.type", type);
        if (type == "constant") value.ic_etadot = new IC::Constant(value.geom, pp, "etadot.ic.constant");
        else if (type == "laminate") value.ic_etadot = new IC::Laminate(value.geom, pp, "etadot.ic.laminate");
        else if (type == "expression") value.ic_etadot = new IC::Expression(value.geom, pp, "etadot.ic.expression");
        else Util::Abort(INFO, "Invalid eta.ic: ", type);
    }
    {
        std::string type = "constant";
        pp.query("Velocity.ic.type", type);
        if (type == "constant") value.ic_Velocity = new IC::Constant(value.geom, pp, "Velocity.ic.constant");
        else if (type == "expression") value.ic_Velocity = new IC::Expression(value.geom, pp, "Velocity.ic.expression");
        else Util::Abort(INFO, "Invalid Velocity.ic: ", type);
    }
    {
        std::string type = "constant";
        pp.query("Pressure.ic.type", type);
        if (type == "constant") value.ic_Pressure = new IC::Constant(value.geom, pp, "Pressure.ic.constant");
        else if (type == "expression") value.ic_Pressure = new IC::Expression(value.geom, pp, "Pressure.ic.expression");
        else Util::Abort(INFO, "Invalid Pressure.ic: ", type);
    }
}


void Hydro::Initialize(int lev)
{
    BL_PROFILE("Integrator::Hydro::Initialize");

    ic_eta->Initialize(lev, eta_mf, 0.0);
    ic_eta->Initialize(lev, eta_old_mf, 0.0);
    etadot_mf[lev]->setVal(0.0);

    ic_Velocity->Initialize(lev, Velocity_mf, 0.0);
    
    ic_Pressure->Initialize(lev, Pressure_mf, 0.0);

    Util::Message(INFO, eta_mf[lev] -> nComp());

    for (amrex::MFIter mfi(*eta_mf[lev], false); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        amrex::Array4<Set::Scalar> const& E_mix = (*EnergyMix_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& rho_mix = (*DensityMix_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& M_mix = (*MomentumMix_mf[lev]).array(mfi);

	amrex::Array4<Set::Scalar> const& etaE_old = (*etaEnergy_old_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const& etaE = (*etaEnergy_mf[lev]).array(mfi);
	
        amrex::Array4<Set::Scalar> const& etarho_old = (*etaDensity_old_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const& etarho = (*etaDensity_mf[lev]).array(mfi);
	
        amrex::Array4<Set::Scalar> const& etaM_old = (*etaMomentum_old_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const& etaM = (*etaMomentum_mf[lev]).array(mfi);

        amrex::Array4<Set::Scalar> const& v = (*Velocity_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& p = (*Pressure_mf[lev]).array(mfi);

        amrex::Array4<const Set::Scalar> const& eta = (*eta_mf[lev]).array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
	  rho_mix(i, j, k) = eta(i, j, k) * rho_fluid + (1.0 - eta(i, j, k)) * rho_solid;
	  etarho(i, j, k) = eta(i, j, k) * rho_fluid;
	  etarho_old(i, j, k) = etarho(i, j, k);

	  E_mix(i, j, k) = eta(i, j, k) * (0.5 * rho_mix(i, j, k) * (v(i, j, k, 0) * v(i, j, k, 0) + v(i, j, k, 1) * v(i, j, k, 1)) + p(i, j, k)/(gamma - 1.0)) + (1.0 - eta(i, j, k)) * E_solid;
    	  etaE(i, j, k) = eta(i, j, k) * E_mix(i, j, k);
	  etaE_old(i, j, k) = etaE(i, j, k);

	  M_mix(i, j, k, 0) = rho_mix(i, j, k) * v(i, j, k, 0);
	  M_mix(i, j, k, 1) = rho_mix(i, j, k) * v(i, j, k, 0);
	  etaM(i, j, k, 0) = M_mix(i, j, k, 0) * eta(i, j, k);
	  etaM(i, j, k, 1) = M_mix(i, j, k, 1) * eta(i, j, k);
	  etaM_old(i, j, k, 0) = etaM(i, j, k, 0);
	  etaM_old(i, j, k, 1) = etaM(i, j, k, 1);
        });
    }

    c_max = 0.0;
    vx_max = 0.0;
    vy_max = 0.0;
}

void Hydro::UpdateEta(Set::Scalar)
{
    //for (int lev = 0; lev <= finest_level; ++lev)
    //{
    //    ic_eta->Initialize(lev, eta_mf, time);
    //    ic_etadot->Initialize(lev, etadot_mf, time);
    //}
}

void Hydro::TimeStepBegin(Set::Scalar, int /*iter*/)
{
    //UpdateEta(time);
}


void Hydro::TimeStepComplete(Set::Scalar, int lev)
{
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
    ic_eta->Initialize(lev, eta_mf, time);
    for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox();
        amrex::Array4<const Set::Scalar> const& eta_new = (*eta_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& eta = (*eta_old_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar>       const& etadot = (*etadot_mf[lev]).array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            etadot(i, j, k) = (eta_new(i, j, k) - eta(i, j, k)) / dt;
        });
    }

    std::swap(etaMomentum_old_mf[lev], etaMomentum_mf[lev]);
    std::swap(etaEnergy_old_mf[lev], etaEnergy_mf[lev]);
    std::swap(etaDensity_old_mf[lev], etaDensity_mf[lev]);

    const Set::Scalar* DX = geom[lev].CellSize();

    

    for (amrex::MFIter mfi(*eta_mf[lev], false); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        amrex::Array4<Set::Scalar> const& E_mix = (*EnergyMix_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& rho_mix = (*DensityMix_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& M_mix = (*MomentumMix_mf[lev]).array(mfi);

        amrex::Array4<const Set::Scalar> const& etaE = (*etaEnergy_old_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& etarho = (*etaDensity_old_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& etaM = (*etaMomentum_old_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& eta = (*eta_old_mf[lev]).array(mfi);

        amrex::Array4<Set::Scalar> const& etaE_new = (*etaEnergy_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& etarho_new = (*etaDensity_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& etaM_new = (*etaMomentum_mf[lev]).array(mfi);

        amrex::Array4<Set::Scalar> const& eta_new = (*eta_mf[lev]).array(mfi);

        amrex::Array4<const Set::Scalar> const& etadot = (*etadot_mf[lev]).array(mfi);

        amrex::Array4<Set::Scalar> const& v = (*Velocity_mf[lev]).array(mfi);

	amrex::Array4<Set::Scalar> const& p = (*Pressure_mf[lev]).array(mfi);

        amrex::Array4<Set::Scalar> const& omega = (*Vorticity_mf[lev]).array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
	    Set::Vector grad_eta = Numeric::Gradient(eta, i, j, k, 0, DX);
            Set::Scalar grad_eta_mag = grad_eta.lpNorm<2>();

	    //Compute New Primitive Variables
	    v(i, j, k, 0)  = M_mix(i, j, k, 0)/rho_mix(i, j, k);
	    v(i, j, k, 1)  = M_mix(i, j, k, 1)/rho_mix(i, j, k);
	    p(i, j, k)     = (E_mix(i, j, k) - 0.5 * rho_mix(i, j, k) * (v(i, j, k, 0) * v(i, j, k, 0) + v(i, j, k, 1) * v(i, j, k, 1))) * (gamma - 1);
	    
            //Godunov flux
            Solver::Local::Riemann::Roe::State state_x(rho_mix(i, j, k), M_mix(i, j, k, 0), M_mix(i, j, k, 1), E_mix(i, j, k), eta(i, j, k));
            Solver::Local::Riemann::Roe::State state_y(rho_mix(i, j, k), M_mix(i, j, k, 1), M_mix(i, j, k, 0), E_mix(i, j, k), eta(i, j, k));
            Solver::Local::Riemann::Roe::State lo_statex(rho_mix(i - 1, j, k), M_mix(i - 1, j, k, 0), M_mix(i - 1, j, k, 1), E_mix(i - 1, j, k), eta(i - 1, j, k));
            Solver::Local::Riemann::Roe::State hi_statex(rho_mix(i + 1, j, k), M_mix(i + 1, j, k, 0), M_mix(i + 1, j, k, 1), E_mix(i + 1, j, k), eta(i + 1, j, k));
            Solver::Local::Riemann::Roe::State lo_statey(rho_mix(i, j - 1, k), M_mix(i, j - 1, k, 1), M_mix(i, j - 1, k, 0), E_mix(i, j - 1, k), eta(i, j - 1, k));
            Solver::Local::Riemann::Roe::State hi_statey(rho_mix(i, j + 1, k), M_mix(i, j + 1, k, 1), M_mix(i, j + 1, k, 0), E_mix(i, j + 1, k), eta(i, j + 1, k));

            Solver::Local::Riemann::Roe::Flux flux_xlo, flux_ylo, flux_xhi, flux_yhi;

            //lo interface fluxes
            flux_xlo = Solver::Local::Riemann::Roe::Solve(lo_statex, state_x, gamma, eta(i, j, k));
            flux_ylo = Solver::Local::Riemann::Roe::Solve(lo_statey, state_y, gamma, eta(i, j, k));

            //hi interface fluxes
            flux_xhi = Solver::Local::Riemann::Roe::Solve(state_x, hi_statex, gamma, eta(i, j, k));
            flux_yhi = Solver::Local::Riemann::Roe::Solve(state_y, hi_statey, gamma, eta(i, j, k));

            //Godunov fluxes
            etaE_new(i, j, k) =
	      etaE(i, j, k)
	      + (flux_xlo.energy - flux_xhi.energy) * dt / DX[0]
	      + (flux_ylo.energy - flux_yhi.energy) * dt / DX[1];

            etarho_new(i, j, k) =
	      etarho(i, j, k)
	      + (flux_xlo.mass - flux_xhi.mass) * dt / DX[0]
	      + (flux_ylo.mass - flux_yhi.mass) * dt / DX[1];

            etaM_new(i, j, k, 0) =
	      etaM(i, j, k, 0)
	      + (flux_xlo.momentum_normal - flux_xhi.momentum_normal) * dt / DX[0]
	      + (flux_ylo.momentum_tangent - flux_yhi.momentum_tangent) * dt / DX[1];

            etaM_new(i, j, k, 1) =
	      etaM(i, j, k, 1)
	      + (flux_xlo.momentum_tangent - flux_xhi.momentum_tangent) * dt / DX[0]
	      + (flux_ylo.momentum_normal - flux_yhi.momentum_normal) * dt / DX[1];

	    //Viscous Terms
            Set::Scalar lap_ux = Numeric::Laplacian(v, i, j, k, 0, DX);
            Set::Scalar lap_uy = Numeric::Laplacian(v, i, j, k, 1, DX);
            Set::Vector grad_ux = Numeric::Gradient(v, i, j, k, 0, DX);
            Set::Vector grad_uy = Numeric::Gradient(v, i, j, k, 1, DX);
            //Set::Scalar div_u      = grad_ux(0) + grad_uy(1);
            //Set::Scalar symgrad_u  = grad_ux(1) + grad_uy(0);
            //Set::Matrix hess_u = Numeric::Hessian(v, i, j, k, 0, DX);

            etaM_new(i, j, k, 0) += mu * dt * eta(i, j, k) * lap_ux;// + mu * hess_u(0)/3.);
            etaM_new(i, j, k, 1) += mu * dt * eta(i, j, k) * lap_uy;// + mu * hess_u(1)/3.);
            //etaE_new(i, j, k)    += 2. * mu * (div_u * div_u + div_u * symgrad_u) - 2./3. * mu * div_u * div_u;

	    //Compute New Mixed Fields
            rho_mix(i, j, k)  = etarho_new(i, j, k) + (1.0 - eta_new(i, j, k)) * rho_solid;
            E_mix(i, j, k)    = etaE_new(i, j, k)   + (1.0 - eta_new(i, j, k)) * E_solid;
            M_mix(i, j, k, 0) = etaM_new(i, j, k, 0);
            M_mix(i, j, k, 1) = etaM_new(i, j, k, 1);

	    //Diffuse Sources
	    std::array<Set::Scalar, 4> source;
            source[0] = rho_interface * (v_x_interface * grad_eta(0) + v_y_interface * grad_eta(1));
            source[1] = (rho_interface * v_x_interface * v_x_interface) * grad_eta_mag;
            source[2] = (rho_interface * v_y_interface * v_y_interface) * grad_eta_mag;
            source[3] = 0.5 * rho_interface * (v_x_interface * v_x_interface * v_x_interface * grad_eta(0) + v_y_interface * v_y_interface * v_y_interface * grad_eta(1));

	    E_mix(i, j, k)    += source[3] * dt + E_mix(i, j, k) * etadot(i, j, k) * dt;
            rho_mix(i, j, k)  += source[0] * dt + rho_mix(i, j, k) * etadot(i, j, k) * dt;
            M_mix(i, j, k, 0) += source[1] * dt + M_mix(i, j, k, 0) * etadot(i, j, k) * dt;
            M_mix(i, j, k, 1) += source[2] * dt + M_mix(i, j, k, 1) * etadot(i, j, k) * dt;

	    //Compute fluid fields
	    etaE_new(i, j, k)    = eta(i, j, k) * E_mix(i, j, k);
            etarho_new(i, j, k)  = eta(i, j, k) * rho_mix(i, j, k);
            etaM_new(i, j, k, 0) = eta(i, j, k) * M_mix(i, j, k, 0);
            etaM_new(i, j, k, 1) = eta(i, j, k) * M_mix(i, j, k, 1);

	    omega(i, j, k) = (grad_uy(0) - grad_ux(1));
	    
        });
    }
}//end Advance

void Hydro::Regrid(int lev, Set::Scalar /* time */)
{
    BL_PROFILE("Integrator::Hydro::Regrid");
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
//      Density_mf[lev] -> FillBoundary();
//      Energy_mf[lev] -> FillBoundary();
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
