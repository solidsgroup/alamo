#include "Hydro.H"
#include "IO/ParmParse.H"
#include "BC/Constant.H"
#include "Numeric/Stencil.H"
#include "IC/Laminate.H"
#include "IC/PSRead.H"
#include "IC/Expression.H"
#include "Solver/Local/Riemann/Roe.H"
#include "Solver/Local/Riemann/HLLC.H"

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
	pp.query("mu", value.cfl);

        pp.query("rho_solid", value.rho_solid);
        pp.query("rho_fluid", value.rho_fluid);
        pp.query("E_solid", value.E_solid);
        pp.query("E_fluid", value.E_fluid);

        pp.query("Mx_init", value.Mx_init);
        pp.query("My_init", value.My_init);

        pp.query("eps", value.eps);

        pp.query("mdot", value.mdot);
        pp.query("Pdot_x", value.Pdot_x);
        pp.query("Pdot_y", value.Pdot_y);
        pp.query("Qdot", value.Qdot);
	pp.query("Ldot_x", value.Ldot_x);
        pp.query("Ldot_y", value.Ldot_y);

        value.bc_eta = new BC::Constant(1, pp, "pf.eta.bc");
        value.bc_rho = new BC::Constant(1, pp, "rho.bc");
        value.bc_E = new BC::Constant(1, pp, "E.bc");
        value.bc_M = new BC::Constant(2, pp, "M.bc");
    }
    // Register FabFields:
    {
        int nghost = 2;
        value.RegisterNewFab(value.eta_mf, value.bc_eta, 1, nghost, "eta", true);

        value.RegisterNewFab(value.etadot_mf, value.bc_eta, 1, nghost, "etadot", true);

        value.RegisterNewFab(value.Density_mf, value.bc_rho, 1, nghost, "Density", true);
        value.RegisterNewFab(value.Density_old_mf, value.bc_rho, 1, nghost, "rho_old", false);

        value.RegisterNewFab(value.Energy_mf, value.bc_E, 1, nghost, "Energy", true);
        value.RegisterNewFab(value.Energy_old_mf, value.bc_E, 1, nghost, "E_old", false);

        value.RegisterNewFab(value.Momentum_mf, value.bc_M, 2, nghost, "Momentum", true);
        value.RegisterNewFab(value.Momentum_old_mf, value.bc_M, 2, nghost, "M_old", false);

        value.RegisterNewFab(value.Velocity_mf, value.bc_M, 2, nghost, "Velocity", true);

        value.RegisterNewFab(value.Pressure_mf, value.bc_E, 1, nghost, "Pressure", true);
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
}


void Hydro::Initialize(int lev)
{
    BL_PROFILE("Integrator::Hydro::Initialize");

    UpdateEta(0.0);

    for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi) {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<const Set::Scalar> const& eta = (*eta_mf[lev]).array(mfi);

        amrex::Array4<Set::Scalar> const& E_new = (*Energy_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& E = (*Energy_old_mf[lev]).array(mfi);

        amrex::Array4<Set::Scalar> const& rho_new = (*Density_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& rho = (*Density_old_mf[lev]).array(mfi);

        amrex::Array4<Set::Scalar> const& M_new = (*Momentum_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& M = (*Momentum_old_mf[lev]).array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
	    const Set::Scalar* DX = geom[lev].CellSize();
            Set::Vector grad_eta = Numeric::Gradient(eta, i, j, k, 0, DX);
            Set::Scalar grad_eta_mag = grad_eta.lpNorm<2>() * eps;
	    
            rho(i, j, k) = rho_solid * (1.0 - eta(i, j, k)) + rho_fluid * eta(i, j, k) + mdot * grad_eta_mag;
            rho_new(i, j, k) = rho(i, j, k);

            M(i, j, k, 0) = 0.0 * (1.0 - eta(i, j, k)) + Mx_init * eta(i, j, k) + Pdot_x * grad_eta[0];
            M_new(i, j, k, 0) = M(i, j, k, 0);
            ///
            M(i, j, k, 1) = 0.0 * (1.0 - eta(i, j, k)) + My_init * eta(i, j, k) + Pdot_y * grad_eta[1];
            M_new(i, j, k, 1) = M(i, j, k, 1);

            E(i, j, k) = E_solid * (1.0 - eta(i, j, k)) + E_fluid * eta(i, j, k) + Qdot * grad_eta_mag;
            E_new(i, j, k) = E(i, j, k);

        });
    }

    c_max = 0.0;
    grad_eta_max = 1.0;
    vx_max = 0.0;
    vy_max = 0.0;
}

void Hydro::UpdateEta(Set::Scalar time)
{
    for (int lev = 0; lev <= finest_level; ++lev)
    {
        ic_eta->Initialize(lev, eta_mf, time);
        ic_etadot->Initialize(lev, etadot_mf, time);
    }
}

void Hydro::TimeStepBegin(Set::Scalar time, int /*iter*/)
{
    UpdateEta(time);
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



void Hydro::Advance(int lev, Set::Scalar, Set::Scalar dt)
{
    std::swap(Momentum_old_mf[lev], Momentum_mf[lev]);
    ///
    std::swap(Energy_old_mf[lev], Energy_mf[lev]);
    ///
    std::swap(Density_old_mf[lev], Density_mf[lev]);

    const Set::Scalar* DX = geom[lev].CellSize();

    for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox();
	
        amrex::Array4<const Set::Scalar> const& E = (*Energy_old_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& rho = (*Density_old_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& M = (*Momentum_old_mf[lev]).array(mfi);

        amrex::Array4<Set::Scalar> const& v = (*Velocity_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& p = (*Pressure_mf[lev]).array(mfi);

        //Compute primitive variables

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            v(i, j, k, 0) = M(i, j, k, 0) / rho(i, j, k);
            v(i, j, k, 1) = 0.0; //M(i, j, k, 1) / rho(i, j, k);

            Set::Scalar ke = 0.5 * rho(i, j, k) * (v(i, j, k, 0) * v(i, j, k, 0) + v(i, j, k, 1) * v(i, j, k, 1));

            p(i, j, k, 0) = (gamma - 1.0) * (E(i, j, k) - ke);

            // Set::Scalar c = sqrt(gamma * p(i, j, k) / rho(i, j, k));

            // if (c > c_max) { c_max = c; }
            // if (v(i, j, k, 0) > vx_max) { vx_max = v(i, j, k, 0); }
            // if (v(i, j, k, 1) > vy_max) { vy_max = v(i, j, k, 1); }
        });
    }
    Pressure_mf[lev]->FillBoundary();
    Velocity_mf[lev]->FillBoundary();

    for (amrex::MFIter mfi(*eta_mf[lev], false); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();

        amrex::Array4<const Set::Scalar> const& eta = (*eta_mf[lev]).array(mfi);

        amrex::Array4<const Set::Scalar> const& E = (*Energy_old_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& rho = (*Density_old_mf[lev]).array(mfi);
        amrex::Array4<const Set::Scalar> const& M = (*Momentum_old_mf[lev]).array(mfi);

        amrex::Array4<Set::Scalar> const& E_new = (*Energy_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& rho_new = (*Density_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& M_new = (*Momentum_mf[lev]).array(mfi);

	amrex::Array4<Set::Scalar> const& v = (*Velocity_mf[lev]).array(mfi);
        amrex::Array4<Set::Scalar> const& p = (*Pressure_mf[lev]).array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
	    //Godunov flux
            Solver::Local::Riemann::HLLC::State state(rho(i, j, k), M(i, j, k, 0), M(i, j, k, 1), E(i, j, k), eta(i, j, k));
            Solver::Local::Riemann::HLLC::State lo_statex(rho(i - 1, j, k), M(i - 1, j, k, 0), M(i - 1, j, k, 1), E(i - 1, j, k), eta(i - 1, j, k) );
            Solver::Local::Riemann::HLLC::State hi_statex(rho(i + 1, j, k), M(i + 1, j, k, 0), M(i + 1, j, k, 1), E(i + 1, j, k), eta(i + 1, j, k) );
            Solver::Local::Riemann::HLLC::State lo_statey(rho(i, j - 1, k), M(i, j - 1, k, 0), M(i, j - 1, k, 1), E(i, j - 1, k), eta(i, j - 1, k) ); 
            Solver::Local::Riemann::HLLC::State hi_statey(rho(i, j + 1, k), M(i, j + 1, k, 0), M(i, j + 1, k, 1), E(i, j + 1, k), eta(i, j + 1, k) );
	    
            Solver::Local::Riemann::HLLC::Flux flux_xlo, flux_ylo, flux_xhi, flux_yhi, flux_test;

            //lo interface fluxes
            flux_xlo = Solver::Local::Riemann::HLLC::Solve(lo_statex, state, gamma);
            flux_ylo = Solver::Local::Riemann::HLLC::Solve(lo_statey, state, gamma);

            //hi interface fluxes
            flux_xhi = Solver::Local::Riemann::HLLC::Solve(state, hi_statex, gamma);
            flux_yhi = Solver::Local::Riemann::HLLC::Solve(state, hi_statey, gamma);


            //Godunov fluxes
            E_new(i, j, k) =
                E(i, j, k)
	      + (flux_xlo.energy - flux_xhi.energy) * dt / DX[0]
	      + (flux_ylo.energy - flux_yhi.energy) * dt / DX[1];

            rho_new(i, j, k) =
                rho(i, j, k)
	      + (flux_xlo.mass - flux_xhi.mass) * dt / DX[0]
	      + (flux_ylo.mass - flux_yhi.mass) * dt / DX[1];

            M_new(i, j, k, 0) =
                M(i, j, k, 0)
	      + (flux_xlo.momentum(0) - flux_xhi.momentum(0)) * dt / DX[0]
	      + (flux_ylo.momentum(1) - flux_yhi.momentum(1)) * dt / DX[1];

            M_new(i, j, k, 1) =
                M(i, j, k, 1)
	      + (flux_xlo.momentum(1) - flux_xhi.momentum(1)) * dt / DX[0]
	      + (flux_ylo.momentum(0) - flux_yhi.momentum(0)) * dt / DX[1];

	    ///////////////////////////
            ///////VISCOUS TERMS///////
            ///////////////////////////
	    
	    Set::Scalar lap_ux     = Numeric::Laplacian(v, i, j, k, 0, DX);
            Set::Scalar lap_uy     = Numeric::Laplacian(v, i, j, k, 1, DX);
	    Set::Vector grad_ux    = Numeric::Gradient(v, i, j, k, 0, DX);
	    Set::Vector grad_uy    = Numeric::Gradient(v, i, j, k, 1, DX);
	    Set::Scalar div_u      = grad_ux(0) + grad_uy(1);
	    Set::Scalar symgrad_u  = grad_ux(1) + grad_uy(0);

            Set::Matrix hess_u = Numeric::Hessian(v, i, j, k, 0, DX);

            M_new(i, j, k, 0) += (mu * lap_ux + mu * hess_u(0)/3.);
            M_new(i, j, k, 1) += (mu * lap_uy + mu * hess_u(1)/3.);

	    E_new(i, j, k)    += 2. * mu * (div_u * div_u + div_u * symgrad_u) - 2./3. * mu * div_u * div_u;

            ///////////////////////////
            //////DIFFUSE SOURCES//////
            ///////////////////////////

	    Set::Vector grad_eta = Numeric::Gradient(eta, i, j, k, 0, DX);
            Set::Scalar grad_eta_mag = grad_eta.lpNorm<2>() * eps;

	    Set::Matrix Ldot_outer_grad_eta;
	    Ldot_outer_grad_eta(0,0) = Ldot_x * grad_eta(0);
	    Ldot_outer_grad_eta(0,1) = Ldot_y * grad_eta(0);
	    Ldot_outer_grad_eta(1,0) = Ldot_x * grad_eta(1);
	    Ldot_outer_grad_eta(1,1) = Ldot_y * grad_eta(1);

	    Set::Vector div_Ldot_outer_grad_eta = Numeric::Divergence(Ldot_outer_grad_eta, i, j, k, DX); //not working

            std::array<Set::Scalar, 3> source;
            source[0] = mdot * grad_eta_mag;
            source[1] = Pdot_x * grad_eta[0] - div_Ldot_outer_grad_eta(0);
            source[2] = Pdot_y * grad_eta[1] - div_Ldot_outer_grad_eta(1);
            source[3] = Qdot * grad_eta_mag  - (div_Ldot_outer_grad_eta(0) * v(i, j, k, 0) + div_Ldot_outer_grad_eta(1) * v(i, j, k, 1));

            E_new(i, j, k) += source[3];
            rho_new(i, j, k) += source[0];
            M_new(i, j, k, 0) += source[1];
            M_new(i, j, k, 1) += source[2];

	    //////////////////////////
            //////SOLID STAND-IN//////
            //////////////////////////

	    E_new(i, j, k) = E_solid * (1 - eta(i,j,k)) + E_new(i, j, k) * eta(i,j,k);
            rho_new(i, j, k) = rho_solid * (1 - eta(i,j,k)) + rho_new(i, j, k) * eta(i,j,k);
            M_new(i, j, k, 0) = M_new(i, j, k, 0) * eta(i,j,k);
            M_new(i, j, k, 1) = M_new(i, j, k, 1) * eta(i,j,k);
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
