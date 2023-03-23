#include "Hydro.H"
#include "IO/ParmParse.H"
#include "BC/Constant.H"
#include "Numeric/Stencil.H"
#include "IC/Laminate.H"
#include "IC/PSRead.H"
#include "IC/Expression.H"
#include "Solver/Local/Riemann_ROE.H"

namespace Integrator
{

    //Hydro(IO::ParmParse &pp) : Hydro() {pp.queryclass(*this);}

   
    void Hydro::Initialize(int lev)
    {
      BL_PROFILE("Integrator::Hydro::Initialize");

      ic_eta -> Initialize(lev, eta_mf);
      ic_eta -> Initialize(lev, eta_old_mf);

      Density_fluid_mf[lev] -> setVal(rho_fluid);
      Density_solid_mf[lev] -> setVal(rho_solid);

      Energy_fluid_mf[lev] -> setVal(E_fluid);
      Energy_solid_mf[lev] -> setVal(E_solid);

      Momentum_solid_mf[lev] -> setVal(0.0);

      for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi){
	const amrex::Box &bx = mfi.tilebox();
	amrex::Array4<Set::Scalar> const &eta = (*eta_mf[lev]).array(mfi);
	/////
	amrex::Array4<Set::Scalar> const &E = (*Energy_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &E_old = (*Energy_old_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &Ef = (*Energy_fluid_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &Es = (*Energy_solid_mf[lev]).array(mfi);
	/////
	amrex::Array4<Set::Scalar> const &rho = (*Density_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &rho_old = (*Density_old_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &rhof = (*Density_fluid_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &rhos = (*Density_solid_mf[lev]).array(mfi);
	/////
	amrex::Array4<Set::Scalar> const &M = (*Momentum_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &M_old = (*Momentum_old_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &Mf = (*Momentum_fluid_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &Ms = (*Momentum_solid_mf[lev]).array(mfi);
	/////
	amrex::Array4<Set::Scalar> const &p = (*Pressure_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &p_old = (*Pressure_old_mf[lev]).array(mfi);
	/////
	amrex::Array4<Set::Scalar> const &v = (*Velocity_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &v_old = (*Velocity_old_mf[lev]).array(mfi);
	
	amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
	{
	  rho(i,j,k) = rhos(i,j,k) * (1 - eta(i,j,k)) + rhof(i,j,k) * eta(i,j,k);
	  rho_old(i,j,k) = rho(i,j,k);

	  Mf(i,j,k,0) = Mx_init;
	  M(i,j,k,0) = Ms(i,j,k,0) * (1 - eta(i,j,k)) + Mf(i,j,k,0) * eta(i,j,k);
	  M_old(i,j,k,0) = M(i,j,k,0);
	  ///
	  Mf(i,j,k,1) = My_init;
	  M(i,j,k,1) = Ms(i,j,k,1) * (1 - eta(i,j,k)) + Mf(i,j,k,1) * eta(i,j,k);
	  M_old(i,j,k,1) = M(i,j,k,1);
	    
	  E(i,j,k) = Es(i,j,k) * (1 - eta(i,j,k)) + Ef(i,j,k) * eta(i,j,k);
	  E_old(i,j,k) = E(i,j,k);
	  
	  p(i,j,k) = (gamma - 1) * rho(i,j,k) * E(i,j,k);
	  p_old(i,j,k) = p(i,j,k);

	  v(i,j,k,0) = M(i,j,k,0) / rho(i,j,k);
	  v_old(i,j,k,0) = v(i,j,k,0);
	  v(i,j,k,1) = M(i,j,k,1) / rho(i,j,k);
	  v_old(i,j,k,1) = v(i,j,k,1);
	});
      }
      
      c_max = 0.0;
      vx_max = 0.0;
      vy_max = 0.0;
    }

    void Hydro::TimeStepBegin(Set::Scalar , int )
    {
      BL_PROFILE("Integrator::Hydro::TimeStepBegin");
    }

  
    void Hydro::TimeStepComplete(Set::Scalar , int lev) 
    {
      const Set::Scalar *DX = geom[lev].CellSize();
      
      amrex::ParallelDescriptor::ReduceRealMax(c_max);
      amrex::ParallelDescriptor::ReduceRealMax(vx_max);
      amrex::ParallelDescriptor::ReduceRealMax(vy_max);

      Set::Scalar new_timestep = cfl / ((c_max + vx_max) / DX[0] + (c_max + vy_max) / DX[1]);
      
      Util::Assert(INFO,TEST(AMREX_SPACEDIM == 2));
      
      SetTimestep(new_timestep);
    }

    

    void Hydro::Advance(int lev, Set::Scalar , Set::Scalar dt)
    {
      std::swap(eta_old_mf[lev], eta_mf[lev]);
      ///
      std::swap(Momentum_old_mf[lev], Momentum_mf[lev]);
      std::swap(Momentum_fluid_old_mf[lev], Momentum_fluid_mf[lev]);
      ///
      std::swap(Velocity_old_mf[lev], Velocity_mf[lev]);
      std::swap(Velocity_fluid_old_mf[lev], Velocity_fluid_mf[lev]);
      ///
      std::swap(Energy_old_mf[lev], Energy_mf[lev]);
      std::swap(Energy_fluid_old_mf[lev], Energy_fluid_mf[lev]);
      ///
      std::swap(Density_old_mf[lev], Density_mf[lev]);
      std::swap(Density_fluid_old_mf[lev], Density_fluid_mf[lev]);
      ///
      std::swap(Pressure_old_mf[lev], Pressure_mf[lev]);
      std::swap(Pressure_fluid_old_mf[lev], Pressure_fluid_mf[lev]);
      

      const Set::Scalar *DX = geom[lev].CellSize();

      for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi)
      {
	const amrex::Box &bx = mfi.tilebox();

	amrex::Array4<Set::Scalar> const &eta = (*eta_mf[lev]).array(mfi);
	///
	amrex::Array4<Set::Scalar> const &E = (*Energy_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &Ef = (*Energy_fluid_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &Es = (*Energy_solid_mf[lev]).array(mfi);
	///
	amrex::Array4<Set::Scalar> const &rho = (*Density_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &rhof = (*Density_fluid_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &rhos = (*Density_solid_mf[lev]).array(mfi);
	///
	amrex::Array4<Set::Scalar> const &M = (*Momentum_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &Mf = (*Momentum_fluid_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &Ms = (*Momentum_solid_mf[lev]).array(mfi);
	///
	amrex::Array4<Set::Scalar> const &v = (*Velocity_mf[lev]).array(mfi);
	///
	amrex::Array4<Set::Scalar> const &p = (*Pressure_mf[lev]).array(mfi);
	
        //Computes Velocity and Pressure over the domain
 
	amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k){
	    v(i, j, k, 0) = M(i, j, k, 0) / rho(i, j, k);
	    v(i, j, k, 1) = M(i, j, k, 1) / rho(i, j, k);

	    Set::Scalar ke = 0.5 * (v(i, j, k, 0) * v(i, j, k, 0) + v(i, j, k, 1) * v(i, j, k, 1));

	    p(i,j,k) = (gamma - 1.0) * rho(i,j,k) * (E(i,j,k) / rho(i,j,k) - ke);

	    //if (p(i,j,k) < 0) {cout << "Negative Pressure";};

	    Set::Scalar c = sqrt(gamma * p(i,j,k) / rho(i,j,k) );

	    if (c > c_max){ c_max = c;}
	    if (v(i, j, k, 0) > vx_max) {vx_max = v(i, j, k, 0);}
	    if (v(i, j, k, 1) > vy_max) {vy_max = v(i, j, k, 1);}

	    Set::Vector drho = Numeric::Gradient(rho, i, j, k, 0, DX);
	    Set::Vector dp   = Numeric::Gradient(p, i, j, k, 0, DX);
	    Set::Vector dvx  = Numeric::Gradient(v, i, j, k, 0, DX);
	    Set::Vector dvy  = Numeric::Gradient(v, i, j, k, 0, DX);

	    Set::Vector drho_n = Numeric::Gradient(rho, i-1, j, k, 0, DX);
	    Set::Vector dvx_n  = Numeric::Gradient(v, i-1, j, k, 0, DX);
	    Set::Vector dvy_n  = Numeric::Gradient(v, i-1, j, k, 1, DX);
	    Set::Vector dp_n   = Numeric::Gradient(p, i-1, j, k, 0, DX);

	    Set::Scalar rho_slope_right, vx_slope_right, vy_slope_right, p_slope_right;
	    Set::Scalar rho_right, vx_right, vy_right, p_right;
	    Set::Scalar rho_slope_left, vx_slope_left, vy_slope_left, p_slope_left;
	    Set::Scalar rho_left, vx_left, vy_left, p_left;
	    std::array<Set::Scalar,4> flux_x, flux_y;

	    // slopes in current cell
	    rho_slope_right = (-v(i,j,k,0) * drho[0] - dvx[0] * rho(i,j,k)) * dt      + (-v(i,j,k,1) * drho[1] - dvy[1] * rho(i,j,k)) * dt;
	    vx_slope_right  = (-v(i,j,k,0) * dvx[0] - dp[0] / rho(i,j,k)) * dt        + (-v(i,j,k,1) * dvy[1]) * dt / DX[1];
	    vy_slope_right  = (-v(i,j,k,0) * dvy[0]) * dt                             + (-v(i,j,k,1) * dvy[1] - dp[1] / rho(i,j,k)) * dt;
	    p_slope_right   = (-v(i,j,k,0) * dp[0] - dvx[0] * gamma * p(i,j,k)) * dt  + (-v(i,j,k,1) * dp[1] - dvy[1] * gamma * p(i,j,k)) * dt;

	    // compute reconstructed states at left interface along x in current cell
	    // left interface: right state
	    rho_right = rho(i,j,k) + 0.5 * rho_slope_right - drho[0] * DX[0];
	    vx_right  = v(i,j,k,0) + 0.5 * vx_slope_right - dvx[0] * DX[0];
	    vy_right  = v(i,j,k,1) + 0.5 * vy_slope_right - dvy[0] * DX[0];
	    p_right   = p(i,j,k)   + 0.5 * p_slope_right - dp[0] * DX[0];

	    //if (p_right < 0) {cout << "Negative Pressure";};

	    // left interface: left state
	    rho_slope_left = (-v(i-1,j,k,0) * drho[0] - dvx[0] * rho(i-1,j,k)) * dt      + (-v(i-1,j,k,1) * drho[1] - dvy[1] * rho(i-1,j,k)) * dt;
	    vx_slope_left  = (-v(i-1,j,k,0) * dvx[0] - dp[0] / rho(i-1,j,k)) * dt        + (-v(i-1,j,k,1) * dvy[1]) * dt;
	    vy_slope_left  = (-v(i-1,j,k,0) * dvy[0]) * dt                               + (-v(i-1,j,k,1) * dvy[1] - dp[1] / rho(i-1,j,k)) * dt;
	    p_slope_left   = (-v(i-1,j,k,0) * dp[0] - dvx[0] * gamma * p(i-1,j,k)) * dt  + (-v(i-1,j,k,1) * dp[1] - dvy[1] * gamma * p(i-1,j,k)) * dt;

	    rho_left = rho(i-1,j,k) + 0.5 * rho_slope_left + drho_n[0] * DX[0];
	    vx_left  = v(i-1,j,k,0) + 0.5 * vx_slope_left + dvx_n[0] * DX[0];
	    vy_left  = v(i-1,j,k,1) + 0.5 * vy_slope_left + dvy_n[0] * DX[0];
	    p_left   = p(i-1,j,k) + 0.5 * p_slope_left + dp_n[0] * DX[0];

	    if (p_left < 0) {cout << "Negative Pressure";};
	    
	    double left_state[5] = {rho_left, vx_left, vy_left, p_left, eta(i-1, j, k)};
	    double right_state[5] = {rho_right, vx_right, vy_right, p_right, eta(i, j, k)};

	    flux_x = Solver::Local::Riemann_ROE(left_state, right_state, eta(i, j, k), gamma);

	    //
	    // left interface along y direction
	    //

	    // compute slopes in left neighbor
	    drho_n = Numeric::Gradient(rho, i, j-1, k, 0, DX);
	    dvx_n  = Numeric::Gradient(v, i, j-1, k, 0, DX);
	    dvy_n  = Numeric::Gradient(v, i, j-1, k, 1, DX);
	    dp_n   = Numeric::Gradient(p, i, j-1, k, 0, DX);

	    // compute reconstructed states at left interface along y in current cell
	    // left interface: right state
	    rho_right = rho(i,j,k) + 0.5 * rho_slope_right - drho[1] * DX[1];
	    vx_right  = v(i,j,k,0) + 0.5 * vx_slope_right - dvx[1] * DX[1];
	    vy_right  = v(i,j,k,1) + 0.5 * vy_slope_right - dvy[1] * DX[1];
	    p_right   = p(i,j,k) + 0.5 * p_slope_right - dp[1] * DX[1];
	    
	    if (p_right < 0) {cout << "Negative Pressure";};
	    
	    // left interface: left state
	    rho_slope_left = (-v(i,j-1,k,0) * drho[0] - dvx[0] * rho(i,j-1,k)) * dt    + (-v(i,j-1,k,1) * drho[1] - dvy[1] * rho(i,j-1,k)) * dt;
	    vx_slope_left  = (-v(i,j-1,k,0) * dvx[0] - dp[0] / rho(i,j-1,k)) * dt      + (-v(i,j-1,k,1) * dvy[1]) * dt;
	    vy_slope_left  = (-v(i,j-1,k,0) * dvy[0]) * dt                             + (-v(i,j-1,k,1) * dvy[1] - dp[1] / rho(i,j-1,k)) * dt;
	    p_slope_left   = (-v(i,j-1,k,0) * dp[0]- dvx[0] * gamma * p(i,j-1,k)) * dt + (-v(i,j-1,k,1) * dp[1] - dvy[1] * gamma * p(i,j-1,k)) * dt;

	    rho_left = rho(i,j-1,k) + 0.5 * rho_slope_left + drho_n[1] * DX[1];
	    vx_left  = v(i,j-1,k,0) + 0.5 * vx_slope_left + dvx_n[1] * DX[1];
	    vy_left  = v(i,j-1,k,1) + 0.5 * vy_slope_left + dvy_n[1] * DX[1];
	    p_left   = p(i,j-1,k) + 0.5 * p_slope_left + dp_n[1] * DX[1];

	    if (p_left < 0) {cout << "Negative Pressure";};

	    // x, y permutations
	    std::swap(vx_left, vy_left);
	    std::swap(vx_right, vy_right);

	    left_state[0] = rho_left;
	    left_state[1] = vx_left;
	    left_state[2] = vy_left;
	    left_state[3] = p_left;
	    left_state[4] = eta(i, j-1, k);

	    right_state[0] = rho_right;
	    right_state[1] = vx_right;
	    right_state[2] = vy_right;
	    right_state[3] = p_right;
	    right_state[4] = eta(i, j-1, k);

	    flux_y = Solver::Local::Riemann_ROE(left_state, right_state, eta(i, j, k), gamma);

	    // swap flux_y components
	    std::swap(flux_y[2], flux_y[3]);

	    // update fluid values - Godunov fluxes
	    Ef(i-1,j,k)   += -flux_x[1] * dt / DX[0];
	    Ef(i,j,k)     += flux_x[1] * dt / DX[0];
	    Ef(i,j-1,k)   += -flux_y[1] * dt / DX[1];
	    Ef(i,j,k)     += flux_y[1] * dt / DX[1];
	    //Ef(i,j,k)     += source[3];
	  
	    rhof(i-1,j,k) += -flux_x[0] * dt / DX[0];
	    rhof(i,j,k)   += flux_x[0] * dt / DX[0];
	    rhof(i,j-1,k) += -flux_y[0] * dt / DX[1];
	    rhof(i,j,k)   += flux_y[0] * dt / DX[1];
	    //rhof(i,j,k)   += source[0];
	  
	    Mf(i-1,j,k,0) += -flux_x[2] * dt / DX[0];
	    Mf(i,j,k,0)   += flux_x[2] * dt / DX[0];
	    Mf(i,j-1,k,0) += -flux_y[2] * dt / DX[1];
	    Mf(i,j,k,0)   += flux_y[2] * dt / DX[1];
	    //Mf(i,j,k,0)   += source[1]; 

	    Mf(i-1,j,k,1) += -flux_x[3] * dt / DX[0];
	    Mf(i,j,k,1)   += flux_x[3] * dt / DX[0];
	    Mf(i,j-1,k,1) += -flux_y[3] * dt / DX[1];
	    Mf(i,j,k,1)   += flux_y[3] * dt / DX[1];
	    //Mf(i,j,k,1)   += source[2];

	    // update fluid values - viscous terms
	    Set::Scalar lap_ux = Numeric::Laplacian(v, i, j, k, 0, DX);
	    Set::Scalar lap_uy = Numeric::Laplacian(v, i, j, k, 1, DX);
	    
	    //Set::Scalar div_ux = Numeric::Divergence(v, i, j, k, 0, DX);
	    //Set::Scalar div_uy = Numeric::Divergence(v, i, j, k, 1, DX);
	    
	    //Set::Vector grad_div_ux = Numeric::Gradient(div_ux, i, j, k, 0, DX);
	    //Set::Vector grad_div_uy = Numeric::Gradient(div_uy, i, j, k, 0, DX);
	    
	    //M(i,j,k,0) += mu * lap_ux; // + mu * grad_div_u/3.;
	    //M(i,j,k,1) += mu * lap_uy; // + mu * grad_div_u/3.;

	    /// Diffuse Interface Source Terms
	    Set::Vector grad_eta = Numeric::Gradient(eta, i, j, k, 0, DX);
	    Set::Scalar grad_eta_mag = grad_eta.lpNorm<2>();
	    
	    // std::array<Set::Scalar,3> source;
	    // source[0] = grad_eta_mag * (rho_solid * std::sqrt(V_cross*V_cross + V_norm*V_norm));
	    // source[1] = grad_eta_mag * (rho_solid * (V_norm*V_norm));
	    // source[2] = grad_eta_mag * (rho_solid * (V_cross*V_cross));
	    // source[3] = grad_eta_mag * (0.5 * rho_solid * (V_cross*V_cross + V_norm*V_norm)*(V_cross*V_cross + V_norm*V_norm)*(V_cross*V_cross + V_norm*V_norm));

	    // update total fields
	    E(i,j,k)   = Es(i,j,k) * (1 - eta(i,j,k)) + Ef(i,j,k) * eta(i,j,k); //+ source[3];
	    rho(i,j,k) = rhos(i,j,k) * (1 - eta(i,j,k)) + rhof(i,j,k) * eta(i,j,k); //+ source[0];
	    M(i,j,k,0)   = Ms(i,j,k,0) * (1 - eta(i,j,k)) + Mf(i,j,k,0) * eta(i,j,k); //+ source[1];
	    M(i,j,k,1)   = Ms(i,j,k,1) * (1 - eta(i,j,k)) + Mf(i,j,k,1) * eta(i,j,k); //+ source[2];

	    //Advance eta

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
  void Hydro::TagCellsForRefinement(int , amrex::TagBoxArray &, Set::Scalar , int )
  {
  //   BL_PROFILE("Integrator::Flame::TagCellsForRefinement");
        
  //   const Set::Scalar *DX = geom[lev].CellSize();
  //   Set::Scalar dr = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));

  //   // Eta criterion for refinement
  //   for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi){
  //     const amrex::Box &bx = mfi.tilebox();
  //     amrex::Array4<char> const &tags = a_tags.array(mfi);
  //     amrex::Array4<const Set::Scalar> const &Eta = (*eta_mf[lev]).array(mfi);

  //     amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k){
  // 		Set::Vector gradeta = Numeric::Gradient(Eta, i, j, k, 0, DX);
  //               if (gradeta.lpNorm<2>() * dr * 2 > m_refinement_criterion) tags(i, j, k) = amrex::TagBox::SET; });
  //   }
    /*
    // Energy criterion
    for (amrex::MFIter mfi(*Energy_mf[lev], true), mfi.isValid(), ++mfi)
    {
      const amrex::Box &bx = mfi.tilebox();
      amrex::Array4<char> const &tags = a_tags.array(mfi);
      amrex::Array4<const Set::Scalar> const &E = (*Energy_mf[lev]).array(mfi);
      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k){
  	Set::Vector grad_E = Numeric::Gradient(E, i, j, k, 0, DX);
  	if(grad_E.lpNorm<2>() * dr > e_refinement_criterion) tags(i,j,k) = amrex::TagBox::SET;
  	});
    }
    
    // Density criterion
    for (amrex::MFIter mfi(*Density_mf[lev] , true), mfi.isValid(), ++mfi)
    {
      const amrex::Box &bx = mfi.tilebox();
      amrex::Array4<char> const &tags = a_tags.array(mfi);
      amrex::Array4<const Set::Scalar> const &rho = (*Density_mf[lev]).array(mfi);
      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k){
  	Set::Vector grad_rho = Numeric::Gradient(rho, i, j, k, 0, DX);
  	if (grad_rho.lpNorm<2>() * dr > r_refinment_criterion) tags(i,j,k) = amrex::TagBox::SET;
      });
    }

    // Momentum criterion
    for (amrex:::MFIter mfi(*Momentum_mf[lev], true), mfi.isValid(), ++mfi)
    {
      const amrex::Box &bx = mfi.tilebox();
      amrex::Array4<char> const &tags = a_tags.array(mfi);
      amrex::Array<const Set::Scalar> const &M = (*Momentum_mf[lev]).array(mfi);
      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k, int m){
  	Set::Vector grad_M = Numeric::Gradient(M, i, j, k, m);
  	if (grad_M.lpNorm<2>() * dr > m_refinement_criterion) tags(i,j,k) = amrex::TagBox::SET;     
      });
    }
    */


  }//end TagCells

  // void Hydro::Integrate(int amrlev, Set::Scalar /*time*/, int /*step*/, const amrex::MFIter &mfi, const amrex::Box &box)
  // {
  //   BL_PROFILE("Hydro::Integrate");
  //   const Set::Scalar *DX = geom[amrlev].CellSize();
  //   Set::Scalar dv = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);
  //   amrex::Array4<amrex::Real> const &eta = (*eta_mf[amrlev]).array(mfi);
  //   amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k){
  //   	volume += eta(i, j, k, 0) * dv;
  //   	Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX);
  //   	Set::Scalar normgrad = grad.lpNorm<2>();
  //   	Set::Scalar da = normgrad * dv;
  //   	area += da;
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
//	amrex::Box bx = mfi.nodaltilebox();
//	amrex::Array4<model_type> const &model = model_mf[lev]->array(mfi);
//
//	amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k){
//	  // TODO
//
//	});
//
//
//     } // end For2
     
//     Util::RealFillBoundary(*model_mf[lev], geom[lev]);
//     amrex::MultiFab::Copy(*psi_mf[lev], *eta_mf[lev], 0, 0, 1, psi_mf[lev]-> nGrow());
     
  //    } //end For1
  //}//end update

  
}//end code

