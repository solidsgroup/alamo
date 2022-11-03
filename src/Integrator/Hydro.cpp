#include "Hydro.H"
#include "IO/ParmParse.H"
#include "BC/Constant.H"
#include "Numeric/Stencil.H"
#include "IC/Laminate.H"
#include "IC/PSRead.H"
#include "IC/Expression.H"

namespace Integrator
{
    Hydro::Hydro(IO::ParmParse &pp) : Hydro() {pp.queryclass(*this);}
  
    void Hydro::Parse(Hydro &Value, IO::ParmParse &pp)
    {
      BL_PROFILE("Integrator::Hydro::Hydro()");
      //General Variables Input Read:
      {
	pp.query("r_refinement_criterion", value.r_refinement_criterion);
	pp.query("e_refinement_criterion", value.e_refinement_criterion);
	pp.query("m_refinement_criterion", value.m_refinement_criterion);
	pp.query("gamma", value.gamma);
	pp.query("cfl", value.cfl);

      }
      // Register FabFields:
      {
	value.RegisterNewFab(value.eta_mf, value.bc_eta, 1, 2, "eta", true);
	value.RegisterNewFab(value.eta_old_mf, value.bc_eta, 1, 2, "eta_old", false);

	value.RegisterNewFab(value.Density_mf, value.bc_rho, 1, 2, "Density", true);
	value.RegisterNewFab(value.Density_old_mf, 1, 2, "rho_old", false);

	value.RegisterNewFab(value.Energy_mf, value.bc_E, 1, 2, "Energy", true);
	value.RegisterNewFab(value.Energy_old_mf, value.bc_E, 1, 2, "E_old", false);

	value.RegisterNewFab(value.Momentum_mf, value.bc_M, 1, 2, "Momentum", true);
	value.RegisterNewFab(value.Momentum_old_mf, value.bc_M, 1, 2, "M_old", false);

	value.RegisterNewFab(value.Velocity_mf, value.bc_M, 1, 2, "Velocity", true);
	value.RegisterNewFab(value.Pressure_mf, value.bc_M, 1, 2, "Pressure", true);

      }
    }


    void Hydro::Initialize(int lev)
    {
      BL_PROFILE("Integrator::Hydro::Initialize");

      ic_eta -> Initialize(lev, eta_mf);
      ic_eta -> Initialize(lev, eta_old_mf);

      Energy_mf[lev] -> setVal(0.0);
      Energy_old_mf[lev] -> setVal(0.0);
      
      Density_mf -> setVal(0.0);
      Density_old_mf -> setVal(0.0);

      Momentum_mf -> setVal(0.0);
      Momentum_old_mf -> setVal(0.0);

      flux_x_mf -> setVal(0.0);
      flux_y_mf -> setVal(0.0);
      
      c_max = 0.0;
      vx_max = 0.0;
      vy_max = 0.0;
    }

    void Hydro::TimeStepBegin(Set::Scalar a_time, int a_iter)
    {
      BL_PROFILE("Integrator::Hydro::TimeStepBegin");
    }

  
    void Hydro::TimeStepEnd(Set::Scalar a_time, int a_iter)
    {
      // Syncronize c_max between processors so that they all have the same minimum value
      amrex::ParallelDescriptor::ReduceRealMax(c_max);
      amrex::ParallelDescriptor::ReduceRealMax(vx_max);
      amrex::ParallelDescriptor::ReduceRealMax(vy_max);

      Set::Scalar new_timestep = cfl / ((c_max + vx_max) / DX[0] + (c_max + vy_max) / DX[1]);
      
      SetTimestep(new_timestep);
    }


    void Hydro::Advance(int lev, Set::Scalar time, Set::Scalar dt)
    {

      std::swap(eta_old_mf[lev], eta_mf[lev]);
      std::swap(Momentum_old_mf[lev], Momentum_mf[lev]);
      std::swap(Energy_old_mf[lev], Energy_mf[lev]);
      std::swap(Density_old_mf[lev], Density_mf[lev]);

      const Set::Scalar *DX = geom[lev].CellSize();

      Set::Scalar rho0, V = 0.5, 1.0 //Diffuse Interface Parameters

      for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi)
      {
	const amrex::Box &bx = mfi.tilebox();

	amrex::Array4<Set::Scalar> const &eta = (*eta_mf[lev]).array(mfi);
	amrex::Array4<const Set::Scalar> const &etaold = (*eta_old_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &E = (*Energy_mf[lev]).array(mfi);
	amrex::Array4<const Set::Scalar> const &Eold = (*Energy_old_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &rho = (*Densitiy_mf[lev]).array(mfi);
	amrex::Array4<const Set::Scalar> const &rhoold = (*Density_old_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &M = (*Momentum_mf[lev]).array(mfi);
	amrex::Array4<const Set::Scalar> const &Mold = (*Momentum_old_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &v = (*Velocity_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &p = (*Pressure_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &flux_x = (*flux_x_mf[lev]).array(mfi);
	amrex::Array4<Set::Scalar> const &flux_y = (*flux_y_mf[lev]).array(mfi);

	
        //Computes Velocity and Pressure over the domain
 
	amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k){
	    v(i, j, k, 0) = M(i, j, k, 0) / rho(i,j,k);
	    v(i, j, k, 1) = M(i, j, k, 1) / rho(i,j,k);
	    v(i, j, k, 2) = M(i, j, k, 2) / rho(i,j,k);

	    Set::Scalar ke = V(i, j, k, 0) * v(i, j, k, 0) + v(i, j, k, 1) * v(i, j, k, 1);

	    p(i,j,k) = (gamma - 1.0) * rho(i,j,k) * (E(i,j,k) / rho(i,j,k) - 0.5 * ke * ke);

	    Set::Scalar c = sqrt(gamma * p(i,j,k) / rho(i,j,k) );

	    if (c > c_max){ c_max = c;}
	    if (v(i, j, k, 0) > vx_max) {vx_max = V(i, j, k, 0);}
	    if (v(i, j, k, 1) > vy_max) {vy_max = V(i, j, k, 1);}
 	});

	
	//this loop will be running the godnov solver over the space
	amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
	{
	  Set::Scalar drhoX, drhoY, dvxX, dvyY, dpX, dpY;
	  Set::Scalar drhoX_n, drhoY_n, dvxX_n, dvyY_n, dpX_n, dpY_n;
	  Set::Scalar rho_slope_right, vx_slope_right, vy_slope_right, p_slope_right;
	  Set::Scalar rho_right, vx_right, vy_right, p_right;
	  Set::Scalar rho_slope_left, vx_slope_left, vy_slope_left, p_slope_left;
	  Set::Scalar rho_left, vx_left, vy_left, p_left;

	  double flux_x[4];
          double flux_y[4];

	  // compute slopes in current cell
          drhoX = Numeric::CSolver(rho, i, j, k, 0, DX[0]);
          drhoY = Numeric::CSolver(rho, i, j, k, 1, DX[1]);

	  dvxX	= Numeric::CSolver(v, i, j, k, 0, 0, DX[0]);
	  dvxY	= Numeric::CSolver(v, i, j, k, 0, 1, DX[1]);

	  dvyX	= Numeric::CSolver(v, i, j, k, 1, 0, DX[0]);
	  dvyY	= Numeric::CSolver(v, i, j, k, 1, 1, DX[1]);

          dpX	= Numeric::CSolver(p, i, j, k, 0, DX[0]);
	  dpY	= Numeric::CSolver(p, i, j, k, 1, DX[1]);
    
	  //
          // left interface along x direction
          //

          // compute slopes in left neighbor
	  drhoX_n = Numeric::LSolver(rho, i, j, k, 0, DX[0]);
          drhoY_n = Numeric::LSolver(rho, i, j, k, 0, DX[1]);

	  dvxX_n  = Numeric::LSolver(v, i, j, k, 0, 0, DX[0]);
	  dvxY_n  = Numeric::LSolver(v, i, j, k, 0, 0, DX[1]);

	  dvyX_n  = Numeric::LSolver(v, i, j, k, 1, 0, DX[0]);
	  dvyY_n  = Numeric::LSolver(v, i, j, k, 1, 0, DX[1]);

          dpX_n	  = Numeric::LSolver(p, i, j, k, 0, DX[0]);
	  dpY_n	  = Numeric::LSolver(p, i, j, k, 0, DX[1]);

	  // slopes in current cell
	  rho_slope_right = (-v(i,j,k,0) * drhoX - dvxX * rho(i,j,k)) * dt / DX[0]          + (-v(i,j,k,1) * drhoY - dvyY * rho(i,j,k)) * dt / DX[1];
	  vx_slope_right  = (-v(i,j,k,0) * dvxX - dpX / rho(i,j,k)) * dt / DX[0]            + (-v(i,j,k,1) * dvyY) * dt / DX[1];
	  vy_slope_right  = (-v(i,j,k,0) * dvyX) * dt / DX[0]                               + (-v(i,j,k,1) * dvyY - dpY / rho(i,j,k)) * dt / DX[1];
	  p_slope_right	  = (-v(i,j,k,0) * dpX - dvxX * gamma * p(i,j,k)) * dt / DX[0] + (-v(i,j,k,1) * dpY - dvyY * gamma * p(i,j,k)) * dt / DX[1];
                    
	  // compute reconstructed states at left interface along x in current cell
	  // left interface: right state
	  rho_right = rho(i,j,k) + 0.5 * rho_slope_right - drhoX;
	  vx_right  = v(i,j,k,0) + 0.5 * vx_slope_right - dvxX;
	  vy_right  = v(i,j,k,1) + 0.5 * vy_slope_right - dvyX;
	  p_right   = p(i,j,k) + 0.5 * p_slope_right - dpX;

	  // left interface: left state
	  rho_slope_left = (-v(i-1,j,k,0) * drhoX - dvxX * rho(i-1,j,k)) * dt / DX[0]          + (-v(i-1,j,k,1) * drhoY - dvyY * rho(i-1,j,k)) * dt / DX[1];
	  vx_slope_left  = (-v(i-1,j,k,0) * dvxX - dpX / rho(i-1,j,k)) * dt / DX[0]            + (-v(i-1,j,k,1) * dvyY) * dt / DX[1];
	  vy_slope_left  = (-v(i-1,j,k,0) * dvyX) * dt / DX[0]                                 + (-v(i-1,j,k,1) * dvyY - dpY / rho(i-1,j,k)) * dt / DX[1];
	  p_slope_left	 = (-v(i-1,j,k,0) * dpX - dvxX * gamma * p(i-1,j,k)) * dt / DX[0] + (-v(i-1,j,k,1) * dpY - dvyY * gamma * p(i-1,j,k)) * dt / DX[1];

	  rho_left = rho(i-1,j,k) + 0.5 * rho_slope_left + drhoX_n;
	  vx_left  = v(i-1,j,k,0) + 0.5 * vx_slope_left + dvxX_n;
	  vy_left  = v(i-1,j,k,1) + 0.5 * vy_slope_left + dvyX_n;
	  p_left   = p(i-1,j,k) + 0.5 * p_slope_left + dpX_n;

	  double left_state[4]  = {rho_left, vx_left, vy_left, p_left, eta(i-1, j, k)};
	  double right_state[4] = {rho_right, vx_right, vy_right, p_right, eta(i, j, k)};

	  flux_x = Solver::Local::Rieman_ROE(left_state, right_state, eta(i, j, k));

	  //
	  // left interface along y direction
	  //

	  // compute slopes in left neighbor
          drhoX_n = Numeric::LSolver(rho, i, j, k, 1, DX[0]);
          drhoY_n = Numeric::LSolver(rho, i, j, k, 1, DX[1]);

	  dvxX_n  = Numeric::LSolver(v, i, j, k, 0, 1, DX[0]);
	  dvxY_n  = Numeric::LSolver(v, i, j, k, 0, 1, DX[1]);

	  dvyX_n  = Numeric::LSolver(v, i, j, k, 1, 1, DX[0]);
	  dvyY_n  = Numeric::LSolver(v, i, j, k, 1, 1, DX[1]);

          dpX_n	  = Numeric::LSolver(p, i, j, k, 1, DX[0]);
	  dpY_n	  = Numeric::LSolver(p, i, j, k, 1, DX[1]);

	  // compute reconstructed states at left interface along y in current cell
          // left interface: right state
	  rho_right = rho(i,j,k) + 0.5 * rho_slope_right - drhoY;
	  vx_right  = v(i,j,k,0) + 0.5 * vx_slope_right - dvxY;
	  vy_right  = v(i,j,k,1) + 0.5 * vy_slope_right - dvyY;
	  p_right   = p(i,j,k) + 0.5 * p_slope_right - dpY;
	    
	  // left interface: left state
	  rho_slope_left = (-v(i,j-1,k,0) * drhoX - dvxX * rho(i,j-1,k)) * dt / DX[0]          + (-v(i,j-1,k,1) * drhoY - dvyY * rho(i,j-1,k)) * dt / DX[1];
	  vx_slope_left  = (-v(i,j-1,k,0) * dvxX - dpX / rho(i,j-1,k)) * dt / DX[0]            + (-v(i,j-1,k,1) * dvyY) * dt / DX[1];
	  vy_slope_left  = (-v(i,j-1,k,0) * dvyX) * dt / DX[0]                                 + (-v(i,j-1,k,1) * dvyY - dpY / rho(i,j-1,k)) * dt / DX[1];
	  p_slope_left	 = (-v(i,j-1,k,0) * dpX - dvxX * gamma * p(i,j-1,k)) * dt / DX[0] + (-v(i,j-1,k,1) * dpY - dvyY * gamma * p(i,j-1,k)) * dt / DX[1];

	  rho_left = rho(i,j-1,k) + 0.5 * rho_slope_left + drhoY_n;
	  vx_left  = v(i,j-1,k,0) + 0.5 * vx_slope_left + dvxY_n;
	  vy_left  = v(i,j-1,k,1) + 0.5 * vy_slope_left + dvyY_n;
	  p_left   = p(i,j-1,k) + 0.5 * p_slope_left + dpY_n;
                
	  // x, y permutations
	  std::swap(vx_left, vy_left);
	  std::swap(vx_right, vy_right);

	  double left_state[4]  = {rho_left, vx_left, vy_left, p_left, eta(i, j-1, k)};
	  double right_state[4] = {rho_right, vx_right, vy_right, p_right, eta(i, j, k)};
                
	  flux_y = Solver::Local::Rieman_ROE(left_state, right_state, eta(i, j, k));
                
	  // swap flux_y components
	  std::swap(flux_y[2], flux_y[3]);
                
	  // update hydro array
	  E(i-1,j,k)   += -flux_x(i, j, k, 0) *	dtdx;
	  E(i,j,k)     += flux_x(i, j, k, 0) * dtdx;
	  E(i,j-1,k)   += -flux_y(i, j, k, 0) * dtdy;
	  E(i,j,k)     += flux_y(i, j, k, 0) * dtdy;
	  
	  rho(i-1,j,k) += -flux_x(i, j, k, 1) * dtdx;
	  rho(i,j,k)   += flux_x(i, j, k, 1) * dtdx;
	  rho(i,j-1,k) += -flux_y(i, j, k, 1) *	dtdy;
	  rho(i,j,k)   += flux_y(i, j, k, 1) * dtdy;
	  
          M(i-1,j,k,0) += -flux_x(i, j, k, 2) *	dtdx;
	  M(i,j,k,0)   += flux_x(i, j, k, 2) *	dtdx;
	  M(i,j-1,k,0) += -flux_y(i, j, k, 2) *	dtdy;
	  M(i,j,k,0)   += flux_y(i, j, k, 2) *	dtdy;

	  M(i-1,j,k,1) += -flux_x(i, j, k, 3) *	dtdx;
	  M(i,j,k,1)   += flux_x(i, j, k, 3) *	dtdx;
	  M(i,j-1,k,1) += -flux_y(i, j, k, 3) *	dtdy;
	  M(i,j,k,1)   += flux_y(i, j, k, 3) *	dtdy;

	  /// Diffuse Interface Source Terms
	  Set::Vector grad_eta = Numeric::Gradient(eta, i, j, k, 0, DX);
	  
	  Set::Matrix source = Set::Matrix::Zero();
	  source[0] = grad_eta * (rho0 * V);
	  source[1] = grad_eta * (rho0 * V**2) 
	  source[2] = grad_eta * (0.5 * rho0 * V**3)
	  
	});
      }      
    }//end Advance

  void Hydro::Regrid(int lev, Set::Scalar /* time */)
  {
    BL_PROFILE("Integrator::Hydro::Regrid");
    if (lev < finest_level) return;

    Util::Message(INFO, "Regridding on level", lev);
  }//end regrid

  void Hydro::TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, Set::Scalar time, int ngrow)
  {
    BL_PROFILE("Integrator::Hydro::TagCellsForRefinement");
    Base::Mechanics<Model::Solid::Affine::Isotropic>::TagCellsForRefinement(lev, a_tags, time, ngrow);

    const Set::Scalar *DX = geom[lev].CellSize();
    Set::Scalar dr = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));

    //Eta Criterion
    for (amrex::MFIter mfi(*eta_mf[lev], true), mfi.isValid(), ++mfi)
    {
      const amrex::Box &bx = mfi.tilebox();
      amrex::Array4<char> const &tags = a_tags.array(mfi);
      amrex::Array4<const Set::Scalar> const &eta = (*eta_mf[lev]).array(mfi);

      amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k){
	Set::Vector grad_eta = Numeric::Gradient(eta, i, j, k, 0, DX);
	if(grad_eta.lpnorm<2>() * dr * 2 > n_refinement_criterion) tags(i,j,k) = amrex::TagBox::SET;
     });
    }
    
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


  }//end TagCells

  void Hydro::Integrate(int amrlev, Set::Scalar /*time*/, int /*step*/, const amrex::MFIter &mfi, const amrex::Box &box)
  {
    BL_PROFILE("Hydro::Integrate");
    const Set::Scalar *DX = geom[amrlev].CellSize();
    Set::Scalar dv = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);
    amrex::Array4<amrex::Real> const &eta = (*eta_mf[amrlev]).array(mfi);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k){
      volume += eta(i, j, k, 0) * dv;
      Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX);
      Set::Scalar normgrad = grad.lpNorm<2>();
      Set::Scalar da = normgrad * dv;
      area += da;
    });
      
  
  }//end Integrate

  void Hydro::UpdateModel(int /*a_step*/)
  {
    for (int lev = 0; lev <= finest_level; ++lev)
    {
      eta_mf[lev] -> FillBoundary();
      Density_mf[lev] -> FillBoundary();
      Energy_mf[lev] -> FillBoundary();
      Momentum[lev] -> FillBoundary();
      
      for (MFIter mfi(*model_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
      {
	amrex::Box bx = mfi.nodaltilebox();
	amrex::Array4<model_type> const &model = model_mf[lev]->array(mfi);

	amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k){
	  // TODO

	});


      } // end For2
      
      Util::RealFillBoundary(*model_mf[lev], geom[lev]);
      amrex::MultiFab::Copy(*psi_mf[lev], *eta_mf[lev], 0, 0, 1, psi_mf[lev]-> nGrow());
      
    } //end For1
  }//end update

  
}//end code

