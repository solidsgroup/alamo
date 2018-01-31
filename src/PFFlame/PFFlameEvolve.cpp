
#include <limits>

#include <AMReX_MultiFabUtil.H>
#include <PFFlame.H>

using namespace amrex;

void
PFFlame::Advance (int lev, Real time, Real dt)
{
  std::swap(phi_old[0][lev], phi_new[0][lev]);
  const Real* dx = geom[lev].CellSize();

  amrex::Array<std::unique_ptr<amrex::MultiFab> > Sborder(number_of_fabs);

  for (int n=0; n<number_of_fabs; n++)
    Sborder[n].reset(new amrex::MultiFab(grids[lev], dmap[lev], number_of_components,1)); 

  FillPatch(lev,t_old[lev],Sborder,0); // TODO Put this up in timestep

  for ( MFIter mfi(*phi_new[0][lev],true); mfi.isValid(); ++mfi )
    {
      const Box& bx = mfi.tilebox();

      amrex::BaseFab<Real> &old_phi = (*Sborder[0])[mfi];
      amrex::BaseFab<Real> &new_phi = (*phi_new[0][lev])[mfi];

      amrex::Array<amrex::Real> Laplacian(number_of_fabs);
      for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
	for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
	  {
	    //
	    // Phase field evolution
	    //

	    amrex::Real eta_lap = 
	      (old_phi(amrex::IntVect(i+1,j),0) - 2.*old_phi(amrex::IntVect(i,j),0) + old_phi(amrex::IntVect(i-1,j),0))/dx[0]/dx[0] +
	      (old_phi(amrex::IntVect(i,j+1),0) - 2.*old_phi(amrex::IntVect(i,j),0) + old_phi(amrex::IntVect(i,j-1),0))/dx[1]/dx[1];
	     
	    amrex::Real a0=w0, a1=0.0, a2= -5*w1 + 16*w12 - 11*w0, a3=14*w1 - 32*w12 + 18*w0, a4=-8*w1 + 16*w12 - 8*w0;
	     
	    amrex::Real oldphi = old_phi(amrex::IntVect(i,j),0);
	    new_phi(amrex::IntVect(i,j),0) = oldphi -
	      M * dt * (a1 + 2*a2*oldphi + 3*a3*oldphi*oldphi + 4*a4*oldphi*oldphi*oldphi
			- kappa*eta_lap);

	    //
	    // Temperature evolution
	    //

	    if (new_phi(amrex::IntVect(i,j),0)>0.01)
	      {
		amrex::Real temperature_delay = 0.05;
		if (time<temperature_delay) continue;

		amrex::Real eta_gradx = (old_phi(amrex::IntVect(i+1,j),0) - old_phi(amrex::IntVect(i-1,j),0))/(2*dx[0]);
		amrex::Real eta_grady = (old_phi(amrex::IntVect(i,j+1),0) - old_phi(amrex::IntVect(i,j-1),0))/(2*dx[1]);
		amrex::Real T_gradx = (old_phi(amrex::IntVect(i+1,j),1) - old_phi(amrex::IntVect(i-1,j),1))/(2*dx[0]);
		amrex::Real T_grady = (old_phi(amrex::IntVect(i,j+1),1) - old_phi(amrex::IntVect(i,j-1),1))/(2*dx[1]);

		amrex::Real eta_grad_mag = sqrt(eta_gradx*eta_gradx + eta_grady*eta_grady);

		amrex::Real T_lap = 
		  (old_phi(amrex::IntVect(i+1,j),1) - 2.*old_phi(amrex::IntVect(i,j),1) + old_phi(amrex::IntVect(i-1,j),1))/dx[0]/dx[0] +
		  (old_phi(amrex::IntVect(i,j+1),1) - 2.*old_phi(amrex::IntVect(i,j),1) + old_phi(amrex::IntVect(i,j-1),1))/dx[1]/dx[1];
	     
		amrex::Real rho = (rho1-rho0)*old_phi(amrex::IntVect(i,j),0) + rho0;
		amrex::Real k   = (k1-k0)*old_phi(amrex::IntVect(i,j),0) + k0;
		amrex::Real cp  = (cp1-cp0)*old_phi(amrex::IntVect(i,j),0) + cp0;

		new_phi(amrex::IntVect(i,j),1) = old_phi(amrex::IntVect(i,j),1) + (dt/rho/cp) *
		  ((k1-k0)*(eta_gradx*T_gradx + eta_grady*T_grady)
		   + k*T_lap
		   + (w1 - w0 - qdotburn)*eta_grad_mag);
	      }
	    else
	      new_phi(amrex::IntVect(i,j),1) = 0;
	  }
    }
}

