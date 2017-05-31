
#include <limits>

#include <AMReX_MultiFabUtil.H>
#include <PFAmr.H>

using namespace amrex;

void
PFAmr::Evolve ()
{
  Real cur_time = t_new[0];
  int last_plot_file_step = 0;

  for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
      if (ParallelDescriptor::IOProcessor()) {
	std::cout << "\nSTEP " << step+1 << " starts ..." << std::endl;
      }
      int lev = 0;
      int iteration = 1;
      timeStep(lev, cur_time, iteration);

      cur_time += dt[0];

      if (ParallelDescriptor::IOProcessor()) {
	std::cout << "STEP " << step+1 << " ends." << " TIME = " << cur_time << " DT = " << dt[0]
		  << std::endl;
      }

      // sync up time
      for (int lev = 0; lev <= finest_level; ++lev) {
	t_new[lev] = cur_time;
      }

      if (plot_int > 0 && (step+1) % plot_int == 0) {
	last_plot_file_step = step+1;
	WritePlotFile();
      }

      if (cur_time >= stop_time - 1.e-6*dt[0]) break;
    }

  if (plot_int > 0 && istep[0] > last_plot_file_step) {
    WritePlotFile();
  }
}

void
PFAmr::timeStep (int lev, Real time, int iteration)
{

  if (regrid_int > 0)  // We may need to regrid
    {
      static Array<int> last_regrid_step(max_level+1, 0);

      // regrid doesn't change the base level, so we don't regrid on max_level
      if (lev < max_level && istep[lev] > last_regrid_step[lev])
	{
          if (istep[lev] % regrid_int == 0)
	    {
	      int old_finest = finest_level; // regrid changes finest_level
	      regrid(lev, time, false);
	      for (int k = lev; k <= finest_level; ++k) {
		last_regrid_step[k] = istep[k];
	      }

	      for (int k = old_finest+1; k <= finest_level; ++k) {
		dt[k] = dt[k-1] / MaxRefRatio(k-1);
	      }
  	    }
  	}
    }

  if (//Verbose() &&
      ParallelDescriptor::IOProcessor()) {
    std::cout << "[Level " << lev 
	      << " step " << istep[lev]+1 << "] ";
    std::cout << "ADVANCE with dt = "
	      << dt[lev]
	      << std::endl;
  }

  Advance(lev, time, dt[lev]);

  ++istep[lev];

  if (//Verbose() &&
      ParallelDescriptor::IOProcessor())
    {
      std::cout << "[Level " << lev
		<< " step " << istep[lev] << "] ";
      std::cout << "Advanced "
		<< CountCells(lev)
		<< " cells"
		<< std::endl;
    }

  if (lev < finest_level)
    {
      // for (int i = 1; i <= nsubsteps[lev+1]; ++i)
      // {
      //    timeStep(lev+1, time+(i-1)*dt[lev+1], i);
      // }
      timeStep(lev+1, time, 1);

      for (int n = 0; n < number_of_grains; n++)
	amrex::average_down(*phi_new[n][lev+1], *phi_new[n][lev],
			    geom[lev+1], geom[lev],
			    0, phi_new[n][lev]->nComp(), refRatio(lev));
    }
}

void
PFAmr::Advance (int lev, Real time, Real dt)
{
  dt = timestep;

  std::swap(phi_old[0][lev], phi_new[0][lev]);
  std::swap(phi_old[1][lev], phi_new[1][lev]);

  const Real* dx = geom[lev].CellSize();
  
  amrex::Array<std::unique_ptr<amrex::MultiFab> > Sborder(number_of_grains);
  for (int n=0; n<number_of_grains; n++)
    Sborder[n].reset(new amrex::MultiFab(grids[lev], dmap[lev], 1,1)); 
  FillPatch(lev,t_old[lev],Sborder,0);

  for ( MFIter mfi(*phi_new[0][lev],true); mfi.isValid(); ++mfi )
    {
      const Box& bx = mfi.tilebox();
    
      amrex::BaseFab<Real> &old_phi1 = (*Sborder[0])[mfi];
      amrex::BaseFab<Real> &new_phi1 = (*phi_new[0][lev])[mfi];

      amrex::BaseFab<Real> &old_phi2 = (*Sborder[1])[mfi];
      amrex::BaseFab<Real> &new_phi2 = (*phi_new[1][lev])[mfi];
      amrex::Real L=1.,mu=10.,gamma=1.,kappa=0.01;
      for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
	for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
	  {
	    // //--- Heat Equation ----
	    // new_phi1(amrex::IntVect(i,j)) = old_phi1(amrex::IntVect(i,j)) +
	    //   dt*(old_phi1(amrex::IntVect(i+1,j)) - 2.*old_phi1(amrex::IntVect(i,j)) + old_phi1(amrex::IntVect(i-1,j)))/dx[0]/dx[0] +
	    //   dt*(old_phi1(amrex::IntVect(i,j+1)) - 2.*old_phi1(amrex::IntVect(i,j)) + old_phi1(amrex::IntVect(i,j-1)))/dx[1]/dx[1];
	    // new_phi2(amrex::IntVect(i,j)) = old_phi2(amrex::IntVect(i,j)) +
	    //   dt*(old_phi2(amrex::IntVect(i+1,j)) - 2.*old_phi2(amrex::IntVect(i,j)) + old_phi2(amrex::IntVect(i-1,j)))/dx[0]/dx[0] +
	    //   dt*(old_phi2(amrex::IntVect(i,j+1)) - 2.*old_phi2(amrex::IntVect(i,j)) + old_phi2(amrex::IntVect(i,j-1)))/dx[1]/dx[1];


	    // --- Allen Cahn ----
	    amrex::Real
	      Lap1 = 
	      (old_phi1(amrex::IntVect(i+1,j)) - 2.*old_phi1(amrex::IntVect(i,j)) + old_phi1(amrex::IntVect(i-1,j)))/dx[0]/dx[0] +
	      (old_phi1(amrex::IntVect(i,j+1)) - 2.*old_phi1(amrex::IntVect(i,j)) + old_phi1(amrex::IntVect(i,j-1)))/dx[1]/dx[1],
	      Lap2 =
	      (old_phi2(amrex::IntVect(i+1,j)) - 2.*old_phi2(amrex::IntVect(i,j)) + old_phi2(amrex::IntVect(i-1,j)))/dx[0]/dx[0] +
	      (old_phi2(amrex::IntVect(i,j+1)) - 2.*old_phi2(amrex::IntVect(i,j)) + old_phi2(amrex::IntVect(i,j-1)))/dx[1]/dx[1];

	    new_phi1(amrex::IntVect(i,j)) = old_phi1(amrex::IntVect(i,j)) -
	      L * dt * ( mu*( old_phi1(amrex::IntVect(i,j))*old_phi1(amrex::IntVect(i,j))*old_phi1(amrex::IntVect(i,j))
	    		      - old_phi1(amrex::IntVect(i,j))
	    		      + 2*gamma*old_phi1(amrex::IntVect(i,j))*old_phi2(amrex::IntVect(i,j))*old_phi2(amrex::IntVect(i,j)))
	    		 - kappa*Lap1);

	    new_phi2(amrex::IntVect(i,j)) = old_phi2(amrex::IntVect(i,j)) -
	      L * dt * ( mu*( old_phi2(amrex::IntVect(i,j))*old_phi2(amrex::IntVect(i,j))*old_phi2(amrex::IntVect(i,j))
	    		      - old_phi2(amrex::IntVect(i,j))
	    		      + 2*gamma*old_phi2(amrex::IntVect(i,j))*old_phi1(amrex::IntVect(i,j))*old_phi1(amrex::IntVect(i,j)))
	    		 - kappa*Lap2);

	  }
    }
}
