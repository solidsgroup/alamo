
#include <limits>

#include <AMReX_MultiFabUtil.H>
#include "GeneralAMRIntegrator.H"

using namespace amrex;

void
GeneralAMRIntegrator::Evolve ()
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
      TimeStep(lev, cur_time, iteration);
      cur_time += dt[0];

      if (ParallelDescriptor::IOProcessor()) {
	std::cout << "STEP " << step+1 << " ends."
		  << " TIME = " << cur_time << " DT = " << dt[0]
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
GeneralAMRIntegrator::TimeStep (int lev, Real time, int iteration)
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

  if (Verbose() && ParallelDescriptor::IOProcessor()) {
    std::cout << "[Level " << lev 
	      << " step " << istep[lev]+1 << "] ";
    std::cout << "ADVANCE with dt = "
	      << dt[lev]
	      << std::endl;
  }
  Advance(lev, time, dt[lev]);
  ++istep[lev];

  if (Verbose() && ParallelDescriptor::IOProcessor())
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
      for (int i = 1; i <= nsubsteps[lev+1]; ++i)
	TimeStep(lev+1, time+(i-1)*dt[lev+1], i);

      for (int n = 0; n < number_of_fabs; n++)
	{
	  amrex::average_down(*(*fab_array[n])[lev+1], *(*fab_array[n])[lev],
			      geom[lev+1], geom[lev],
			      0, (*fab_array[n])[lev]->nComp(), refRatio(lev));
	}
    }
}

// void
// GeneralAMRIntegrator::Advance (int lev, Real time, Real dt)
// {
//   // TODO - replace this
//   // std::swap(phi_old[0][lev], phi_new[0][lev]);
//   const Real* dx = geom[lev].CellSize();

//   amrex::Array<std::unique_ptr<amrex::MultiFab> > Sborder(number_of_fabs);

//   // TODO - fix this
//   for (int n=0; n<number_of_fabs; n++)
//     Sborder[n].reset(new amrex::MultiFab(grids[lev], dmap[lev], number_of_grains+2,nghost)); 

//   FillPatch(lev,t_old[lev],Sborder,0);

//   // for ( MFIter mfi(*phi_new[0][lev],true); mfi.isValid(); ++mfi )
//   //   {
//   //     const Box& bx = mfi.tilebox();

//   //     amrex::BaseFab<Real> &old_phi = (*Sborder[0])[mfi];
//   //     amrex::BaseFab<Real> &new_phi = (*phi_new[0][lev])[mfi];

//   //     amrex::Array<amrex::Real> Laplacian(number_of_fabs);
//   //     for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
//   // 	for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
//   // 	  {
//   // 	    new_phi(amrex::IntVect(i,j)) = i*j*dx[0]*dx[1]; // here's where we put the numerics
//   // 	  }
//   //   }
// }

