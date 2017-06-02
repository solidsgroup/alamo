#include <stdlib.h>
#include <time.h>
#include <AMReX_MultiFabUtil.H>
#include <PFAmr.H>

using namespace amrex;

void
PFAmr::InitData ()
{
    if (restart_chkfile.empty())
    {
	const Real time = 0.0;
	InitFromScratch(time);
	for (int lev = finest_level-1; lev >= 0; --lev)
	  {
	    amrex::average_down(*phi_new[0][lev+1], *phi_new[0][lev],
				geom[lev+1], geom[lev],
				0, phi_new[0][lev]->nComp(), refRatio(lev));
	  }

	if (plot_int > 0) {
	    WritePlotFile();
	}
    }
    else
    {
	InitFromCheckpoint();
    }
}

void PFAmr::MakeNewLevelFromScratch (int lev, Real t, const BoxArray& ba,
				      const DistributionMapping& dm)
{
  const int ncomp = number_of_grains;
  const int nghost = 1;

  phi_new[0][lev].reset(new MultiFab(ba, dm, ncomp, nghost));
  phi_old[0][lev].reset(new MultiFab(ba, dm, ncomp, nghost));

  t_new[lev] = t;
  t_old[lev] = t - 1.e200;

  const Real* dx = geom[lev].CellSize();
  const Real* prob_lo = geom[lev].ProbLo();
  Real cur_time = t_new[lev];
  for (MFIter mfi(*phi_new[0][lev],true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();

      amrex::BaseFab<Real> &phi_box = (*phi_new[0][lev])[mfi];

      //double offset_x = 0.25, offset_y=0.25;
      double offset_x = 0., offset_y=0.;
      for (int i = box.loVect()[0]; i<=box.hiVect()[0]; i++)
	for (int j = box.loVect()[1]; j<=box.hiVect()[1]; j++)
	  {
	    amrex::Real x = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];
	    amrex::Real y = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];
	    amrex::Real r = sqrt((x-offset_x)*(x-offset_x) + (y-offset_y)*(y-offset_y));
	    // // circular distribution
	    // if (r<0.5) phi_box(amrex::IntVect(i,j),0) =  1.;
	    // else phi_box(amrex::IntVect(i,j),0) =  0. + exp(-((r-0.5)*(r-0.5))/0.001);
	    // if (x<0)
	    //   {
	    // 	phi_box(amrex::IntVect(i,j),1) =  1. - phi_box(amrex::IntVect(i,j),0 );
	    // 	phi_box(amrex::IntVect(i,j),2) =  0;
	    //   }
	    // else
	    //   {
	    // 	phi_box(amrex::IntVect(i,j),2) =  1. - phi_box(amrex::IntVect(i,j),0 );
	    // 	phi_box(amrex::IntVect(i,j),1) =  0;
	    //   }


	    // random distribution
	    phi_box(amrex::IntVect(i,j)) = (amrex::Real)rand()/(amrex::Real)RAND_MAX;
	    amrex::Real sum = 0;
	    for (int n = 0; n < number_of_grains; n++)
	      {
		amrex::Real rand_num = (amrex::Real)rand()/(amrex::Real)RAND_MAX;
	     	phi_box(amrex::IntVect(i,j),n) = rand_num;
		sum += rand_num;
	      }
	    for (int n = 0; n < number_of_grains; n++)
	      phi_box(amrex::IntVect(i,j),n) /= sum;
	    
	      
	  }
    }
}
