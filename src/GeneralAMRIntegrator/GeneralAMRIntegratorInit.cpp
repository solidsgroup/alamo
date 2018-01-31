#include <stdlib.h>
#include <time.h>
#include <AMReX_MultiFabUtil.H>
#include "GeneralAMRIntegrator.H"

#include "GeneralAMRIntegratorBC.H"

using namespace amrex;

void
GeneralAMRIntegrator::InitData ()
{
    if (restart_chkfile.empty())
    {
	const Real time = 0.0;
	InitFromScratch(time);
	for (int lev = finest_level-1; lev >= 0; --lev)
	  {
	    for (int n = 0; n < number_of_fabs; n++)
	      amrex::average_down(*(*fab_array[n])[lev+1], *(*fab_array[n])[lev],
				geom[lev+1], geom[lev],
				  0, (*fab_array[n])[lev]->nComp(), refRatio(lev));
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

void GeneralAMRIntegrator::MakeNewLevelFromScratch (int lev, Real t, const BoxArray& ba,
						    const DistributionMapping& dm)
{
  for (int n = 0 ; n < number_of_fabs; n++)
    (*fab_array[n])[lev].reset(new MultiFab(ba, dm, ncomp_array[n], nghost_array[n]));
  
  // TODO - update this
  // t_new[lev] = t;
  // t_old[lev] = t - 1.e200;

  // const amrex::Real width = geom[lev].ProbHi()[0] - geom[lev].ProbHi()[1];
  // for (MFIter mfi(*phi_new[0][lev],true); mfi.isValid(); ++mfi)
  //   {
  //     const Box& box = mfi.tilebox();

  //     amrex::BaseFab<Real> &phi_box = (*phi_new[0][lev])[mfi];
  //     amrex::BaseFab<Real> &phi_box_old = (*phi_old[0][lev])[mfi];

  //     for (int i = box.loVect()[0]-nghost; i<=box.hiVect()[0]+nghost; i++) 
  //      	for (int j = box.loVect()[1]-nghost; j<=box.hiVect()[1]+nghost; j++)
  // 	    {

  // 	    }
  //   }
  
  // physbc.SetLevel(lev);
  // physbc.FillBoundary(*phi_new[0][lev],0,0,t);
  // physbc.FillBoundary(*phi_old[0][lev],0,0,t);

}
