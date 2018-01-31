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
  
  // TODO - still don't understand this?
  t_new[lev] = t;
  t_old[lev] = t - 1.e200;

  Initialize(lev);
  
  for (int n = 0 ; n < number_of_fabs; n++)
    {
      physbc_array[n]->SetLevel(lev);
      physbc_array[n]->FillBoundary(*(*fab_array[n])[lev],0,0,t);
    }
}

