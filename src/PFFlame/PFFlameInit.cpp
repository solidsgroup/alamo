#include <stdlib.h>
#include <time.h>
#include <AMReX_MultiFabUtil.H>
#include <PFFlame.H>
#include <PFFlameBC.H>

using namespace amrex;

// void
// PFFlame::InitData ()
// {
//     if (restart_chkfile.empty())
//     {
// 	const Real time = 0.0;
// 	InitFromScratch(time);
// 	for (int lev = finest_level-1; lev >= 0; --lev)
// 	  {
// 	    amrex::average_down(*phi_new[0][lev+1], *phi_new[0][lev],
// 				geom[lev+1], geom[lev],
// 				0, phi_new[0][lev]->nComp(), refRatio(lev));
// 	  }

// 	if (plot_int > 0) {
// 	    WritePlotFile();
// 	}
//     }
//     else
//     {
// 	InitFromCheckpoint();
//     }
// }

 void PFFlame::MakeNewLevelFromScratch (int lev, Real t, const BoxArray& ba,
 				      const DistributionMapping& dm)
 {
   const int nghost = 1;

   phi_new[0][lev].reset(new MultiFab(ba, dm, number_of_components, nghost));
   phi_old[0][lev].reset(new MultiFab(ba, dm, number_of_components, nghost));


   t_new[lev] = t;
   t_old[lev] = t - 1.e200;

   const Real* dx = geom[lev].CellSize();
   const Real* prob_lo = geom[lev].ProbLo();
   Real cur_time = t_new[lev];
   for (MFIter mfi(*phi_new[0][lev],true); mfi.isValid(); ++mfi)
     {
       const Box& box = mfi.tilebox();

       amrex::BaseFab<Real> &phi_box = (*phi_new[0][lev])[mfi];

       double offset_x = 0., offset_y=0.;
       for (int i = box.loVect()[0]-nghost; i<=box.hiVect()[0]+nghost; i++) // todo
        	for (int j = box.loVect()[1]-nghost; j<=box.hiVect()[1]+nghost; j++)
 	  {
 	    amrex::Real x = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];

 	    //phi_box(amrex::IntVect(i,j),0) =  1.;

 	    // if (x > -0.5 & x < 0.5)
 	    //   {
 	    //  	phi_box(amrex::IntVect(i,j),0) =  1.;
 	    //  	phi_box(amrex::IntVect(i,j),1) =  0.;
 	    //   }
 	    // else
 	    //   {
 	    // 	phi_box(amrex::IntVect(i,j),0) =  0.;
 	    //  	phi_box(amrex::IntVect(i,j),1) =  0.;
 	    //   }
	    
 	    // DEBUG;
 	    // if (x < -0.75)
 	    //   {
 	    // 	DEBUG;
 	    //   	phi_box(amrex::IntVect(i,j),0) =  0.;
 	    // 	DEBUG;
 	    //   	//phi_box(amrex::IntVect(i,j),1) =  1.;
 	    //   }
 	    // else
 	    //   {
 	    // 	DEBUG;
 	    phi_box(amrex::IntVect(i,j),0) =  1.0;

 	    phi_box(amrex::IntVect(i,j),1) =  0.0;
 	    // 	DEBUG;
 	    //   	//phi_box(amrex::IntVect(i,j),1) =  0.;
 	    //   }
 	  }
     }

   PFFlamePhysBC physbc(geom[lev]);
   physbc.FillBoundary(*phi_new[0][lev],0,0,t);
   physbc.FillBoundary(*phi_old[0][lev],0,0,t);
 }
