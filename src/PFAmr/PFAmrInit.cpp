#include <AMReX.H>
#if AMREX_SPACEDIM==2
#include <stdlib.h>
#include <time.h>
#include <AMReX_MultiFabUtil.H>
#include "PFAmr.H"

#include "PFAmrBC.H"

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


#if BL_SPACEDIM==2
#define INTVECT amrex::IntVect(i,j)
#elif BL_SPACEDIM==3
#define INTVECT amrex::IntVect(i,j,k)
#endif
void PFAmr::MakeNewLevelFromScratch (int lev, Real t, const BoxArray& ba,
				      const DistributionMapping& dm)
{
  phi_new[0][lev].reset(new MultiFab(ba, dm, number_of_grains+2, nghost));
  phi_old[0][lev].reset(new MultiFab(ba, dm, number_of_grains+2, nghost));


  t_new[lev] = t;
  t_old[lev] = t - 1.e200;

  //const Real* dx = geom[lev].CellSize();
  const amrex::Real width = geom[lev].ProbHi()[0] - geom[lev].ProbHi()[1];
  //Real cur_time = t_new[lev];
  for (MFIter mfi(*phi_new[0][lev],true); mfi.isValid(); ++mfi)
    {
      const Box& box = mfi.tilebox();

      amrex::BaseFab<Real> &phi_box = (*phi_new[0][lev])[mfi];
      amrex::BaseFab<Real> &phi_box_old = (*phi_old[0][lev])[mfi];

      for (int i = box.loVect()[0]-nghost; i<=box.hiVect()[0]+nghost; i++) 
       	for (int j = box.loVect()[1]-nghost; j<=box.hiVect()[1]+nghost; j++)
#if BL_SPACEDIM==3
	  for (int k = box.loVect()[2]-nghost; k<=box.hiVect()[2]+nghost; k++)
#endif
	    {
	      amrex::Real x = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];
	      amrex::Real y = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];
#if BL_SPACEDIM==3
	      amrex::Real z = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];
#endif

	      if (anisotropy)
		{
		  if (ic_type == "circle")
		    {
		      //
		      // circle IC
		      //
		      phi_box_old(INTVECT,0) = 0.; // good practice to initialize all new memory
		      phi_box_old(INTVECT,1) = 0.; // good practice to initialize all new memory

		      if (sqrt(x*x+y*y)<=0.5)

			{
			  phi_box(INTVECT,0) = 1.;     
			  phi_box(INTVECT,1) = 0.;     
			}
		      else
			{
			  phi_box(INTVECT,0) = 0.;     
			  phi_box(INTVECT,1) = 1.;     
			}
		    }
		  else if (ic_type == "perturbed_bar")
		    {
		      //
		      // perturbed bar IC
		      //
		      phi_box_old(INTVECT,0) = 0.; // good practice to initialize all new memory
		      phi_box_old(INTVECT,1) = 0.; // good practice to initialize all new memory

		      amrex::Real pi = 3.14159265359;
		      amrex::Real bdry  =
			(+ 0.1 * sin(x*(2*pi)/width)
			+ 0.1 * cos(2.0*x*(2*pi)/width)
			+ 0.1 * sin(3.0*x*(2*pi)/width)
			+ 0.1 * cos(4.0*x*(2*pi)/width)
			+ 0.1 * sin(5.0*x*(2*pi)/width)
			+ 0.15 * cos(6.0*x*(2*pi)/width)
			+ 0.15 * sin(7.0*x*(2*pi)/width)
			+ 0.15 * cos(8.0*x*(2*pi)/width)
			+ 0.15 * sin(9.0*x*(2*pi)/width)
			 + 0.15 * cos(10.0*x*(2*pi)/width))*0.1
			;

		      if (y < bdry)

			{
			  phi_box(INTVECT,0) = 1.;     
			  phi_box(INTVECT,1) = 0.;     
			}
		      else
			{
			  phi_box(INTVECT,0) = 0.;     
			  phi_box(INTVECT,1) = 1.;     
			}
		    }
		  else if (ic_type == "random")
		    {
		      phi_box(INTVECT,0) = ((amrex::Real)rand())/((amrex::Real)RAND_MAX);
		      phi_box(INTVECT,1) = 1.0 - phi_box(INTVECT,0);
		    }


		  phi_box(INTVECT,number_of_grains) = 0.;
		  phi_box(INTVECT,number_of_grains+1) = 0.;
		}	      
	      else
		{
		  //
		  // voronoi tesselation IC
		  //
		  amrex::Real len_x = (geom[0].ProbHi(0)-geom[0].ProbLo(0));
		  amrex::Real len_y = (geom[0].ProbHi(1)-geom[0].ProbLo(1));
#if BL_SPACEDIM==3
		  amrex::Real len_z = (geom[0].ProbHi(2)-geom[0].ProbLo(2));
#endif
		  amrex::Real min_distance = std::numeric_limits<amrex::Real>::infinity();
		  int min_grain_id = -1;
		  for (int n = 0; n<number_of_grains; n++)
		    {
		      phi_box(INTVECT,n) = 0.;     // initialize
		      phi_box_old(INTVECT,n) = 0.; // good practice to initialize all new memory
		      for (amrex::Real offset_x = -len_x; offset_x <= len_x; offset_x += len_x)
			for (amrex::Real offset_y = -len_y; offset_y <= len_y; offset_y += len_y)
#if BL_SPACEDIM==3
			  for (amrex::Real offset_z = -len_z; offset_z <= len_z; offset_z += len_z)
#endif
			    {
#if BL_SPACEDIM==2
			      amrex::Real d = sqrt((x-voronoi_x[n]-offset_x)*(x-voronoi_x[n]-offset_x) + (y-voronoi_y[n]-offset_y)*(y-voronoi_y[n]-offset_y));
#elif BL_SPACEDIM==3
			      amrex::Real d = sqrt((x-voronoi_x[n]-offset_x)*(x-voronoi_x[n]-offset_x) + (y-voronoi_y[n]-offset_y)*(y-voronoi_y[n]-offset_y) + (z-voronoi_z[n]-offset_z)*(z-voronoi_z[n]-offset_z));
#endif
			      if (d<min_distance )  {min_distance = d;  min_grain_id = n;}
			    }
		    }
		  phi_box(INTVECT,min_grain_id) = 1.;

		  phi_box(INTVECT,number_of_grains) = (amrex::Real)min_grain_id;
		  phi_box(INTVECT,number_of_grains+1) = 0;
		}	      

	      phi_box_old(INTVECT,number_of_grains) = 0;   // Good practice to initialize
	      phi_box_old(INTVECT,number_of_grains+1) = 0; // all newly created data

	    }
    }

  
  //PFAmrPhysBC physbc(geom,bc_lo,bc_hi);
  physbc.SetLevel(lev);
  physbc.FillBoundary(*phi_new[0][lev],0,0,t);
  physbc.FillBoundary(*phi_old[0][lev],0,0,t);
}
#endif
