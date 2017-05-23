#include <iostream>
#include <fstream>
#include <iomanip>

#include <AMReX_ParmParse.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>
#include <AMReX_VisMF.H>
#include <AMReX_BoxIterator.H>

#include <AMReX_MultiFabUtil.H>
#include "MyAmr.H"
#include "writePlotFile.H"

#include "AmrAdv.H"

using namespace amrex;

void advance (MultiFab& old_phi, MultiFab& new_phi,
	      Real dt, const Geometry& geom)
{
  // Fill the ghost cells of each grid from the other grids
  // includes periodic domain boundaries
  old_phi.FillBoundary(geom.periodicity());

  int Ncomp = old_phi.nComp();
  int ng_p = old_phi.nGrow();

  const Real* dx = geom.CellSize();

  // Update
  for ( MFIter mfi(old_phi); mfi.isValid(); ++mfi )
    {
      const Box& bx = mfi.validbox();
    
      amrex::BaseFab<Real> &old_phi_box = old_phi[mfi];
      amrex::BaseFab<Real> &new_phi_box = new_phi[mfi];

      for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
 	for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
	  new_phi_box(amrex::IntVect(i,j)) = old_phi_box(amrex::IntVect(i,j)) +
	    dt*(old_phi_box(amrex::IntVect(i+1,j)) - 2.*old_phi_box(amrex::IntVect(i,j)) + old_phi_box(amrex::IntVect(i-1,j)))/dx[0]/dx[0] +
	    dt*(old_phi_box(amrex::IntVect(i,j+1)) - 2.*old_phi_box(amrex::IntVect(i,j)) + old_phi_box(amrex::IntVect(i,j-1)))/dx[1]/dx[1];
    }
}

int main (int argc, char* argv[])
{
  amrex::Initialize(argc,argv);




  AmrAdv amradv;
  amradv.InitData();

  amradv.Evolve();

  amrex::Finalize();
  exit(0);








  Real strt_time = ParallelDescriptor::second();

  std::cout << std::setprecision(15);

  int n_cell, max_grid_size, nsteps, plot_int, is_periodic[BL_SPACEDIM];

  //MyAmr myamr(n_cell, max_grid_size, plot_int, nsteps);

  //
  // READ INPUT PARAMETERS
  //

  // {
  //   ParmParse pp;
  //   pp.get("n_cell",n_cell); // number of cells on each edge
  //   pp.get("max_grid_size",max_grid_size);
  //   plot_int = 1; pp.query("plot_int",plot_int);
  //   nsteps = 0; pp.query("nsteps",nsteps);
  // }


  //
  // CONSTRUCT GEOMETRY
  //

  BoxArray ba;
  Geometry geom;
  RealBox real_box;
  {
    IntVect dom_lo(IntVect(D_DECL(0,0,0)));
    IntVect dom_hi(IntVect(D_DECL(n_cell-1, n_cell-1, n_cell-1)));
    Box domain(dom_lo, dom_hi);

    // Initialize the boxarray "ba" from the single box "bx"
    ba.define(domain);
    // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
    ba.maxSize(max_grid_size);

    // Define box over domain [-1,1]x[-1,1]
    for (int n = 0; n < BL_SPACEDIM; n++) {
      real_box.setLo(n,-1.0);
      real_box.setHi(n, 1.0);
    }

    // Cartesian coordinates
    int coord = 0;
	
    // Periodic boundary conditions
    int is_periodic[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++) {
      is_periodic[i] = 1;
    }

    // This defines a Geometry object
    geom.define(domain,&real_box,coord,is_periodic);
  }


  exit(0);


  //
  // SET UP AMR
  //


  /// BEGIN - AMR
  int max_level = 2;
  amrex::Array<int> n_cells{D_DECL(n_cell,n_cell,n_cell)};
  MyAmr amr(&real_box,max_level,n_cells);
  amr.MakeNewGrids();
  /// END - AMR
    std::cout << __LINE__ << std::endl;

  // define dx[]
  const Real* dx = geom.CellSize();

  
  int Nghost = 1;  // number of ghost cells for each array 

  
  int Ncomp  = 1;  // number of components for each array

  
  Real time = 0.0; // starting time in the simulation
  
    std::cout << __LINE__ << std::endl;
  DistributionMapping dm(ba);

  Array<std::unique_ptr<MultiFab> > phi(2);
  //phi[0].reset(new MultiFab(ba, dm, Ncomp, Nghost));
  //phi[1].reset(new MultiFab(ba, dm, Ncomp, Nghost));

  std::cout << __LINE__ << std::endl;
  phi[0].reset(new MultiFab(amr.boxArray(0), amr.DistributionMap(0), Ncomp, Nghost));
  //phi[1].reset(new MultiFab(amr.boxArray(0), amr.DistributionMap(0), Ncomp, Nghost));

  std::cout << __LINE__ << std::endl;
  phi[0]->setVal(0.0);
  //phi[1]->setVal(0.0);

    std::cout << __LINE__ << std::endl;


  //
  // INITIALIZE PHI
  //

  int init_index = 0;
  // for ( MFIter mfi(*phi[init_index]); mfi.isValid(); ++mfi )
  //   {
  //     const Box& bx = mfi.validbox();
  //     amrex::BaseFab<Real> &phi_box = (*phi[init_index])[mfi];
  //     for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
  // 	for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
  // 	  {
  // 	    amrex::Real x = geom.ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom.CellSize()[0];
  // 	    amrex::Real y = geom.ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom.CellSize()[1];
  // 	    // gaussian
  // 	    phi_box(amrex::IntVect(i,j)) =  1. + exp(-( (x-0.25)*(x-0.25) + (y-0.25)*(y-0.25))/0.01);
  // 	    // multiple
  // 	    //phi_box(amrex::IntVect(i,j)) =  x*y;
  // 	  }
  //   }


  //
  // INITIAL MESH REFINEMENT
  //
  MultiFab mf(amr.boxArray(0), amr.DistributionMap(0), Ncomp, Nghost);
  mf.setVal(0.);

  IntVect ref_ratio = IntVect::TheUnitVector();
  for (int lev = 1; lev <= amr.finestLevel(); ++lev)
    {
      MultiFab fmf(amr.boxArray(lev), amr.DistributionMap(lev), Ncomp, Nghost);
      fmf.setVal(static_cast<Real>(lev));
      //ref_ratio *= amr.refRatio(lev-1);
      // amrex::average_down(fmf,    // fine fab
      // 			  *phi[0],        // course fab
      // 			  0,          // int scomp
      // 			  1,          // int ncomp
      //			  ref_ratio); // int rr
      //amrex::average_down(*phi[1],fmf, 0, 1, ref_ratio);
      if (lev==2)
	amrex::average_down(fmf,mf, 0, Ncomp, ref_ratio);
      //const std::string& pltfile = amrex::Concatenate("pltamr_",lev,5);
      //writePlotFile(pltfile, fmf, amr.Geom()[lev], time);

    }

  writePlotFile("pltamr_00000", mf, amr.Geom()[0], time);


  // compute the time step
  Real dt = 0.9*dx[0]*dx[0] / (2.0*BL_SPACEDIM);

  // // Write a plotfile of the initial data if plot_int > 0 (plot_int was defined in the inputs file)
  // if (plot_int > 0)
  //   {
  //     int n = 0;
  //     //const std::string& pltfile = amrex::Concatenate("plt",n,5);
  //     //writePlotFile(pltfile, *phi[init_index], geom, time);

  //     const std::string& pltfile0 = amrex::Concatenate("plt0_",n,5);
  //     MultiFab fmf0(amr.boxArray(0), amr.DistributionMap(0), Ncomp, Nghost);
  //     writePlotFile(pltfile0, fmf0, amr.Geom()[0], time);

  //     const std::string& pltfile1 = amrex::Concatenate("plt1_",n,5);
  //     MultiFab fmf1(amr.boxArray(1), amr.DistributionMap(1), Ncomp, Nghost);
  //     writePlotFile(pltfile1, fmf1, amr.Geom()[1], time);

  //     const std::string& pltfile2 = amrex::Concatenate("plt2_",n,5);
  //     MultiFab fmf2(amr.boxArray(2), amr.DistributionMap(2), Ncomp, Nghost);
  //     writePlotFile(pltfile2, fmf2, amr.Geom()[2], time);

  //   }


  exit(0);

  //
  // ITERATE
  // 

  int old_index = init_index;
  for (int n = 1; n <= nsteps; n++, old_index = 1 - old_index)
    {
      // swap index
      int new_index = 1 - old_index;

      // integrate over one timestep
      advance(*phi[old_index], *phi[new_index], dt, geom); 
      time = time + dt;

      // print notification
      if (ParallelDescriptor::IOProcessor())
	std::cout << "Advanced step " << n << std::endl;

      // print data
      if (plot_int > 0 && n%plot_int == 0)
	{
	  const std::string& pltfile = amrex::Concatenate("plt",n,5);
	  writePlotFile(pltfile, *phi[new_index], geom, time);
	}
    }



  //
  // FINISH
  //

  Real stop_time = ParallelDescriptor::second() - strt_time;
  const int IOProc = ParallelDescriptor::IOProcessorNumber();
  ParallelDescriptor::ReduceRealMax(stop_time,IOProc);

  if (ParallelDescriptor::IOProcessor()) {
    std::cout << "Run time = " << stop_time << std::endl;
  }
  amrex::Finalize();
  return 0;
}
