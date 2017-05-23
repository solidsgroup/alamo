#ifndef MY_AMR_H_
#define MY_AMR_H_

#include <AMReX_AmrMesh.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Geometry.H>

class MyAmr
  : public amrex::AmrMesh
{
public:
  using amrex::AmrMesh::AmrMesh;  // inherit AmrMesh constructors

  // MyAmr(int _n_cell, int _max_grid_size, int _plot_int, int _nsteps) :
  //   n_cell(_n_cell), max_grid_size(_max_grid_size), plot_int(_plot_int), nsteps(_nsteps)
  // {
  //   Initialize();
    
  //   amrex::BoxArray ba; 
  //   // Geometry geom; // < -- defined in AmrMesh
  //   {
  //     amrex::IntVect dom_lo(amrex::IntVect(D_DECL(0,0,0)));
  //     amrex::IntVect dom_hi(amrex::IntVect(D_DECL(n_cell-1, n_cell-1, n_cell-1)));
  //     amrex::Box domain(dom_lo, dom_hi);

  //     // Initialize the boxarray "ba" from the single box "bx"
  //     ba.define(domain);
  //     // Break up boxarray "ba" into chunks no larger than "max_grid_size" along a direction
  //     ba.maxSize(max_grid_size);

  //     // Define box over domain [-1,1]x[-1,1]
  //     for (int n = 0; n < BL_SPACEDIM; n++) {
  //   	real_box.setLo(n,-1.0);
  //   	real_box.setHi(n, 1.0);
  //     }

  //     // Cartesian coordinates
  //     int coord = 0;
    	
  //     // Periodic boundary conditions
  //     int is_periodic[BL_SPACEDIM];
  //     for (int i = 0; i < BL_SPACEDIM; i++) {
  //   	is_periodic[i] = 1;
  //     }

  //     // This defines a Geometry object
  //     //geom.define(domain,&real_box,coord,is_periodic);
  //     amrex::Geometry::Setup(&real_box,coord);
  //     grids[0] = ba;
  //   }
    

  // }

  virtual ~MyAmr() {}

  //! Manually tag.  Note that tags is built on level lev grids coarsened by bf_lev.
  virtual void ManualTagsPlacement (int lev, amrex::TagBoxArray& tags,
				    const amrex::Array<amrex::IntVect>& bf_lev) override;

  // We could also override 
  //! Tag cells for refinement.  TagBoxArray tags is built on level lev grids.
  // virtual void ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow) override;

private:
  int n_cell, max_grid_size,plot_int,nsteps; // input parameters

  amrex::RealBox real_box;

  amrex::Array<std::unique_ptr<amrex::MultiFab> > phi_new;
  amrex::Array<std::unique_ptr<amrex::MultiFab> > phi_old;
};

#endif
