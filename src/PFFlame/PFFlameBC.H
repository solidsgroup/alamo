#ifndef PFFlameBC_H_
#define PFFlameBC_H_

#include <AMReX_PhysBCFunct.H>

// In this tutorial, we assume the boundaries are all periodic.  
// Thus there are no physical boundaries. And FillBoundary does nothing.

class PFFlamePhysBC
  : public amrex::PhysBCFunctBase
{
public:
  PFFlamePhysBC (amrex::Geometry &_geom) : geom(_geom)  {}
  virtual ~PFFlamePhysBC () {}
  virtual void FillBoundary (amrex::MultiFab& mf, int, int, amrex::Real time) override
  {
    amrex::Box domain(geom.Domain());

    mf.FillBoundary(geom.periodicity());

    for (amrex::MFIter mfi(mf,true); mfi.isValid(); ++mfi)
      {
	const amrex::Box& box = mfi.tilebox();

	amrex::BaseFab<amrex::Real> &phi_box = mf[mfi];

	for (int i = box.loVect()[0] - mf.nGrow(); i<=box.hiVect()[0] + mf.nGrow(); i++)
	  for (int j = box.loVect()[1] - mf.nGrow(); j<=box.hiVect()[1] + mf.nGrow(); j++)
	    {
	      if (i < domain.loVect()[0]) // Left boundary
		{
		  // if (j > (domain.hiVect()[1]-domain.loVect()[1])/4 && j <  3*(domain.hiVect()[1]-domain.loVect()[1])/4)
		  //   {
		      phi_box(amrex::IntVect(i,j),0) =  0.;
		      phi_box(amrex::IntVect(i,j),1) =  1.;
		  //   }
		  // else
		  //   {
		  //     phi_box(amrex::IntVect(i,j),0) =  phi_box(amrex::IntVect(i+1,j),0);
		  //     phi_box(amrex::IntVect(i,j),1) =  phi_box(amrex::IntVect(i+1,j),1);
		  //   }
		}
	      if (i > domain.hiVect()[0]) // Right boundary
		{
		  phi_box(amrex::IntVect(i,j),0) =  1.;
		  phi_box(amrex::IntVect(i,j),1) =  0.;
		}
	      if (j < domain.loVect()[1]) // Bottom boundary
		{
		  // Nothing
		}
	      if (j > domain.hiVect()[1]) // Top boundary
		{
		  // Nothing
		}

	    }
      }






  }
private:
  amrex::Geometry &geom;
};

#endif
