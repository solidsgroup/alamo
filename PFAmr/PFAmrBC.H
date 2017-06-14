#ifndef PFAmrBC_H_
#define PFAmrBC_H_

#include <AMReX_PhysBCFunct.H>

// In this tutorial, we assume the boundaries are all periodic.  
// Thus there are no physical boundaries. And FillBoundary does nothing.

class PFAmrPhysBC
  : public amrex::PhysBCFunctBase
{
public:
  PFAmrPhysBC (amrex::Geometry &_geom) : geom(_geom)  {}
  virtual ~PFAmrPhysBC () {}
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
		  // Nothing
		}
	      if (i > domain.hiVect()[0]) // Right boundary
		{
		  // Nothing
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
