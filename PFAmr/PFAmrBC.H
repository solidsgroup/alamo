#ifndef PFAmrBC_H_
#define PFAmrBC_H_

#include <AMReX_PhysBCFunct.H>

// In this tutorial, we assume the boundaries are all periodic.  
// Thus there are no physical boundaries. 

class PFAmrPhysBC
  : public amrex::PhysBCFunctBase
{
public:
  PFAmrPhysBC (amrex::Geometry &_geom, int *bc_lo, int *bc_hi) : geom(_geom)  {}
  virtual ~PFAmrPhysBC () {}
  virtual void FillBoundary (amrex::MultiFab& mf, int, int, amrex::Real time) override
  {
    amrex::Box domain(geom.Domain());

    mf.FillBoundary(geom.periodicity());

    for (amrex::MFIter mfi(mf,true); mfi.isValid(); ++mfi)
      {
	const amrex::Box& box = mfi.tilebox();

	amrex::BaseFab<amrex::Real> &phi_box = mf[mfi];

	bool print = true;
	for (int i = box.loVect()[0] - mf.nGrow(); i<=box.hiVect()[0] + mf.nGrow(); i++)
	  {
	    for (int j = box.loVect()[1] - mf.nGrow(); j<=box.hiVect()[1] + mf.nGrow(); j++)
	      for (int n = 0; n < mf.nComp(); n++)
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
		      if (n<2)
			{
			  if (j == domain.loVect()[1]-1) phi_box(amrex::IntVect(i,j),n) = (amrex::Real)(1-n);//phi_box(amrex::IntVect(i,j+1),n); //1.0;//
			  if (j == domain.loVect()[1]-2) phi_box(amrex::IntVect(i,j),n) = (amrex::Real)(1-n);//phi_box(amrex::IntVect(i,j+3),n); //1.0;//
			}
		    }
		  if (j > domain.hiVect()[1]) // Top boundary
		    {
		      if (n<2)
			{
			  if (j == domain.hiVect()[1]+1) phi_box(amrex::IntVect(i,j),n) = (amrex::Real)(n);//phi_box(amrex::IntVect(i,j-1),n); //0.0;//
			  if (j == domain.hiVect()[1]+2) phi_box(amrex::IntVect(i,j),n) = (amrex::Real)(n);//phi_box(amrex::IntVect(i,j-3),n); //0.0;//
			}
		    }
		}
	  }
      }
  }
private:
  amrex::Geometry &geom;
};

#endif
