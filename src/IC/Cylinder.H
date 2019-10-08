#ifndef HEATCONDUCTION_INITIALCONDITION_CYLINDER_H_
#define HEATCONDUCTION_INITIALCONDITION_CYLINDER_H_

#include <AMReX_ParallelDescriptor.H>

#include "IC/IC.H"

namespace IC
{
class Cylinder : public IC
{
public:
  Cylinder (amrex::Vector<amrex::Geometry> &_geom) : IC(_geom) {}
  /// Constructor requires temperature inputs for either side of the boundary.
  /// The radius is defaulted to 1, and the cylinder is centered on the axis.
  Cylinder (amrex::Vector<amrex::Geometry> &_geom, Set::Scalar _Tin, Set::Scalar _Tout) : IC(_geom) 
  {
	  Define(_Tin,_Tout);
  }

void Define (Set::Scalar _Tin, Set::Scalar _Tout)
  {
	Tin = _Tin;
	Tout = _Tout;
  }

virtual void Add(const int lev,
		amrex::Vector<amrex::MultiFab * > &field)
  {
    const amrex::Real L1 = geom[lev].ProbHi()[0] - geom[lev].ProbLo()[0];
    const amrex::Real L2 = geom[lev].ProbHi()[1] - geom[lev].ProbLo()[1];

    for (amrex::MFIter mfi(*field[lev],true); mfi.isValid(); ++mfi)
      {
	const amrex::Box& box = mfi.tilebox();

	amrex::FArrayBox &field_box = (*field[lev])[mfi];

	for (int i = box.loVect()[0]-field[lev]->nGrow(); i<=box.hiVect()[0]+field[lev]->nGrow(); i++) 
	  for (int j = box.loVect()[1]-field[lev]->nGrow(); j<=box.hiVect()[1]+field[lev]->nGrow(); j++)
#if BL_SPACEDIM==3
	    for (int k = box.loVect()[2]-field[lev]->nGrow(); k<=box.hiVect()[2]+field[lev]->nGrow(); k++)
#endif
	      {
		amrex::Real x = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];
		amrex::Real y = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];

		amrex::Real eX = (x - (geom[lev].ProbLo()[0] + 0.5*L1))/(0.25*L1);
		amrex::Real eY = (y - (geom[lev].ProbLo()[1] + 0.5*L2))/(0.25*L2);


		if (eX*eX + eY*eY < 1.0)
		  field_box(amrex::IntVect(AMREX_D_DECL(i,j,k))) = Tin;     
		else
		  field_box(amrex::IntVect(AMREX_D_DECL(i,j,k))) = Tout;     
	      }
      }

  };
  
private:
	Set::Scalar Tin, Tout;
};
}
#endif
