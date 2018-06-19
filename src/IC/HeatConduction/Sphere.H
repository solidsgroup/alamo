#ifndef PerturbedInterface_H_
#define PerturbedInterface_H_

#include "IC/IC.H"

/// \class Cylinder
/// \brief Prescribe a cylindrical (3D) / circular (2D) region at a certain temperature.
///
namespace HeatConduction
{
namespace IC
{
class Sphere : public IC
{
public:
  PerturbedInterface (amrex::Array<amrex::Geometry> &_geom) :
    IC::IC(_geom)
  {
      amrex::ParmParse pp("ic");
      pp.queryarr("Tin",Tin);
      pp.queryarr("Tout",Tout);
  }
  
void Initialize(const int lev,
		amrex::Array<std::unique_ptr<amrex::MultiFab> > &field)
  {
    const amrex::Real L1 = geom[lev].ProbHi()[0] - geom[lev].ProbLo()[0];
    const amrex::Real L2 = geom[lev].ProbHi()[1] - geom[lev].ProbLo()[1];

    for (amrex::MFIter mfi(*field[lev],true); mfi.isValid(); ++mfi)
      {
	const amrex::Box& box = mfi.tilebox();

	amrex::BaseFab<Real> &field_box = (*field[lev])[mfi];

	for (int i = box.loVect()[0]-field[lev]->nGrow(); i<=box.hiVect()[0]+field[lev]->nGrow(); i++) 
	  for (int j = box.loVect()[1]-field[lev]->nGrow(); j<=box.hiVect()[1]+field[lev]->nGrow(); j++)
#if BL_SPACEDIM==3
	    for (int k = box.loVect()[2]-field[lev]->nGrow(); k<=box.hiVect()[2]+field[lev]->nGrow(); k++)
#endif
	      {
		amrex::real

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
  amrex::Real Tin, Tout;
};
}
}
#endif
