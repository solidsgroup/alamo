#ifndef IC_CIRCLE_H_
#define IC_CIRCLE_H_

#include "IC/IC.H"
#include "Util/Util.H"

/// \class Circle
/// \brief Initialize a circular inclusion
namespace IC
{
class Circle : public IC
{
public:
	Circle (amrex::Vector<amrex::Geometry> &_geom) :
		IC(_geom)
	{
		amrex::ParmParse pp("ic");
		pp.queryarr("center",center);
		pp.query("radius",radius);
	}
  
	void Add(const int lev,
			amrex::Vector<amrex::MultiFab * > &field)
	{
		for (amrex::MFIter mfi(*field[lev],true); mfi.isValid(); ++mfi)
		{
			const amrex::Box& box = mfi.tilebox();

			amrex::BaseFab<amrex::Real> &field_box = (*field[lev])[mfi];

			AMREX_D_TERM(for (int i = box.loVect()[0]-field[lev]->nGrow(); i<=box.hiVect()[0]+field[lev]->nGrow(); i++),
				     for (int j = box.loVect()[1]-field[lev]->nGrow(); j<=box.hiVect()[1]+field[lev]->nGrow(); j++),
				     for (int k = box.loVect()[2]-field[lev]->nGrow(); k<=box.hiVect()[2]+field[lev]->nGrow(); k++))
					{
						amrex::Real x = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];
						amrex::Real y = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];

						if ( ((x-center[0])*(x-center[0]) + (y-center[1])*(y-center[1])) < radius*radius)
						{
							field_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),0) = 1.;     
							field_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),1) = 0.;     
						}
						else
						{
							field_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),0) = 0.;     
							field_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),1) = 1.;     
						}
					}
		}

	};
  
private:
	amrex::Vector<amrex::Real> center;
	amrex::Real radius;
};
}
#endif
