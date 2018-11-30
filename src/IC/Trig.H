#ifndef IC_TRIG_H_
#define IC_TRIG_H_

#include "AMReX_Vector.H"
#include "IC/IC.H"
#include "Util/Util.H"
#include "Set/Set.H"

namespace IC
{
class Trig : public IC
{
public:
	Trig (amrex::Vector<amrex::Geometry> &_geom) :
		IC(_geom)
	{}

	void SetComp(int a_comp) {comp = a_comp;}

	void Initialize(const int lev,
			amrex::Vector<amrex::MultiFab *> &field)
	{
		for (amrex::MFIter mfi(*field[lev],true); mfi.isValid(); ++mfi)
		{
			const amrex::Box& box = mfi.tilebox();

			amrex::BaseFab<amrex::Real> &field_box = (*field[lev])[mfi];

			AMREX_D_TERM(for (int i = box.loVect()[0]-field[lev]->nGrow(); i<=box.hiVect()[0]+field[lev]->nGrow(); i++),
				     for (int j = box.loVect()[1]-field[lev]->nGrow(); j<=box.hiVect()[1]+field[lev]->nGrow(); j++),
				     for (int k = box.loVect()[2]-field[lev]->nGrow(); k<=box.hiVect()[2]+field[lev]->nGrow(); k++))
			{
				// Cell-centered version
				// AMREX_D_TERM(amrex::Real x1 = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];,
				// 	     amrex::Real x2 = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];,
				// 	     amrex::Real x3 = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];);

				AMREX_D_TERM(amrex::Real x1 = geom[lev].ProbLo()[0] + ((amrex::Real)(i)) * geom[lev].CellSize()[0];,
					     amrex::Real x2 = geom[lev].ProbLo()[1] + ((amrex::Real)(j)) * geom[lev].CellSize()[1];,
					     amrex::Real x3 = geom[lev].ProbLo()[2] + ((amrex::Real)(k)) * geom[lev].CellSize()[2];);

				Set::Vector x(AMREX_D_DECL(x1,x2,x3));

				field_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),comp) = std::sin(2*3.14159*x1);// * std::sin(2*3.14159*x2) ;

			}
		}

	};


	void Initialize(const int lev,
			amrex::Vector<std::unique_ptr<amrex::MultiFab> > &field)
	{
		amrex::Vector<amrex::MultiFab *> pfield = amrex::GetVecOfPtrs(field);
		Initialize(lev,pfield);
	};
	void Initialize(const int lev,
			amrex::Vector<amrex::MultiFab> &field)
	{
		amrex::Vector<amrex::MultiFab *> pfield = amrex::GetVecOfPtrs(field);
		Initialize(lev,pfield);
	};

private:
	int comp;
};
}
#endif
