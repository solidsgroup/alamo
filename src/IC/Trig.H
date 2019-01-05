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
	Trig (amrex::Vector<amrex::Geometry> &_geom,
	      Set::Scalar _alpha = 1.0,
	      AMREX_D_DECL(int _n1 = 1, int _n2 = 1, int _n3 = 1)) :
	      IC(_geom),
	      alpha(_alpha),
	      AMREX_D_DECL(n1(_n1),n2(_n2),n3(_n3))
	{}

	virtual void SetComp(int a_comp) {comp = a_comp;}

	void Initialize(const int lev,
			amrex::Vector<amrex::MultiFab *> &field)
	{
		const amrex::Real AMREX_D_DECL(L1 = geom[lev].ProbLo()[0] - geom[lev].ProbHi()[0],
					       L2 = geom[lev].ProbLo()[1] - geom[lev].ProbHi()[1],
					       L3 = geom[lev].ProbLo()[2] - geom[lev].ProbHi()[2]);
				   

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

				field_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),comp)
					= alpha* AMREX_D_TERM(std::sin(n1*Set::Constant::Pi*x1 / L1),
							      * std::sin(n2*Set::Constant::Pi*x2 / L2),
							      * std::sin(n3*Set::Constant::Pi*x3 / L3));

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
	Set::Scalar alpha;
	int AMREX_D_DECL(n1, n2, n3);
};
}
#endif
