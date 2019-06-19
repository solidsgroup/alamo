#ifndef IC_RANDOM_H_
#define IC_RANDOM_H_

#include "AMReX_Vector.H"

#include "IC/IC.H"
#include "Util/Util.H"
#include "Set/Set.H"

namespace IC
{
/// \brief Set each point to a random value.
class Random : public IC
{
public:
	Random (amrex::Vector<amrex::Geometry> &_geom, Set::Scalar a_mult = 1.0) :
		IC(_geom), mult(a_mult)
	{}
  
	void Add(const int lev,
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
				field_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),comp) += mult * Util::Random();
			}
		}

	};
	using IC::Add;
private:
	Set::Scalar mult;

};
}
#endif
