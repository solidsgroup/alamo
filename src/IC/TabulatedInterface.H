#ifndef IC_TABULATEDINTERFACE_H_
#define IC_TABULATEDINTERFACE_H_

#include "IC/IC.H"
#include "Util/Util.H"
#include "Numeric/Interpolator/Linear.H"

/// \class TabulatedInterface
/// \brief Initialize a perturbed interface using a linear interpolation
///
namespace IC
{
class TabulatedInterface : public IC
{
public:
	TabulatedInterface (amrex::Vector<amrex::Geometry> &_geom) :
		IC(_geom)
	{
		amrex::ParmParse pp("ic");
		pp.queryarr("xs",xs);
		pp.queryarr("ys",ys);

		if (xs.size() != ys.size())
			Util::Abort("Number of x points must equal number of y points");

	}
  
	void Initialize(const int lev,
			amrex::Vector<std::unique_ptr<amrex::MultiFab> > &field)
	{
		//const amrex::Real width = geom[lev].ProbLo()[0] - geom[lev].ProbHi()[0];

		Numeric::Interpolator::Linear<Set::Scalar> f(ys,xs);

		for (amrex::MFIter mfi(*field[lev],true); mfi.isValid(); ++mfi)
		{
			const amrex::Box& box = mfi.tilebox();

			amrex::BaseFab<amrex::Real> &field_box = (*field[lev])[mfi];

			for (int i = box.loVect()[0]-field[lev]->nGrow(); i<=box.hiVect()[0]+field[lev]->nGrow(); i++) 
				for (int j = box.loVect()[1]-field[lev]->nGrow(); j<=box.hiVect()[1]+field[lev]->nGrow(); j++)
#if BL_SPACEDIM==3
					for (int k = box.loVect()[2]-field[lev]->nGrow(); k<=box.hiVect()[2]+field[lev]->nGrow(); k++)
#endif
					{
						amrex::Real x = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];
						amrex::Real y = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];


						if (y > f(x))
						{
							field_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),0) = 0.;     
							field_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),1) = 1.;     
						}
						else
						{
							field_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),0) = 1.;     
							field_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),1) = 0.;     
						}
					}
		}

	};
  
private:
	amrex::Vector<Set::Scalar> xs;
	amrex::Vector<Set::Scalar> ys;
};
}
#endif
