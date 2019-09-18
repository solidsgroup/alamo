#ifndef IC_AFFINE_H_
#define IC_AFFINE_H_

#include "AMReX_Vector.H"

#include "IC/IC.H"
#include "Util/Util.H"
#include "Set/Set.H"

namespace IC
{
/// \brief Initialize with the function \f$\alpha[(\mathbf{x}-\mathbf{b})\cdot\mathbf{n}]^m\f$
class Affine : public IC
{
public:
	Affine (amrex::Vector<amrex::Geometry> &_geom) :
		IC(_geom)
	{
		amrex::ParmParse pp("ic");
		amrex::Vector<amrex::Real> ntmp; pp.queryarr("n",ntmp);
		for (int i = 0; i < AMREX_SPACEDIM; i++) n(i) = ntmp[i];
		Util::Warning(INFO,"This function is not complete");
		pp.query("alpha",alpha);
	}
  
	Affine (amrex::Vector<amrex::Geometry> &_geom,
		Set::Vector a_n,
		Set::Scalar a_alpha,
		Set::Vector a_b = Set::Vector::Zero(),
		bool a_halfspace = false,
		Set::Scalar a_m = 1.0) :
		IC(_geom), n(a_n), alpha(a_alpha), b(a_b), halfspace(a_halfspace), m(a_m)
	{}

	void SetHalfspace(bool a_halfspace) {halfspace = a_halfspace;}

	void Add(const int lev,
		 amrex::Vector<amrex::MultiFab *> &field)
	{
		bool cellcentered = (field[0]->boxArray().ixType() == amrex::IndexType(amrex::IntVect::TheCellVector()));

		for (amrex::MFIter mfi(*field[lev],true); mfi.isValid(); ++mfi)
		{
			const amrex::Box& box = mfi.tilebox();

			amrex::BaseFab<amrex::Real> &field_box = (*field[lev])[mfi];

			AMREX_D_TERM(for (int i = box.loVect()[0]-field[lev]->nGrow(); i<=box.hiVect()[0]+field[lev]->nGrow(); i++),
				     for (int j = box.loVect()[1]-field[lev]->nGrow(); j<=box.hiVect()[1]+field[lev]->nGrow(); j++),
				     for (int k = box.loVect()[2]-field[lev]->nGrow(); k<=box.hiVect()[2]+field[lev]->nGrow(); k++))
			{
				amrex::Real AMREX_D_DECL(x1,x2,x3);
				
				if (cellcentered)
				{
					AMREX_D_TERM(x1 = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];,
						     x2 = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];,
						     x3 = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];);
				}
				{
					AMREX_D_TERM(x1 = geom[lev].ProbLo()[0] + ((amrex::Real)(i)) * geom[lev].CellSize()[0];,
						     x2 = geom[lev].ProbLo()[1] + ((amrex::Real)(j)) * geom[lev].CellSize()[1];,
						     x3 = geom[lev].ProbLo()[2] + ((amrex::Real)(k)) * geom[lev].CellSize()[2];);
				}

				Set::Vector x(AMREX_D_DECL(x1,x2,x3));

				field_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),comp) = (x - b).dot(n);
				Set::Scalar val = (x - b).dot(n);
				if (halfspace)
				{
					if (val >= 0)
						field_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),comp) = alpha*std::pow(val,m);
					else
						field_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),comp) = 0.0;
				}
			}
		}

	};

private:
	Set::Vector n;
	amrex::Real alpha;
	Set::Vector b;
	bool halfspace = false;
	amrex::Real m = 1.0;
};
}
#endif
