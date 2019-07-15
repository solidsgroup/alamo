#ifndef IC_TRIG2_H_
#define IC_TRIG2_H_

#include <complex>

#include "AMReX_Vector.H"
#include "IC/IC.H"
#include "Util/Util.H"
#include "Set/Set.H"

namespace IC
{
/// \brief Initialize using a trigonometric series
class Trig2 : public IC
{
public:
	Trig2 (amrex::Vector<amrex::Geometry> &_geom,
	      Set::Scalar _alpha = 1.0,
	      AMREX_D_DECL(Set::Scalar _theta1 = 0,
			   Set::Scalar _theta2 = 0,
			   Set::Scalar _theta3 = 0),
	      AMREX_D_DECL(int _n1 = 0,
			   int _n2 = 0,
			   int _n3 = 0),
	      int _dim = AMREX_SPACEDIM) :
	      IC(_geom)
	{
		Define(_alpha,
		       AMREX_D_DECL(_theta1,_theta2,_theta3),
		       AMREX_D_DECL(_n1,_n2,_n3),
		       _dim);
	}


	void Define(Set::Scalar _alpha = 1.0,
		    AMREX_D_DECL(Set::Scalar _theta1 = 0,
				 Set::Scalar _theta2 = 0,
				 Set::Scalar _theta3 = 0),
		    AMREX_D_DECL(int _n1 = 0,
				 int _n2 = 0,
				 int _n3 = 0),
		    int _dim = AMREX_SPACEDIM)
	{
		alpha = _alpha;
		AMREX_D_TERM(theta1 = _theta1;, theta2 = _theta2;, theta3 = _theta3;);
		AMREX_D_TERM(n1 = _n1;, n2 = _n2;, n3 = _n3;);
		dim = _dim;
	}


	void Add(const int lev,
		 amrex::Vector<amrex::MultiFab *> &field)
	{
		bool cellcentered = (field[0]->boxArray().ixType() == amrex::IndexType(amrex::IntVect::TheCellVector()));

		const std::complex<Set::Scalar> I(0.0,1.0);

		const amrex::Real AMREX_D_DECL(L1 = geom[lev].ProbHi()[0] - geom[lev].ProbLo()[0],
					       L2 = geom[lev].ProbHi()[1] - geom[lev].ProbLo()[1],
					       L3 = geom[lev].ProbHi()[2] - geom[lev].ProbLo()[2]);
		for (amrex::MFIter mfi(*field[lev],true); mfi.isValid(); ++mfi)
		{
			const amrex::Box& box = mfi.validbox();
			amrex::BaseFab<amrex::Real> &field_box = (*field[lev])[mfi];

			AMREX_D_TERM(for (int i = box.loVect()[0] - field[lev]->nGrow(); i<=box.hiVect()[0] + field[lev]->nGrow(); i++),
				     for (int j = box.loVect()[1] - field[lev]->nGrow(); j<=box.hiVect()[1] + field[lev]->nGrow(); j++),
				     for (int k = box.loVect()[2] - field[lev]->nGrow(); k<=box.hiVect()[2] + field[lev]->nGrow(); k++))
			{

				Set::Scalar AMREX_D_DECL(x1,x2,x3);
				if (cellcentered)
				{
					AMREX_D_TERM(x1 = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];,
						     x2 = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];,
						     x3 = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];);
				}
				else
				{
					AMREX_D_TERM(x1 = geom[lev].ProbLo()[0] + ((amrex::Real)(i)) * geom[lev].CellSize()[0];,
						     x2 = geom[lev].ProbLo()[1] + ((amrex::Real)(j)) * geom[lev].CellSize()[1];,
						     x3 = geom[lev].ProbLo()[2] + ((amrex::Real)(k)) * geom[lev].CellSize()[2];);
				}

				
				Set::Scalar trigfn = 1.0;
#if AMREX_SPACEDIM > 0
				if (dim > 0) trigfn *= std::exp(I*(theta1 + ((amrex::Real)n1)*x1*Set::Constant::Pi / L1)).real();

#endif
#if AMREX_SPACEDIM > 1
 				if (dim > 1) trigfn *= std::exp(I*(theta2 + ((amrex::Real)n2)*x2*Set::Constant::Pi / L2)).real();
#endif
#if AMREX_SPACEDIM > 2
 				if (dim > 2) trigfn *= std::exp(I*(theta3 + ((amrex::Real)n3)*x3*Set::Constant::Pi / L3)).real();
#endif
				field_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),comp) += alpha * trigfn;

			}
		}
	}
	using IC::Add;

private:
	int dim = AMREX_SPACEDIM;
	Set::Scalar alpha;
	int AMREX_D_DECL(n1, n2, n3);
	Set::Scalar AMREX_D_DECL(theta1, theta2, theta3);
};
}
#endif
