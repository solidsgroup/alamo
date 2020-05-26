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
	enum Type {Partition, Values};
	TabulatedInterface (amrex::Vector<amrex::Geometry> &_geom) :
		IC(_geom) {}
	
	TabulatedInterface (amrex::Vector<amrex::Geometry> &_geom,
						std::vector<Set::Scalar> a_xs, std::vector<Set::Scalar> a_ys) :
	IC(_geom)
	{
		Define(a_xs,a_ys);
	}

	TabulatedInterface (amrex::Vector<amrex::Geometry> &_geom,
						std::vector<Set::Scalar> a_xs, std::vector<Set::Scalar> a_ys,
						Set::Scalar a_alpha1, Set::Scalar a_alpha2, Type a_type = Type::Values) :
	IC(_geom)
	{
		Define(a_xs,a_ys,a_type,a_alpha1,a_alpha2);
	}

	void Define(std::vector<Set::Scalar> a_xs, std::vector<Set::Scalar> a_ys, 
				Type a_type = Type::Partition, Set::Scalar a_alpha1 = 1.0, Set::Scalar a_alpha2 = 1.0)
	{
		type = a_type;
		alpha1 = a_alpha1;
		alpha2 = a_alpha2;
		xs = a_xs;
		ys = a_ys;
	}
  
	void Add(const int lev,
		 amrex::Vector<amrex::MultiFab * > &a_field)
	{
		Numeric::Interpolator::Linear<Set::Scalar> f(ys,xs);

		for (amrex::MFIter mfi(*a_field[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
		{
			amrex::Box bx = mfi.tilebox();
			bx.grow(a_field[lev]->nGrow());
			amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
				amrex::Real x = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];
				amrex::Real y = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];
				if (y > f(x))
				{
					if (type == Type::Partition)
					{
						field(i,j,k,0) += 0.;     
						field(i,j,k,1) += alpha2;     
					}
					else
						field(i,j,k,0) += alpha2;					
				}
				else
				{
					if (type == Type::Partition)
					{
						field(i,j,k,0) += alpha1;     
						field(i,j,k,1) += 0.;     
					}
					else
						field(i,j,k,0) += alpha1;
				}
			});
		}
	};
  
private:
	Type type;
	std::vector<Set::Scalar> xs;
	std::vector<Set::Scalar> ys;
	Set::Scalar alpha1, alpha2;

public:
	static void Parse(TabulatedInterface & value, amrex::ParmParse & pp)
	{
		pp.queryarr("xs",value.xs);
		pp.queryarr("ys",value.ys);
		value.Define(value.xs,value.ys);
	}
};
}
#endif
