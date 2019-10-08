#ifndef IC_ELLIPSOID_H_
#define IC_ELLIPSOID_H_

#include "IC/IC.H"
#include "Util/Util.H"
#include "IO/ParmParse.H"

/// \class Ellipsoid
/// \brief Initialize an ellipsoidal inclusion
namespace IC
{
class Ellipsoid : public IC
{
public:
	enum Mollifier {Dirac, Gaussian};

	Ellipsoid (amrex::Vector<amrex::Geometry> &_geom) :
		IC(_geom)
	{
		IO::ParmParse pp("ic");
		
		if(pp.contains("center")) pp.queryarr("center",center);
		
		if (pp.contains("A")) pp.queryarr("A",A);
		if (pp.contains("radius"))
		{ 
			Set::Vector a_radius; 
			pp.queryarr("radius",a_radius);  
			AMREX_D_TERM( A(0,0) = 1./a_radius[0];, A(1,1) = 1./a_radius[1];, A(2,2) = 1./a_radius[2]; );
		}

		pp.query("eps", eps);
		if(eps <= 0.) eps = 1.e-5;

		std::string mollifier;
        pp.query("mollifier",mollifier);
        if(mollifier == "Dirac" || mollifier == "dirac") 
            moll = Mollifier::Dirac;
        else 
            moll = Mollifier::Gaussian;
	}
  
	void Add(const int lev,
			amrex::Vector<amrex::MultiFab * > &a_field)
	{
		bool cellcentered = (a_field[0]->boxArray().ixType() == amrex::IndexType(amrex::IntVect::TheCellVector()));

		for (amrex::MFIter mfi(*a_field[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
		{
			amrex::Box bx = mfi.tilebox();
			bx.grow(a_field[lev]->nGrow());

			amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {

					Set::Vector x;
					// NODE
					if (a_field[0]->ixType() == amrex::IndexType::TheNodeType())
					{
						AMREX_D_TERM(x(0) = geom[lev].ProbLo()[0] + ((amrex::Real)(i)) * geom[lev].CellSize()[0];,
									x(1) = geom[lev].ProbLo()[1] + ((amrex::Real)(j)) * geom[lev].CellSize()[1];,
									x(2) = geom[lev].ProbLo()[2] + ((amrex::Real)(k)) * geom[lev].CellSize()[2];);
					}
					else if (a_field[0]->ixType() == amrex::IndexType::TheCellType())
					{
						AMREX_D_TERM(x(0) = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];,
									x(1) = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];,
									x(2) = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];);
					}

					Set::Scalar norm = (A*(x-center)).lpNorm<2>();
                	field(i,j,k) = 0.5 - 0.5*std::erf(((x-center).transpose() * A * (x-center) - 1.0) / eps / norm);
				});
		}
		a_field[lev]->FillBoundary();
	};
  
private:
	int dim = AMREX_SPACEDIM;
	Set::Vector center = Set::Vector::Zero();
	Set::Scalar eps = 0.01, eps_sq = eps*eps;
	Set::Scalar window = 6.0*eps;
	Set::Matrix A = Set::Matrix::Identity();
	Mollifier moll = Mollifier::Dirac;
};
}
#endif
