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

	Ellipsoid (amrex::Vector<amrex::Geometry> &_geom) : IC(_geom) {}
  
	void Add(const int &lev,Set::Field<Set::Scalar> &a_field)
	{
		for (amrex::MFIter mfi(*a_field[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
		{
			amrex::Box bx = mfi.tilebox();
			bx.grow(a_field[lev]->nGrow());
			amrex::IndexType type = a_field[lev]->ixType();

			amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {

					Set::Vector x;
					// NODE
					if (type == amrex::IndexType::TheNodeType())
					{
						AMREX_D_TERM(x(0) = geom[lev].ProbLo()[0] + ((amrex::Real)(i)) * geom[lev].CellSize()[0];,
									x(1) = geom[lev].ProbLo()[1] + ((amrex::Real)(j)) * geom[lev].CellSize()[1];,
									x(2) = geom[lev].ProbLo()[2] + ((amrex::Real)(k)) * geom[lev].CellSize()[2];);
					}
					else if (type == amrex::IndexType::TheCellType())
					{
						AMREX_D_TERM(x(0) = geom[lev].ProbLo()[0] + ((amrex::Real)(i) + 0.5) * geom[lev].CellSize()[0];,
									x(1) = geom[lev].ProbLo()[1] + ((amrex::Real)(j) + 0.5) * geom[lev].CellSize()[1];,
									x(2) = geom[lev].ProbLo()[2] + ((amrex::Real)(k) + 0.5) * geom[lev].CellSize()[2];);
					}

					Set::Scalar min_value = 100.0;
					for (int i=0 ; i < center.size(); i++)
					{
						Set::Scalar value = 0.;
						Set::Scalar norm = (A[i]*(x-center[i])).lpNorm<2>();
						value = 0.5 + 0.5*std::erf(((x-center[i]).transpose() * A[i] * (x-center[i]) - 1.0) / eps[i] / norm);
						min_value = value < min_value ? value : min_value;
					}
					field(i,j,k) = min_value;
					if (field(i,j,k) < 0.) field(i,j,k) = 0.;
                	if (field(i,j,k) > 1.) field(i,j,k) = 1.;
				});
		}
		a_field[lev]->FillBoundary();
	};
  
private:
	int dim = AMREX_SPACEDIM;
	amrex::Vector<Set::Vector> center; 
	amrex::Vector<Set::Scalar> eps;
	amrex::Vector<Set::Matrix> A;
	Mollifier moll = Mollifier::Dirac;

public:
    static void Parse(Ellipsoid & value, IO::ParmParse & pp)
    {
		amrex::Vector<Set::Scalar> center;
        if(pp.contains("center")) pp.queryarr("center",center);

		if(center.size() < AMREX_SPACEDIM) value.center.push_back(Set::Vector::Zero());
		else
		{
			for (int i = 0; i<center.size(); i+=AMREX_SPACEDIM)
				value.center.push_back(Set::Vector(AMREX_D_DECL(center[i],center[i+1],center[i+2])));
		}
		
		amrex::Vector<Set::Scalar> A;
		if (pp.contains("A"))
		{
			pp.queryarr("A", A);
			if(A.size() < AMREX_SPACEDIM*AMREX_SPACEDIM) value.A.push_back(Set::Matrix::Identity());
			else
			{
				int i=0, j=0;
				Set::Matrix temp = Set::Matrix::Zero();
				for (i=0; i<A.size(); i+=AMREX_SPACEDIM)
				{
					for (j=0; j<AMREX_SPACEDIM; j++)
						temp(i,j) = A[i+j];
					
					if((i+j) % (AMREX_SPACEDIM*AMREX_SPACEDIM) == 0) 
					{
						value.A.push_back(temp);
						temp = Set::Matrix::Zero();
					}
				}
			}
		}
		if (pp.contains("radius"))
		{ 
			amrex::Vector<Set::Scalar> a_radius; 
			pp.queryarr("radius",a_radius);  

			if (a_radius.size() < AMREX_SPACEDIM) value.A.push_back(Set::Matrix::Identity());
			else
			{
				for (int i = 0; i<a_radius.size(); i+=AMREX_SPACEDIM)
				{
					Set::Matrix temp = Set::Matrix::Zero();
					AMREX_D_TERM( temp(0,0) = 1./(a_radius[i]*a_radius[i]) ;, temp(1,1) = 1./(a_radius[i+1]*a_radius[i+1]) ;, temp(2,2) = 1./(a_radius[i+2]*a_radius[i+2]) ; );
					value.A.push_back(temp);
				}
			}
		}

		amrex::Vector<Set::Scalar> eps;
		pp.queryarr("eps", eps);

		if(eps.size() < 1) value.eps.push_back(1.e-5);
		for (int i =0; i < value.A.size(); i++)
		{
			if (eps[i] <= 0.0) eps[i] = 1.e-5;
			value.eps.push_back(eps[i]);
		}
		//if(value.eps <= 0.) value.eps = 1.e-5;

		std::string mollifier;
        pp.query("mollifier",mollifier);
        if(mollifier == "Dirac" || mollifier == "dirac") 
            value.moll = Mollifier::Dirac;
        else 
            value.moll = Mollifier::Gaussian;
    }
};
}
#endif
