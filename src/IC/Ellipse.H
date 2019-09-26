#ifndef IC_ELLIPSE_H_
#define IC_ELLIPSE_H_

#include "Set/Set.H"
#include "IC/IC.H"

namespace IC
{
class Ellipse : public IC
{
public:
    enum Mollifier {Dirac, Gaussian};

	Ellipse (amrex::Vector<amrex::Geometry> &_geom) : IC(_geom) 
    {
        Define(Set::Matrix::Identity(), Set::Vector::Zero(), Mollifier::Dirac, 0.0);
    }

	Ellipse (amrex::Vector<amrex::Geometry> &_geom, Set::Matrix a_A, Set::Vector a_x0, Mollifier a_moll, Set::Scalar a_eps) : IC(_geom) 
    {
        Define(a_A,a_x0,a_moll,a_eps);
    }

	Ellipse (amrex::Vector<amrex::Geometry> &_geom, Set::Vector a_a, Set::Vector a_x0, Mollifier a_moll, Set::Scalar a_eps) : IC(_geom) 
    {
        Set::Matrix A = Set::Matrix::Zero();
        for (int d = 0; d < AMREX_SPACEDIM; d++) A(d,d) = 1./a_a(d)/a_a(d);
        Define(A,a_x0,a_moll,a_eps);
    }

  	void Define (Set::Matrix a_A, Set::Vector a_x0, Mollifier a_moll, Set::Scalar a_eps)
	{
        moll = a_moll;
        A = a_A;
        x0 = a_x0; 
        eps = a_eps;
        eps_sq = eps*eps;
        window = 6.0*eps;
  	}

	
	void Add(const int lev, amrex::Vector<amrex::MultiFab * > &a_field)
	{
        Set::Vector DX(geom[lev].CellSize());

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


                Set::Scalar norm = (A*(x-x0)).lpNorm<2>();
                field(i,j,k) = 0.5 - 0.5*std::erf(((x-x0).transpose() * A * (x-x0) - 1.0) / eps / norm);

//                if (moll == Mollifier::Dirac)
//                {
//                    if ((x-x0).transpose() * A * (x-x0) < 1.0) field(i,j,k) += 1.0;
//                }
//                else if (moll == Mollifier::Gaussian)
//                {
//                    Set::Scalar conv = 0.0;
//                    Set::Vector _x;
//                    AMREX_D_TERM(for (_x(0) = x(0) - window; _x(0) <= x(0) + window; _x(0) += DX(0)),
//                                 for (_x(1) = x(1) - window; _x(1) <= x(1) + window; _x(1) += DX(1)),
//                                 for (_x(2) = x(2) - window; _x(2) <= x(2) + window; _x(2) += DX(2)))
//                        {
//#if AMREX_SPACEDIM==2
//                            if ((_x-x0).transpose() * A * (_x-x0) < 1.0)
//                            {
//                                Set::Scalar mollifier = std::exp( - (_x - x).squaredNorm() / eps_sq) / (eps_sq * Set::Constant::Pi);
//                                conv += mollifier * DX(0) * DX(1);
//                            }
//#elif AMREX_SPACEDIM==3
//                            if ((_x-x0).transpose() * A * (_x-x0) < 1.0)
//                            {
//                                Set::Scalar mollifier = std::exp( - (_x - x).squaredNorm() / eps_sq) / pow(eps_sq * Set::Constant::Pi,1.5);
//                                conv += mollifier * DX(0) * DX(1) * DX(2);
//                            }
//#endif
//                        }
//                    field(i,j,k) += conv;
//                }

			});
		}
        a_field[lev]->FillBoundary();
	}
	
private:
    Set::Vector x0;
    Set::Scalar eps = 0.01, eps_sq = eps*eps;
    Set::Scalar window = 6.0*eps;
    Set::Matrix A = Set::Matrix::Identity();    
    Mollifier moll = Mollifier::Dirac;
};
}
#endif
