#ifndef IC_NOTCH_H_
#define IC_NOTCH_H_

#include "Set/Set.H"
#include "IC/IC.H"

// Note: right now this is meant for 2D. We need to rethink this implementation for 3D.
namespace IC
{
class Notch : public IC
{
public:
    enum Mollifier {Dirac, Gaussian};

	Notch (amrex::Vector<amrex::Geometry> &_geom) : IC(_geom) 
    {
        nCenter = Set::Vector::Zero();
        nOrientation = Set::Vector::Random();
        nThickness = 0.01;
        nLength = 0.1;
        moll = Mollifier::Dirac;
        eps = 1e-5;

        amrex::ParmParse pp("ic");

        amrex::Vector<Set::Scalar> center;
		pp.queryarr("center",center);
        if(center.size()>=AMREX_SPACEDIM){AMREX_D_TERM(nCenter(0) = center[0];,nCenter(1) = center[1];,nCenter(2) = center[2];);}

        amrex::Vector<Set::Scalar> orientation;
		pp.queryarr("orientation",orientation);
        if(orientation.size()>=AMREX_SPACEDIM){AMREX_D_TERM(nOrientation(0) = orientation[0];,nOrientation(1) = orientation[1];,nOrientation(2) = orientation[2];);}
        if(nOrientation.lpNorm<2>() <= 0.) nOrientation = Set::Vector::Random();

        pp.query("thickness", nThickness);
        pp.query("length", nLength);
        pp.query("eps", eps);

        if(nThickness <= 0.) nThickness = 0.01;
        if(nLength <= 0.) nLength = 0.1;
        if(eps <= 0.) eps = 1.e-5;

        std::string mollifier;
        pp.query("mollifier",mollifier);
        if(mollifier == "Dirac" || mollifier == "dirac") 
            moll = Mollifier::Dirac;
        else 
            moll = Mollifier::Gaussian;

        nOrientation = nOrientation / nOrientation.lpNorm<2>();
        eps_sq = eps*eps;
        window = 6.0*eps;

        if(nOrientation(0) != 0.)
        {
            AMREX_D_TERM(nNormal(0) = 1.;, nNormal(1) = 1.;, nNormal(2) = 1.;);
            nNormal(0) = -(AMREX_D_TERM(0.,+nOrientation(1),+nOrientation(2)))/nOrientation(0);
            nNormal = nNormal/nNormal.lpNorm<2>();
        }
        else if(nOrientation(1) != 0.)
        {
            AMREX_D_TERM(nNormal(0) = 1.;, nNormal(1) = 1.;, nNormal(2) = 1.;);
            nNormal(1) = -(AMREX_D_TERM(nOrientation(0), + 0.0, + nOrientation(2)))/nOrientation(1);
            nNormal = nNormal/nNormal.lpNorm<2>();
        }
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

                if(std::abs((x-nCenter).transpose()*nOrientation) <= nLength)
                {
                    Set::Scalar t = std::abs((x-nCenter).transpose()*nNormal) - (nThickness/2.0);
                    field(i,j,k) = 0.5 + 0.5*std::erf(t / eps);
                }

                Set::Vector nLeft =  nCenter - 0.5*nLength*nOrientation;
                Set::Vector nRight = nCenter + 0.5*nLength*nOrientation;
                Set::Scalar distLeft = (x-nLeft).transpose()*nOrientation;
                Set::Scalar distRight = (x-nRight).transpose()*nOrientation;

                if(distLeft < 0.)
                    field(i,j,k) = 0.5 + 0.5*std::erf( ((x-nLeft).lpNorm<2>() - nThickness/2.0) / eps );
                if(distRight > 0.)
                    field(i,j,k) = 0.5 + 0.5*std::erf( ((x-nRight).lpNorm<2>() - nThickness/2.0) / eps );
            

//                field(i,j,k) = 0.5 - 0.5*std::erf()
//                Set::Scalar norm = (A*(x-x0)).lpNorm<2>();
//                field(i,j,k) = 0.5 - 0.5*std::erf(((x-x0).transpose() * A * (x-x0) - 1.0) / eps / norm);
//
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
    Set::Vector nCenter, nOrientation, nNormal;
    Set::Scalar eps = 1.e-5, eps_sq = eps*eps, nThickness, nLength;
    Set::Scalar window = 6.0*eps;
    Mollifier moll = Mollifier::Dirac;
    bool nRad = true;
};
}
#endif