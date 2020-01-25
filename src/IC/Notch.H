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

                if(std::abs((x-nCenter).transpose()*nOrientation) <= 0.5*nLength)
                {
                    Set::Scalar t = std::abs((x-nCenter).transpose()*nNormal) - (nThickness/2.0);
                    field(i,j,k) = 0.5 + 0.5*std::erf(t / eps);
                }
                else
                {
                    Set::Vector nLeft =  nCenter - 0.5*nLength*nOrientation;
                    Set::Vector nRight = nCenter + 0.5*nLength*nOrientation;
                    Set::Scalar distLeft = (x-nLeft).transpose()*nOrientation;
                    Set::Scalar distRight = (x-nRight).transpose()*nOrientation;

                    if(distLeft < 0.)
                        field(i,j,k) = 0.5 + 0.5*std::erf( ((x-nLeft).lpNorm<2>() - nThickness/2.0) / eps );
                    if(distRight > 0.)
                        field(i,j,k) = 0.5 + 0.5*std::erf( ((x-nRight).lpNorm<2>() - nThickness/2.0) / eps );
                }
            
                if (field(i,j,k) < 0.) field(i,j,k) = 0.;
                if (field(i,j,k) > 1.) field(i,j,k) = 1.;
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

public:
    static void Parse(Notch & value, IO::ParmParse & pp)
    {
        amrex::Vector<Set::Scalar> center;
		pp.queryarr("center",center);
        if(center.size()>=AMREX_SPACEDIM)
        {
            AMREX_D_TERM(value.nCenter(0) = center[0];,value.nCenter(1) = center[1];,value.nCenter(2) = center[2];);
        }
        
        amrex::Vector<Set::Scalar> orientation;
		pp.queryarr("orientation",orientation);

        if(orientation.size()>=AMREX_SPACEDIM)
        {
            AMREX_D_TERM(value.nOrientation(0) = orientation[0];,value.nOrientation(1) = orientation[1];,value.nOrientation(2) = orientation[2];);
        }

        if(value.nOrientation.lpNorm<2>() <= 0.) value.nOrientation = Set::Vector::Random();

        pp.query("thickness", value.nThickness);
        pp.query("length", value.nLength);
        pp.query("eps", value.eps);

        if(value.nThickness <= 0.) value.nThickness = 0.01;
        if(value.nLength <= 0.) value.nLength = 0.1;
        if(value.eps <= 0.) value.eps = 1.e-5;

        std::string mollifier;
        pp.query("mollifier",mollifier);
        if(mollifier == "Dirac" || mollifier == "dirac") 
            value.moll = Mollifier::Dirac;
        else 
            value.moll = Mollifier::Gaussian;

        value.nOrientation = value.nOrientation / value.nOrientation.lpNorm<2>();
        value.eps_sq = value.eps*value.eps;
        value.window = 6.0*value.eps;

        if(value.nOrientation(0) != 0.)
        {
            AMREX_D_TERM(value.nNormal(0) = 1.;, value.nNormal(1) = 1.;, value.nNormal(2) = 1.;);
            value.nNormal(0) = -(AMREX_D_TERM(0.,+value.nOrientation(1),+value.nOrientation(2)))/value.nOrientation(0);
            value.nNormal = value.nNormal/value.nNormal.lpNorm<2>();
        }
        else if(value.nOrientation(1) != 0.)
        {
            AMREX_D_TERM(value.nNormal(0) = 1.;, value.nNormal(1) = 1.;, value.nNormal(2) = 1.;);
            value.nNormal(1) = -(AMREX_D_TERM(value.nOrientation(0), + 0.0, + value.nOrientation(2)))/value.nOrientation(1);
            value.nNormal = value.nNormal/value.nNormal.lpNorm<2>();
        }
        //Util::Message(INFO,"nOrientation = (", nNormal(0),",",nNormal(1),")");
    }
};
}
#endif