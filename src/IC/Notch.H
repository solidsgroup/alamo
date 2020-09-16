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
    	nCenter.resize(1); nCenter[0] = Set::Vector::Zero();
        nOrientation.resize(1); nOrientation[0] = Set::Vector::Random();
        nThickness.resize(1); nThickness[0] = 0.01;
        nLength.resize(1); nLength[0] = 0.1;
        moll = Mollifier::Dirac;
        eps = 1e-5;
    }
	
	void Add(const int &lev, Set::Field<Set::Scalar> &a_field)
	{
        Set::Vector DX(geom[lev].CellSize());

		for (amrex::MFIter mfi(*a_field[lev],amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
		{
			amrex::Box bx = mfi.tilebox();
		    bx.grow(a_field[lev]->nGrow());
			amrex::Array4<Set::Scalar> const& field = a_field[lev]->array(mfi);
			amrex::IndexType type = a_field[lev]->ixType();
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
                
                for (int i = 0; i < nCenter.size(); i++)
				{
					Set::Scalar value = 0;
	                if(std::abs((x-nCenter[i]).transpose()*nOrientation[i]) <= 0.5*nLength[i])
 	               {
 	                   Set::Scalar t = std::abs((x-nCenter[i]).transpose()*nNormal[i]) - (nThickness[i]/2.0);
 	                   value = 0.5 + 0.5*std::erf(t / eps);
 	               }
 	               else
 	               {
 	                   Set::Vector nLeft =  nCenter[i] - 0.5*nLength[i]*nOrientation[i];
 	                   Set::Vector nRight = nCenter[i] + 0.5*nLength[i]*nOrientation[i];
 	                   Set::Scalar distLeft = (x-nLeft).transpose()*nOrientation[i];
 	                   Set::Scalar distRight = (x-nRight).transpose()*nOrientation[i];

 	                   if(distLeft < 0.)
 	                       value = 0.5 + 0.5*std::erf( ((x-nLeft).lpNorm<2>() - nThickness[i]/2.0) / eps );
 	                   if(distRight > 0.)
 	                       value = 0.5 + 0.5*std::erf( ((x-nRight).lpNorm<2>() - nThickness[i]/2.0) / eps );
 	               }
 					min_value = value < min_value ? value : min_value;
            	}
            	field(i,j,k) = min_value;
            
                if (field(i,j,k) < 0.) field(i,j,k) = 0.;
                if (field(i,j,k) > 1.) field(i,j,k) = 1.;
			});
		}
        a_field[lev]->FillBoundary();
	}
	
private:
    amrex::Vector<Set::Vector> nCenter, nOrientation, nNormal;
    Set::Scalar eps = 1.e-5, eps_sq = eps*eps;
	amrex::Vector<Set::Scalar> nThickness, nLength;
    Set::Scalar window = 6.0*eps;
    Mollifier moll = Mollifier::Dirac;
    bool nRad = true;

public:
    static void Parse(Notch & value, IO::ParmParse & pp)
    {
    	value.nCenter.clear();
    	value.nOrientation.clear();
    	value.nNormal.clear();
    	value.nThickness.clear();
    	value.nLength.clear();
    
        amrex::Vector<Set::Scalar> center;
		pp.queryarr("center",center);
		
        if(center.size() >= AMREX_SPACEDIM)
        {
        	for (int i = 0; i<center.size(); i+=AMREX_SPACEDIM)
        		value.nCenter.push_back(Set::Vector(AMREX_D_DECL(center[i], center[i+1], center[i+2])));
            //AMREX_D_TERM(value.nCenter(0) = center[0];,value.nCenter(1) = center[1];,value.nCenter(2) = center[2];);
        }
        else
        	Util::Abort(INFO, "Insufficient values in center");
        
        amrex::Vector<Set::Scalar> orientation;
		pp.queryarr("orientation",orientation);

        if(orientation.size()>=AMREX_SPACEDIM && orientation.size() == center.size())
        {
        	for (int i=0; i < orientation.size(); i+=AMREX_SPACEDIM)
        		value.nOrientation.push_back(Set::Vector(AMREX_D_DECL(orientation[i], orientation[i+1], orientation[i+2])));
            //AMREX_D_TERM(value.nOrientation(0) = orientation[0];,value.nOrientation(1) = orientation[1];,value.nOrientation(2) = orientation[2];);
        }
        else
        	Util::Abort(INFO, "Insufficient values in orientation");
        
        for (int i =0; i < value.nOrientation.size(); i++)
        	if(value.nOrientation[i].lpNorm<2>() <= 0.) value.nOrientation[i] = Set::Vector::Random();

        pp.queryarr("thickness", value.nThickness);
        pp.queryarr("length", value.nLength);
        pp.query("eps", value.eps);

		if(value.nThickness.size() != value.nCenter.size()) Util::Abort(INFO, "Inconsistent size of thickness and centers");
		for (int i=0; i<value.nThickness.size(); i++)
			if(value.nThickness[i] <= 0.) value.nThickness[i] = 0.01;
        
        if(value.nLength.size() != value.nCenter.size()) Util::Abort(INFO, "Inconsistent size of length and centers");
        for (int i=0; i<value.nLength.size(); i++)
			if(value.nLength[i] <= 0.) value.nLength[i] = 0.1;

        if(value.eps <= 0.) value.eps = 1.e-5;

        std::string mollifier;
        pp.query("mollifier",mollifier);
        if(mollifier == "Dirac" || mollifier == "dirac") 
            value.moll = Mollifier::Dirac;
        else 
            value.moll = Mollifier::Gaussian;

		for (int i=0; i<value.nOrientation.size(); i++)
			value.nOrientation[i] = value.nOrientation[i] / value.nOrientation[i].lpNorm<2>();
		
		value.nNormal.resize(value.nOrientation.size());

        value.eps_sq = value.eps*value.eps;
        value.window = 6.0*value.eps;

		for (int i = 0; i<value.nOrientation.size(); i++)
		{
			value.nNormal[i] = Set::Vector::Zero();
	        if(value.nOrientation[i](0) != 0.)
 	       {
 	           AMREX_D_TERM(value.nNormal[i](0) = 1.;, value.nNormal[i](1) = 1.;, value.nNormal[i](2) = 1.;);
 	           value.nNormal[i](0) = -(AMREX_D_TERM(0.,+value.nOrientation[i](1),+value.nOrientation[i](2)))/value.nOrientation[i](0);
 	           value.nNormal[i] = value.nNormal[i]/value.nNormal[i].lpNorm<2>();
 	       }
 	       else if(value.nOrientation[i](1) != 0.)
 	       {
 	           AMREX_D_TERM(value.nNormal[i](0) = 1.;, value.nNormal[i](1) = 1.;, value.nNormal[i](2) = 1.;);
 	           value.nNormal[i](1) = -(AMREX_D_TERM(value.nOrientation[i](0), + 0.0, + value.nOrientation[i](2)))/value.nOrientation[i](1);
 	           value.nNormal[i] = value.nNormal[i]/value.nNormal[i].lpNorm<2>();
 	       }
 	       //Util::Message(INFO,"nOrientation = (", nNormal(0),",",nNormal(1),")");
 		}
    }
};
}
#endif