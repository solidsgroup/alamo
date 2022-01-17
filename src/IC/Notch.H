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
    enum Mollifier {Dirac, Gaussian, Cosine};

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
        Set::Scalar pi = std::atan(1.0)*4;

        for (amrex::MFIter mfi(*a_field[lev],true); mfi.isValid(); ++mfi)
        {
            //amrex::Box bx = mfi.grownnodaltilebox();
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
                
                Set::Scalar min_value = field(i,j,k);
                
                for (int m = 0; m < nCenter.size(); m++)
                {
                    Set::Scalar value = 1.0;
                    if(std::abs((x-nCenter[m]).transpose()*nOrientation[m]) <= 0.5*nLength[m])
                    {
                        Set::Scalar t = std::abs((x-nCenter[m]).transpose()*nNormal[m]) - (nThickness[m]/2.0);
                        if (moll == Mollifier::Gaussian)
                            value = 0.5 + 0.5*std::erf(t / eps);
                            // value = std::erf(t/eps);
                        else if (moll == Mollifier::Dirac)
                        {
                            if (t < 0.0) value = 0.0;
                            else value = 1.0;
                        }
                        else if (moll == Mollifier::Cosine)
                        {
                            if (t < 0.0) value = 0.0;
                            else if (t < eps) value = 0.5 - 0.5*std::cos(pi*t/eps);
                            else value = 1.0;
                        }
                    }
                    else
                    {
                        Set::Vector nLeft =  nCenter[m] - 0.5*nLength[m]*nOrientation[m];
                        Set::Vector nRight = nCenter[m] + 0.5*nLength[m]*nOrientation[m];
                        Set::Vector nLeft2 = nLeft, nRight2 = nRight;
                        Set::Vector correction = nOrientation[m]*(std::sqrt(nRadius[m]*nRadius[m] - nThickness[m]*nThickness[m]/4.0)); 
                        if (nRadius[m] > 0.5*nThickness[m])
                        {
                            nLeft2 = nLeft + correction;
                            nRight2 = nRight - correction;
                        }
                        
                        Set::Scalar cosR = (nRadius[m] <= 0.0) ? 1.0 : std::sqrt(nRadius[m]*nRadius[m] - nThickness[m]*nThickness[m]/4.0)/nRadius[m];
                        Set::Scalar cosLeft = (x-nLeft2).dot(nOrientation[m]) / ((x-nLeft2).lpNorm<2>());
                        Set::Scalar sinLeft = (x-nLeft2).dot(nNormal[m]) / ((x-nLeft2).lpNorm<2>());
                        Set::Scalar cosRight = (x-nRight2).dot(nOrientation[m]) / ((x-nRight2).lpNorm<2>());
                        Set::Scalar sinRight = (x-nRight2).dot(nNormal[m]) / ((x-nRight2).lpNorm<2>());
                        
                        Set::Scalar distLeft = (x-nLeft).dot(nOrientation[m]);
                        Set::Scalar distRight = (x-nRight).dot(nOrientation[m]);
                        

                        if(distLeft <= 0.)
                        {
                            Set::Scalar t = (x-nLeft2).lpNorm<2>() - nRadius[m];
                            if (moll == Mollifier::Gaussian) 
                            {
                                value = 0.5 + 0.5*std::erf( ((x-nLeft).lpNorm<2>() - nThickness[m]/2.0) / eps );
                            }
                            else if (moll == Mollifier::Dirac)
                            {
                                if (t < 0.0) value = 0.0;
                                else value = 1.0;
                            }
                            else if (moll == Mollifier::Cosine)
                            {
                                Set::Vector nTopLeft = nLeft + 0.5*nThickness[m]*nNormal[m];
                                Set::Vector nBotLeft = nLeft - 0.5*nThickness[m]*nNormal[m];

                                if (nRadius[m] <= 0.0)
                                {
                                    Set::Scalar t2 = -(x-nTopLeft).dot(nOrientation[m]);
                                    
                                    Set::Scalar angTop = (x-nTopLeft).dot(nNormal[m]);
                                    Set::Scalar angBot = (x-nBotLeft).dot(nNormal[m]);

                                    if(angTop <= 0.0 && angBot >= 0.0)
                                    {
                                        if (t2 < 0.0) value = 0.0;
                                        else if (t2 < eps) value = 0.5 - 0.5*std::cos(pi*t2/eps);
                                        else value = 1.0;
                                    }
                                    else if(angTop > 0.0  && angBot > 0.0)
                                    {
                                        Set::Scalar t3 = (x-nTopLeft).lpNorm<2>();
                                        if (t3 < 0.0) value = 0.0;
                                        else if (t3 < eps) value = 0.5 - 0.5*std::cos(pi*t3/eps);
                                        else value = 1.0;
                                    }
                                    else
                                    {
                                        Set::Scalar t3 = (x-nBotLeft).lpNorm<2>();
                                        if (t3 < 0.0) value = 0.0;
                                        else if (t3 < eps) value = 0.5 - 0.5*std::cos(pi*t3/eps);
                                        else value = 1.0;
                                    }
                                }

                                else 
                                {
                                    if (std::abs(cosLeft) > cosR)
                                    {
                                        if (t < 0.0) value = 0.0;
                                        else if (t < eps) value = 0.5 - 0.5*std::cos(pi*t/eps);
                                        else value = 1.0;
                                    }

                                    else
                                    {
                                        if (sinLeft >= 0)
                                        {
                                            Set::Scalar t2 = (x-nTopLeft).lpNorm<2>();
                                            if (t2 < 0.0) value = 0.0;
                                            else if (t2 < eps) value = 0.5 - 0.5*std::cos(pi*t2/eps);
                                            else value = 1.0;
                                        }
                                        else
                                        {
                                            Set::Scalar t2 = (x-nBotLeft).lpNorm<2>();
                                            if (t2 < 0.0) value = 0.0;
                                            else if (t2 < eps) value = 0.5 - 0.5*std::cos(pi*t2/eps);
                                            else value = 1.0;
                                        }
                                    }
                                }
                            }
                        }
                        if(distRight >= 0.)
                        {
                            Set::Scalar t = (x-nRight2).lpNorm<2>() - nRadius[m];
                            if (moll == Mollifier::Gaussian) 
                            {
                                value = 0.5 + 0.5*std::erf( ((x-nRight).lpNorm<2>() - nThickness[m]/2.0) / eps );
                            }
                            else if (moll == Mollifier::Dirac)
                            {
                                if (t < 0.0) value = 0.0;
                                else value = 1.0;
                            }
                            else if (moll == Mollifier::Cosine)
                            {
                                Set::Vector nTopRight = nRight + 0.5*nThickness[m]*nNormal[m];
                                Set::Vector nBotRight = nRight - 0.5*nThickness[m]*nNormal[m];
                                
                                if (nRadius[m] <= 0.0)
                                {
                                    Set::Scalar t2 = (x-nTopRight).dot(nOrientation[m]);
                                    Set::Scalar angTop = (x-nTopRight).dot(nNormal[m]);
                                    Set::Scalar angBot = (x-nBotRight).dot(nNormal[m]);

                                    if(angTop <= 0.0 && angBot >= 0.0)
                                    {
                                        if (t2 < 0.0) value = 0.0;
                                        else if (t2 < eps) value = 0.5 - 0.5*std::cos(pi*t2/eps);
                                        else value = 1.0;
                                    }
                                    else if(angTop > 0.0  && angBot > 0.0)
                                    {
                                        Set::Scalar t3 = (x-nTopRight).lpNorm<2>();
                                        if (t3 < 0.0) value = 0.0;
                                        else if (t3 < eps) value = 0.5 - 0.5*std::cos(pi*t3/eps);
                                        else value = 1.0;
                                    }
                                    else
                                    {
                                        Set::Scalar t3 = (x-nBotRight).lpNorm<2>();
                                        if (t3 < 0.0) value = 0.0;
                                        else if (t3 < eps) value = 0.5 - 0.5*std::cos(pi*t3/eps);
                                        else value = 1.0;
                                    }
                                }

                                else 
                                {
                                    if (std::abs(cosRight) > cosR)
                                    {
                                        if (t < 0.0) value = 0.0;
                                        else if (t < eps) value = 0.5 - 0.5*std::cos(pi*t/eps);
                                        else value = 1.0;
                                    }

                                    else
                                    {
                                        if (sinRight >= 0)
                                        {
                                            Set::Scalar t2 = (x-nTopRight).lpNorm<2>();
                                            if (t2 < 0.0) value = 0.0;
                                            else if (t2 < eps) value = 0.5 - 0.5*std::cos(pi*t2/eps);
                                            else value = 1.0;
                                        }
                                        else
                                        {
                                            Set::Scalar t2 = (x-nBotRight).lpNorm<2>();
                                            if (t2 < 0.0) value = 0.0;
                                            else if (t2 < eps) value = 0.5 - 0.5*std::cos(pi*t2/eps);
                                            else value = 1.0;
                                        }
                                    }
                                }
                            }
                        }
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
    amrex::Vector<Set::Scalar> nThickness, nLength, nRadius;
    Mollifier moll = Mollifier::Dirac;

public:
    static void Parse(Notch & value, IO::ParmParse & pp)
    {
        value.nCenter.clear();
        value.nOrientation.clear();
        value.nNormal.clear();
        value.nRadius.clear();
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
        pp.queryarr("radius", value.nRadius);
        pp.query("eps", value.eps);
        
        if(value.nThickness.size() == 0)
        {
            value.nThickness.resize(value.nCenter.size());
            for (int i=0; i<value.nThickness.size(); i++) value.nThickness[i] = 0.0;
        }
        else if(value.nThickness.size() != value.nCenter.size()) Util::Abort(INFO, "Inconsistent size of thickness and centers");
        else
        {
            for (int i=0; i<value.nThickness.size(); i++)
                if(value.nThickness[i] <= 0.) value.nThickness[i] = 0.01;
        }
        
        if(value.nLength.size() != value.nCenter.size()) Util::Abort(INFO, "Inconsistent size of length and centers");
        for (int i=0; i<value.nLength.size(); i++)
            if(value.nLength[i] <= 0.) value.nLength[i] = 0.1;

        if (value.nRadius.size() == 0)
        {
            value.nRadius.resize(value.nThickness.size());
            for (int i=0; i<value.nRadius.size(); i++)
                value.nRadius[i] = value.nThickness[i]/2.0;
        }
        else if(value.nRadius.size() != value.nCenter.size()) Util::Abort(INFO, "Inconsistent size of radius and centers");
        else 
        {
            for (int i=0; i<value.nRadius.size(); i++)
            {
                if(value.nRadius[i] <= 0.0) value.nRadius[i] = 0.0;
                else if(value.nRadius[i] <= value.nThickness[i]) value.nRadius[i] = value.nThickness[i];
            }
        }
        //Util::Message(INFO, "value.nRadius.size() = ", value.nRadius.size(), ". value.nRadius[0] = ", value.nRadius[0]);
        if(value.eps <= 0.) value.eps = 1.e-5;

        std::string mollifier;
        pp.query("mollifier",mollifier);
        if(mollifier == "Dirac" || mollifier == "dirac") 
            value.moll = Mollifier::Dirac;
        else if (mollifier == "gauss" || mollifier == "Gauss")
        {
            value.moll = Mollifier::Gaussian;
            for (int i=0; i<value.nRadius.size(); i++)
                value.nRadius[i] = 0.5*value.nThickness[i];
        }
        else if (mollifier == "erf" || mollifier == "Erf")
        {
            value.moll = Mollifier::Gaussian;
            for (int i=0; i<value.nRadius.size(); i++)
                value.nRadius[i] = 0.5*value.nThickness[i];
        }
        else if (mollifier == "cos" || mollifier == "Cosine")
            value.moll = Mollifier::Cosine;
        else
            value.moll = Mollifier::Dirac;

        for (int i=0; i<value.nOrientation.size(); i++)
            value.nOrientation[i] = value.nOrientation[i] / value.nOrientation[i].lpNorm<2>();
        
        value.nNormal.resize(value.nOrientation.size());

        value.eps_sq = value.eps*value.eps;
        
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
