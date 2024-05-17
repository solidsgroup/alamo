//
// This class implements basic hexagonal elasticity.
//

#ifndef MODEL_SOLID_LINEAR_HEXAGONAL_H_
#define MODEL_SOLID_LINEAR_HEXAGONAL_H_

#include "Model/Solid/Solid.H"
#include "IO/ParmParse.H"

namespace Model
{
namespace Solid
{
namespace Linear
{
class Hexagonal : public Solid<Set::Sym::MajorMinor>
{
public:

    Hexagonal() {};
    Hexagonal(Solid<Set::Sym::MajorMinor> base) : Solid<Set::Sym::MajorMinor>(base) {};
    virtual ~Hexagonal() {};

    void
    Define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C13, Set::Scalar C33, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2)
    {
        Eigen::Matrix3d m;
        m = Eigen::AngleAxisd(phi2, Eigen::Vector3d::UnitX()) *
            Eigen::AngleAxisd(Phi,  Eigen::Vector3d::UnitZ()) *
            Eigen::AngleAxisd(phi1, Eigen::Vector3d::UnitX());
        Define(C11,C12,C13,C33,C44,m);
    }
    void
    Define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C13, Set::Scalar C33, Set::Scalar C44, Eigen::Matrix3d R)
    {
    
        amrex::Real Ctmp[3][3][3][3];
        ddw = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor>::Zero();

        for(int i = 0; i < 3; i++) 
            for(int j = 0; j < 3; j++) 
                for(int k = 0; k < 3; k++) 
                    for(int l = 0; l < 3; l++)
                    {
                        if(i == j && j == k && k == l)
                        {
                            if      (i==0) Ctmp[i][j][k][l] = C11;
                            else if (i==1) Ctmp[i][j][k][l] = C11;
                            else if (i==2) Ctmp[i][j][k][l] = C33;
                        }  
                        else if (i==k && j==l) 
                        {
                            if ((i==0 && j==1) || (i==1 && j==0))
                                Ctmp[i][j][k][l] = 0.5*(C11-C12);
                            else
                                Ctmp[i][j][k][l] = C44;
                        }
                        else if (i==j && k==l) 
                        {
                            if      ((i==0 && k==1) || (i==1 && k==0)) Ctmp[i][j][k][l] = C12;
                            else if ((i==0 && k==2) || (i==2 && k==0)) Ctmp[i][j][k][l] = C13;
                            else if ((i==1 && k==2) || (i==2 && k==1)) Ctmp[i][j][k][l] = C13;
                        }
                        else Ctmp[i][j][k][l] = 0.0;
                    }
        for(int p = 0; p < AMREX_SPACEDIM; p++) 
            for(int q = 0; q < AMREX_SPACEDIM; q++) 
                for(int s = 0; s < AMREX_SPACEDIM; s++) 
                    for(int t = 0; t < AMREX_SPACEDIM; t++)
                    {
                        ddw(p,q,s,t) = 0.0;
                        for(int i = 0; i < 3; i++) 
                            for(int j = 0; j < 3; j++) 
                                for(int k = 0; k < 3; k++) 
                                    for(int l = 0; l < 3; l++) 
                                        ddw(p,q,s,t) += R(p,i)*R(s,k)*Ctmp[i][j][k][l]*R(q,j)*R(t,l);
                    }
    }
    Set::Scalar W(const Set::Matrix & gradu) const override
    {
        return ( 0.5 * gradu.transpose() * (ddw*gradu) ).trace();
    }
    Set::Matrix DW(const Set::Matrix & gradu) const override
    {
        return ddw*gradu;
    }
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor> DDW(const Set::Matrix & /*gradu*/) const override
    {
        return ddw;
    }
    virtual void Print(std::ostream &out) const override 
    {
        out << ddw;
    }

public:
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor> ddw;
    static const KinematicVariable kinvar = KinematicVariable::gradu;

    AMREX_FORCE_INLINE
    static Hexagonal Combine(const std::vector<Hexagonal> &models, const std::vector<Set::Scalar> &eta)
    {
        Hexagonal ret;
        ret.ddw = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor>::Zero();
        Set::Scalar etasum = 0.;
        for (unsigned int n = 0 ; n < models.size(); n++) etasum += eta[n];
        for (unsigned int n = 0 ; n < models.size(); n++)
        {
            ret.ddw += models[n].ddw * (eta[n] / etasum);
        }
        return ret;
    }

    static Hexagonal Zero()
    {
        Hexagonal ret;
        ret.Define(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
        return ret;
        
    }
    static Hexagonal Random()
    {
        return Random(Util::Random(), Util::Random(), Util::Random(),Util::Random(),Util::Random());
    }
    static Hexagonal Random(Set::Scalar C11, Set::Scalar C12, Set::Scalar C13, Set::Scalar C33, Set::Scalar C44)
    {
        Hexagonal ret;
        Set::Scalar phi1 = 2.0*Set::Constant::Pi * Util::Random();
        Set::Scalar Phi  = 2.0*Set::Constant::Pi * Util::Random();
        Set::Scalar phi2 = 2.0*Set::Constant::Pi * Util::Random();
        ret.Define(C11,C12,C13,C33,C44,phi1,Phi,phi2);
        return ret;
    }

    static void Parse(Hexagonal & value, IO::ParmParse & pp)
    {
        Set::Scalar C11 = NAN, C12 = NAN, C13 = NAN, C33 = NAN, C44 = NAN;
        pp.query_required("C11",C11); // Elastic constant
        pp.query_required("C12",C12); // Elastic constant
        pp.query_required("C13",C13); // Elastic constant
        pp.query_required("C33",C33); // Elastic constant
        pp.query_required("C44",C44); // Elastic constant

        if (pp.contains("random"))
        {
            value = Hexagonal::Random(C11,C12,C13,C33,C44);
            return;
        }

        Set::Scalar phi1 = NAN, Phi = NAN, phi2 = NAN;
        pp.query_default("phi1",phi1,0.0);  // Bunge Euler angle :math:`\phi_1`
        pp.query_default("Phi", Phi, 0.0);  // Bunge Euler angle :math:`\Phi`
        pp.query_default("phi2",phi2,0.0);  // Bunge Euler angle :math:`\phi_2`
        value.Define(C11,C12,C13,C33,C44,phi1,Phi,phi2);
    }

    #define OP_CLASS Hexagonal
    #define OP_VARS X(ddw)
    #include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"

}
}
}

#endif
