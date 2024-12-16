//
// This class implements basic Transvesely Isotropic elasticity .
// For a discussion on Transeversely Isotropics, please see this link <http://solidmechanics.org/text/Chapter3_2/Chapter3_2.htm#Sect3_2_16>_.
//

#ifndef MODEL_SOLID_LINEAR_TRANSVERSELYISOTROPIC_H_
#define MODEL_SOLID_LINEAR_TRANSVERSELYISOTROPIC_H_

#include "Model/Solid/Solid.H"
//# include "Model/Solid/Linear/TransverselyIsotropic.H"
#include "IO/ParmParse.H"

namespace Model
{
namespace Solid
{
namespace Linear
{
class TransverselyIsotropic : public Solid<Set::Sym::MajorMinor>
{
public:

    TransverselyIsotropic() {};
    TransverselyIsotropic(Solid<Set::Sym::MajorMinor> base) : Solid<Set::Sym::MajorMinor>(base) {};
    virtual ~TransverselyIsotropic() {};
    
    void
    Define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C13, Set::Scalar C33, Set::Scalar C44, Set::Scalar phi1, Set::Scalar Phi, Set::Scalar phi2)
    {
        ddw = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor>::Transverse(C11,C12,C13,C33,C44,phi1,Phi,phi2);
    }

    void
    Define(Set::Scalar C11, Set::Scalar C12, Set::Scalar C13, Set::Scalar C33, Set::Scalar C44, 
                                            Eigen::Matrix3d R = Eigen::Matrix3d::Identity())
    {
        ddw = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor>::Transverse(C11,C12,C13,C33,C44,R);
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
    static TransverselyIsotropic Combine(const std::vector<TransverselyIsotropic> &models, const std::vector<Set::Scalar> &eta)
    {
        TransverselyIsotropic ret;
        ret.ddw = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor>::Zero();
        Set::Scalar etasum = 0.;
        for (unsigned int n = 0 ; n < models.size(); n++) etasum += eta[n];
        for (unsigned int n = 0 ; n < models.size(); n++)
        {
            ret.ddw += models[n].ddw * (eta[n] / etasum);
        }
        return ret;
    }

    static TransverselyIsotropic Zero()
    {
        TransverselyIsotropic ret;
        ret.Define(0.0,0.0,0.0,0.0,0.0);
        return ret;
        
    }
    static TransverselyIsotropic Random()
    {
        return Random(Util::Random(), Util::Random(), Util::Random(), Util::Random(), Util::Random());
    }
    static TransverselyIsotropic Random(Set::Scalar C11, Set::Scalar C12, Set::Scalar C13, Set::Scalar C33, Set::Scalar C44)
    {
        TransverselyIsotropic ret;
        Set::Scalar phi1 = 2.0*Set::Constant::Pi * Util::Random();
        Set::Scalar Phi  = 2.0*Set::Constant::Pi * Util::Random();
        Set::Scalar phi2 = 2.0*Set::Constant::Pi * Util::Random();
        ret.Define(C11,C12,C13,C33,C44,phi1,Phi,phi2);
        return ret;
    }

    static void Parse(TransverselyIsotropic & value, IO::ParmParse & pp)
    {
        Set::Scalar C11 = NAN, C12 = NAN, C13 = NAN, C33 = NAN , C44 = NAN; 
        
        pp_query_default("C11",C11,1.68); // Elastic constant
        pp_query_default("C12",C12,1.21); // Elastic constant
        pp_query_default("C13",C13,0.75); // Elastic constant
        pp_query_default("C33",C33,1.68); // Elastic constant
        pp_query_default("C44",C44,1.68); // Elastic constant

        if (pp.contains("random"))
        {
            value = TransverselyIsotropic::Random(C11,C12,C13,C33,C44);
            return;
        }

        Set::Scalar phi1 = NAN, Phi = NAN, phi2 = NAN;
        pp_query_default("phi1",phi1,0.0);  // Bunge Euler angle :math:\phi_1
        pp_query_default("Phi", Phi, 0.0);  // Bunge Euler angle :math:\Phi
        pp_query_default("phi2",phi2,0.0);  // Bunge Euler angle :math:\phi_2
        value.Define(C11,C12,C13,C33,C44,phi1,Phi,phi2);
    }

    #define OP_CLASS TransverselyIsotropic
    #define OP_VARS X(ddw)
    #include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"

}
}
}

#endif
