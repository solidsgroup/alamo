#ifndef MODEL_SOLID_ELASTIC_PSEUDOLINEARCUBIC_H_
#define MODEL_SOLID_ELASTIC_PSEUDOLINEARCUBIC_H_

#include "IO/ParmParse.H"
#include "Model/Solid/Solid.H"

namespace Model
{
namespace Solid
{
namespace Elastic
{
class PseudoLinearCubic : public Solid<Set::Sym::Major>
{
public:
    PseudoLinearCubic() {};
    PseudoLinearCubic(Solid<Set::Sym::Major> base) : Solid<Set::Sym::Major>(base) {};
    virtual ~PseudoLinearCubic() {};

    Set::Scalar W(const Set::Matrix & F) const override
    {
        Set::Matrix E = 0.5 * (F.transpose() * F - Set::Matrix::Identity());
        return 0.5 * (E * (CC*E).transpose()).trace();
    }
    Set::Matrix DW(const Set::Matrix & F) const override
    {
        Set::Matrix C = F.transpose() * F;
        Set::Matrix dw = Set::Matrix::Zero();
        for (int i = 0; i < AMREX_SPACEDIM; i++)
        for (int J = 0; J < AMREX_SPACEDIM; J++)
        for (int A = 0; A < AMREX_SPACEDIM; A++)
        for (int R = 0; R < AMREX_SPACEDIM; R++)
        for (int S = 0; S < AMREX_SPACEDIM; S++)
        {
            dw(i,J) += 0.5 * F(i,A)*CC(J,A,R,S) * (C(R,S) - ((R==S) ? 1.0 : 0.0));
        }
        return dw;
    }
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Major> DDW(const Set::Matrix & F) const override
    {
        Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Major> ddw;
        Set::Matrix C = F.transpose() * F;
        Set::Matrix I = Set::Matrix::Identity();
        
        for (int i = 0; i < AMREX_SPACEDIM; i++)
        for (int J = 0; J < AMREX_SPACEDIM; J++)
        for (int k = 0; k < AMREX_SPACEDIM; k++)
        for (int L = 0; L < AMREX_SPACEDIM; L++)
        {
            ddw(i,J,k,L) = 0.0;
            for (int A = 0; A < AMREX_SPACEDIM; A++)
                for (int B = 0; B < AMREX_SPACEDIM; B++)
                {
                    ddw(i,J,k,L) += 0.5*I(i,k)*CC(J,L,A,B) * (C(A,B) - I(A,B));
                    ddw(i,J,k,L) += F(i,A)*F(k,B)*CC(A,J,B,L);
                }
        }

        return ddw;
    }
    virtual void Print(std::ostream &out) const override 
    {
        out << "CC = " << CC;
    }
    
public:
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor> CC;
    static constexpr KinematicVariable kinvar = KinematicVariable::F;

public:
    static PseudoLinearCubic Zero()
    {
        PseudoLinearCubic ret;
        ret.ddw = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor>::Zero();
        return ret;
    }
    static PseudoLinearCubic Random()
    {
        return Random(Util::Random(), Util::Random(), Util::Random());
    }
    static PseudoLinearCubic Random(Set::Scalar C11, Set::Scalar C12, Set::Scalar C44)
    {
        PseudoLinearCubic ret;
        Set::Scalar phi1 = 2.0*Set::Constant::Pi * Util::Random();
        Set::Scalar Phi  = 2.0*Set::Constant::Pi * Util::Random();
        Set::Scalar phi2 = 2.0*Set::Constant::Pi * Util::Random();
        ret.CC = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor>::Cubic(C11,C12,C44,phi1,Phi,phi2);
        return ret;
    }
    static void Parse(PseudoLinearCubic & value, IO::ParmParse & pp)
    {
        Set::Scalar C11 = 1.68, C12 = 1.21, C44 = 0.75;
        pp.query("C11",C11); // Elastic constant (default: 1.68)
        pp.query("C12",C12); // Elastic constant (default: 1.21)
        pp.query("C44",C44); // Elastic constant (default: 0.75)

        if (pp.contains("random"))
        {
            value = PseudoLinearCubic::Random(C11,C12,C44);
            return;
        }

        Set::Scalar phi1 = 0.0, Phi = 0.0, phi2 = 0.0;
        pp.query("phi1",phi1);  // Bunge Euler angle :math:`\phi_1`
        pp.query("Phi",Phi);    // Bunge Euler angle :math:`\Phi`
        pp.query("phi2",phi2);  // Bunge Euler angle :math:`\phi_2`

        value.CC = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor>::Cubic(C11,C12,C44,phi1,Phi,phi2);
    }

#define OP_CLASS PseudoLinearCubic
#define OP_VARS X(CC)
#include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"

}
}
}

#endif
