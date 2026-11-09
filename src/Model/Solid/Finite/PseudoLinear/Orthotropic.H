#ifndef MODEL_SOLID_FINITE_PSEUDOLINEAR_ORTHOTROPIC_H_
#define MODEL_SOLID_FINITE_PSEUDOLINEAR_ORTHOTROPIC_H_

#include "IO/ParmParse.H"
#include "Model/Solid/Solid.H"

namespace Model
{
namespace Solid
{
namespace Finite
{
namespace PseudoLinear
{
/// St. Venant--Kirchhoff finite elasticity with orthorhombic material symmetry.
/// The crystal axes are rotated into the laboratory frame by q. Specifying only
/// C11, C12, and C44 recovers the historical cubic behavior.
class Orthotropic : public Solid<Set::Sym::Major>
{
public:
    Orthotropic() = default;
    ~Orthotropic() = default;

    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor> Elasticity() const
    {
        Set::Quaternion normalized_q = q;
        normalized_q.normalize();
        const Set::Matrix3d R = normalized_q.toRotationMatrix();
        Set::Matrix4<AMREX_SPACEDIM,Set::Sym::MajorMinor> CC;

        auto normal = [=](int a, int c) -> Set::Scalar
        {
            if (a == 0 && c == 0) return C11;
            if (a == 1 && c == 1) return C22;
            if (a == 2 && c == 2) return C33;
            if ((a == 0 && c == 1) || (a == 1 && c == 0)) return C12;
            if ((a == 0 && c == 2) || (a == 2 && c == 0)) return C13;
            return C23;
        };
        auto shear = [=](int a, int b) -> Set::Scalar
        {
            if ((a == 1 && b == 2) || (a == 2 && b == 1)) return C44;
            if ((a == 0 && b == 2) || (a == 2 && b == 0)) return C55;
            return C66;
        };
        auto material = [&](int a, int b, int c, int d) -> Set::Scalar
        {
            if (a == b && c == d) return normal(a,c);
            if (a != b && c != d &&
                ((a == c && b == d) || (a == d && b == c)))
                return shear(a,b);
            return 0.0;
        };

        for (int i = 0; i < AMREX_SPACEDIM; ++i)
        for (int j = 0; j < AMREX_SPACEDIM; ++j)
        for (int k = 0; k < AMREX_SPACEDIM; ++k)
        for (int l = 0; l < AMREX_SPACEDIM; ++l)
        {
            CC(i,j,k,l) = 0.0;
            for (int a = 0; a < 3; ++a)
            for (int b = 0; b < 3; ++b)
            for (int c = 0; c < 3; ++c)
            for (int d = 0; d < 3; ++d)
                CC(i,j,k,l) += R(i,a)*R(j,b)*R(k,c)*R(l,d)*
                               material(a,b,c,d);
        }
        return CC;
    }

    Set::Scalar W(const Set::Matrix& F) const
    {
        const auto CC = Elasticity();
        const Set::Matrix E = 0.5*(F.transpose()*F - Set::Matrix::Identity());
        return 0.5*(E*(CC*E).transpose()).trace();
    }

    Set::Matrix DW(const Set::Matrix& F) const
    {
        const auto CC = Elasticity();
        const Set::Matrix C = F.transpose()*F;
        Set::Matrix dw = Set::Matrix::Zero();
        for (int i = 0; i < AMREX_SPACEDIM; ++i)
        for (int J = 0; J < AMREX_SPACEDIM; ++J)
        for (int A = 0; A < AMREX_SPACEDIM; ++A)
        for (int R = 0; R < AMREX_SPACEDIM; ++R)
        for (int S = 0; S < AMREX_SPACEDIM; ++S)
            dw(i,J) += 0.5*F(i,A)*CC(J,A,R,S)*
                       (C(R,S) - (R == S ? 1.0 : 0.0));
        return dw;
    }

    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Major> DDW(const Set::Matrix& F) const
    {
        const auto CC = Elasticity();
        Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Major> ddw;
        const Set::Matrix C = F.transpose()*F;
        const Set::Matrix I = Set::Matrix::Identity();
        for (int i = 0; i < AMREX_SPACEDIM; ++i)
        for (int J = 0; J < AMREX_SPACEDIM; ++J)
        for (int k = 0; k < AMREX_SPACEDIM; ++k)
        for (int L = 0; L < AMREX_SPACEDIM; ++L)
        {
            ddw(i,J,k,L) = 0.0;
            for (int A = 0; A < AMREX_SPACEDIM; ++A)
            for (int B = 0; B < AMREX_SPACEDIM; ++B)
            {
                ddw(i,J,k,L) += 0.5*I(i,k)*CC(J,L,A,B)*(C(A,B)-I(A,B));
                ddw(i,J,k,L) += F(i,A)*F(k,B)*CC(A,J,B,L);
            }
        }
        return ddw;
    }

    void Print(std::ostream& out) const
    {
        out << "C11=" << C11 << " C22=" << C22 << " C33=" << C33
            << " C12=" << C12 << " C13=" << C13 << " C23=" << C23
            << " C44=" << C44 << " C55=" << C55 << " C66=" << C66;
    }

public:
    static constexpr KinematicVariable kinvar = KinematicVariable::F;
    Set::Scalar C11=1.68, C22=1.68, C33=1.68;
    Set::Scalar C12=1.21, C13=1.21, C23=1.21;
    Set::Scalar C44=0.75, C55=0.75, C66=0.75;
    mutable Set::Quaternion q = Set::Quaternion(1.0,0.0,0.0,0.0);

    static Orthotropic Zero()
    {
        Orthotropic ret;
        ret.C11=ret.C22=ret.C33=0.0;
        ret.C12=ret.C13=ret.C23=0.0;
        ret.C44=ret.C55=ret.C66=0.0;
        ret.q = Set::Quaternion(0.0,0.0,0.0,0.0);
        return ret;
    }

    static void Parse(Orthotropic& value, IO::ParmParse& pp)
    {
        pp_query_default("C11",value.C11,1.68);
        pp_query_default("C12",value.C12,1.21);
        pp_query_default("C44",value.C44,0.75);
        pp_query_default("C22",value.C22,value.C11);
        pp_query_default("C33",value.C33,value.C11);
        pp_query_default("C13",value.C13,value.C12);
        pp_query_default("C23",value.C23,value.C12);
        pp_query_default("C55",value.C55,value.C44);
        pp_query_default("C66",value.C66,value.C44);

        if (pp.contains("random"))
        {
            value.q = Set::Quaternion::UnitRandom();
            return;
        }

        Set::Scalar phi1=0.0, Phi=0.0, phi2=0.0;
        pp_query_default("phi1",phi1,0.0);
        pp_query_default("Phi",Phi,0.0);
        pp_query_default("phi2",phi2,0.0);
        Eigen::Matrix3d R;
        R = Eigen::AngleAxisd(phi2,Eigen::Vector3d::UnitX()) *
            Eigen::AngleAxisd(Phi, Eigen::Vector3d::UnitZ()) *
            Eigen::AngleAxisd(phi1,Eigen::Vector3d::UnitX());
        value.q = Set::Quaternion(R);
    }

#define OP_CLASS Orthotropic
#define OP_VARS X(C11) X(C22) X(C33) X(C12) X(C13) X(C23) X(C44) X(C55) X(C66) X(q)
#include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"

}
}
}
}

#endif
