#ifndef MODEL_SOLID_ELASTIC_NEOHOOKEAN_H_
#define MODEL_SOLID_ELASTIC_NEOHOOKEAN_H_

#include "IO/ParmParse.H"
#include "Model/Solid/Solid.H"

namespace Model
{
namespace Solid
{
namespace Elastic
{
class NeoHookean : public Solid<Set::Sym::Major>
{
public:
    NeoHookean() {};
    NeoHookean(Solid<Set::Sym::Major> base) : Solid<Set::Sym::Major>(base) {};
    virtual ~NeoHookean() {};

    Set::Scalar W(const Set::Matrix & F) const override
    {
        Set::Scalar J = F.determinant();
        Set::Scalar J23 = std::pow(fabs(J),2./3.);
        Set::Scalar w = 0.0;
        w += 0.5 * mu * ((F*F.transpose()).trace() / J23 - 3.);
        w += 0.5 * kappa * (J - 1.0) * (J - 1.0);
        return w;
    }
    Set::Matrix DW(const Set::Matrix & F) const override
    {
        Set::Scalar J = F.determinant();
        Set::Scalar J23 = std::pow(fabs(J),2./3.);
        Set::Matrix FinvT = F.inverse().transpose();

        Set::Matrix dw = Set::Matrix::Zero();

        dw += mu * (F/J23 - (F*F.transpose()).trace()*FinvT / (3.*J23));
        dw += kappa*(J-1)*J*FinvT;

        return dw;
    }
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Major> DDW(const Set::Matrix & F) const override
    {
        Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Major> ddw;
        Set::Scalar J = F.determinant();
        Set::Scalar J23 = std::pow(fabs(J),2./3.);
        Set::Matrix FinvT = F.inverse().transpose();
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    for (int l = 0; l < 3; l++)
                    {
                        ddw(i,j,k,l) = 0.0;

                        Set::Scalar t1 = 0.0, t2 = 0.0;

                        if (i==k && j==l) t1 += 1.0;
                        t1 -= (2./3.) * F(i,j)*FinvT(k,l);
                        t1 -= (2./3.) * FinvT(i,j)*F(k,l);
                        t1 += (2./9.) * (F*F.transpose()).trace() * FinvT(i,j) * FinvT(k,l);
                        t1 += (1./3.) * (F*F.transpose()).trace() * FinvT(i,l) * FinvT(k,j);

                        t2 += (2.*J - 1.) * FinvT(i,j)*FinvT(k,l);
                        t2 += (1. - J) * FinvT(i,l) * FinvT(k,j);

                        ddw(i,j,k,l) = (mu/J23)*t1 + kappa*J*t2;
                    }
        return ddw;
    }
    virtual void Print(std::ostream &out) const override 
    {
        out << "mu = " << mu << " kappa = " << kappa;
    }
    
public:
    Set::Scalar mu = NAN, kappa = NAN;
    static constexpr KinematicVariable kinvar = KinematicVariable::F;

public:
    static NeoHookean Zero()
    {
        NeoHookean ret;
        ret.mu = 0.0;
        ret.kappa = 0.0;
        return ret;
    }
    static NeoHookean Random()
    {
        NeoHookean ret;
        ret.mu = Util::Random();
        ret.kappa = Util::Random();
        return ret;
    }
    static void Parse(NeoHookean & value, IO::ParmParse & pp)
    {
        pp.query("mu",value.mu);
        pp.query("kappa",value.kappa);
    }

    AMREX_FORCE_INLINE
    void operator += (const NeoHookean &rhs)
    {
        kappa += rhs.kappa;
        mu += rhs.mu;
    }    
    AMREX_FORCE_INLINE
        virtual NeoHookean operator * (const Set::Scalar alpha) const
    {
        NeoHookean ret;
        ret.kappa = kappa*alpha;
        ret.mu = kappa*mu;
        return ret;
    }    
    friend NeoHookean operator * (const Set::Scalar alpha, const NeoHookean b);
    friend NeoHookean operator + (const NeoHookean a, const NeoHookean b);
    friend NeoHookean operator - (const NeoHookean a, const NeoHookean b);
};
AMREX_FORCE_INLINE
NeoHookean operator * (const Set::Scalar alpha, const NeoHookean b)
{
    NeoHookean ret;
    ret.mu = alpha * b.mu;
    ret.kappa = alpha * b.kappa;
    return ret;
}    
AMREX_FORCE_INLINE
NeoHookean operator + (const NeoHookean a, const NeoHookean b)
{
    NeoHookean ret;
    ret.mu = a.mu + b.mu;
    ret.kappa = a.kappa + b.kappa;
    return ret;
}    
AMREX_FORCE_INLINE
NeoHookean operator - (const NeoHookean a, const NeoHookean b)
{
    NeoHookean ret;
    ret.mu = a.mu - b.mu;
    ret.kappa = a.kappa - b.kappa;
    return ret;
}    




}
}
}

#endif
