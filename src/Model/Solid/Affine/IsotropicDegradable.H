#ifndef MODEL_SOLID_AFFINE_ISOTROPICDEGRADABLE_H_
#define MODEL_SOLID_AFFINE_ISOTROPICDEGRADABLE_H_

#include "AMReX.H"
#include "Model/Solid/Solid.H"
#include "Model/Solid/Affine/Isotropic.H"

namespace Model
{
namespace Solid
{
namespace Affine
{
class IsotropicDegradable : public Isotropic
{
public:
    IsotropicDegradable() {};
    IsotropicDegradable(Solid<Set::Sym::Isotropic> base) : Isotropic(base) {};
    IsotropicDegradable(Isotropic base) : Isotropic(base) {}
    IsotropicDegradable(Set::Scalar a_mu, Set::Scalar a_lambda)
    {
        Define(a_mu, a_lambda);
    }

    void Define(Set::Scalar a_mu, Set::Scalar a_lambda)
    {
        Define(a_mu, a_mu, a_lambda, a_lambda);
    }
    
    void Define(Set::Scalar a_mu, Set::Scalar a_mu0, Set::Scalar a_lambda, Set::Scalar a_lambda0, Set::Matrix a_F0=Set::Matrix::Zero())
    {
        m_mu = a_mu0;
        mu0 = a_mu0; lambda0 = a_lambda0; 
        Isotropic::Define(a_mu, a_lambda, a_F0);
    }

    void DegradeModulus(const Set::Scalar eta)
    {
        if(eta < 0.0 || eta > 1.0) Util::Abort(INFO, "Invalid value of eta", eta);
        if(std::isnan(lambda0)) Util::Abort(INFO, "Nans in lambda0");
        if(std::isnan(mu0)) Util::Abort(INFO, "Nans in mu0");
        Isotropic::Define( mu0*(1.-eta), lambda0*(1.-eta), _F0 );
        m_mu = mu0*(1-eta);
    }

    void DegradeModulus(const amrex::Vector<Set::Scalar> &eta)
    {
        if(eta.size() == 0) return;
        else if(eta.size() == 1)
        {
            if(std::isnan(eta[0])) Util::Abort(INFO,"Nans in eta");
            if(std::isinf(eta[0])) Util::Abort(INFO,"Infs in eta");
            if(eta[0] < 0.0 || eta[0] > 1.0) Util::Abort(INFO,"Invalid value of eta",eta[0]);
            Isotropic::Define( mu0*(1.-eta[0]), lambda0*(1.-eta[0]), _F0 );
            m_mu = mu0*(1-eta[0]);
        }
        else if(eta.size() == 2)
        {
            if(std::isnan(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Nans in eta");
            if(std::isinf(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Infs in eta");
            if(eta[0] < 0.0 || eta[0] > 1.0 || eta[1] < 0.0 || eta[1] > 1.0) Util::Abort(INFO,"Invalid value of eta");
            Isotropic::Define( mu0*(1.-eta[0]), lambda0*(1.-eta[1]), _F0 );
            m_mu = mu0*(1-eta[0]);
        }
        else
            return;
    }
    AMREX_FORCE_INLINE
    void SetF0(Set::Matrix a_F0) {_F0 = a_F0; Isotropic::SetF0(a_F0); }
    Set::Scalar GetMu () {return m_mu;}

public:
    Set::Scalar mu0 = NAN, lambda0 = NAN;
    Set::Scalar m_mu;
    Set::Matrix _F0;

public:
    static IsotropicDegradable Random()
    {
        IsotropicDegradable ret;
        ret.Define(Util::Random(), Util::Random(), Util::Random(), Util::Random(), Set::Matrix::Random());
        return ret;
    }

    static void Parse(IsotropicDegradable & value, IO::ParmParse & pp)
    {
        Set::Scalar mu, lambda;
        if (pp.contains("lame") && pp.contains("shear"))
        {
            pp.query("lame",lambda);
            pp.query("shear",mu);
        }
        else if (pp.contains("E") && pp.contains("nu"))
        {
            Set::Scalar E, nu;
            pp.query("E",E);
            pp.query("nu",nu);
            lambda = E * nu / (1.0 + nu) / (1.0 - 2.0*nu);
            mu = E / 2.0 / (1.0 + nu);
        }
        value.Define(mu,mu,lambda,lambda,Set::Matrix::Zero());
    }
    AMREX_FORCE_INLINE
    void operator += (const IsotropicDegradable &rhs)
    {
        ddw += rhs.ddw;
        F0 += rhs.F0;
        _F0 += rhs._F0;
        mu0 += rhs.mu0;
        m_mu += rhs.m_mu;
        lambda0 += rhs.lambda0;
    }    
    AMREX_FORCE_INLINE
        IsotropicDegradable operator * (const Set::Scalar alpha) const
    {
        IsotropicDegradable ret;
        ret.ddw = ddw*alpha;
        ret.F0 = F0*alpha;
        ret._F0 = _F0*alpha;
        ret.mu0 = mu0*alpha;
        ret.m_mu = m_mu * alpha;
        ret.lambda0 = lambda0*alpha;
        return ret;
    }    
    friend IsotropicDegradable operator * (const Set::Scalar alpha, const IsotropicDegradable b);
    friend IsotropicDegradable operator + (const IsotropicDegradable a, const IsotropicDegradable b);
    friend IsotropicDegradable operator - (const IsotropicDegradable a, const IsotropicDegradable b);
};
AMREX_FORCE_INLINE
IsotropicDegradable operator * (const Set::Scalar alpha, const IsotropicDegradable b)
{
    IsotropicDegradable ret;
    ret.ddw = b.ddw*alpha;
    ret.F0 = b.F0*alpha;
    ret._F0 = b._F0*alpha;
    ret.mu0 = b.mu0*alpha;
    ret.m_mu = b.m_mu * alpha;
    ret.lambda0 = b.lambda0*alpha;
    return ret;
}    
AMREX_FORCE_INLINE
IsotropicDegradable operator + (const IsotropicDegradable a, const IsotropicDegradable b)
{
    IsotropicDegradable ret;
    ret.ddw = a.ddw + b.ddw;
    ret.F0 = a.F0 + b.F0;
    ret._F0 = a._F0 + b._F0;
    ret.mu0 = a.mu0 + b.mu0;
    ret.m_mu = a.m_mu + b.m_mu;
    ret.lambda0 = a.lambda0 + b.lambda0;
    return ret;
}    
AMREX_FORCE_INLINE
IsotropicDegradable operator - (const IsotropicDegradable a, const IsotropicDegradable b)
{
    IsotropicDegradable ret;
    ret.ddw = a.ddw - b.ddw;
    ret.F0 = a.F0 - b.F0;
    ret._F0 = a._F0 - b._F0;
    ret.mu0 = a.mu0 - b.mu0;
    ret.m_mu = a.m_mu - b.m_mu;
    ret.lambda0 = a.lambda0 - b.lambda0;
    return ret;
} 

}
}
}
#endif
