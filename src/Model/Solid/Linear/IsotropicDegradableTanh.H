#ifndef MODEL_SOLID_LINEAR_ISOTROPICDEGRADABLETANH_H_
#define MODEL_SOLID_LINEAR_ISOTROPICDEGRADABLETANH_H_

#include "Model/Solid/Solid.H"
#include "IO/ParmParse.H"
#include "Model/Solid/Linear/Isotropic.H"

namespace Model
{
namespace Solid
{
namespace Linear
{
class IsotropicDegradableTanh : public Isotropic
{
public:
    IsotropicDegradableTanh() {};
    IsotropicDegradableTanh(Solid<Set::Sym::Isotropic> base) : Isotropic(base) {};
    IsotropicDegradableTanh(Isotropic base) : Isotropic(base) {};
    IsotropicDegradableTanh(Set::Scalar _lambda, Set::Scalar _mu)
    {
        Set::Scalar E1 = _mu*(3.*_lambda + 2.*_mu)/(_lambda + _mu);
        Set::Scalar nu = _lambda/(2.*(_lambda+_mu));
        Define(E1,0.2*E1,313.,20., nu, 298.);
    };
    IsotropicDegradableTanh(Set::Scalar _E1, Set::Scalar _E2, Set::Scalar _Tg, Set::Scalar _Ts, Set::Scalar _nu, Set::Scalar _temp=298.0)
    {
        Define(_E1,_E2,_Tg,_Ts,_nu,_temp);
    };
    void Define (Set::Scalar _E1, Set::Scalar _E2, Set::Scalar _Tg, Set::Scalar _Ts, Set::Scalar _nu, Set::Scalar _temp)
    {
        m_E10 = _E1; 
        m_E20 = _E2; 
        m_Tg0 = _Tg; 
        m_Ts0 = _Ts; 
        m_nu = _nu; 
        m_temp = _temp;
        
        Set::Scalar E = 0.5*(_E1 + _E2) - 0.5*(_E1 - _E2)*tanh((_temp - _Tg)/_Ts);
        Set::Scalar lambda = E*_nu/((1.+_nu)*(1.-2.*_nu));
        Set::Scalar mu = E/(2.*(1.+_nu));
        Isotropic::Define(lambda, mu);
    }
    void DegradeModulus(const Set::Scalar eta)
    {
        if(eta < 0.0 || eta > 1.0) Util::Abort(INFO, "Invalid value of eta", eta);
        Set::Scalar E1 = m_E10*(1.-eta);

        Set::Scalar E = 0.5*(E1 + m_E20) - 0.5*(E1-m_E20)*tanh((m_temp - m_Tg0)/m_Ts0);
        Set::Scalar lambda = E*m_nu/((1.+m_nu)*(1.-2.*m_nu));
        Set::Scalar mu = E/(2.*(1.+m_nu));
        Isotropic::Define(lambda,mu);
    }
    void DegradeModulus(const Set::Scalar eta1, const Set::Scalar eta2, const Set::Scalar eta3)
    {
        if(eta1 < 0.0 || eta1 > 1.0) Util::Abort(INFO,"Invalid value of eta1");
        if(eta2 < 0.0 || eta2 > 1.0) Util::Abort(INFO,"Invalid value of eta2");
        if(eta3 < 0.0 || eta3 > 1.0) Util::Abort(INFO,"Invalid value of eta3");
        
        Set::Scalar E1 = m_E10*(1.-eta1), Tg = m_Tg0*(1.-eta2), Ts = m_Ts0*(1.+eta3);
        Set::Scalar E = 0.5*(E1 + m_E20) - 0.5*(E1-m_E20)*tanh((m_temp - Tg)/Ts);
        Set::Scalar lambda = E*m_nu/((1.+m_nu)*(1.-2.*m_nu));
        Set::Scalar mu = E/(2.*(1.+m_nu));
        Isotropic::Define(lambda,mu);
    }
    void DegradeModulus(const amrex::Vector<Set::Scalar> &eta)
    {
        if (eta.size() == 0) return;
        else if (eta.size() == 1) DegradeModulus(eta[0]);
        else if (eta.size() == 2) DegradeModulus(eta[0],eta[1],eta[1]);
        else if (eta.size() == 3) DegradeModulus(eta[0],eta[1],eta[2]);
        else return;
    }
public:
    Set::Scalar m_E10 = NAN, m_E20 = NAN, m_Tg0 = NAN, m_Ts0 = NAN, m_nu = NAN, m_temp = NAN;

public:
    static void Parse(IsotropicDegradableTanh & value, IO::ParmParse & pp)
    {
        Set::Scalar E1, E2, Tg, Ts, nu, tmp=298.;
        pp.query("E1", E1);
        pp.query("E2", E2);
        pp.query("Tg", Tg);
        pp.query("Ts", Ts);
        pp.query("nu", nu);
        pp.query("temp", tmp);
        value.Define(E1,E2,Tg,Ts,nu,tmp);
    }
    AMREX_FORCE_INLINE
    void operator += (const IsotropicDegradableTanh &rhs)
    {
        ddw += rhs.ddw;
        m_E10 += rhs.m_E10;
        m_E20 += rhs.m_E20;
        m_Tg0 += rhs.m_Tg0;
        m_Ts0 += rhs.m_Ts0;
        m_nu += rhs.m_nu;
        m_temp += rhs.m_temp;
    }    
    AMREX_FORCE_INLINE
        IsotropicDegradableTanh operator * (const Set::Scalar alpha) const
    {
        IsotropicDegradableTanh ret;
        ret.ddw = ddw*alpha;
        ret.m_E10 = m_E10*alpha;
        ret.m_E20 = m_E20*alpha;
        ret.m_Tg0 = m_Tg0*alpha;
        ret.m_Ts0 = m_Ts0*alpha;
        ret.m_nu = m_nu*alpha;
        ret.m_temp = m_temp*alpha;
        return ret;
    }    
    friend IsotropicDegradableTanh operator * (const Set::Scalar alpha, const IsotropicDegradableTanh b);
    friend IsotropicDegradableTanh operator + (const IsotropicDegradableTanh a, const IsotropicDegradableTanh b);
    friend IsotropicDegradableTanh operator - (const IsotropicDegradableTanh a, const IsotropicDegradableTanh b);
};
AMREX_FORCE_INLINE
IsotropicDegradableTanh operator * (const Set::Scalar alpha, const IsotropicDegradableTanh b)
{
    IsotropicDegradableTanh ret;
    ret.ddw = b.ddw*alpha;
    ret.m_E10 = b.m_E10*alpha;
    ret.m_E20 = b.m_E20*alpha;
    ret.m_Tg0 = b.m_Tg0*alpha;
    ret.m_Ts0 = b.m_Ts0*alpha;
    ret.m_nu = b.m_nu*alpha;
    ret.m_temp = b.m_temp*alpha;
    return ret;
}    
AMREX_FORCE_INLINE
IsotropicDegradableTanh operator + (const IsotropicDegradableTanh a, const IsotropicDegradableTanh b)
{
    IsotropicDegradableTanh ret;
    ret.ddw = a.ddw + b.ddw;
    ret.m_E10 = a.m_E10 + b.m_E10;
    ret.m_E20 = a.m_E20 + b.m_E20;
    ret.m_Tg0 = a.m_Tg0 + b.m_Tg0;
    ret.m_Ts0 = a.m_Ts0 + b.m_Ts0;
    ret.m_nu = a.m_nu + b.m_nu;
    ret.m_temp = a.m_temp + b.m_temp;
    return ret;
}    
AMREX_FORCE_INLINE
IsotropicDegradableTanh operator - (const IsotropicDegradableTanh a, const IsotropicDegradableTanh b)
{
    IsotropicDegradableTanh ret;
    ret.ddw = a.ddw - b.ddw;
    ret.m_E10 = a.m_E10 - b.m_E10;
    ret.m_E20 = a.m_E20 - b.m_E20;
    ret.m_Tg0 = a.m_Tg0 - b.m_Tg0;
    ret.m_Ts0 = a.m_Ts0 - b.m_Ts0;
    ret.m_nu = a.m_nu - b.m_nu;
    ret.m_temp = a.m_temp - b.m_temp;
    return ret;
} 

}
}
}
#endif
