#ifndef MODEL_SOLID_LINEAR_ISOTROPIC_H_
#define MODEL_SOLID_LINEAR_ISOTROPIC_H_

#include "Model/Solid/Solid.H"
#include "IO/ParmParse.H"

namespace Model
{
namespace Solid
{
namespace Linear
{
class Isotropic : public Solid<Set::Sym::Isotropic>
{
public:

	Isotropic() {};
    Isotropic(Solid<Set::Sym::Isotropic> base) : Solid<Set::Sym::Isotropic>(base) {};
    Isotropic(Set::Scalar a_mu, Set::Scalar a_lambda) 
    {
        Define(a_mu,a_lambda);
    };
	virtual ~Isotropic() {};

	void Define(Set::Scalar a_mu, Set::Scalar a_lambda)
	{
		ddw = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic>(a_lambda,a_mu);
	}

    Set::Scalar W(const Set::Matrix & gradu) const override
    {
        return ( 0.5 * gradu.transpose() * (ddw*gradu) ).trace();
    }
    Set::Matrix DW(const Set::Matrix & gradu) const override
    {
        return ddw*gradu;
    }
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic> DDW(const Set::Matrix & /*gradu*/) const override
    {
        return ddw;
    }
    virtual void Print(std::ostream &out) const override 
    {
        out << ddw;
    }

public:
    static const KinematicVariable kinvar = KinematicVariable::gradu;

public:
    static Isotropic Random()
    {
        Isotropic ret;
        ret.Define(Util::Random(),Util::Random());
        return ret;
    }
    static void Parse(Isotropic & value, IO::ParmParse & pp)
    {
        Set::Scalar mu, lambda;
        if (pp.contains("lame") && pp.contains("shear"))
        {
            pp.query("lame",lambda);
            pp.query("shear",mu);
        }
        if (pp.contains("lambda") && pp.contains("mu"))
        {
            pp.query("lambda",lambda);
            pp.query("mu",mu);
        }
        else if (pp.contains("E") && pp.contains("nu"))
        {
            Set::Scalar E, nu;
            pp.query("E",E);
            pp.query("nu",nu);
            lambda = E * nu / (1.0 + nu) / (1.0 - 2.0*nu);
            mu = E / 2.0 / (1.0 + nu);
        }
        value.Define(mu,lambda);
    }
    AMREX_FORCE_INLINE
	void operator += (const Isotropic &rhs)
    {
        ddw += rhs.ddw;
	}    
    AMREX_FORCE_INLINE
   	Isotropic operator * (const Set::Scalar alpha) const
	{
        Isotropic ret;
        ret.ddw = ddw*alpha;
        return ret;
	}    
    friend Isotropic operator * (const Set::Scalar alpha, const Isotropic b);
    friend Isotropic operator + (const Isotropic a, const Isotropic b);
    friend Isotropic operator - (const Isotropic a, const Isotropic b);
};
AMREX_FORCE_INLINE
Isotropic operator * (const Set::Scalar alpha, const Isotropic b)
{
    Isotropic ret;
    ret.ddw = b.ddw*alpha;
    return ret;
}    
AMREX_FORCE_INLINE
Isotropic operator + (const Isotropic a, const Isotropic b)
{
    Isotropic ret;
    ret.ddw = a.ddw + b.ddw;
    return ret;
}    
AMREX_FORCE_INLINE
Isotropic operator - (const Isotropic a, const Isotropic b)
{
    Isotropic ret;
    ret.ddw = a.ddw - b.ddw;
    return ret;
}    



}
}
}

#endif

