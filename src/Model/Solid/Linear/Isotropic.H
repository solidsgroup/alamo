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
		//mu = a_mu; lambda = a_lambda;
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
};


}
}
}

#endif

