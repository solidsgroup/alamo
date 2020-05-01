#ifndef MODEL_SOLID_AFFINE_ISOTROPIC_H_
#define MODEL_SOLID_AFFINE_ISOTROPIC_H_

#include "AMReX.H"
#include "Model/Solid/Solid.H"
#include "IO/ParmParse.H"

namespace Model
{
namespace Solid
{
namespace Affine
{
class Isotropic : public Solid<Set::Sym::Isotropic>
{
public:

	Isotropic() {};
    Isotropic(Solid<Set::Sym::Isotropic> base) : Solid<Set::Sym::Isotropic>(base) {};
    Isotropic(Set::Scalar a_mu, Set::Scalar a_lambda, Set::Matrix a_F0=Set::Matrix::Zero()) 
    {
        Define(a_mu,a_lambda,a_F0);
    };
	virtual ~Isotropic() {};

	void Define(Set::Scalar a_mu, Set::Scalar a_lambda, Set::Matrix a_F0)
	{
		F0 = a_F0;
        ddw = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic>(a_lambda,a_mu);
	}

    Set::Scalar W(const Set::Matrix & F) const override
    {
        return ((F-F0).transpose() * (ddw*((F-F0)))).trace();
    }
    Set::Matrix DW(const Set::Matrix & F) const override
    {
        return ddw*(F-F0);
    }
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic> DDW(const Set::Matrix & /*F*/) const override
    {
        return ddw;
    }
	
public:
    Set::Matrix F0;
    KinematicVariable kinvar = KinematicVariable::gradu;

public:
    static Isotropic Random()
    {
        Isotropic ret;
        Set::Scalar mu = Util::Random();
        Set::Scalar lambda = Util::Random();
        Set::Matrix F0 = Set::Matrix::Random();
        ret.Define(mu,lambda,F0);
        return ret;
    }
    static void Parse(Isotropic & value, IO::ParmParse & pp)
    {
        Set::Scalar mu, lambda;
        Set::Matrix F0 = Set::Matrix::Zero();
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
        if (pp.contains("F0"))
        {
            pp.queryarr("F0",F0);
        }
        value.Define(mu,lambda,F0);
    }
};



//template<Set::Sym sym>
//std::ostream& operator<< (std::ostream& os, const LinearElastic<sym>& b);

}
}
}

#endif

