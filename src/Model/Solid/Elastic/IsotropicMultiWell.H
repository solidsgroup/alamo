#ifndef MODEL_SOLID_ELASTIC_ISOTROPICMULTIWELL_H_
#define MODEL_SOLID_ELASTIC_ISOTROPICMULTIWELL_H_

#include "Model/Solid/Solid.H"

namespace Model
{
namespace Solid
{
namespace Elastic
{
class IsotropicMultiWell : public Solid<Set::Sym::Isotropic>
{
public:

	IsotropicMultiWell() {};
    IsotropicMultiWell(Solid<Set::Sym::Isotropic> base) : Solid<Set::Sym::Isotropic>(base) {};
    IsotropicMultiWell(Set::Scalar a_mu, Set::Scalar a_lambda, Set::Matrix a_F0=Set::Matrix::Zero()) 
    {
        Util::Message(INFO);
        Define(a_mu,a_lambda,a_F0);
    };
	virtual ~IsotropicMultiWell() {};

	void Define(Set::Scalar a_mu, Set::Scalar a_lambda, Set::Matrix a_F0)
	{
		mu = a_mu; lambda = a_lambda; F0 = a_F0;
        ddw = Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic>(lambda,mu);
	}

    Set::Scalar W(Set::Matrix F)
    {
        return ((F-F0).transpose() * (ddw*((F-F0)))).trace();
    }
    Set::Matrix DW(Set::Matrix F)
    {
        return ddw*(F-F0);
    }
    Set::Matrix4<AMREX_SPACEDIM,Set::Sym::Isotropic> DDW(Set::Matrix /*F*/)
    {
        return ddw;
    }
	
public:
    Set::Scalar mu = NAN, lambda = NAN;
    Set::Matrix F0;
    KinematicVariable kinvar = KinematicVariable::gradu;

public:
    static IsotropicMultiWell Random()
    {
        IsotropicMultiWell ret;
        ret.mu = Util::Random();
        ret.lambda = Util::Random();
        return ret;
    }
};



//template<Set::Sym sym>
//std::ostream& operator<< (std::ostream& os, const LinearElastic<sym>& b);

}
}
}

#endif

