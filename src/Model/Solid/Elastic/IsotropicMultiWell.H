#ifndef MODEL_SOLID_ELASTIC_ISOTROPICMULTIWELL_H_
#define MODEL_SOLID_ELASTIC_ISOTROPICMULTIWELL_H_

#include "Model/Solid/Elastic/Elastic.H"

namespace Model
{
namespace Solid
{
namespace Elastic
{
class IsotropicMultiWell : public Elastic<Set::Sym::Isotropic>
{
public:
	IsotropicMultiWell() {};
    IsotropicMultiWell(Set::Scalar a_mu, Set::Scalar a_lambda, Set::Matrix a_F0=Set::Matrix::Zero()) 
    {
        Util::Message(INFO);
        Define(a_mu,a_lambda,a_F0);
    };
	virtual ~IsotropicMultiWell() {};

	void Define(Set::Scalar a_mu, Set::Scalar a_lambda, Set::Matrix a_F0)
	{
		mu = a_mu; lambda = a_lambda; F0 = a_F0;
        ddw = Set::Matrix4<3,Set::Sym::Isotropic>(lambda,mu);
	}

    Set::Scalar W(Set::Matrix F)
    {
        return ((F-F0).transpose() * (ddw*((F-F0)))).trace();
    }
    Set::Matrix DW(Set::Matrix F)
    {
        //Util::Message(INFO,F);
        //Util::Message(INFO,F0);
        //Util::Message(INFO,ddw);
        //Util::Message(INFO,ddw*(F-F0));
        return ddw*(F-F0);
    }
    Set::Matrix4<3,Set::Sym::Isotropic> DDW(Set::Matrix F)
    {
        return ddw;
    }
	
public:
    Set::Scalar mu = NAN, lambda = NAN;
    Set::Matrix F0;

	IsotropicMultiWell operator + (const IsotropicMultiWell &rhs) const
	{
        IsotropicMultiWell ret;
        ret.mu    = mu + rhs.mu;
        ret.lambda = lambda + rhs.lambda;
        ret.w     = w + rhs.w;
        ret.dw    = dw + rhs.dw;
        ret.divdw = divdw + rhs.divdw;
        ret.ddw   = ddw; ret.ddw += rhs.ddw;
        return ret;
	}

	void operator += (const IsotropicMultiWell &rhs)
	{
        mu += rhs.mu;
        lambda += rhs.lambda;
        w += rhs.w;
        dw += rhs.dw;
        divdw += rhs.divdw;
        ddw += rhs.ddw;
	}

	IsotropicMultiWell operator - (const IsotropicMultiWell &rhs) const
	{
        IsotropicMultiWell ret;
        ret.mu = mu - rhs.mu;
        ret.lambda = lambda - rhs.lambda;
        ret.w = w - rhs.w;
        ret.dw = dw - rhs.dw;
        ret.divdw = divdw - rhs.divdw;
        ret.ddw = ddw; ret.ddw -= rhs.ddw;
        return ret;
	}

	IsotropicMultiWell & operator = (const IsotropicMultiWell &rhs)
	{
		mu    = rhs.mu; 
        lambda = rhs.lambda;
        w     = rhs.w;
        dw    = rhs.dw;
        divdw = rhs.divdw;
        ddw   = rhs.ddw;
		return *this;
	}

	virtual IsotropicMultiWell operator / (const Set::Scalar alpha) const
	{
        IsotropicMultiWell ret;
        ret.mu = mu/alpha;
        ret.lambda = lambda/alpha;
        ret.w     = w;     ret.w     /= alpha;
        ret.dw    = dw;    ret.dw    /= alpha;
        ret.divdw = divdw; ret.divdw /= alpha;
        ret.ddw   = ddw;   ret.ddw   /= alpha;
        return ret;
	}

	virtual IsotropicMultiWell operator * (const Set::Scalar alpha) const
	{
        IsotropicMultiWell ret;
        ret.mu = mu*alpha;
        ret.lambda = lambda*alpha;
        ret.w     = w;     ret.w     *= alpha;
        ret.dw    = dw;    ret.dw    *= alpha;
        ret.divdw = divdw; ret.divdw *= alpha;
        ret.ddw   = ddw;   ret.ddw   *= alpha;
        return ret;
	}

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

