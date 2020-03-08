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
        mu0 = a_mu0; lambda0 = a_lambda0; 
        mu = a_mu0; lambda = a_lambda0;
        F0 = a_F0;
        Isotropic::Define(a_mu, a_lambda, a_F0);
    }

    void DegradeModulus(const Set::Scalar eta)
	{
		if(eta < 0.0 || eta > 1.0) Util::Abort(INFO, "Invalid value of eta", eta);
		if(std::isnan(lambda0)) Util::Abort(INFO, "Nans in lambda0");
		if(std::isnan(mu0)) Util::Abort(INFO, "Nans in mu0");
        mu = mu0*(1.-eta); lambda = lambda0*(1.-eta);
		Isotropic::Define( mu, lambda );
	}

    void DegradeModulus(const amrex::Vector<Set::Scalar> &eta)
	{
		if(eta.size() == 0) return;
		else if(eta.size() == 1)
		{
			if(std::isnan(eta[0])) Util::Abort(INFO,"Nans in eta");
			if(std::isinf(eta[0])) Util::Abort(INFO,"Infs in eta");
			if(eta[0] < 0.0 || eta[0] > 1.0) Util::Abort(INFO,"Invalid value of eta",eta[0]);
            mu = mu0*(1.-eta[0]); lambda = lambda0*(1.-eta[0]);
			Isotropic::Define( mu, lambda );
		}
		else if(eta.size() == 2)
		{
			if(std::isnan(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Nans in eta");
			if(std::isinf(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Infs in eta");
			if(eta[0] < 0.0 || eta[0] > 1.0 || eta[1] < 0.0 || eta[1] > 1.0) Util::Abort(INFO,"Invalid value of eta");
            mu = mu0*(1.-eta[0]); lambda = lambda0*(1.-eta[1]);
			Isotropic::Define( mu, lambda );
		}
		else
			return;
	}
    void UpdateF0(Set::Matrix a_F0) { F0 = a_F0; Isotropic::UpdateF0(a_F0); }
private:
    Set::Scalar mu0 = NAN, lambda0 = NAN;
    Set::Matrix F0;

protected:
    Set::Scalar mu = NAN, lambda = NAN;

public:
    static IsotropicDegradable Random()
    {
        IsotropicDegradable ret;
        Set::Scalar mu = Util::Random();
        Set::Scalar lambda = Util::Random();
        Set::Matrix F0 = Set::Matrix::Random();
        ret.Define(mu, mu, lambda, lambda, F0);
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
};
}
}
}
#endif