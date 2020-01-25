#ifndef MODEL_SOLID_LINEAR_ISOTROPICDEGRADABLEPLASTIC_H_
#define MODEL_SOLID_LINEAR_ISOTROPICDEGRADABLEPLASTIC_H_

#include "AMReX.H"
#include "Model/Solid/Solid.H"
#include "Model/Solid/Affine/Isotropic.H"

namespace Model
{
namespace Solid
{
namespace Affine
{
class IsotropicDegradablePlastic : public Isotropic
{
public:
    IsotropicDegradablePlastic() {};
    IsotropicDegradablePlastic(Solid<Set::Sym::Isotropic> base) : Isotropic(base) {};
    IsotropicDegradablePlastic(Isotropic base) : Isotropic(base) {}
    IsotropicDegradablePlastic(Set::Scalar a_mu, Set::Scalar a_lambda, Set::Scalar a_yield, Set::Scalar a_hardening)
    {
        Define(a_mu, a_mu, a_lambda, a_lambda, a_yield, a_yield, a_hardening, a_hardening);
    }
    
    void Define(Set::Scalar a_mu, Set::Scalar a_mu0, Set::Scalar a_lambda, Set::Scalar a_lambda0, Set::Scalar a_yield, Set::Scalar a_yield0, Set::Scalar a_hardening, Set::Scalar a_hardening0, Set::Matrix F0=Set::Matrix::Zero())
    {
        mu0 = a_mu0; lambda0 = a_lambda0; 
        yield_strength = a_yield;
        yield_strength0 = a_yield0;
        hardening_modulus0 = a_hardening0;
        hardening_modulus = a_hardening;
        Isotropic::Define(a_mu, a_lambda, F0);
    }

    void DegradeModulus(const Set::Scalar eta)
	{
		if(eta < 0.0 || eta > 1.0) Util::Abort(INFO, "Invalid value of eta", eta);
		if(std::isnan(lambda0)) Util::Abort(INFO, "Nans in lambda0");
		if(std::isnan(mu0)) Util::Abort(INFO, "Nans in mu0");

		Isotropic::Define( mu0*(1.-eta), lambda0*(1.-eta) );
	}

    void DegradeModulus(const amrex::Vector<Set::Scalar> &eta)
	{
		if(eta.size() == 0) return;
		else if(eta.size() == 1)
		{
			if(std::isnan(eta[0])) Util::Abort(INFO,"Nans in eta");
			if(std::isinf(eta[0])) Util::Abort(INFO,"Infs in eta");
			if(eta[0] < 0.0 || eta[0] > 1.0) Util::Abort(INFO,"Invalid value of eta",eta[0]);
			Isotropic::Define( mu0*(1.-eta[0]), lambda0*(1.-eta[0]) );
		}
		else if(eta.size() == 2)
		{
			if(std::isnan(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Nans in eta");
			if(std::isinf(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Infs in eta");
			if(eta[0] < 0.0 || eta[0] > 1.0 || eta[1] < 0.0 || eta[1] > 1.0) Util::Abort(INFO,"Invalid value of eta");
			Isotropic::Define( mu0*(1.-eta[0]), lambda0*(1.-eta[1]) );
		}
		else
			return;
	}

    void DegradeYieldSurface(const Set::Scalar eta)
	{
		if(eta < 0.0 || eta > 1.0) Util::Abort(INFO, "Invalid value of eta", eta);
		if(std::isnan(eta)) Util::Abort(INFO, "Nans in eta");
        if(std::isinf(eta)) Util::Abort(INFO, "Infs in eta");
		yield_strength = yield_strength0*(1. - eta);
        hardening_modulus = hardening_modulus0*(1.-eta);
	}

    void DegradeYieldSurface(const amrex::Vector<Set::Scalar> &eta)
    {
        if (eta.size() == 0) return;
        else if (eta.size() == 1)
        {
            if(std::isnan(eta[0])) Util::Abort(INFO,"Nans in eta");
			if(std::isinf(eta[0])) Util::Abort(INFO,"Infs in eta");
			if(eta[0] < 0.0 || eta[0] > 1.0) Util::Abort(INFO,"Invalid value of eta",eta[0]);
			yield_strength = yield_strength0*(1. - eta[0]);
            hardening_modulus = hardening_modulus0*(1.-eta[0]);
        }
        else if (eta.size() == 2)
        {
            if(std::isnan(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Nans in eta");
			if(std::isinf(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Infs in eta");
			if(eta[0] < 0.0 || eta[0] > 1.0 || eta[1] < 0.0 || eta[1] > 1.0) Util::Abort(INFO,"Invalid value of eta");
            yield_strength = yield_strength0*(1. - eta[0]);
            hardening_modulus = hardening_modulus0*(1.-eta[1]);
        }
        else return;
    }

    Set::Scalar YieldSurface (Set::Scalar alpha)
    {
        return yield_strength + hardening_modulus*alpha;
    }

    Set::Scalar OriginalYieldSurface (Set::Scalar alpha)
    {
        return yield_strength0 + hardening_modulus0*alpha;
    }

private:
    Set::Scalar mu0 = NAN, lambda0 = NAN;
    Set::Scalar yield_strength0 = NAN, yield_strength = NAN;
    Set::Scalar hardening_modulus0 = NAN, hardening_modulus = NAN;

public:
    static IsotropicDegradablePlastic Random()
    {
        IsotropicDegradablePlastic ret;
        Set::Scalar mu = Util::Random();
        Set::Scalar lambda = Util::Random();
        Set::Matrix F0 = Set::Matrix::Random();
        Set::Scalar yield_strength = Util::Random();
        Set::Scalar hardening_modulus = Util::Random();

        ret.Define(mu, mu, lambda, lambda, yield_strength, yield_strength, hardening_modulus, hardening_modulus, F0);
        return ret;
    }

    static void Parse(IsotropicDegradablePlastic & value, IO::ParmParse & pp)
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

        Set::Scalar yield = 1.0, hardening = 1.0;
        pp.query("yield", yield);
        pp.query("hardening", hardening);
        value.Define(mu,mu,lambda,lambda,yield,yield,hardening,hardening);
    }
};
}
}
}
#endif