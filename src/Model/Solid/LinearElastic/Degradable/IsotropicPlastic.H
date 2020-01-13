#ifndef MODEL_SOLID_LINEARELASTIC_DEGRADABLE_ISOTROPICPLASTIC_H_
#define MODEL_SOLID_LINEARELASTIC_DEGRADABLE_ISOTROPICPLASTIC_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>

#include "Set/Set.H"
#include "Model/Solid/Solid.H"
#include "Util/Util.H"

namespace Model
{
namespace Solid
{
namespace LinearElastic
{
namespace Degradable
{
class IsotropicPlastic : public Model::Solid::LinearElastic::LinearElastic<Set::Sym::Isotropic>
{
public:
	IsotropicPlastic() { } ;
	IsotropicPlastic(Set::Scalar _mu, Set::Scalar _lambda, Set::Scalar _yield, Set::Scalar _hardening, Set::Matrix _epsp=Set::Matrix::Zero()) 
	{
		mu = _mu; mu0 = _mu;
		lambda = _lambda; lambda0 = _lambda;
		epsp = _epsp;
        yield_strength = _yield; yield_strength0 = _yield;
        hardening_modulus = _hardening; hardening_modulus0 = _hardening;
		
		for (int p = 0; p < AMREX_SPACEDIM; p++) 
			gradEpsp[p] = Set::Matrix::Zero();
    } ;
	~IsotropicPlastic() {} ;

    virtual Set::Matrix operator () (Set::Matrix &eps,bool a_homogeneous=true) const
	{
		if (!a_homogeneous) eps -= epsp;
		return mu*(eps + eps.transpose()) + lambda*Set::Matrix::Identity()*eps.trace();
	};
	virtual Set::Vector operator () (std::array<Set::Matrix,AMREX_SPACEDIM> &gradeps,bool a_homogeneous=true)
	{
		Set::Vector f = Set::Vector::Zero();
		if (!a_homogeneous) for (int i = 0; i < AMREX_SPACEDIM; i++) gradeps[i] -= gradEpsp[i];

		for (int i = 0; i < AMREX_SPACEDIM; i++)
			for (int j=0; j < AMREX_SPACEDIM; j++)
				f(i) += mu*(gradeps[i](j,j) + gradeps[j](i,j))  + lambda*gradeps[j](j,i);
		return f;
	}

	void Randomize()
	{
		mu = Util::Random();
		lambda = Util::Random();
		mu0 = mu; lambda0 = lambda;
        yield_strength = Util::Random();
        hardening_modulus = Util::Random();
        yield_strength0 = yield_strength; hardening_modulus0 = hardening_modulus;

		amrex::Vector<Set::Scalar> eta(2);
		eta[0] = Util::Random();
		eta[1] = Util::Random();
		DegradeModulus(eta);
		
	};

	IsotropicPlastic operator + (const IsotropicPlastic &rhs) const
	{
		IsotropicPlastic ret(mu+rhs.mu, lambda+rhs.lambda, yield_strength+rhs.yield_strength, hardening_modulus-rhs.hardening_modulus, epsp+rhs.epsp);
		for (int i = 0; i < AMREX_SPACEDIM; i++) ret.gradEpsp[i] = gradEpsp[i] + rhs.gradEpsp[i];
		return ret;
	}

	void operator += (const IsotropicPlastic &rhs)
	{
		mu+=rhs.mu; lambda+=rhs.lambda;
		epsp += rhs.epsp;
		mu0 += rhs.mu0; lambda0 += rhs.lambda0;
        yield_strength += rhs.yield_strength; yield_strength0 += rhs.yield_strength0;
        hardening_modulus += rhs.hardening_modulus; hardening_modulus0 += rhs.hardening_modulus0;

		for (int i = 0; i < AMREX_SPACEDIM; i++) gradEpsp[i] += rhs.gradEpsp[i];	
	}

	IsotropicPlastic & operator = (const IsotropicPlastic &rhs)
	{
		mu=rhs.mu; lambda=rhs.lambda; epsp = rhs.epsp;
		mu0 = rhs.mu0; lambda0 = rhs.lambda0;
        yield_strength = rhs.yield_strength; yield_strength0 = rhs.yield_strength0;
        hardening_modulus = rhs.hardening_modulus; hardening_modulus0 = rhs.hardening_modulus0;

		for (int i = 0; i < AMREX_SPACEDIM; i++) gradEpsp[i] = rhs.gradEpsp[i];
		return *this;
	}

	IsotropicPlastic operator - (const IsotropicPlastic &rhs) const
	{
		IsotropicPlastic ret(mu-rhs.mu, lambda-rhs.lambda, yield_strength-rhs.yield_strength, hardening_modulus-rhs.hardening_modulus, epsp-rhs.epsp);
		for (int i = 0; i < AMREX_SPACEDIM; i++) ret.gradEpsp[i] = gradEpsp[i] - rhs.gradEpsp[i];
		return ret;
	}

	virtual IsotropicPlastic operator / (const Set::Scalar alpha) const
	{
		IsotropicPlastic ret(mu/alpha, lambda/alpha, yield_strength/alpha, hardening_modulus/alpha, epsp/alpha);
		for (int i = 0; i < AMREX_SPACEDIM; i++) ret.gradEpsp[i] = gradEpsp[i]/alpha;
		return ret;
	}

	virtual IsotropicPlastic operator * (const Set::Scalar alpha) const
	{
		IsotropicPlastic ret(mu*alpha, lambda*alpha, yield_strength*alpha, hardening_modulus*alpha, epsp*alpha);
		for (int i = 0; i < AMREX_SPACEDIM; i++) ret.gradEpsp[i] = gradEpsp[i]*alpha;
		return ret;
	}

	void DegradeModulus(const Set::Scalar eta)
	{
		if(eta < 0.0 || eta > 1.0) Util::Abort(INFO, "Invalid value of eta", eta);
		if(std::isnan(lambda0)) Util::Abort(INFO, "Nans in lambda0");
		if(std::isnan(mu0)) Util::Abort(INFO, "Nans in mu0");

		lambda = lambda0*(1.0 - eta);
		mu = mu0*(1.0-eta);
		
	}

	void DegradeModulus(const amrex::Vector<Set::Scalar> &eta)
	{
		if(eta.size() == 0) return;
		else if(eta.size() == 1)
		{
			if(std::isnan(eta[0])) Util::Abort(INFO,"Nans in eta");
			if(std::isinf(eta[0])) Util::Abort(INFO,"Infs in eta");
			if(eta[0] < 0.0 || eta[0] > 1.0) Util::Abort(INFO,"Invalid value of eta",eta[0]);
			lambda = lambda0*(1.0-eta[0]);
			mu = mu0*(1.0-eta[0]);
		}
		else if(eta.size() == 2)
		{
			if(std::isnan(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Nans in eta");
			if(std::isinf(eta[0]) || std::isnan(eta[1])) Util::Abort(INFO,"Infs in eta");
			if(eta[0] < 0.0 || eta[0] > 1.0 || eta[1] < 0.0 || eta[1] > 1.0) Util::Abort(INFO,"Invalid value of eta");
			lambda = lambda0*(1.0-eta[0]);
			mu = mu0*(1.0-eta[1]);
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

	Set::Scalar GetMu() {return mu0;}
	Set::Scalar GetLambda() {return lambda0;}

	virtual void Print (std::ostream& os) const
	{
		os <<    "mu = " << mu <<
			" lambda = " << lambda <<
			" mu0=" << mu0 <<
			" lambda0=" << lambda0 << std::endl;
	}


private:
	Set::Scalar mu = NAN, lambda = NAN, mu0 = NAN, lambda0 = NAN;
    Set::Scalar yield_strength0 = NAN, yield_strength = NAN;
    Set::Scalar hardening_modulus0 = NAN, hardening_modulus = NAN;
public:
	Set::Matrix epsp = Set::Matrix::Zero();
	std::array<Set::Matrix,AMREX_SPACEDIM> gradEpsp;
};

}
}
}
}
#endif


