

#ifndef MODEL_SOLID_LINEARELASTIC_DEGRADABLE_ISOTROPIC2_H_
#define MODEL_SOLID_LINEARELASTIC_DEGRADABLE_ISOTROPIC2_H_

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
class Isotropic2 : public Model::Solid::LinearElastic::LinearElastic
{

public:
	static Set::Scalar E10_iso, E20_iso, Tg0_iso, Ts0_iso, temp_iso;
	static Set::Scalar nu_iso;

	Isotropic2() {};

	Isotropic2(Set::Scalar _lambda, Set::Scalar _mu) : lambda(_lambda), mu(_mu) {};

	Isotropic2(Set::Scalar _E1, Set::Scalar _E2, Set::Scalar _Tg, Set::Scalar _Ts, Set::Scalar _nu, Set::Scalar _temp=298.0)
	{
		E10_iso = _E1; E20_iso = _E2; Tg0_iso = _Tg; Ts0_iso = _Ts; temp_iso = _temp; nu_iso = _nu;
		Set::Scalar E = 0.5*(E10_iso + E20_iso) - 0.5*(E10_iso-E20_iso)*tanh((temp_iso - Tg0_iso)/Ts0_iso);
		lambda = E*nu_iso/((1.+nu_iso)*(1.-2.*nu_iso));
		mu = E/(2.*(1.+nu_iso));
	}

	~Isotropic2() {} ;
	virtual Set::Matrix operator () (Set::Matrix &eps) const
	{
		return mu*(eps + eps.transpose()) + lambda*Set::Matrix::Identity()*eps.trace();
	};
	virtual Set::Vector operator () (std::array<Set::Matrix,AMREX_SPACEDIM> &gradeps)
	{
		Set::Vector f = Set::Vector::Zero();
		for (int i = 0; i < AMREX_SPACEDIM; i++)
			for (int j=0; j < AMREX_SPACEDIM; j++)
				f(i) += mu*(gradeps[i](j,j) + gradeps[j](i,j))  + lambda*gradeps[j](j,i);
		return f;
	}

	void Randomize()
	{
		lambda = Util::Random();
		mu = Util::Random();
		Set::Scalar eta = Util::Random();
		DegradeModulus(eta);
		
	};

	Isotropic2 operator + (const Isotropic2 &rhs) const
	{return Isotropic2(mu+rhs.mu, lambda+rhs.lambda);}

	void operator += (const Isotropic2 &rhs)
	{lambda+=rhs.lambda; mu+=rhs.mu;}

	Isotropic2 operator - (const Isotropic2 &rhs) const
	{return Isotropic2(lambda-rhs.lambda, mu-rhs.mu);}

	Isotropic2 operator * (const Isotropic2 &rhs) const
	{return Isotropic2(lambda*rhs.lambda, mu*rhs.mu);}

	Isotropic2 operator / (const Isotropic2 &rhs) const
	{return Isotropic2(lambda/rhs.lambda, mu/rhs.mu);}

	virtual Isotropic2 operator / (const Set::Scalar alpha) const
	{return Isotropic2(lambda/alpha, mu/alpha);}

	virtual Isotropic2 operator * (const Set::Scalar alpha) const
	{return Isotropic2(lambda*alpha, mu*alpha);}

	void DegradeModulus(const Set::Scalar eta)
	{
		if(eta < 0.0 || eta > 1.0) Util::Abort(INFO, "Invalid value of eta", eta);
		Set::Scalar E1 = E10_iso*(1.-eta); //Tg = E20*(1.-eta2); Ts = Ts0*(1.+eta3);
		Set::Scalar E = 0.5*(E1 + E20_iso) - 0.5*(E1-E20_iso)*tanh((temp_iso - Tg0_iso)/Ts0_iso);
		lambda = E*nu_iso/((1.+nu_iso)*(1.-2.*nu_iso));
		mu = E/(2.*(1.+nu_iso));
	}

	void DegradeModulus(const Set::Scalar eta1, const Set::Scalar eta2, const Set::Scalar eta3)
	{
		if(eta1 < 0.0 || eta1 > 1.0) Util::Abort(INFO,"Invalid value of eta1");
		if(eta2 < 0.0 || eta2 > 1.0) Util::Abort(INFO,"Invalid value of eta2");
		if(eta3 < 0.0 || eta3 > 1.0) Util::Abort(INFO,"Invalid value of eta3");
		
		Set::Scalar E1 = E10_iso*(1.-eta1), Tg = Tg0_iso*(1.-eta2), Ts = Ts0_iso*(1.+eta3);
		Set::Scalar E = 0.5*(E1 + E20_iso) - 0.5*(E1-E20_iso)*tanh((temp_iso - Tg)/Ts);
		lambda = E*nu_iso/((1.+nu_iso)*(1.-2.*nu_iso));
		mu = E/(2.*(1.+nu_iso));
		if(std::isnan(lambda)) 
		{
			Util::Message(INFO, "E10_iso = ", E10_iso, ". E20_iso = ", E20_iso, ". Tg0_iso = ", Tg0_iso, ". Ts0_iso = ", Ts0_iso);
			Util::Message(INFO, "nu_iso = ", nu_iso, ". temp_iso = ", temp_iso);
			Util::Abort(INFO, "Nan in lambda, E = ",E, ", nu = ", mu);
		}
		if(std::isnan(mu)) Util::Abort(INFO, "Nan in mu");
	}

	virtual void Print (std::ostream& os) const
	{
		os <<    "mu = " << mu <<
			" lambda = " << lambda << std::endl;
	}


private:
	
	Set::Scalar lambda = NAN, mu = NAN;
};

}
}
}
}
#endif


