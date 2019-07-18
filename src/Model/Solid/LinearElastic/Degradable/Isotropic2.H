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
	Isotropic2() { } ;
	Isotropic2(Set::Scalar _lambda, Set::Scalar _mu)
	{
		Set::Scalar lambda = _lambda;
		Set::Scalar mu = _mu;
		Set::Scalar E = mu*(3.*lambda + 2*mu)/(lambda + mu);
		nu = lambda/(2*(lambda+mu));
		Tg = 319.0; //some default value
		Ts = 17.0; //some default value
		temp = 298.0;
		E1 = E/0.928;
		E2 = 0.1*E1;
	}
	Isotropic2(Set::Scalar _E1, Set::Scalar _E2, Set::Scalar _Tg, Set::Scalar _Ts, Set::Scalar _nu, Set::Scalar _temp=298.0)
	{
		E10 = _E1; E20 = _E2; Tg0 = _Tg; Ts0 = _Ts; temp = _temp; nu = _nu;

		E1 = _E1; E2 = _E2; Tg = _Tg; Ts = _Ts;

		//E = 0.5*(E1 + E2) - 0.5*(E1-E2)*tanh((temp - Tg)/Ts);
		//lambda = E*nu/((1.+nu)*(1.-2.*nu));
		//mu = E/(2.*(1.+nu));
	}
	~Isotropic2() {} ;
	virtual Set::Matrix operator () (Set::Matrix &eps) const
	{
		Set::Scalar E = 0.5*(E1 + E2) - 0.5*(E1-E2)*tanh((temp - Tg)/Ts);
		Set::Scalar lambda = E*nu/((1.+nu)*(1.-2.*nu));
		Set::Scalar mu = E/(2.*(1.+nu));
		return mu*(eps + eps.transpose()) + lambda*Set::Matrix::Identity()*eps.trace();
	};
	virtual Set::Vector operator () (std::array<Set::Matrix,AMREX_SPACEDIM> &gradeps)
	{
		Set::Scalar E = 0.5*(E1 + E2) - 0.5*(E1-E2)*tanh((temp - Tg)/Ts);
		Set::Scalar lambda = E*nu/((1.+nu)*(1.-2.*nu));
		Set::Scalar mu = E/(2.*(1.+nu));

		Set::Vector f = Set::Vector::Zero();
		for (int i = 0; i < AMREX_SPACEDIM; i++)
			for (int j=0; j < AMREX_SPACEDIM; j++)
				f(i) += mu*(gradeps[i](j,j) + gradeps[j](i,j))  + lambda*gradeps[j](j,i);
		return f;
	}

	void Randomize()
	{
		E10 = Util::Random();
		E1 = E10;
		E2 = Util::Random();
		temp = Util::Random();
		Tg = Util::Random();
		Ts = Util::Random();
		Set::Scalar eta = Util::Random();
		DegradeModulus(eta);
		
	};

	Isotropic2 operator + (const Isotropic2 &rhs) const
	{return Isotropic2(E10+rhs.E10, E20+rhs.E20, Tg, Ts, nu, temp);}

	void operator += (const Isotropic2 &rhs)
	{E1+=rhs.E1; E2+=rhs.E2; E10+=rhs.E10; E20+=rhs.E20;}

	Isotropic2 operator - (const Isotropic2 &rhs) const
	{return Isotropic2(E10-rhs.E10, E20-rhs.E20, Tg, Ts, nu, temp);}

	Isotropic2 operator * (const Isotropic2 &rhs) const
	{return Isotropic2(E10*rhs.E10, E20*rhs.E20, Tg, Ts, nu, temp);}

	Isotropic2 operator / (const Isotropic2 &rhs) const
	{return Isotropic2(E10/rhs.E10, E20/rhs.E20, Tg, Ts, nu, temp);}

	virtual Isotropic2 operator / (const Set::Scalar alpha) const
	{return Isotropic2(E10/alpha, E20/alpha, Tg, Ts, nu, temp);}

	virtual Isotropic2 operator * (const Set::Scalar alpha) const
	{return Isotropic2(E10*alpha, E20*alpha, Tg, Ts, nu, temp);}

	void DegradeModulus(const Set::Scalar eta)
	{
		if(eta < 0.0 || eta > 1.0) Util::Abort(INFO, "Invalid value of eta", eta);
		E1 = E10*(1.-eta); //Tg = E20*(1.-eta2); Ts = Ts0*(1.+eta3);
		//E = 0.5*(E1 + E2) - 0.5*(E1-E2)*tanh((temp - Tg)/Ts);
		//lambda = E*nu/((1.+nu)*(1.-2.*nu));
		//mu = E/(2.*(1.+nu));
	}

	void DegradeModulus(const Set::Scalar eta1, const Set::Scalar eta2, const Set::Scalar eta3)
	{
		if(eta1 < 0.0 || eta1 > 1.0) Util::Abort(INFO,"Invalid value of eta1");
		if(eta2 < 0.0 || eta2 > 1.0) Util::Abort(INFO,"Invalid value of eta2");
		if(eta3 < 0.0 || eta3 > 1.0) Util::Abort(INFO,"Invalid value of eta3");
		
		E1 = E10*(1.-eta1); Tg = Tg0*(1.-eta2); Ts = Ts0*(1.+eta3);
		//E = 0.5*(E1 + E2) - 0.5*(E1-E2)*tanh((temp - Tg)/Ts);
		//lambda = E*nu/((1.+nu)*(1.-2.*nu));
		//mu = E/(2.*(1.+nu));
	}

	virtual void Print (std::ostream& os) const
	{
		os <<
			" Tg = " << Tg <<
			" Ts = " << Ts <<
			" nu = " << nu <<
			" E10 = " << E10 <<
			" E20 = " << E10 <<
			" Tg0 = " << Tg0 <<
			" Ts0 = " << Ts0;
	}


private:
	Set::Scalar E10 = NAN, E20 = NAN, Tg0 = NAN, Ts0 = NAN, temp = 298.0;
	Set::Scalar E1 = NAN, E2 = NAN, Tg = NAN, Ts = NAN;
	Set::Scalar nu = NAN;
};

}
}
}
}
#endif


