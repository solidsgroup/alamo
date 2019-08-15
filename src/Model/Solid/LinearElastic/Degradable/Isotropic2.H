

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

	Isotropic2(Set::Scalar _E1, Set::Scalar _E2, Set::Scalar _Tg, Set::Scalar _Ts) : E1(_E1), E2(_E2), Tg(_Tg), Ts(_Ts) {};
	
	Isotropic2(Set::Scalar _lambda, Set::Scalar _mu)
	{
		Set::Scalar lambda = _lambda;
		Set::Scalar mu = _mu;
		Set::Scalar E = mu*(3.*lambda + 2*mu)/(lambda + mu);
		Set::Scalar nu = lambda/(2*(lambda+mu));
		Tg = 319.0; //some default value
		Ts = 17.0; //some default value
		E1 = E/0.928;
		E2 = 0.1*E1;
	}
	Isotropic2(Set::Scalar _E1, Set::Scalar _E2, Set::Scalar _Tg, Set::Scalar _Ts, Set::Scalar _nu, Set::Scalar _temp=298.0)
	{
		E10_iso = _E1; E20_iso = _E2; Tg0_iso = _Tg; Ts0_iso = _Ts; temp_iso = _temp; nu_iso = _nu;

		E1 = _E1; E2 = _E2; Tg = _Tg; Ts = _Ts;

		Util::Message(INFO, "E10_iso = ", E10_iso,
						". E20_iso = ", E20_iso, 
						". Tg0_iso = ", Tg0_iso, 
						". Ts0_iso = ", Ts0_iso, 
						". temp_iso = ", temp_iso,
						". nu_iso = ", nu_iso);

		//E = 0.5*(E1 + E2) - 0.5*(E1-E2)*tanh((temp - Tg)/Ts);
		//lambda = E*nu/((1.+nu)*(1.-2.*nu));
		//mu = E/(2.*(1.+nu));
	}

	~Isotropic2() {} ;
	virtual Set::Matrix operator () (Set::Matrix &eps) const
	{
		Set::Scalar E = 0.5*(E1 + E2) - 0.5*(E1-E2)*tanh((temp_iso - Tg)/Ts);
		Set::Scalar lambda = E*nu_iso/((1.+nu_iso)*(1.-2.*nu_iso));
		Set::Scalar mu = E/(2.*(1.+nu_iso));
		return mu*(eps + eps.transpose()) + lambda*Set::Matrix::Identity()*eps.trace();
	};
	virtual Set::Vector operator () (std::array<Set::Matrix,AMREX_SPACEDIM> &gradeps)
	{
		Set::Scalar E = 0.5*(E1 + E2) - 0.5*(E1-E2)*tanh((temp_iso - Tg)/Ts);
		Set::Scalar lambda = E*nu_iso/((1.+nu_iso)*(1.-2.*nu_iso));
		Set::Scalar mu = E/(2.*(1.+nu_iso));

		Set::Vector f = Set::Vector::Zero();
		for (int i = 0; i < AMREX_SPACEDIM; i++)
			for (int j=0; j < AMREX_SPACEDIM; j++)
				f(i) += mu*(gradeps[i](j,j) + gradeps[j](i,j))  + lambda*gradeps[j](j,i);
		return f;
	}

	void Randomize()
	{
		E1 = Util::Random();
		E2 = Util::Random();
		//temp = Util::Random();
		Tg = Util::Random();
		Ts = Util::Random();
		Set::Scalar eta = Util::Random();
		DegradeModulus(eta);
		
	};

	Isotropic2 operator + (const Isotropic2 &rhs) const
	{return Isotropic2(E1+rhs.E1, E2+rhs.E2, Tg, Ts);}

	void operator += (const Isotropic2 &rhs)
	{E1+=rhs.E1; E2+=rhs.E2;}

	Isotropic2 operator - (const Isotropic2 &rhs) const
	{return Isotropic2(E1-rhs.E1, E2-rhs.E2, Tg, Ts);}

	Isotropic2 operator * (const Isotropic2 &rhs) const
	{return Isotropic2(E1*rhs.E1, E2*rhs.E2, Tg, Ts);}

	Isotropic2 operator / (const Isotropic2 &rhs) const
	{return Isotropic2(E1/rhs.E1, E2/rhs.E2, Tg, Ts);}

	virtual Isotropic2 operator / (const Set::Scalar alpha) const
	{return Isotropic2(E1/alpha, E2/alpha, Tg, Ts);}

	virtual Isotropic2 operator * (const Set::Scalar alpha) const
	{return Isotropic2(E1*alpha, E2*alpha, Tg, Ts);}

	void DegradeModulus(const Set::Scalar eta)
	{
		if(eta < 0.0 || eta > 1.0) Util::Abort(INFO, "Invalid value of eta", eta);
		E1 = E10_iso*(1.-eta); //Tg = E20*(1.-eta2); Ts = Ts0*(1.+eta3);
		//E = 0.5*(E1 + E2) - 0.5*(E1-E2)*tanh((temp - Tg)/Ts);
		//lambda = E*nu/((1.+nu)*(1.-2.*nu));
		//mu = E/(2.*(1.+nu));
	}

	void DegradeModulus(const Set::Scalar eta1, const Set::Scalar eta2, const Set::Scalar eta3)
	{
		if(eta1 < 0.0 || eta1 > 1.0) Util::Abort(INFO,"Invalid value of eta1");
		if(eta2 < 0.0 || eta2 > 1.0) Util::Abort(INFO,"Invalid value of eta2");
		if(eta3 < 0.0 || eta3 > 1.0) Util::Abort(INFO,"Invalid value of eta3");
		
		E1 = E10_iso*(1.-eta1); Tg = Tg0_iso*(1.-eta2); Ts = Ts0_iso*(1.+eta3);
		//E = 0.5*(E1 + E2) - 0.5*(E1-E2)*tanh((temp - Tg)/Ts);
		//lambda = E*nu/((1.+nu)*(1.-2.*nu));
		//mu = E/(2.*(1.+nu));
	}

	virtual void Print (std::ostream& os) const
	{
		os <<
			" Tg = " << Tg <<
			" Ts = " << Ts <<
			" nu = " << nu_iso <<
			" E1 = " << E1 <<
			" E2 = " << E2;
	}


private:
	
	Set::Scalar E1 = NAN, E2 = NAN, Tg = NAN, Ts = NAN;
};

}
}
}
}
#endif


