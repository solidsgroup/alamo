#ifndef MODEL_SOLID_LINEARELASTIC_DEGRADABLE_ISOTROPIC_H_
#define MODEL_SOLID_LINEARELASTIC_DEGRADABLE_ISOTROPIC_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>

#include "Set/Set.H"
#include "Model/Solid/Solid.H"
#include "Model/Solid/LinearElastic/Isotropic.H"
#include "Util/Util.H"

namespace Model
{
namespace Solid
{
namespace LinearElastic
{
namespace Degradable
{
class Isotropic : public Model::Solid::LinearElastic::LinearElastic<Set::Sym::Isotropic>
{
public:
	Isotropic() { } ;
	Isotropic(Set::Scalar _mu, Set::Scalar _lambda) : mu(_mu) , lambda(_lambda), mu0(_mu), lambda0(_lambda) {} ;
	~Isotropic() {} ;
	virtual Set::Matrix operator () (Set::Matrix &eps,bool =true) const
	{
		return mu*(eps + eps.transpose()) + lambda*Set::Matrix::Identity()*eps.trace();
	};
	virtual Set::Vector operator () (Set::Matrix3 &gradeps,bool =true)
	{
		Set::Vector f = Set::Vector::Zero();
		for (int i = 0; i < AMREX_SPACEDIM; i++)
			for (int j=0; j < AMREX_SPACEDIM; j++)
				f(i) += mu*(gradeps(i,j,j) + gradeps(j,i,j))  + lambda*gradeps(j,j,i);
		return f;
	}

	void Randomize()
	{
		mu = Util::Random();
		lambda = Util::Random();
		mu0 = mu;
		lambda0 = lambda;
		amrex::Vector<Set::Scalar> eta(2);
		eta[0] = Util::Random();
		eta[1] = Util::Random();
		DegradeModulus(eta);
		
	};

	Isotropic operator + (const Isotropic &rhs) const
	{return Isotropic(mu+rhs.mu, lambda+rhs.lambda);}

	void operator += (const Isotropic &rhs)
	{mu+=rhs.mu; lambda+=rhs.lambda;}

	Isotropic operator - (const Isotropic &rhs) const
	{return Isotropic(mu-rhs.mu, lambda-rhs.lambda);}

	Isotropic operator * (const Isotropic &rhs) const
	{return Isotropic(mu*rhs.mu, lambda*rhs.lambda);}

	Isotropic operator / (const Isotropic &rhs) const
	{return Isotropic(mu/rhs.mu, lambda/rhs.lambda);}

	virtual Isotropic operator / (const Set::Scalar alpha) const
	{return Isotropic(mu/alpha, lambda/alpha);}

	virtual Isotropic operator * (const Set::Scalar alpha) const
	{return Isotropic(mu*alpha, lambda*alpha);}

	void DegradeModulus(const Set::Scalar eta)
	{
		if(eta < 0.0 || eta > 1.0) Util::Abort(INFO, "Invalid value of eta", eta);

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

	virtual void Print (std::ostream& os) const
	{
		os <<    "mu = " << mu <<
			" lambda = " << lambda <<
			" mu0=" << mu0 <<
			" lambda0=" << lambda0 << std::endl;
	}


private:
	Set::Scalar mu = NAN, lambda = NAN, mu0 = NAN, lambda0 = NAN;
};

}
}
}
}
#endif
