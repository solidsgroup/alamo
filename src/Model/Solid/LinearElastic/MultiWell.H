#ifndef MODEL_SOLID_LINEARELASTIC_MULTIWELL_H_
#define MODEL_SOLID_LINEARELASTIC_MULTIWELL_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>

#include "Util/Util.H"
#include "Set/Set.H"
#include "Model/Solid/LinearElastic/LinearElastic.H"

namespace Model
{
namespace Solid
{
namespace LinearElastic
{
class Multiwell : public Model::Solid::LinearElastic::LinearElastic<Set::Sym::Isotropic>
{
public:
	Multiwell() {} ;
	Multiwell(Set::Scalar _mu, Set::Scalar _lambda, Set::Matrix _Fgb=Set::Matrix::Zero()) {Define(_mu,_lambda,_Fgb);} ;
	~Multiwell() {} ;
	void Define(Set::Scalar _mu, Set::Scalar _lambda,Set::Matrix _Fgb)
	{
		mu = _mu; lambda = _lambda; Fgb = _Fgb;
		for (int p = 0; p < AMREX_SPACEDIM; p++) gradFgb[p] = Set::Matrix::Zero();
	}

	AMREX_FORCE_INLINE
	virtual Set::Matrix operator () (Set::Matrix &gradu,bool a_homogeneous=true) const
	{
		Set::Matrix eps = gradu;
		if (!a_homogeneous) eps -= Fgb;
		return mu*(eps + eps.transpose()) + lambda*Set::Matrix::Identity()*eps.trace();
	};

	AMREX_FORCE_INLINE
	virtual Set::Vector operator () (std::array<Set::Matrix,AMREX_SPACEDIM> &gradgradu,bool a_homogeneous=true)
	{
		std::array<Set::Matrix,AMREX_SPACEDIM> gradeps;
		for (int i = 0; i < AMREX_SPACEDIM; i++) gradeps[i] = gradgradu[i];
		if (!a_homogeneous) for (int i = 0; i < AMREX_SPACEDIM; i++) gradeps[i] -= gradFgb[i];
		Set::Vector f = Set::Vector::Zero();
		for (int i = 0; i < AMREX_SPACEDIM; i++)
			for (int j=0; j < AMREX_SPACEDIM; j++)
				{
					f(i) += mu*(gradeps[i](j,j) + gradeps[j](i,j))  + lambda*gradeps[j](j,i);
				}
		return f;
	}

	void Randomize()
	{
		mu = Util::Random();
		lambda = Util::Random();
	};

	virtual void Print (std::ostream& os) const
	{
		os <<    "mu = " << mu <<
			" lambda = " << lambda << std::endl;
	}

	Multiwell operator + (const Multiwell &rhs) const
	{
		Multiwell ret(mu+rhs.mu, lambda+rhs.lambda,Fgb + rhs.Fgb);
		for (int i = 0; i < AMREX_SPACEDIM; i++) ret.gradFgb[i] = gradFgb[i] + rhs.gradFgb[i];
		return ret;
	}

	Multiwell & operator = (const Multiwell &rhs)
	{
		mu=rhs.mu; lambda=rhs.lambda; Fgb = rhs.Fgb;
		for (int i = 0; i < AMREX_SPACEDIM; i++) gradFgb[i] = rhs.gradFgb[i];
		return *this;
	}

	void operator += (const Multiwell &rhs)
	{
		mu+=rhs.mu; lambda+=rhs.lambda; Fgb += rhs.Fgb;
		for (int i = 0; i < AMREX_SPACEDIM; i++) gradFgb[i] += rhs.gradFgb[i];
	}

	Multiwell operator - (const Multiwell &rhs) const
	{
		Multiwell ret(mu-rhs.mu, lambda-rhs.lambda,Fgb-rhs.Fgb);
		for (int i = 0; i < AMREX_SPACEDIM; i++) ret.gradFgb[i] = gradFgb[i] - rhs.gradFgb[i];
		return ret;
	}

	//Multiwell operator * (const Multiwell &rhs) const
	//{return Multiwell(mu*rhs.mu, lambda*rhs.lambda,Fgb);}

	//Multiwell operator / (const Multiwell &rhs) const
	//{return Multiwell(mu/rhs.mu, lambda/rhs.lambda);}

	virtual Multiwell operator / (const Set::Scalar alpha) const
	{
		Multiwell ret(mu/alpha, lambda/alpha,Fgb/alpha);
		for (int i = 0; i < AMREX_SPACEDIM; i++) ret.gradFgb[i] = gradFgb[i]/alpha;
		return ret;
	}

	virtual Multiwell operator * (const Set::Scalar alpha) const
	{
		Multiwell ret(mu*alpha, lambda*alpha,Fgb*alpha);
		for (int i = 0; i < AMREX_SPACEDIM; i++) ret.gradFgb[i] = gradFgb[i]*alpha;
		return ret;
	}
	
private:
	Set::Scalar mu; ///< First Lame parameter \f$\mu\f$
	Set::Scalar lambda; ///< Second Lame parameter \f$\lambda\f$
public:
	Set::Matrix Fgb = Set::Matrix::Zero();
	std::array<Set::Matrix,AMREX_SPACEDIM> gradFgb;
};

}
}
}

#endif


