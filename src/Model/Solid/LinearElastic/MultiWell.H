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
	Multiwell(Set::Scalar _mu, Set::Scalar _lambda) : mu(_mu) , lambda(_lambda) {} ;
	~Multiwell() {} ;

#ifndef DOXYGEN_SHOULD_SKIP_THIS
	AMREX_FORCE_INLINE
#endif
	virtual Set::Matrix operator () (Set::Matrix &gradu) const
	{

		return mu*(gradu + gradu.transpose()) + lambda*Set::Matrix::Identity()*gradu.trace();
	};

#ifndef DOXYGEN_SHOULD_SKIP_THIS
	AMREX_FORCE_INLINE
#endif
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
		mu = Util::Random();
		lambda = Util::Random();
	};

	virtual void Print (std::ostream& os) const
	{
		os <<    "mu = " << mu <<
			" lambda = " << lambda << std::endl;
	}

	Multiwell operator + (const Multiwell &rhs) const
	{return Multiwell(mu+rhs.mu, lambda+rhs.lambda);}

	void operator += (const Multiwell &rhs)
	{mu+=rhs.mu; lambda+=rhs.lambda;}

	Multiwell operator - (const Multiwell &rhs) const
	{return Multiwell(mu-rhs.mu, lambda-rhs.lambda);}

	Multiwell operator * (const Multiwell &rhs) const
	{return Multiwell(mu*rhs.mu, lambda*rhs.lambda);}

	Multiwell operator / (const Multiwell &rhs) const
	{return Multiwell(mu/rhs.mu, lambda/rhs.lambda);}

	virtual Multiwell operator / (const Set::Scalar alpha) const
	{return Multiwell(mu/alpha, lambda/alpha);}

	virtual Multiwell operator * (const Set::Scalar alpha) const
	{return Multiwell(mu*alpha, lambda*alpha);}

private:
	Set::Scalar mu; ///< First Lame parameter \f$\mu\f$
	Set::Scalar lambda; ///< Second Lame parameter \f$\lambda\f$
	Set::Matrix Fgb = Set::Matrix::Identity();
};

}
}
}

#endif


