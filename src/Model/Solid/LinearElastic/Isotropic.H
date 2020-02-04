#ifndef MODEL_SOLID_LINEARELASTIC_ISOTROPIC_H_
#define MODEL_SOLID_LINEARELASTIC_ISOTROPIC_H_

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
///
/// Isotropic implements the following elastic model
/// \f[ \sigma_{ij} = 2\mu\varepsilon_{ij} + \lambda\varepsilon_{kk}\delta_{ij}  \f]
///
class Isotropic : public Model::Solid::LinearElastic::LinearElastic<Set::Sym::Isotropic>
{
public:
	Isotropic() {} ;
	Isotropic(Set::Scalar _mu, Set::Scalar _lambda) : mu(_mu) , lambda(_lambda) {} ;
	~Isotropic() {} ;

	void Define(Set::Scalar a_mu, Set::Scalar a_lambda)
	{
		mu = a_mu;
		lambda = a_lambda;
	}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
	AMREX_FORCE_INLINE
#endif
	virtual Set::Matrix operator () (Set::Matrix &eps,bool =true) const
	{
		return mu*(eps + eps.transpose()) + lambda*Set::Matrix::Identity()*eps.trace();
	};

#ifndef DOXYGEN_SHOULD_SKIP_THIS
	AMREX_FORCE_INLINE
#endif
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
	};

	virtual void Print (std::ostream& os) const
	{
		os <<    "mu = " << mu <<
			" lambda = " << lambda << std::endl;
	}

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

private:
	Set::Scalar mu; ///< First Lame parameter \f$\mu\f$
	Set::Scalar lambda; ///< Second Lame parameter \f$\lambda\f$
};

}
}
}

#endif


