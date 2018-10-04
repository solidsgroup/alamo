#ifndef MODEL_SOLID_ELASTIC_ISOTROPIC_H_
#define MODEL_SOLID_ELASTIC_ISOTROPIC_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>

#include "Util/Util.H"
#include "Set/Set.H"
#include "Model/Solid/Elastic/Elastic.H"

namespace Model
{
namespace Solid
{
namespace Elastic
{
namespace Isotropic
{
class Isotropic : public Model::Solid::Elastic::Elastic
{
public:
	Isotropic() {} ;
	Isotropic(Set::Scalar _mu, Set::Scalar _lambda) : mu(_mu) , lambda(_lambda) {} ;
	virtual Set::Matrix operator () (Set::Matrix &eps)
	{
		return 2*mu*eps + lambda*Set::Matrix::Identity()*eps.trace();
	};
	virtual Set::Vector operator () (std::array<Set::Matrix,AMREX_SPACEDIM> &gradeps)
	{
		Set::Vector f = Set::Vector::Zero();
		for (int i = 0; i < AMREX_SPACEDIM; i++)
			for (int j=0; j < AMREX_SPACEDIM; j++)
				f(i) += 2*mu*gradeps[i](j,j)  + lambda*gradeps[j](j,i);
		return f;
	}

	Isotropic operator + (const Isotropic &rhs)
	{return Isotropic(mu+rhs.mu, lambda+rhs.lambda);}

	Isotropic operator - (const Isotropic &rhs)
	{return Isotropic(mu-rhs.mu, lambda-rhs.lambda);}

	virtual Isotropic operator / (const Set::Scalar alpha)
	{return Isotropic(mu/alpha, lambda/alpha);}

private:
	Set::Scalar mu, lambda;
};

}
}
}
}

#endif


