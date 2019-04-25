#ifndef MODEL_SOLID_LINEARELASTIC_SIMPLE_H_
#define MODEL_SOLID_LINEARELASTIC_SIMPLE_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>

#include "Util/Util.H"
#include "Set/Set.H"
#include "Model/Solid/LinearElastic/Simple.H"

namespace Model
{
namespace Solid
{
namespace LinearElastic
{
///
/// Isotropic implements the following elastic model
/// \f[ \mathbb{C}_{ijkl} = \alpha\delta_{ik}\delta_{jl}  \f]
/// corresponding to the elastic model
/// \f[ \sigma_{ij} = \alpha\delta_{ik}\delta_{jl}u_{k,l} = \alpha u_{i,j} \f]
/// \f[ f_{i} = \mathbb{C}_{ijkl}u_{k,jl} = \alpha\delta_{ik}\delta_{jl}u_{k,jl} = \alpha u_{i,jj} \f]
///
class Simple : public Model::Solid::LinearElastic::LinearElastic
{
public:
	Simple() {} ;
	Simple(Set::Scalar _alpha) : alpha(_alpha)  {} ;
	virtual Set::Matrix operator () (Set::Matrix &gradu) const
	{
		return alpha * 0.5*(gradu + gradu.transpose());
	};
	virtual Set::Vector operator () (std::array<Set::Matrix,AMREX_SPACEDIM> &gradgradu)
	{
		Set::Vector f = Set::Vector::Zero();
		for (int i = 0; i < AMREX_SPACEDIM; i++)
			for (int j=0; j < AMREX_SPACEDIM; j++)
				f(i) += alpha * 0.5*(gradgradu[i](j,j) + gradgradu[j](i,j));
		return f;
	}

	void Randomize()
	{
		alpha = Util::Random();
	};

	Simple operator + (const Simple &rhs) const
	{return Simple(alpha+rhs.alpha);}

	void operator += (const Simple &rhs)
	{alpha+=rhs.alpha;}

	Simple operator - (const Simple &rhs) const
	{return Simple(alpha-rhs.alpha);}

	Simple operator * (const Simple &rhs) const
	{return Simple(alpha*rhs.alpha);}

	Simple operator / (const Simple &rhs) const
	{return Simple(alpha/rhs.alpha);}

	virtual Simple operator / (const Set::Scalar beta) const
	{return Simple(alpha/beta);}

	virtual Simple operator * (const Set::Scalar beta) const
	{return Simple(alpha*beta);}

private:
	Set::Scalar alpha; ///< Coefficient
};

}
}
}

#endif


