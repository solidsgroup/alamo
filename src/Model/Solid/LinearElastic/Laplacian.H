#ifndef MODEL_SOLID_LINEARELASTIC_LAPLACIAN_H_
#define MODEL_SOLID_LINEARELASTIC_LAPLACIAN_H_

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
/// \f[ \mathbb{C}_{ijkl} = \alpha\delta_{ik}\delta_{jl}  \f]
/// corresponding to the elastic model
/// \f[ \sigma_{ij} = \alpha\delta_{ik}\delta_{jl}u_{k,l} = \alpha u_{i,j} \f]
/// \f[ f_{i} = \mathbb{C}_{ijkl}u_{k,jl} = \alpha\delta_{ik}\delta_{jl}u_{k,jl} = \alpha u_{i,jj} \f]
///
class Laplacian : public Model::Solid::LinearElastic::LinearElastic<Set::Sym::Isotropic>
{
public:
	Laplacian() {} ;
	Laplacian(Set::Scalar _alpha) : alpha(_alpha)  {} ;
	~Laplacian() {} ;
	AMREX_FORCE_INLINE
	virtual Set::Matrix operator () (Set::Matrix &gradu,bool =true) const
	{
		//return alpha * 0.5*(gradu + gradu.transpose());
		return alpha * gradu;
	};
	AMREX_FORCE_INLINE
	virtual Set::Vector operator () (std::array<Set::Matrix,AMREX_SPACEDIM> &gradgradu,bool =true)
	{
		Set::Vector f = Set::Vector::Zero();
		for (int i = 0; i < AMREX_SPACEDIM; i++)
			for (int j=0; j < AMREX_SPACEDIM; j++)
				f(i) += alpha *gradgradu[i](j,j);
		return f;
	}

	void Randomize()
	{
		alpha = Util::Random();
	};

	Laplacian operator + (const Laplacian &rhs) const
	{return Laplacian(alpha+rhs.alpha);}

	void operator += (const Laplacian &rhs)
	{alpha+=rhs.alpha;}

	Laplacian operator - (const Laplacian &rhs) const
	{return Laplacian(alpha-rhs.alpha);}

	Laplacian operator * (const Laplacian &rhs) const
	{return Laplacian(alpha*rhs.alpha);}

	Laplacian operator / (const Laplacian &rhs) const
	{return Laplacian(alpha/rhs.alpha);}

	virtual Laplacian operator / (const Set::Scalar beta) const
	{return Laplacian(alpha/beta);}

	virtual Laplacian operator * (const Set::Scalar beta) const
	{return Laplacian(alpha*beta);}

public:
	Set::Scalar alpha = NAN; ///< Coefficient
};

}
}
}

#endif


