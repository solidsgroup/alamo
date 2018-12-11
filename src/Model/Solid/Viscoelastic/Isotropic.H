#ifndef MODEL_SOLID_VISCOELASTIC_ISOTROPIC_H_
#define MODEL_SOLID_VISCOELASTIC_ISOTROPIC_H_

#include <AMReX.H>
#include <AMReX_REAL.H>
#include <eigen3/Eigen/Core>

#include "Util/Util.H"
#include "Set/Set.H"
#include "Model/Solid/Viscoelastic/Viscoelastic.H"

namespace Model
{
namespace Solid
{
namespace Viscoelastic
{
///
/// Isotropic implements the following elastic model
/// \f[ \sigma_{ij} = 2\mu\varepsilon_{ij} + \lambda\varepsilon_{kk}\delta_{ij}  \f]
///
class Isotropic : public Model::Solid::Viscoelastic::Viscoelastic
{
public:
	Isotropic() {} ;
	Isotropic(Set::Scalar _nu, amrex::Vector<Set::Scalar> _modulus, amrex::Vector<Set::Scalar> _relaxation) : 
		nu(_nu) , youngs_modulus(_modulus), relaxation_times(_relaxation) {} ;
	virtual Set::Matrix operator () (Set::Matrix &eps, Set::Scalar time) const
	{
		Set::Scalar e = youngs_modulus[youngs_modulus.size()-1];
		for(int i = 0; i < relaxation_times.size(); i++)
			e += youngs_modulus[i]*std::exp(-time/relaxation_times[i]);
		
		Set::Scalar lambda = e*nu/((1.0+nu)*(1-2*nu));
		Set::Scalar mu = e/(2.0*(1.0+nu));

		return mu*(eps + eps.transpose()) + lambda*Set::Matrix::Identity()*eps.trace();
	};
	virtual Set::Vector operator () (std::array<Set::Matrix,AMREX_SPACEDIM> &gradeps, Set::Scalar time)
	{
		Set::Scalar e = youngs_modulus[youngs_modulus.size()-1];
		for(int i = 0; i < relaxation_times.size(); i++)
			e += youngs_modulus[i]*std::exp(-time/relaxation_times[i]);
		
		Set::Scalar lambda = e*nu/((1.0+nu)*(1-2*nu));
		Set::Scalar mu = e/(2.0*(1.0+nu));

		Set::Vector f = Set::Vector::Zero();
		for (int i = 0; i < AMREX_SPACEDIM; i++)
			for (int j=0; j < AMREX_SPACEDIM; j++)
				f(i) += mu*(gradeps[i](j,j) + gradeps[j](i,j))  + lambda*gradeps[j](j,i);
		return f;
	}

	Isotropic operator + (const Isotropic &rhs) const
	{return Isotropic(nu+rhs.nu, youngs_modulus+rhs.youngs_modulus, relaxation_times+rhs.relaxation_times);}

	void operator += (const Isotropic &rhs)
	{nu+=rhs.nu; youngs_modulus+=rhs.youngs_modulus; relaxation_times+=rhs.relaxation_times;}

	Isotropic operator - (const Isotropic &rhs) const
	{return Isotropic(nu-rhs.nu, youngs_modulus-rhs.youngs_modulus, relaxation_times-rhs.relaxation_times);}

	Isotropic operator * (const Isotropic &rhs) const
	{return Isotropic(nu*rhs.nu, youngs_modulus*rhs.youngs_modulus, relaxation_times*rhs.relaxation_times);}

	Isotropic operator / (const Isotropic &rhs) const
	{return Isotropic(nu/rhs.nu, youngs_modulus/rhs.youngs_modulus, relaxation_times/rhs.relaxation_times);}

	virtual Isotropic operator / (const Set::Scalar alpha) const
	{return Isotropic(nu/alpha, youngs_modulus/alpha);}

	virtual Isotropic operator * (const Set::Scalar alpha) const
	{return Isotropic(u*alpha, youngs_modulus*alpha);}

private:
	Set::Scalar nu; ///< Poisson Ratio
	//Set::Scalar lambda; ///< Second Lame parameter \f$\lambda\f$
	amrex::Vector<Set::Scala> youngs_modulus; ///< Prony Series elastic modulus
	amrex::Vector<Set::Scalar> relaxation_times; ///< Prony Series relaxation times
};

}
}
}

#endif


