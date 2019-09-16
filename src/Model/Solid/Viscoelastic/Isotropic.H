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
	~Isotropic() {} ;
	Isotropic(Set::Scalar _nu, int _number_of_terms, std::array<Set::Scalar,9> _modulus, std::array<Set::Scalar,8> _relaxation) : 
		nu(_nu) , number_of_terms(_number_of_terms) 
	{
		if(number_of_terms < 1)
			Util::Abort(INFO, "invalid value of number of terms");

		if(number_of_terms > 8)
			Util::Abort(INFO,"number of terms exceeded 8");

		for (int i = 0; i < number_of_terms; i++)
		{
			youngs_modulus[i] = _modulus[i];
			relaxation_times[i] = _relaxation[i];
		}
		youngs_modulus[number_of_terms] = _modulus[number_of_terms];

		Set::Scalar e = youngs_modulus[number_of_terms];
		for(int i = 0; i < number_of_terms; i++)
			e += youngs_modulus[i];
		
		lambda = e*nu/((1.0+nu)*(1-2*nu));
		mu = e/(2.0*(1.0+nu));
	}

	Isotropic(Set::Scalar _lambda, Set::Scalar _mu): lambda(_lambda), mu(_mu){};

	void SetTime(Set::Scalar time)
	{
		Set::Scalar e = youngs_modulus[number_of_terms];
		for(int i = 0; i < number_of_terms; i++)
			e += youngs_modulus[i]*std::exp(-time/relaxation_times[i]);
		
		lambda = e*nu/((1.0+nu)*(1-2*nu));
		mu = e/(2.0*(1.0+nu));
	}
	
	virtual Set::Matrix operator () (Set::Matrix &eps,bool =true) const
	{
		return mu*(eps + eps.transpose()) + lambda*Set::Matrix::Identity()*eps.trace();
	};
	virtual Set::Vector operator () (std::array<Set::Matrix,AMREX_SPACEDIM> &gradeps,bool =true)
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
	Set::Scalar nu; ///< Poisson Ratio
	//Set::Scalar lambda; ///< Second Lame parameter \f$\lambda\f$
	std::array<Set::Scalar,9> youngs_modulus; ///< Prony Series elastic modulus
	std::array<Set::Scalar,8> relaxation_times; ///< Prony Series relaxation times
	Set::Scalar lambda;
	Set::Scalar mu;
	int number_of_terms = 8;
};

}
}
}

#endif


