#ifndef MODEL_INTERFACE_GB_SH_H
#define MODEL_INTERFACE_GB_SH_H

#include <iostream>
#include <fstream>

#include "AMReX.H"
#include "GB.H"
#include "Set/Set.H"
#include "Util/Util.H"

#include <boost/math/special_functions/spherical_harmonic.hpp>

#define PI 3.14159265 

namespace Model
{
namespace Interface
{
namespace GB
{
///
/// A 2D interface model class. Under construction
///
/// Here are some useful trig identities for reference:
///
///     \f{align}{\sin(2 \arcsin x) &= 2 x \sqrt{1-x^2} &
///               \cos(2 \arcsin x) &= 1 - 2 x^2 \\
///               \sin(2 \arccos x) &= 2 x \sqrt{1-x^2} &
///               \cos(2 \arccos x) &= 2 x^2 - 1 \\
///               \sin(2 \arctan x) &= \frac{2 x}{1 + x^2}  &
///               \cos(2 \arctan x) &= \frac{1-x^2}{1+x^2}	\f}
/// Specializing for tangent:
///     \f{align}{\sin(2 \arctan y/x) &= \frac{2 x y}{x^2 + y^2}	 &
///               \cos(2 \arctan y/x) &= \frac{x^2-y^2}{x^2+y^2}  \f}
///
class SH : public GB
{
public:
	SH()
	{
	};
	SH(amrex::Real _theta0, amrex::Real _sigma0, amrex::Real _sigma1)
	{
		Define(_theta0,_sigma0,_sigma1);
	};
	void Define(amrex::Real _theta0, amrex::Real _sigma0, amrex::Real _sigma1)
	{
		theta0 = _theta0;
		sigma0 = _sigma0;
		sigma1 = _sigma1;
	};
	void Randomize()
	{
		theta0 = Util::Random()*Set::Constant::Pi;
		sigma0 = Util::Random();
		sigma1 = Util::Random();
	};
	amrex::Real W(Set::Scalar a_theta, Set::Scalar a_phi) const
	{
		//int n=3, m=2;
		//return boost::math::spherical_harmonic_r(n,m,a_theta,a_phi);
		return sigma0 - sigma1*cos(2*a_phi)*cos(2*a_phi) * sin(2*a_theta)*sin(2*a_theta);
	};
	std::array<Set::Scalar,2> DW(Set::Scalar a_theta, Set::Scalar a_phi) const
	{
		return {- sigma1 * 4.0 * cos(2*a_phi) * cos(2*a_phi) * sin(2*a_theta) * cos(2*a_theta),
				+ sigma1 * 4.0 * sin(2*a_phi) * cos(2*a_phi) * sin(2*a_theta) * sin(2*a_theta)} ;
	};

	Set::Scalar W(Set::Vector a_n) const
	{
		Set::Vector n = a_n / a_n.lpNorm<2>();
		Set::Scalar theta = std::acos(n(2));
		Set::Scalar phi   = std::atan2(n(1),n(0));
		return W(theta,phi);


		return sigma0 - sigma1*(
			((n(1)*n(1) - n(0)*n(0))*(n(1)*n(1) - n(0)*n(0)) / (n(1)*n(1) + n(0)*n(0)) / (n(1)*n(1) + n(0)*n(0)))
			*
			4.*n(2)*n(2)*(1. - n(2)*n(2)));
	};
	Set::Vector DW(Set::Vector a_n) const
	{
		Set::Vector n = a_n / a_n.lpNorm<2>();
		Set::Vector dw;
		Set::Scalar theta = std::acos(n(2));
		Set::Scalar phi   = std::atan2(n(1),n(0));
		std::array<Set::Scalar,2> dw1 = DW(theta,phi);

		dw(0) = -dw1[1] / (n(0)*n(0) + n(1)*n(1));
		dw(1) = +dw1[1] * n(0) * n(0) / (n(0)*n(0) + n(1)*n(1));

		return dw;
	};


	amrex::Real W(amrex::Real theta)
	{
		Util::Message(INFO);
		
		return NAN;
	};
	amrex::Real DW(amrex::Real theta)
	{
		return NAN;
	};
	amrex::Real DDW(amrex::Real theta)
	{
		return NAN;
	};
  
	//virtual amrex::Real operator()(amrex::Real theta)
	//{
	//amrex::Real d_gbenergy;
	//d_gbenergy= 
	//return d_gbenergy;    
	//};
  
private:
	amrex::Real theta0 = NAN, sigma0 = NAN, sigma1 = NAN;
};
}
}
}
#endif
