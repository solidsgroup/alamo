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
	AMREX_FORCE_INLINE
	amrex::Real W(Set::Scalar a_theta, Set::Scalar a_phi) const
	{
		Util::Message(INFO);
		//int n=3, m=2;
		//return boost::math::spherical_harmonic_r(n,m,a_theta,a_phi);
		return sigma0 - sigma1*cos(2*a_phi)*cos(2*a_phi) * sin(2*a_theta)*sin(2*a_theta);
	};
	AMREX_FORCE_INLINE
	Set::Scalar W(Set::Vector n) const
	{
		// Handy trig identities:
		// sin(2 arcsin(x)) = 2 x sqrt(1-x^2) ## confirmed
		// sin(2 arccos(x)) = 2 x sqrt(1-x^2) ## confirmed
		// sin(2 arctan(x)) = 2 x / (1 + x^2) ## confirmed
		// cos(2 arcsin(x)) = 1 - 2 x^2       ## confirmed
		// cos(2 arccos(x)) = 2 x^2 - 1       ## confirmed
		// cos(2 arctan(x)) = (1-x^2)/(1+x^2) ## confirmed
		
		// tangent extensions:
		// sin(2 arctan(y/x)) = 2 x y / (x^2 + y^2)
		// cos(2 arctan(y/x)) = (x^2-y^2)/(x^2+y^2)

		// Set::Scalar theta = std::acos(normal(2));
		// Set::Scalar phi   = std::atan2(normal(1),normal(0));

		// sin(2 arccos(x)) = 2 x sqrt(1-x^2) ## confirmed
		return sigma0 - sigma1*(
			// (cos(2phi))^2 = cos(2 atan(N1/N0)) = (N1^2-N0^2)^2/(N1^2+N0^2)^2
			((n(1)*n(1) - n(0)*n(0))*(n(1)*n(1) - n(0)*n(0)) / (n(1)*n(1) + n(0)*n(0)) / (n(1)*n(1) + n(0)*n(0)))
			*
			// sin(2theta)^2 = sin(2 acos(N2))^2 = (2*N2*sqrt(1-N2^2))^2 = 4 N2^2 (1-N2^2)
			4.*n(2)*n(2)*(1. - n(2)*n(2)));
	};
	AMREX_FORCE_INLINE
	Set::Vector DW(Set::Vector n) const
	{
		return Set::Vector::Zero();
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
