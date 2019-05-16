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
		int n=3, m=2;
		//return boost::math::spherical_harmonic_r(n,m,a_theta,a_phi);
		return cos(2*a_phi)*cos(2*a_phi) * sin(2*a_theta)*sin(2*a_theta);
	};
	amrex::Real W(amrex::Real theta)
	{
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
