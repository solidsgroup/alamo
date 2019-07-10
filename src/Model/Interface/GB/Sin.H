#ifndef MODEL_INTERFACE_GB_SIN_H
#define MODEL_INTERFACE_GB_SIN_H

#include <iostream>
#include <fstream>

#include "AMReX.H"
#include "GB.H"
#include "Set/Set.H"
#include "Util/Util.H"

namespace Model
{
namespace Interface
{
namespace GB
{
class Sin : public GB
{
public:
	Sin()
	{
	};
	Sin(amrex::Real _theta0, amrex::Real _sigma0, amrex::Real _sigma1)
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
	amrex::Real W(amrex::Real theta)
	{
		//theta = theta*PI/180;
		//theta0 = theta0*PI/180;
		//return sin(2.0*(theta-theta0))*sin(2.0*(theta-theta0))*sigma1+sigma0;
		return sigma0 + 0.5*sigma1*(1.0 - cos(4.0*(theta-theta0)));
	};
	amrex::Real DW(amrex::Real theta)
	{
		//theta = theta*PI/180;
		//theta0 = theta0*PI/180;
		//return 4.0*sigma1*sin(2.0*(theta-theta0))*cos(2.0*(theta-theta0));
		return 2.0*sigma1*sin(4.0*(theta-theta0));
	};
	amrex::Real DDW(amrex::Real theta)
	{
		//theta = theta*PI/180;
		//    theta0 = theta0*PI/180;
		return 8.0*sigma1*cos(4.0*(theta-theta0));
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
