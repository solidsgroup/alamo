#ifndef MODEL_INTERFACE_GB_ABSSIN_H
#define MODEL_INTERFACE_GB_ABSSIN_H

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
class AbsSin : public GB
{
public:
AbsSin()
{};
AbsSin(amrex::Real _theta0, amrex::Real _sigma0, amrex::Real _sigma1)
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
	//sigma(theta)=sigma0+sigma1*fabs(sin(n*(theta-theta0)))
	//n=2:  
	return sigma0+sigma1*fabs(sin(2*(theta-theta0)));  
  
};
amrex::Real DW(amrex::Real theta)
{
	//sigma'(theta)=(n*sigma1*sin(n*(theta-theta0))*cos(n*(theta-theta0)))/fabs(sin(n*(theta-theta0)))
	//n=2:
	return (2*sigma1*sin(2*(theta-theta0))*cos(2*(theta-theta0)))/fabs(sin(2*(theta-theta0)));
};
amrex::Real DDW(amrex::Real theta)
{
    //sigma''(theta)=-n^2*sigma1*fabs(sin(n*(theta-theta0)))
    //n=2:
    return -4*sigma1*fabs(sin(2*(theta-theta0)));
};
 
private:
	amrex::Real theta0 = NAN, sigma0 = NAN, sigma1 = NAN;

public:	
	static void Parse(AbsSin & value, amrex::ParmParse & pp)
	{
		pp.query("theta0",value.theta0);
		value.theta0 *= 0.01745329251; // convert degrees into radians
		pp.query("sigma0",value.sigma0);
		pp.query("sigma1",value.sigma1);
	}

};
}
}
}

#endif
