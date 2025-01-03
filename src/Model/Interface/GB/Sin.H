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
    static constexpr const char* name = "sin";

    Sin() {};
    Sin(IO::ParmParse &pp) {pp_queryclass(*this);};
    Sin(IO::ParmParse &pp,std::string name) {pp_queryclass(name,*this);};
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
    Set::Scalar W(const Set::Scalar theta) const
    {
        return sigma0 + 0.5*sigma1*(1.0 - cos(n*(theta-theta0)));
    };
    Set::Scalar DW(const Set::Scalar theta) const
    {
        return 2.0*sigma1*sin(n*(theta-theta0));
    };
    Set::Scalar DDW(const Set::Scalar theta) const
    {
        return 8.0*sigma1*cos(n*(theta-theta0));
    };
  
private:
    Set::Scalar theta0 = NAN, sigma0 = NAN, sigma1 = NAN, n = 4.0;

public:    
    static void Parse(Sin & value, amrex::ParmParse & pp)
    {
        pp_query("theta0",value.theta0); // Theta offset (degrees)
        value.theta0 *= 0.01745329251; // convert degrees into radians
        pp_query("sigma0",value.sigma0); // Minimum energy
        pp_query("sigma1",value.sigma1); // Energy multiplier
        pp_query("n",value.n); // Frequency number (integer)
    }
    
};
}
}
}
#endif
