#ifndef MODEL_INTERFACE_GB_SH_H
#define MODEL_INTERFACE_GB_SH_H

#include <iostream>
#include <fstream>

#include "AMReX.H"
#include "GB.H"
#include "Set/Set.H"
#include "Util/Util.H"

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
///               \cos(2 \arcsin x) &= 1 - 2 x^2 \\ %
///               \sin(2 \arccos x) &= 2 x \sqrt{1-x^2} &
///               \cos(2 \arccos x) &= 2 x^2 - 1 \\ %
///               \sin(2 \arctan x) &= \frac{2 x}{1 + x^2}  &
///               \cos(2 \arctan x) &= \frac{1-x^2}{1+x^2}    \f}
/// Specializing for tangent:
///     \f{align}{\sin(2 \arctan y/x) &= \frac{2 x y}{x^2 + y^2}     &
///               \cos(2 \arctan y/x) &= \frac{x^2-y^2}{x^2+y^2}  \f}
///
class SH : public GB
{
public:
    SH(){};
    SH(IO::ParmParse &pp){pp_queryclass(*this);};
    SH(IO::ParmParse &pp,std::string name){pp_queryclass(name,*this);};
    SH(const amrex::Real a_theta0, const amrex::Real a_phi0, const amrex::Real a_sigma0, const amrex::Real a_sigma1)
    {
        Define(a_theta0,a_phi0,a_sigma0,a_sigma1);
    };
    void Define(const amrex::Real a_theta0, const amrex::Real a_phi0, const amrex::Real a_sigma0, const amrex::Real a_sigma1)
    {
        theta0 = a_theta0;
        phi0   = a_phi0;
        sigma0 = a_sigma0;
        sigma1 = a_sigma1;
    };
    void Randomize()
    {
        theta0 = Util::Random()*Set::Constant::Pi;
        sigma0 = Util::Random();
        sigma1 = Util::Random();
    };
    Set::Scalar W(const Set::Scalar a_theta, const Set::Scalar a_phi) const
    {
        return sigma0 + sigma1*(1.0 - (cos(2*(a_phi))*cos(2*(a_phi)) * sin(2*(a_theta))*sin(2*(a_theta))));
    };
    std::array<Set::Scalar,2> DW(const Set::Scalar a_theta, const Set::Scalar a_phi) const
    {
        return {- sigma1 * 4.0 * cos(2*a_phi) * cos(2*a_phi) * sin(2*a_theta) * cos(2*a_theta),
                + sigma1 * 4.0 * sin(2*a_phi) * cos(2*a_phi) * sin(2*a_theta) * sin(2*a_theta)} ;
    };

    virtual Set::Scalar W(const Set::Vector &a_n) const override
    {
        Set::Vector n = a_n / a_n.lpNorm<2>();
        Set::Scalar theta = std::acos(n(0));
        Set::Scalar phi   = std::atan2(n(2),n(1));
        return W(theta,phi);
    };
    virtual Set::Scalar DW(const Set::Vector &a_n, const Set::Vector &t_n) const override
    {
        const Set::Scalar alpha = 1E-4;
        return (W(a_n + alpha*t_n) - W(a_n - alpha*t_n)) / 2.0 / alpha;
    };
    virtual Set::Scalar DDW(const Set::Vector &a_n, const Set::Vector &t_n) const override
    {
        const Set::Scalar alpha = 1E-4;
        return (W(a_n + alpha*t_n) - 2.0*W(a_n) + W(a_n - alpha*t_n)) / alpha / alpha;
    };


    Set::Scalar W(const Set::Scalar ) const override   {return NAN;};
    Set::Scalar DW(const Set::Scalar ) const override  {return NAN;};
    Set::Scalar DDW(const Set::Scalar ) const override {return NAN;};
  
private:
    Set::Scalar theta0 = NAN, phi0 = NAN, sigma0 = NAN, sigma1 = NAN;
    
public:
    static void Parse(SH & value, amrex::ParmParse & pp)
    {
        pp_query("theta0",value.theta0); // Theta offset (degrees)
        value.theta0 *= 0.01745329251; // convert degrees into radians
        pp_query("phi0",value.phi0); // Phi offset (radians)
        value.phi0 *= 0.01745329251;   // convert degrees into radians
        pp_query("sigma0",value.sigma0); // Minimum energy value
        pp_query("sigma1",value.sigma1); // Energy multiplier

        std::string reg_str = "wilhelm";
        pp_query("regularization",reg_str); // Type of regularization to use: {wilhelm,k23}
        if (reg_str == "wilhelm") value.regularization = Regularization::Wilhelm;
        else if (reg_str == "k23") value.regularization = Regularization::K23;
        else Util::Abort(INFO,"Invalid regularization string ",reg_str);
    }

};
}
}
}
#endif
