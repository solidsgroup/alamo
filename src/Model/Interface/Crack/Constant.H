#ifndef MODEL_INTERFACE_CRACK_CONSTANT_H
#define MODEL_INTERFACE_CRACK_CONSTANT_H

#include <iostream>
#include <fstream>

#include "AMReX.H"
#include "Crack.H"
#include "Set/Set.H"
#include "Util/Util.H"

namespace Model
{
namespace Interface
{
namespace Crack
{
class Constant : public Crack
{
public:
    Constant()
    {       
    };
    Constant(Set::Scalar _Gc, Set::Scalar _zeta, Set::Scalar _mob)
    {
        Gc = _Gc;
        zeta = _zeta;
        mobility = _mob;
    };
    Set::Scalar w_phi(Set::Scalar c)
    {
        //return 1. - g_phi(c);
        return (1. - c)*(1. - c);
    };
    
    Set::Scalar Dw_phi(Set::Scalar c)
    {
        //return -1.*Dg_phi(c);
        return -2.0*(1.-c);
    };

    Set::Scalar g_phi(Set::Scalar c)
    {
        //return 4.*c*c*c - 3.*c*c*c*c;
        return c*c*c*c;
    };

    Set::Scalar Dg_phi(Set::Scalar c)
    {
        //return 12.*c*c*(1.-c);
        //return 2.*c;
        return 4.*c*c*c;
    };

    //Set::Scalar g_phi(Set::Scalar c) { return 1. - w_phi(c); }
    //Set::Scalar Dg_phi(Set::Scalar c) { return -Dw_phi(c);}
    
    Set::Scalar Epc(Set::Scalar /*theta*/)
    {
        return Gc/(4.0*zeta);
    };
    Set::Scalar kappa(Set::Scalar /*theta*/)
    {
        return 2.0*Gc*zeta;
    };

private:
    amrex::Real Gc = NAN, zeta = NAN, mobility = NAN;
};
}
}
}
#endif