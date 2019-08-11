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
    Constant(Set::Scalar _Gc, Set::Scalar _zeta)
    {
        Gc = _Gc;
        zeta = _zeta;
    };
    Set::Scalar w_phi(Set::Scalar c)
    {
        return 1. - 4.*std::pow(c,3) + 3.*std::pow(c,4);
    };
    
    Set::Scalar Dw_phi(Set::Scalar c)
    {
        return -12.*std::pow(c,2) + 12.*std::pow(c,3);
    };

    //Set::Scalar g_phi(Set::Scalar c) { return 1. - w_phi(c); }
    //Set::Scalar Dg_phi(Set::Scalar c) { return -Dw_phi(c);}
    
    Set::Scalar Epc(Set::Scalar /*theta*/)
    {
        return Gc/(4.*C_phi*zeta);
    };
    Set::Scalar kappa(Set::Scalar /*theta*/)
    {
        return Gc*zeta/(2*C_phi);
    };

private:
    amrex::Real Gc = NAN, zeta = NAN, C_phi = 0.71658;
};
}
}
}
#endif