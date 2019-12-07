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
        mult_Gc = 1.0;
        mult_Lap = 1.0;
    };
    Constant(Set::Scalar _Gc, Set::Scalar _zeta, Set::Scalar _mob,PhiType a_type)
    {
        Gc = _Gc;
        zeta = _zeta;
        mobility = _mob;
        mult_Gc = 1.0;
        mult_Lap = 1.0;
        
        SetPhiType(a_type);
    };
    Constant(Set::Scalar _Gc, Set::Scalar _zeta, Set::Scalar _mob, Set::Scalar _mult_Gc, Set::Scalar _mult_Lap)
    {
        Gc = _Gc;
        zeta = _zeta;
        mobility = _mob;
        mult_Gc = _mult_Gc;
        mult_Lap = _mult_Lap;
    };

    Constant(Set::Scalar _Gc, Set::Scalar _zeta, Set::Scalar _mob, Set::Scalar _mult_Gc, Set::Scalar _mult_Lap, PhiType a_type)
    {
        Gc = _Gc;
        zeta = _zeta;
        mobility = _mob;
        mult_Gc = _mult_Gc;
        mult_Lap = _mult_Lap;

        SetPhiType(a_type);
    };
    /*Set::Scalar w_phi(Set::Scalar c)
    {
        return 1. - g_phi(c);
	    //return (1.-c)*(1.-c)*c*c;
        //return (1. - c)*(1. - c)*(1.-c)*(1.-c);
    };
    
    Set::Scalar Dw_phi(Set::Scalar c)
    {
        return -1.*Dg_phi(c);
        //return 2.0*c*(1.-c)*(1.-2.*c);
	    //return -4.*(1.-c)*(1.-c)*(1.-c);
    };

    Set::Scalar g_phi(Set::Scalar c)
    {
        return 4.*c*c*c - 3.*c*c*c*c;
        //return c*c;
	    //return c*c*c*c;
    };

    Set::Scalar Dg_phi(Set::Scalar c)
    {
        return 12.*c*c*(1.-c);
        //return 4.*c*c*c;
        //return 2.*c;
    };*/

    //Set::Scalar g_phi(Set::Scalar c) { return 1. - w_phi(c); }    
    //Set::Scalar Dg_phi(Set::Scalar c) { return -Dw_phi(c);}
    
    Set::Scalar Epc(Set::Scalar /*theta*/)
    {
        return mult_Gc*Gc/(4.0*zeta);
    };
    Set::Scalar kappa(Set::Scalar /*theta*/)
    {
        return mult_Lap*2.0*Gc*zeta;
    };

private:
    amrex::Real Gc = NAN, zeta = NAN, mobility = NAN;
    amrex::Real mult_Gc = NAN, mult_Lap = NAN;
};
}
}
}
#endif
