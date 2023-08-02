//
// This model implements the base class for Crack models.
// The class provides various options for degradation functions :math:`g(c)` and :math:`w(c)`.
// See `this link <https://www.sciencedirect.com/science/article/pii/S0927025615003602?via%3Dihub>`_
// for more information on degradation functions.
//
//
// Class methods:
// 
// #. :code:`Crack()`: 
//    Basic constructor. Does nothing, and leaves all values initiated as NAN.
// #. :code:`virtual Set::Scalar w_phi(Set::Scalar c, Set::Scalar /*p=0.*/)`
//    Returns the value for :math:`w(c)`
// #. :code:`virtual Set::Scalar g_phi(Set::Scalar c, Set::Scalar p=0.)`
//    Returns the value of :math:`g(c)` given the value of crack field :math:`c` and plastic strain :math:`p`.
// #. :code:`virtual Set::Scalar Dw_phi(Set::Scalar c, Set::Scalar /*p=0.*/)`
//    Returns the value of the first derivative of :math:`w(c)`
// #. :code:`virtual Set::Scalar Dg_phi(Set::Scalar c, Set::Scalar p=0.)`
//    Returns the value of the first derivative of :math:`g(c)`
// #. :code:`virtual Set::Scalar Gc(Set::Scalar theta) = 0;`
//    Virtual function to return the value of fracture energy :math:`G(c)`
// #. :code:`virtual Set::Scalar DGc(Set::Scalar theta) = 0;`
//    Virtual function to return the value of the derivative of fracture energy :math:`G(c)` with angle :math:`\theta`
// #. :code:`virtual Set::Scalar DDGc(Set::Scalar theta) = 0;`
//    Virtual function to return the value of the second derivative of fracture energy :math:`G(c)` with angle :math:`\theta`
// #. :code:`virtual Set::Scalar Zeta(Set::Scalar theta) = 0;`
//    Virtual function to return the crack length scale :math:`\xi`
// #. :code:`virtual Set::Scalar Mobility(Set::Scalar theta) = 0;`
//    Virtual function to return the crack mobility :math:`M`
// #. :code:`virtual Set::Scalar DrivingForceThreshold(Set::Scalar theta) = 0;`
//    Virtual function to return the driving force threshold.
// #. :code:`void ExportToFile(std::string filename, amrex::Real dTheta)`
//    Writes :math:`G_c(\theta)` as a function of :math:`\theta` to a file.
// #. :code:`void SetGType(const GType a_type)`
//    Sets the type of the degradation function :math:`g(c)`
// #. :code:`void SetWType(const WType a_type)`
//    Sets the type of the interpolation function :math:`w(c)`
// #. :code:`void SetDuctileExponent(const Set::Scalar m)`
//    Sets exponent :math:`m` for functions involving plastic strains.
//
#ifndef MODEL_INTERFACE_CRACK_H
#define MODEL_INTERFACE_CRACK_H

#include <AMReX.H>
#include <AMReX_AmrCore.H>

#include <iostream>
#include <fstream>

namespace Model
{
namespace Interface
{
namespace Crack
{
class Crack
{
    public:
    enum GType {GSQUARE, GMULTIWELL, GPHI4C3, GPHI3C2, GCUBICM, GSQUAREP, GSQUAREPM};
    enum WType {WSQUARE, WMULTIWELL, WMULTIWELL2, WPHI4C3};
    
    Crack() {};
    
    AMREX_FORCE_INLINE
    virtual Set::Scalar w_phi(Set::Scalar c, Set::Scalar /*p=0.*/)
    {
        switch(w_type)
        {
            case WSQUARE:        return (1.-c)*(1.-c);
            case WMULTIWELL:        return (1.-c)*(1.-c)*c*c;
            case WMULTIWELL2:    return (1.+c)*(1.+c)*(1.-c)*(1.-c);
            case WPHI4C3:        return 1. - 4.*c*c*c + 3.*c*c*c*c;
            
            default:            return (1.-c)*(1.-c);
        }
    }
    AMREX_FORCE_INLINE
    virtual Set::Scalar g_phi(Set::Scalar c, Set::Scalar p=0.)
    {
        switch(g_type)
        {
            case GSQUARE:        return c*c;
            case GMULTIWELL:     return (2.-c)*(2.-c)*c*c;
            case GPHI4C3:         return 4.*c*c*c - 3.*c*c*c*c;
            case GPHI3C2:       return 3*c*c - 2*c*c*c;
            case GSQUAREP:        return std::pow(c,2.*p);
            case GSQUAREPM:        return std::pow(c,2.*(std::pow(p,m_d_exponent)));
            case GCUBICM:        return m_d_exponent*(c*c*c - c*c) + 3.*c*c - 2.*c*c*c;

            default:             return c*c;
        }
    }
    AMREX_FORCE_INLINE
    virtual Set::Scalar Dw_phi(Set::Scalar c, Set::Scalar /*p=0.*/)
    {
        switch(w_type)
        {
            case WSQUARE:        return -2.*(1.-c);
            case WMULTIWELL:        return 4.*c*c*c - 6.*c*c + 2.*c;
            case WMULTIWELL2:    return 4.*c*c*c - 4.*c;
            case WPHI4C3:        return 12.*(c-1.)*c*c;
            
            default:            return -2.*(1.-c);
        }
    }
    AMREX_FORCE_INLINE
    virtual Set::Scalar Dg_phi(Set::Scalar c, Set::Scalar p=0.)
    {
        switch(g_type)
        {
            case GSQUARE:        return 2.*c;
            case GMULTIWELL:     return 4.*c*c*c - 12.*c*c + 8.*c;
            case GPHI4C3:         return 12.*(1.-c)*c*c;
            case GPHI3C2:       return 6*c - 6*c*c;
            case GSQUAREP:        return 2.*p*std::pow(c,2*p -1.);
            case GSQUAREPM:        return 2.*std::pow(p,m_d_exponent)*(std::pow(c, 2*std::pow(p,m_d_exponent)-1));
            case GCUBICM:        return m_d_exponent*(3.*c*c - 2.*c) + 6.*c - 6.*c*c;

            default:             return 2*c;
        }
    }

    virtual Set::Scalar Gc(Set::Scalar theta) = 0;
    virtual Set::Scalar DGc(Set::Scalar theta) = 0;
    virtual Set::Scalar DDGc(Set::Scalar theta) = 0;
    virtual Set::Scalar Zeta(Set::Scalar theta) = 0;
    virtual Set::Scalar Mobility(Set::Scalar theta) = 0;
    virtual Set::Scalar DrivingForceThreshold(Set::Scalar theta) = 0;

    void ExportToFile(std::string filename, amrex::Real dTheta)
    {
        std::ofstream outFile;
        outFile.open(filename);
    
        for(amrex::Real theta=0; theta<2*pi; theta=theta+dTheta)
            {
                outFile << theta << " " << Gc(theta) << std::endl;   
            }
        outFile.close();
    
    }
    void SetGType(const GType a_type)
    {
        g_type = a_type;
    }

    void SetWType(const WType a_type)
    {
        w_type = a_type;
    }
    
    void SetDuctileExponent(const Set::Scalar m)
    {
        m_d_exponent = m;
    }

protected:
    static constexpr amrex::Real pi = 3.14159265359;
    GType g_type = GType::GSQUARE;
    WType w_type = WType::WSQUARE;
    Set::Scalar m_d_exponent = 1.;
};
}
}
}

#endif
