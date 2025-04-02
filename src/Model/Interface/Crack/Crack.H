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
    enum GType  {GSQUARE, GMULTIWELL, GPHI4C3, GCUBICM, GSQUAREP, GSQUAREPM, GWULINEAR, GAICUBIC};
    enum G2Type {G2WULINEAR};
    enum WType  {WSQUARE, WMULTIWELL, WMULTIWELL2, WPHI4C3, WU};
    
    Crack() {};
    
    // Geometric function involved in representing sharp crack as a band.
    // Often this function is given in terms of a variable d = 1-c.
    // Here we are implementing w as a function of c.
    AMREX_FORCE_INLINE
    virtual Set::Scalar w_phi(Set::Scalar c, Set::Scalar /*p=0.*/)
    {
        switch(w_type)
        {
            case WSQUARE:       return (1.-c)*(1.-c);
            case WMULTIWELL:    return (1.-c)*(1.-c)*c*c;
            case WMULTIWELL2:   return (1.+c)*(1.+c)*(1.-c)*(1.-c);
            case WPHI4C3:       return 1. - 4.*c*c*c + 3.*c*c*c*c;
            case WU:            return 2.0*(1.-c) - (1.-c)*(1.-c);
            
            default:            return (1.-c)*(1.-c);
        }
    }
    AMREX_FORCE_INLINE
    virtual Set::Scalar g_phi(Set::Scalar c, Set::Scalar p=0.)
    {
        switch(g_type)
        {
            case GSQUARE:       return c*c;
            case GMULTIWELL:    return (2.-c)*(2.-c)*c*c;
            case GPHI4C3:       return 4.*c*c*c - 3.*c*c*c*c;
            case GSQUAREP:      return std::pow(c,2.*p);
            case GSQUAREPM:     return std::pow(c,2.*(std::pow(p,m_d_exponent)));
            case GCUBICM:       return m_d_exponent*(c*c*c - c*c) + 3.*c*c - 2.*c*c*c;
            case GWULINEAR:     {
                                    // this is for a linear softening curve for the CZM model
                                    // Obtain from Wu 2018 CMAME
                                    // w(c) = c^p / (c^p + Q(c))
                                    // where, Q(c) = a1 (1-c) + a1 a2 (1-c)^2 + a1 a2 a3 (1-c)^3
                                    // for linear law with p=2, a1 = 4 * l_w / (pi * b)
                                    // a2 = -1/2, a3 = 0
                                    Set::Scalar d = 1.0-c;
                                    return (c*c) / (c*c +   ((2.0 * m_pf_czm_a0 * d) - (m_pf_czm_a0 * d * d)));
                                    // if (c < 1.e-5) c = 1.e-5;
                                    // double phi = m_pf_czm_a0 * m_pf_czm_order* std::sqrt(w_phi(c,0.));
                                    // phi *= std::sqrt(1 - std::pow(c, 2*m_pf_czm_order));
                                    // phi /= std::pow(c, m_pf_czm_order+1);
                                    // return 1.0 / (1.0 + phi);
                                }
            case GAICUBIC:      {
                                    Set::Scalar d = 1.0 - c;
                                    return 1 - (m_pf_czm_a0 * d * d) + ((m_pf_czm_a0 - 1) * d * d * d);
                                }
            default:            return c*c;
        }
    }

    // The derivative of the geometric function is often written in terms of w'(d).
    // Here, we are taking w'(c). It's better to just take the derivative of the function
    // we implemented in w_phi(c) with respect to c.
    AMREX_FORCE_INLINE
    virtual Set::Scalar Dw_phi(Set::Scalar c, Set::Scalar /*p=0.*/)
    {
        switch(w_type)
        {
            case WSQUARE:           return -2.*(1.-c);
            case WMULTIWELL:        return 4.*c*c*c - 6.*c*c + 2.*c;
            case WMULTIWELL2:       return 4.*c*c*c - 4.*c;
            case WPHI4C3:           return 12.*(c-1.)*c*c;
            case WU:                return -2.*c;
            
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
            case GSQUAREP:        return 2.*p*std::pow(c,2*p -1.);
            case GSQUAREPM:        return 2.*std::pow(p,m_d_exponent)*(std::pow(c, 2*std::pow(p,m_d_exponent)-1));
            case GCUBICM:        return m_d_exponent*(3.*c*c - 2.*c) + 6.*c - 6.*c*c;
            case GWULINEAR:      {
                                    // Derivative of c^2 / (c^2 + Q(c))
                                    // 2c / (c^2 + Q(c)) - c^2 (2c + Q'(c)) / (c^2+Q(c))^2
                                    // ( 2c (c^2 + Q(c)) - c^2 (2c + Q'(c)) ) / (c^2 + Q(c))^2
                                    // ( 2c Q(c) - c^2 Q'(c) ) / (c^2 + Q(c))^2
                                    // Q(d) = a1 d + a1 a2 d^2 + a1 a2 a3 d^3
                                    // In this case, a = 2*a0, a2 = -0.5, a3 = 0
                                    // Q(d) = 2 a0 d - a0 d^2
                                    // Q'(d) = 2 a0 - 2 a0 d
                                    // Q'(c) - -Q'(d) = -2 a0 + 2 a0 d
                                    Set::Scalar d = 1.0 -c;
                                    Set::Scalar dQ_c = (-2.0*m_pf_czm_a0) + (2.0*m_pf_czm_a0*d);
                                    Set::Scalar Q_c = (2.0*m_pf_czm_a0*d) - (m_pf_czm_a0*d*d);
                                    return (2.0*c*Q_c - c*c*dQ_c) / ((c*c + Q_c)*(c*c + Q_c));
                                }

            case GAICUBIC:      {
                                    Set::Scalar d = 1.0 - c;
                                    return (-1.0) * (- (2.0 * m_pf_czm_a0 * d) + (3.0 * (m_pf_czm_a0 - 1.0)* d * d));
                                }

            default:             return 2*c;
        }
    }

    AMREX_FORCE_INLINE
    virtual Set::Scalar Dg2_phi(Set::Scalar c)
    {
        switch (g2_type)
        {
            case G2WULINEAR:  {
                                if (c < 1.e-5) c = 1.e-5;
                                // Set::Scalar mu = m_pf_czm_a0 * w_phi(c,0.) / (std::pow(c, 2*m_pf_czm_order));
                                Set::Scalar Dmu = m_pf_czm_a0 * ( (2.*m_pf_czm_order*w_phi(c,0.)) + (c * -1. * Dw_phi(c,0.)) ) / (std::pow(c, 1 + 2*m_pf_czm_order));
                                return Dmu * g_phi(c,0.) * g_phi(c,0.);
                            }

            default:    Util::Abort(INFO, "Dissipation function type must be specified");
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

    void SetG2Type (const G2Type a_type)
    {

    }
    
    void SetDuctileExponent(const Set::Scalar m)
    {
        m_d_exponent = m;
    }

    void SetPFCZMConstants(const Set::Scalar a_pf_czm_order, const Set::Scalar a_pf_czm_a0)
    {
        m_pf_czm_order = a_pf_czm_order;
        m_pf_czm_a0 = a_pf_czm_a0;
    }

protected:
    static constexpr Set::Scalar pi = 3.14159265359;
    GType g_type = GType::GSQUARE;
    G2Type g2_type = G2Type::G2WULINEAR;
    WType w_type = WType::WSQUARE;
    Set::Scalar m_d_exponent = 1.;
    Set::Scalar m_pf_czm_order = 1.;
    Set::Scalar m_pf_czm_a0 = 1.0;
};
}
}
}

#endif
