//
// This SCP model implements the method described in the following paper
//
// .. bibliography::
//    :list: none
//    :filter: False
// 
//    meier2024diffuse
//
// This method imitates inheritance :code:`Propellant::Propellant` using static polymorphism
// 

#ifndef MODEL_PROPELLANT_FULLFEEDBACK_H
#define MODEL_PROPELLANT_FULLFEEDBACK_H

#include "Propellant.H"

namespace Model
{
namespace Propellant
{
class FullFeedback
{
public:
    static constexpr const char* name = "fullfeedback";

    FullFeedback() = default;
    FullFeedback(IO::ParmParse &pp, std::string name) 
    {
        pp_queryclass(name,*this);
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    void set_pressure(const Set::Scalar P) 
    {
        this->P = P;
        
        zeta_2 = 0.000045 - P * 6.42e-6;
        if (arrhenius_dependency == 1) zeta_1 = zeta_2;
        else zeta_1 = zeta_0;

        // Pressure-dependent post calculation variables
        k1 = a1 * P + b1 - zeta_1 / zeta;
        k2 = a2 * P + b2 - zeta_1 / zeta;
        k3 = 4.0 * log((c1 * P * P + a3 * P + b3) - k1 / 2.0 - k2 / 2.0);
        return;
    }
    
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar get_qdot(const Set::Scalar mdot, const Set::Scalar phi) 
    {

        Set::Scalar qflux = k1 * phi +
            k2 * (1.0 - phi) +
            (zeta_1 / zeta) * exp(k3 * phi * (1.0 - phi));
        Set::Scalar mlocal = (mlocal_ap) * phi + (mlocal_htpb) * (1.0 - phi) + mlocal_comb * phi * (1.0 - phi);
        Set::Scalar mdota = fabs(mdot);
        Set::Scalar mbase = tanh(4.0 * mdota / (mlocal));

        return mbase * qflux;
    }


    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar get_K (const Set::Scalar phi)
    {
        return k_ap * phi + k_htpb * (1.0 - phi);
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar get_rho (const Set::Scalar phi)
    {
        return rho_ap * phi + rho_htpb * (1.0 - phi); 
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar get_cp (const Set::Scalar phi)
    {
        return cp_ap * phi + cp_htpb * (1.0 - phi);
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar get_L (Set::Scalar phi, Set::Scalar T) 
    {
        Set::Scalar L = NAN;
        if (mob_ap == 1) L = m_ap * P * exp(-E_ap / T) * phi;
        else L = m_ap * exp(-E_ap / T) * phi;
        L += m_htpb * exp(-E_htpb / T) * (1.0 - phi);
        if (T <= bound) L = 0.0;
        return L;
    }


    static void Parse(FullFeedback& value, IO::ParmParse& pp)
    {
        // AP/HTPB interface length
        pp_query_default("phi.zeta", value.zeta, 1.0e-5); 
        // Reference interface length for heat integration
        pp_query_default("phi.zeta_0", value.zeta_0, 1.0e-5); 
        
        // Surrogate heat flux model paramater - AP
        pp_query_required("a1", value.a1); 
        // Surrogate heat flux model paramater - HTPB
        pp_query_required("a2", value.a2); 
        // Surrogate heat flux model paramater - Total
        pp_query_required("a3", value.a3); 
        // Surrogate heat flux model paramater - AP
        pp_query_required("b1", value.b1); 
        // Surrogate heat flux model paramater - HTPB
        pp_query_required("b2", value.b2); 
        // Surrogate heat flux model paramater - Total
        pp_query_required("b3", value.b3); 
        // Surrogate heat flux model paramater - Total
        pp_query_required("c1", value.c1); 

        // Whether to use pressure to determined the reference Zeta 
        pp_query_default("pressure.dependency", value.arrhenius_dependency, true); 

        // AP mass flux reference value 
        pp_query_default("mlocal_ap", value.mlocal_ap, 0.0); 
        // HTPB mass flux reference value 
        pp_query_default("mlocal_comb", value.mlocal_comb, 0.0); 

        Set::Scalar massfraction = 0.8;
        value.mlocal_htpb = 685000.0 - 850e3 * massfraction; // Replicating but needs to be revised!

        // AP Density
        pp_query_required("rho_ap", value.rho_ap);
        // HTPB Density
        pp_query_required("rho_htpb", value.rho_htpb);
        // AP Thermal Conductivity
        pp_query_required("k_ap", value.k_ap);
        // HTPB Thermal Conductivity
        pp_query_required("k_htpb", value.k_htpb);
        // AP Specific Heat
        pp_query_required("cp_ap", value.cp_ap);
        //HTPB Specific Heat
        pp_query_required("cp_htpb", value.cp_htpb); 

        // AP Pre-exponential factor for Arrhenius Law
        pp_query_required("m_ap", value.m_ap); 
        // HTPB Pre-exponential factor for Arrhenius Law
        pp_query_required("m_htpb", value.m_htpb); 
        // AP Activation Energy for Arrhenius Law
        pp_query_required("E_ap", value.E_ap); 
        // HTPB Activation Energy for Arrhenius Law
        pp_query_required("E_htpb", value.E_htpb); 
        // Whether to include pressure to the arrhenius law [??]
        pp_query_default("mob_ap", value.mob_ap,0); 
        // HTPB Activation Energy for Arrhenius Law
        pp_query_default("bound", value.bound,0.0); 

    }

    // IO Variables
    Set::Scalar zeta = NAN, zeta_0 = NAN; 
    Set::Scalar a1 = NAN, a2 = NAN, a3 = NAN;
    Set::Scalar b1 = NAN, b2 = NAN, b3 = NAN;
    Set::Scalar c1 = NAN;
    //Set::Scalar h1 = NAN, h2 = NAN;
    bool arrhenius_dependency = false;
    Set::Scalar mlocal_ap = NAN, mlocal_comb = NAN, mlocal_htpb = NAN;

    //Thermal properties
    Set::Scalar k_ap = NAN, k_htpb = NAN;
    Set::Scalar rho_ap = NAN, rho_htpb = NAN;
    Set::Scalar cp_ap = NAN, cp_htpb = NAN;

    // Regression variables
    int mob_ap = -1;
    Set::Scalar m_ap = NAN, m_htpb = NAN, m_comb = NAN;
    Set::Scalar E_ap = NAN, E_htpb = NAN;
    Set::Scalar bound = NAN;

    // Calculated variables
    Set::Scalar P = NAN;
    Set::Scalar zeta_1 = NAN, zeta_2 = NAN;
    Set::Scalar k1=NAN, k2=NAN, k3=NAN, k4=NAN;
};

}
}



#endif


