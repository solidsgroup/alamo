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

    void set_pressure(const Set::Scalar a_P) 
    {
        this->Pbar = a_P/Pref;
        
        zeta_2 = zeta_fit_a + Pbar * zeta_fit_b;

        if (arrhenius_dependency == 1) zeta_1 = zeta_2; // Reaction rate equation (Arrhenius Eqn)
        else zeta_1 = zeta_0;

        // Pressure-dependent post calculation variables
        k1 = a1 * Pbar + b1 - zeta_1 / zeta;
        k2 = a2 * Pbar + b2 - zeta_1 / zeta;
        k3 = 4.0 * log((c1*Pbar*Pbar + a3*Pbar + b3) - k1 / 2.0 - k2 / 2.0);

        return;
    }
    
    AMREX_FORCE_INLINE
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
        if (mob_ap) L = m_ap * Pbar * exp(-E_ap / T) * phi;
        else L = m_ap * exp(-E_ap / T) * phi;
        L += m_htpb * exp(-E_htpb / T) * (1.0 - phi);
        if (T <= bound) L = 0.0;
        return L;
    }


    static void Parse(FullFeedback& value, IO::ParmParse& pp)
    {
        // AP/HTPB interface length
        pp.query_default("phi.zeta", value.zeta, "10.0_um", Unit::Length()); 
        // Reference interface length for heat integration
        pp.query_default("phi.zeta_0", value.zeta_0, "10.0_um", Unit::Length()); 
        
        // Constant offset for $\zeta$ adjustment with pressure
        pp.query_default("phi.zeta_fit_a",value.zeta_fit_a, "45.0_um", Unit::Length());
        // Scaling factors for $\zeta$ adjustment with pressure
        pp.query_default("phi.zeta_fit_b",value.zeta_fit_b, "-6.42_um", Unit::Length());
        
        // Surrogate heat flux model paramater - AP
        pp_query_required("a1", value.a1, Unit::Less()); 
        // Surrogate heat flux model paramater - HTPB
        pp_query_required("a2", value.a2,Unit::Less()); 
        // Surrogate heat flux model paramater - Total
        pp_query_required("a3", value.a3,Unit::Less()); 
        // Surrogate heat flux model paramater - AP
        pp_query_required("b1", value.b1, Unit::Less()); 
        // Surrogate heat flux model paramater - HTPB
        pp_query_required("b2", value.b2, Unit::Less()); 
        // Surrogate heat flux model paramater - Total
        pp_query_required("b3", value.b3,Unit::Less()); 
        // Surrogate heat flux model paramater - Total
        pp_query_required("c1", value.c1, Unit::Less()); 

        // Whether to use pressure to determined the reference Zeta 
        pp_query_default("pressure.dependency", value.arrhenius_dependency, true); 

        // Reference pressure for pressure nondimensionalization
        pp_query_default("Pref",value.Pref,"1_MPa",Unit::Pressure());
        
        // AP volumetric mass flux reference value 
        pp.query_default("mlocal_ap", value.mlocal_ap, "0.0", Unit::Density()/Unit::Time()); 
        // combined volumetric mass flux reference value 
        pp_query_default("mlocal_comb", value.mlocal_comb, "0.0", Unit::Density()/Unit::Time()); 
        // HTPB volumetric mass flux reference value
        pp.query_default("mlocal_htpb",value.mlocal_htpb, "5000_kg/m^3/s", Unit::Density()/Unit::Time());

        // Set::Scalar massfraction = 0.8;
        // value.mlocal_htpb = 685000.0 - 850e3 * massfraction; // Replicating but needs to be revised!
        // Util::Message(INFO,value.mlocal_htpb);

        // AP Density
        pp_query_required("rho_ap", value.rho_ap, Unit::Density());
        // HTPB Density
        pp_query_required("rho_htpb", value.rho_htpb, Unit::Density());
        // AP Thermal Conductivity
        pp_query_required("k_ap", value.k_ap, Unit::ThermalConductivity());
        // HTPB Thermal Conductivity
        pp_query_required("k_htpb", value.k_htpb, Unit::ThermalConductivity());
        // AP Specific Heat
        pp_query_required("cp_ap", value.cp_ap,Unit::SpecificHeatCapacity());
        //HTPB Specific Heat
        pp_query_required("cp_htpb", value.cp_htpb,Unit::SpecificHeatCapacity()); 

        // AP Pre-exponential factor for Arrhenius Law
        pp_query_required("m_ap", value.m_ap,Unit::Velocity()); 
        // HTPB Pre-exponential factor for Arrhenius Law
        pp_query_required("m_htpb", value.m_htpb,Unit::Velocity()); 
        // COMBINED AP Activation Energy, divided by Rg, for Arrhenius Law (has units of temperature)
        pp_query_required("E_ap", value.E_ap,Unit::Temperature()); 
        // COMBINED HTPB Activation Energy, divided by Rg, for Arrhenius Law(has units of temperature)
        pp_query_required("E_htpb", value.E_htpb,Unit::Temperature()); 
        // Whether to include pressure to the arrhenius law [??]n
        pp_query_default("mob_ap", value.mob_ap,false); 
        // HTPB Activation Energy for Arrhenius Law
        pp_query_default("bound", value.bound,"0.0", Unit::Temperature()); 
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
    bool mob_ap = false;
    Set::Scalar m_ap = NAN, m_htpb = NAN, m_comb = NAN;
    Set::Scalar E_ap = NAN, E_htpb = NAN;
    Set::Scalar bound = NAN;
    Set::Scalar Pref = NAN;

    // Calculated variables
    Set::Scalar Pbar = NAN;
    Set::Scalar zeta_1 = NAN, zeta_2 = NAN;
    Set::Scalar zeta_fit_a = NAN, zeta_fit_b = NAN;
    Set::Scalar k1=NAN, k2=NAN, k3=NAN, k4=NAN;
};

}
}



#endif
