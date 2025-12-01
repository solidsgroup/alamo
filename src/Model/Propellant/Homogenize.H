//
// 
//
//
//


#ifndef MODEL_PROPELLANT_HOMOGENIZE_H
#define MODEL_PROPELLANT_HOMOGENIZE_H

#include "IO/ParmParse.H"
#include "Propellant.H"
#include "Util/Util.H"

namespace Model
{
namespace Propellant
{
class Homogenize
{
public:
    static constexpr const char* name = "homogenize";

    Homogenize() = default;
    Homogenize(IO::ParmParse &pp, std::string name) 
    {
        pp_queryclass(name,*this);
    }

    Set::Scalar Pbar = NAN;
    void set_pressure(const Set::Scalar P) 
    {
        this->Pbar = P/P_reference;
        k4 = h1 * P + h2;
    }
    
    virtual Set::Scalar get_qdot(const Set::Scalar mdot, const Set::Scalar phi) 
    {
        Set::Scalar qdot = NAN;

        Set::Scalar qflux = k4 * phi;
        Set::Scalar mdota = fabs(mdot);
        Set::Scalar mbase = tanh(1250 * mdota);

        qdot = mbase * qflux;

        if (qdot != qdot)
        {
            Util::Message(INFO, "mdot = ", mdot);
            Util::Message(INFO, "phi = ", phi);
            Util::Message(INFO, "k4 = ", k4);
            Util::Message(INFO, "mdota = ", mdota);
            Util::Message(INFO, "mlocal = ", mlocal);
        }

        return qdot;
    }

    Set::Scalar get_K (const Set::Scalar phi)
    {
        return k_prop * phi + dispersion1 * (1. - phi);
    }

    Set::Scalar get_rho (const Set::Scalar phi)
    {
        return rho_prop * phi + dispersion2 * (1. - phi); 
    }

    Set::Scalar get_cp (const Set::Scalar phi)
    {
        return cp_prop * phi + dispersion3 * (1. - phi);
    }

    Set::Scalar get_L (Set::Scalar phi, Set::Scalar T) 
    {
        Set::Scalar L = NAN;
        if (mob_prop == 1) L = m_prop * pow(Pbar, pressure_exponent) * exp(-E_prop / T) * phi;
        else L = m_prop * exp(-E_prop / T) * phi;
        if (T <= bound) L = 0.0;
        
        if (L != L)
        {
            Util::Message(INFO,phi);
            Util::Message(INFO,Unit::Temperature(T));
            Util::Message(INFO,Pbar);
            Util::Message(INFO,Unit::Temperature(E_prop));
            Util::Message(INFO,phi);
            Util::Message(INFO,Unit::Velocity(m_prop));
            Util::Abort(INFO,L / Unit::Pressure() / Unit::Time());
        }
            
        return L;
    }


    static void Parse(Homogenize& value, IO::ParmParse& pp)
    {
        
        // K; dispersion variables are use to create an inert region for the void grain case. 
        // An inert region is one that dissipates energy fast enough to remove regression and 
        // thermal strain effects.
        pp.query_required("dispersion1", value.dispersion1, Unit::ThermalConductivity());
        // rho; dispersion variables are use to create an inert region for the void grain case. 
        // An inert region is one that dissipates energy fast enough to remove regression and 
        // thermal strain effects.
        pp.query_required("dispersion2", value.dispersion2, Unit::Density());
        // cp; dispersion variables are use to create an inert region for the void grain case. 
        // An inert region is one that dissipates energy fast enough to remove regression and 
        // thermal strain effects.
        pp.query_required("dispersion3", value.dispersion3, Unit::SpecificHeatCapacity()); 


        // Prop Density
        pp.query_required("rho_prop", value.rho_prop, Unit::Density());
        // Prop Thermal Conductivity
        pp.query_required("k_prop", value.k_prop, Unit::ThermalConductivity());
        // Prop Specific Heat
        pp.query_required("cp_prop", value.cp_prop,Unit::SpecificHeatCapacity());


        // Surgate heat flux model paramater - Homogenized
        pp.query_default("h1", value.h1, "1.81_J/m/K",Unit::HeatFlux()); 
        // Surgate heat flux model paramater - Homogenized
        pp.query_default("h2", value.h2, "1.34_J/m/K",Unit::HeatFlux()); 
        // Reference pressure
        pp.query_default("P_reference", value.P_reference, "1_MPa",Unit::Pressure()); 

        // AP mass flux reference value 
        pp.query_default("mlocal_prop", value.mlocal, "0.0_kg/m^2/s", 
                         Unit::Mass()/Unit::Area()/Unit::Time()); 
        
        // Prop Pre-exponential factor for Arrhenius Law
        pp.query_required("m_prop", value.m_prop, Unit::Velocity()); 
        // Prop Activation Energy for Arrhenius Law
        pp.query_required("E_prop", value.E_prop, Unit::Temperature()); 
        // Whether to include pressure to the arrhenius law [??]
        pp.query_default("mob_prop", value.mob_prop,false); 
        // HTPB Activation Energy for Arrhenius Law
        pp.query_default("bound", value.bound,"0_K", Unit::Temperature()); 

        pp.query_default("pressure_exponent", value.pressure_exponent,"0.5_1", Unit::Less());
    }

    // IO variables
    Set::Scalar massfraction = NAN;
    Set::Scalar dispersion1 = NAN, dispersion2 = NAN, dispersion3 = NAN;
    Set::Scalar k_prop = NAN;
    Set::Scalar rho_prop = NAN;
    Set::Scalar cp_prop = NAN;
    Set::Scalar h1 = NAN, h2=NAN;
    Set::Scalar mlocal = NAN;
    // Regression variables
    bool mob_prop = -1;
    Set::Scalar m_prop = NAN;
    Set::Scalar E_prop = NAN;
    Set::Scalar bound = NAN;
    Set::Scalar P_reference = NAN;
    Set::Scalar pressure_exponent = NAN;

    // Calculated variables
    Set::Scalar k4 = NAN;
    Set::Scalar mlocal_prop = NAN;
};

}
}



#endif


