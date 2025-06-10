#ifndef MODEL_PROPELLANT_HOMOGENIZE_H
#define MODEL_PROPELLANT_HOMOGENIZE_H

#include "Model/Propellant/Propellant.H"

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

    void set_pressure(const Set::Scalar P) 
    {
        // TODO
    }
    
    virtual Set::Scalar get_qdot(const Set::Scalar mdot, const Set::Scalar phi) 
    {
        return NAN;
    }

    Set::Scalar get_K (const Set::Scalar phi)
    {
        return (k_ap * massfraction + k_htpb * (1.0 - massfraction)) * phi + dispersion1 * (1. - phi);
    }

    Set::Scalar get_rho (const Set::Scalar phi)
    {
        return (rho_ap * massfraction + rho_htpb * (1.0 - massfraction)) * phi + dispersion2 * (1. - phi); 
    }

    Set::Scalar get_cp (const Set::Scalar phi)
    {
        return (cp_ap * massfraction + cp_htpb * (1.0 - massfraction)) * phi + dispersion3 * (1. - phi);
    }

    Set::Scalar get_L (Set::Scalar phi, Set::Scalar T) 
    {
        return NAN;
    }


    static void Parse(Homogenize& value, IO::ParmParse& pp)
    {
        
        // K; dispersion variables are use to create an inert region for the void grain case. 
        // An inert region is one that dissipates energy fast enough to remove regression and 
        // thermal strain effects.
        pp_query("dispersion1", value.dispersion1);
        // rho; dispersion variables are use to create an inert region for the void grain case. 
        // An inert region is one that dissipates energy fast enough to remove regression and 
        // thermal strain effects.
        pp_query("dispersion2", value.dispersion2);
        // cp; dispersion variables are use to create an inert region for the void grain case. 
        // An inert region is one that dissipates energy fast enough to remove regression and 
        // thermal strain effects.
        pp_query("dispersion3", value.dispersion3); 


        // AP Density
        pp_query_required("thermal.rho_ap", value.rho_ap);
        // HTPB Density
        pp_query_required("thermal.rho_htpb", value.rho_htpb);
        // AP Thermal Conductivity
        pp_query_required("thermal.k_ap", value.k_ap);
        // HTPB Thermal Conductivity
        pp_query_required("thermal.k_htpb", value.k_htpb);
        // AP Specific Heat
        pp_query_required("thermal.cp_ap", value.cp_ap);
        //HTPB Specific Heat
        pp_query_required("thermal.cp_htpb", value.cp_htpb); 

    }

    Set::Scalar massfraction = NAN;
    Set::Scalar dispersion1 = NAN, dispersion2 = NAN, dispersion3 = NAN;

    Set::Scalar k_ap = NAN, k_htpb = NAN;
    Set::Scalar rho_ap = NAN, rho_htpb = NAN;
    Set::Scalar cp_ap = NAN, cp_htpb = NAN;

    // TODO
};

}
}



#endif


