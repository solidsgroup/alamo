#ifndef MODEL_SURFACECOMBUSTION_SURFACECOMBUSTION_H
#define MODEL_SURFACECOMBUSTION_SURFACECOMBUSTION_H

#include "IO/ParmParse.H"

namespace Model
{
namespace SurfaceCombusion
{
class SurfaceCombustion
{
public:
    virtual void set_pressure(const Set::Scalar P){};

    virtual Set::Scalar get_qdot(const Set::Scalar mdot, const Set::Scalar phi) {};
};

class Mesoscale : public SurfaceCombustion
{
public:
    static constexpr const char* name = "mesoscale";

    Mesoscale() = default;
    Mesoscale(IO::ParmParse &pp, std::string name) 
    {
        pp_queryclass(name,*this);
    }


    virtual void set_pressure(const Set::Scalar P) override
    {
        zeta_2 = 0.000045 - P * 6.42e-6;
        if (arrhenius_dependency == 1) zeta_1 = zeta_2;
        else zeta_1 = zeta_0;

        // Pressure-dependent post calculation variables
        k1 = a1 * P + b1 - zeta_1 / zeta;
        k2 = a2 * P + b2 - zeta_1 / zeta;
        k3 = 4.0 * log((c1 * P * P + a3 * P + b3) - k1 / 2.0 - k2 / 2.0);
        k4 = h1 * P + h2;
        return;
    }
    
    virtual Set::Scalar get_qdot(const Set::Scalar mdot, const Set::Scalar phi) override
    {

        Set::Scalar qflux = k1 * phi +
            k2 * (1.0 - phi) +
            (zeta_1 / zeta) * exp(k3 * phi * (1.0 - phi));
        Set::Scalar mlocal = (mlocal_ap) * phi + (mlocal_htpb) * (1.0 - phi) + mlocal_comb * phi * (1.0 - phi);
        Set::Scalar mdota = fabs(mdot);
        Set::Scalar mbase = tanh(4.0 * mdota / (mlocal));
        return mbase * qflux;
        //heatflux(i, j, k) = (thermal.hc * mbase * qflux + laser(i, j, k)) / K;
    }

    static void Parse(Mesoscale& value, IO::ParmParse& pp)
    {
        // AP/HTPB interface length
        pp_query_default("phi.zeta", value.zeta, 1.0e-5); 
        // Reference interface length for heat integration
        pp_query_default("phi.zeta_0", value.zeta_0, 1.0e-5); 
        
        pp_query_required("a1", value.a1); // Surrogate heat flux model paramater - AP
        pp_query_required("a2", value.a2); // Surrogate heat flux model paramater - HTPB
        pp_query_required("a3", value.a3); // Surrogate heat flux model paramater - Total
        pp_query_required("b1", value.b1); // Surrogate heat flux model paramater - AP
        pp_query_required("b2", value.b2); // Surrogate heat flux model paramater - HTPB
        pp_query_required("b3", value.b3); // Surrogate heat flux model paramater - Total
        pp_query_required("c1", value.c1); // Surrogate heat flux model paramater - Total

        // Whether to use pressure to determined the reference Zeta 
        pp_query_default("pressure.dependency", value.arrhenius_dependency, true); 
        // Surgate heat flux model parameter - Homogenized
        pp_query_default("pressure.h1", value.h1, 1.81); 
        // Surgate heat flux model parameter - Homogenized
        pp_query_default("pressure.h2", value.h2, 1.34); 

        // AP mass flux reference value 
        pp_query_default("mlocal_ap", value.mlocal_ap, 0.0); 
        // HTPB mass flux reference value 
        pp_query_default("mlocal_htpb", value.mlocal_htpb, 0.0); 

        Set::Scalar massfraction = 0.8;
        value.mlocal_htpb = 685000.0 - 850e3 * massfraction; // Replicating but needs to be revised!
    }

    // IO Variables
    Set::Scalar zeta = NAN, zeta_0 = NAN; 
    Set::Scalar a1 = NAN, a2 = NAN, a3 = NAN;
    Set::Scalar b1 = NAN, b2 = NAN, b3 = NAN;
    Set::Scalar c1 = NAN;
    Set::Scalar h1 = NAN, h2 = NAN;
    bool arrhenius_dependency = false;
    Set::Scalar mlocal_ap = NAN, mlocal_comb = NAN, mlocal_htpb = NAN;

    // Calculated variables
    Set::Scalar zeta_1 = NAN, zeta_2 = NAN;
    Set::Scalar k1=NAN, k2=NAN, k3=NAN, k4=NAN;
};

}
}



#endif


