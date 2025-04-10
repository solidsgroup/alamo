#ifndef MODEL_REGRESSION_POWERLAW_H
#define MODEL_REGRESSION_POWERLAW_H

#include "Regression.H"
#include "IO/ParmParse.H"

namespace Model
{
namespace Regression
{

class PowerLaw : public Regression
{
public:
    static constexpr const char* name = "powerlaw";
    
    PowerLaw() {} ;
    PowerLaw(IO::ParmParse &pp) {pp_queryclass(*this);} 
    PowerLaw(IO::ParmParse &pp, std::string name) {pp_queryclass(name, *this);} 
    virtual ~PowerLaw(){};

    virtual Set::Scalar
    rate(const Set::Scalar ap_frac) override final
    {
        Set::Scalar L = 0.0;
        if (homogeneous)
        {
            Util::Abort(INFO,"I don't this this is actually still used?");
            // Set::Scalar angle = acos(n(0)) * 180. / Set::Constant::Pi;
            // if (angle > 90) angle = angle - 90.0;
            // if (angle > 180) angle = angle - 180.0;
            // if (angle > 270) angle = angle - 270.0;
            // L = a_fit + b_fit * exp(-c_fit * angle);
        }
        else
        {
            L = fmod_ap * ap_frac
                + fmod_htpb * (1.0 - ap_frac)
                + 4.0 * fmod_comb * ap_frac * (1.0 - ap_frac);
        }
        return L;
    }

    static void
    Parse(PowerLaw &value, IO::ParmParse &pp)
    {
        // Fixed pressure
        pp_query_required("pressure", value.pressure);
        // AP power pressure law parameter (r*P^n)
        pp_query_required("r_ap", value.r_ap);
        // HTPB power pressure law parameter (r*P^n)
        pp_query_required("r_htpb", value.r_htpb);
        // AP/HTPB power pressure law parameter (r*P^n)
        pp_query_required("r_comb", value.r_comb);
        // AP power pressure law parameter (r*P^n)
        pp_query_required("n_ap", value.n_ap);
        // HTPB power pressure law parameter (r*P^n)
        pp_query_required("n_htpb", value.n_htpb);
        // AP/HTPB power pressure law parameter (r*P^n)
        pp_query_required("n_comb", value.n_comb); 
        // Whether to initialize Phi with homogenized properties
        pp_query_default("homogeneousSystem", value.homogeneous, false); 


        value.fmod_ap   = value.r_ap   * pow(value.pressure, value.n_ap);
        value.fmod_htpb = value.r_htpb * pow(value.pressure, value.n_htpb);
        value.fmod_comb = value.r_comb * pow(value.pressure, value.n_comb);

        value.a_fit = -1.16582   * sin(value.pressure) - 0.681788 * cos(value.pressure) + 3.3563;
        value.b_fit = -0.708225  * sin(value.pressure) + 0.548067 * cos(value.pressure) + 1.55985;
        value.c_fit = -0.0130849 * sin(value.pressure) - 0.03597  * cos(value.pressure) + 0.00725694;
    }

private:
    // Internal variables
    Set::Scalar pressure = NAN;
    Set::Scalar fmod_ap = NAN, fmod_htpb = NAN, fmod_comb = NAN;
    Set::Scalar a_fit = NAN, b_fit = NAN, c_fit = NAN;

    // IO Params
    Set::Scalar r_ap = NAN, r_htpb = NAN, r_comb = NAN;
    Set::Scalar n_ap = NAN, n_htpb = NAN, n_comb = NAN;
    bool homogeneous = false;
};
}
}

#endif
