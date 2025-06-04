#ifndef MODEL_REGRESSION_REGRESSION_H
#define MODEL_REGRESSION_REGRESSION_H

#include "Set/Base.H"
#include "IO/ParmParse.H"


namespace Model
{
namespace Regression
{

class Regression
{
public:
    virtual Set::Scalar L() = 0;

    bool needs_normal = false;


    Set::Scalar P = NAN;
    bool needs_pressure = false;
    virtual void set_pressure(const Set::Scalar P) 
    {
        this->P = P;
    }

    Set::Scalar temperature = NAN;
    bool needs_temperature = false;
    virtual void set_temperature(const Set::Scalar temperature) 
    {
        this->temperature = temperature;
    }

    Set::Scalar phi1 = NAN;
    bool needs_phi = false;
    virtual void set_phi(const Set::Scalar phi1)
    {
        this->phi1 = phi1;
    }

    Set::Scalar deltaw = NAN;
    bool needs_deltaw = false;
    virtual void set_deltaw(const Set::Scalar deltaw)
    {
        this->deltaw = deltaw;
    }
};


class Arrhenius : public Regression
{
public:
    static constexpr const char* name = "arrhenius";

    Arrhenius() = default;
    Arrhenius(IO::ParmParse &pp, std::string name) 
    {
        pp_queryclass(name,*this);
    }

    virtual Set::Scalar L() override
    {
        Set::Scalar L = NAN;
        if (mob_ap == 1) L = m_ap * P * exp(-E_ap / temperature) * phi1;
        else L = m_ap * exp(-E_ap / temperature) * phi1;
        L += m_htpb * exp(-E_htpb / temperature) * (1.0 - phi1);
        if (temperature <= bound) L = 0.0;
        return L;
    }


    int mob_ap = NAN;
    Set::Scalar m_ap = NAN, m_htpb = NAN, m_comb = NAN;
    Set::Scalar E_ap = NAN, E_htpb = NAN;
    Set::Scalar bound = NAN;

    static void Parse(Arrhenius& value, IO::ParmParse& pp)
    {
        // AP Pre-exponential factor for Arrhenius Law
        pp_query_required("m_ap", value.m_ap); 
        // HTPB Pre-exponential factor for Arrhenius Law
        pp_query_required("m_htpb", value.m_htpb); 
        // 
        pp_query_required("m_comb", value.m_comb); 

        // AP Activation Energy for Arrhenius Law
        pp_query_required("E_ap", value.E_ap); 
        // HTPB Activation Energy for Arrhenius Law
        pp_query_required("E_htpb", value.E_htpb); 
        
        // HTPB Activation Energy for Arrhenius Law
        pp_query_default("bound", value.bound,1E100); 
    }
};


class PowerLaw : public Regression
{
public:
    static constexpr const char* name = "powerlaw";
    
    PowerLaw() = default;
    PowerLaw(IO::ParmParse &pp, std::string name) 
    {
        pp_queryclass(name,*this);
    }

    virtual Set::Scalar L() override
    {
        Set::Scalar L = NAN;
        // if (homogeneousSystem == 1) {
        //     Set::Vector grad_eta = Numeric::Gradient(eta, i, j, k, 0, DX);
        //     Set::Scalar angle = acos(grad_eta[0] / grad_eta.lpNorm<2>()) * 180 / 3.1415;
        //     if (angle > 90) angle = angle - 90.0;
        //     if (angle > 180) angle = angle - 180.0;
        //     if (angle > 270) angle = angle - 270.0;
        //     L = pressure.power.a_fit + pressure.power.b_fit * exp(-pressure.power.c_fit * angle);
        // }
        // else 
        {
            Set::Scalar fs_actual = fmod_ap * phi1
                + fmod_htpb * (1.0 - phi1)
                + 4.0 * fmod_comb * phi1 * (1.0 - phi1);
            L = fs_actual / gamma / deltaw /*(pf.w1 - pf.w0)*/;
        }

        return L;
    }

    virtual void set_pressure(const Set::Scalar P) override
    {
        // for the pressure power law option
        fmod_ap = r_ap * pow(P, n_ap);
        fmod_htpb = r_htpb * pow(P, n_htpb);
        fmod_comb = r_comb * pow(P, n_comb);

        // a_fit = -1.16582 * sin(P) - 0.681788 * cos(P) + 3.3563;
        // b_fit = -0.708225 * sin(P) + 0.548067 * cos(P) + 1.55985;
        // c_fit = -0.0130849 * sin(P) - 0.03597 * cos(P) + 0.00725694;
    }

    static void Parse(PowerLaw& value, IO::ParmParse& pp)
    {
        Util::Message(INFO,pp.getPrefix());
        pp_query_required("gamma", value.gamma);   // Scaling factor for mobility
        pp_query_required("r_ap", value.r_ap);     // AP power pressure law parameter (r*P^n)
        pp_query_required("r_htpb", value.r_htpb); // HTPB power pressure law parameter (r*P^n)
        pp_query_required("r_comb", value.r_comb); // AP/HTPB power pressure law parameter (r*P^n)
        pp_query_required("n_ap", value.n_ap);     // AP power pressure law parameter (r*P^n)
        pp_query_required("n_htpb", value.n_htpb); // HTPB power pressure law parameter (r*P^n)
        pp_query_required("n_comb", value.n_comb); // AP/HTPB power pressure law parameter (r*P^n)

        value.needs_pressure = true;
        value.needs_phi = true;
        value.needs_deltaw = true;
    }

    Set::Scalar gamma = NAN;
    Set::Scalar r_ap = NAN, r_htpb = NAN, r_comb = NAN;
    Set::Scalar n_ap = NAN, n_htpb = NAN, n_comb = NAN;
    Set::Scalar fmod_ap = NAN, fmod_htpb = NAN, fmod_comb = NAN;
    //Set::Scalar a_fit = NAN, b_fit = NAN, c_fit = NAN;
};


}
}




#endif
