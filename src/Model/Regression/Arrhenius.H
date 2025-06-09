#ifndef MODEL_REGRESSION_ARRHENIUS_H
#define MODEL_REGRESSION_ARRHENIUS_H

#include "Regression.H"

namespace Model
{
namespace Regression
{
class Arrhenius : public Regression
{
public:
    static constexpr const char* name = "arrhenius";

    Arrhenius() = default;
    Arrhenius(IO::ParmParse &pp, std::string name) 
    {
        pp_queryclass(name,*this);
    }

    virtual Set::Scalar operator() (Set::Scalar phi, Set::Scalar T) override
    {
        Set::Scalar L = NAN;
        if (mob_ap == 1) L = m_ap * P * exp(-E_ap / T) * phi;
        else L = m_ap * exp(-E_ap / T) * phi;
        L += m_htpb * exp(-E_htpb / T) * (1.0 - phi);
        if (T <= bound) L = 0.0;
        return L;
    }


    int mob_ap = -1;
    Set::Scalar m_ap = NAN, m_htpb = NAN, m_comb = NAN;
    Set::Scalar E_ap = NAN, E_htpb = NAN;
    Set::Scalar bound = NAN;

    static void Parse(Arrhenius& value, IO::ParmParse& pp)
    {
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
};
}
}

#endif
