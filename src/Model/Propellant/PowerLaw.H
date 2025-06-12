//
// This SCP model implements the method described in the following paper
//
// .. bibliography::
//    :list: none
//    :filter: False
// 
//    kanagarajan2022diffuse
//
// This method imitates inheritance :code:`Propellant::Propellant` using static polymorphism
// 

#ifndef MODEL_PROPELLANT_POWERLAW_H
#define MODEL_PROPELLANT_POWERLAW_H

#include "Propellant.H"

namespace Model
{
namespace Propellant
{
class PowerLaw
{
public:
    static constexpr const char* name = "powerlaw";
    
    PowerLaw() = default;
    PowerLaw(IO::ParmParse &pp, std::string name) 
    {
        pp_queryclass(name,*this);
    }


    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    void set_pressure(const Set::Scalar P) 
    {
        // for the pressure power law option
        fmod_ap = r_ap * pow(P, n_ap);
        fmod_htpb = r_htpb * pow(P, n_htpb);
        fmod_comb = r_comb * pow(P, n_comb);

        // a_fit = -1.16582 * sin(P) - 0.681788 * cos(P) + 3.3563;
        // b_fit = -0.708225 * sin(P) + 0.548067 * cos(P) + 1.55985;
        // c_fit = -0.0130849 * sin(P) - 0.03597 * cos(P) + 0.00725694;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar get_K (const Set::Scalar /*phi*/) {return NAN;}
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar get_rho (const Set::Scalar /*phi*/) {return NAN;}
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar get_cp (const Set::Scalar /*phi*/) {return NAN;}
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar get_qdot (const Set::Scalar /*mdot*/, const Set::Scalar /*phi*/) {return NAN;}

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar get_L (Set::Scalar phi, Set::Scalar /*T*/) 
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
            Set::Scalar fs_actual = fmod_ap * phi
                + fmod_htpb * (1.0 - phi)
                + 4.0 * fmod_comb * phi * (1.0 - phi);
            L = fs_actual / gamma / deltaw /*(pf.w1 - pf.w0)*/;
        }

        return L;
    }

    static void Parse(PowerLaw& value, IO::ParmParse& pp)
    {
        Util::Message(INFO,pp.getPrefix());
        // Scaling factor for mobility
        pp_query_required("gamma", value.gamma);   
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
        // jump in chemical potential
        pp_query_default("deltaw", value.deltaw, 1.0); 
    }

    Set::Scalar gamma = NAN;
    Set::Scalar r_ap = NAN, r_htpb = NAN, r_comb = NAN;
    Set::Scalar n_ap = NAN, n_htpb = NAN, n_comb = NAN;
    Set::Scalar fmod_ap = NAN, fmod_htpb = NAN, fmod_comb = NAN;
    Set::Scalar deltaw = NAN;
};
}
}

#endif
