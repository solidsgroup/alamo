#ifndef SOLVER_LOCAL_RIEMANN_HLLE_H
#define SOLVER_LOCAL_RIEMANN_HLLE_H

#include "IO/ParmParse.H"
#include "Solver/Local/Riemann/Riemann.H"

/// A bunch of solvers
namespace Solver
{
/// Local solvers
namespace Local
{

namespace Riemann
{

class HLLE : public Riemann
{
public:


    static constexpr const char* name = "hlle";
    HLLE (IO::ParmParse &pp, std::string name) 
    {pp_queryclass(name,*this);}
    HLLE (IO::ParmParse &pp) 
    {pp_queryclass(*this);}
    HLLE () 
    {
        IO::ParmParse pp;
        pp_queryclass(*this);
    }

    int lowmach = 0;
    Set::Scalar cutoff = NAN;

    static void Parse(HLLE & value, IO::ParmParse & pp)
    {
        // toggle to apply a low mach correction approximation
        pp.query_default("lowmach",value.lowmach,false);
        // cutoff value if using the low mach approximation
        pp.query_default("cutoff",value.cutoff,0.1);
    }

    virtual Flux Solve(State lo, State hi, Set::Scalar p_ref, Set::Scalar small, std::vector<double> &species_mw, std::vector<double> &species_cp) override
    {
        size_t nspecies = lo.rho.size();
        // Grab and label the input states
        std::vector<double> rho_L = lo.rho       ,  rho_R = hi.rho;
        Set::Scalar Mn_L          = lo.M_normal  ,  Mn_R  = hi.M_normal  ;
        //Set::Scalar Mt_L          = lo.M_tangent ,  Mt_R  = hi.M_tangent ;
        Set::Scalar E_L           = lo.E         ,  E_R   = hi.E         ;

        // Ensure no negative densities
        double rho_sum_L = 0.0;
        double rho_sum_R = 0.0;
        for (size_t n=0; n<nspecies; ++n) {
            rho_L[n] = std::max(0.0,rho_L[n]);
            rho_R[n] = std::max(0.0,rho_R[n]);
            rho_sum_L += rho_L[n];
            rho_sum_R += rho_R[n];
        }
        double cp_L = 0.0;
        double cv_L = 0.0;
        double gamma_L = 0.0;
        double cp_R = 0.0;
        double cv_R = 0.0;
        double gamma_R = 0.0;
        for (size_t n=0; n<nspecies; ++n) {
            cp_L += rho_L[n]/rho_sum_L * species_cp[n];
            cv_L += rho_L[n]/rho_sum_L * (species_cp[n] - 8314.46/species_mw[n]);
            cp_R += rho_R[n]/rho_sum_R * species_cp[n];
            cv_R += rho_R[n]/rho_sum_R * (species_cp[n] - 8314.46/species_mw[n]);
        }
        gamma_L = cp_L/cv_L;
        gamma_R = cp_R/cv_R;

        // Calculate primitives
        Set::Scalar u_L   = Mn_L/(rho_sum_L+small),  u_R   = Mn_R/(rho_sum_R+small);
        Set::Scalar v_L   = Mn_L/(rho_sum_L+small),  v_R   = Mn_R/(rho_sum_R+small);
        Set::Scalar p_L  = (gamma_L - 1.0) * ( E_L - 0.5*rho_sum_L*(u_L*u_L + v_L*v_L) ) + p_ref;
        Set::Scalar p_R  = (gamma_R - 1.0) * ( E_R - 0.5*rho_sum_R*(u_R*u_R + v_R*v_R) ) + p_ref;
        Set::Scalar c_L = std::sqrt(gamma_L * p_L / (rho_sum_L + small));
        Set::Scalar c_R = std::sqrt(gamma_R * p_R / (rho_sum_R + small));

        // Max wave speeds
        Set::Scalar S_L = std::min(u_L-c_L, u_R - c_R);
        Set::Scalar S_R = std::max(u_L+c_L, u_R + c_R);

        
        if (lowmach)
        {
            Set::Scalar u_RL = 0.5 * (u_L + u_R);
            Set::Scalar a_RL = 0.5 * (c_L + c_R);
            Set::Scalar Mach = std::abs(u_RL) / (a_RL + small);
            Set::Scalar psi = std::max(Mach, cutoff); // Mach cutoff
            S_L *= psi;
            S_R *= psi;
        }

        std::vector<double> rhoYu_L(nspecies);
        std::vector<double> rhoYu_R(nspecies);
        for (size_t n=0; n<nspecies; ++n) {
            rhoYu_L[n] = u_L*lo.rho[n];
            rhoYu_R[n] = u_R*hi.rho[n];
        }

        Flux F_L (
            rhoYu_L,
            lo.M_normal * u_L + p_L,
            lo.M_tangent * u_L,
            u_L * (lo.E + p_L)
            );
        Flux F_R (
            rhoYu_R,
            hi.M_normal * u_R + p_R,
            hi.M_tangent * u_R,
            u_R * (hi.E + p_R)
            );

        if (S_L >= 0.0)
        {
            return F_L;

        }
        else if (S_R <= 0.0)
        {
            return F_R;
        }
        else
        {
            State dState = hi - lo;
            Set::Scalar dS = S_R - S_L;
            return (S_R*F_L - S_L*F_R + S_L*S_R*Flux(dState)) / dS;
        }
    }
};
}
}
}

#endif
