//
// HLLC flux calculation based on equations 10.67 - 10.73 of
// _Riemann solvers and numerical methods for fluid dynamics_ by E.F.Toro.
//
//

#ifndef SOLVER_LOCAL_RIEMANN_HLLC_H
#define SOLVER_LOCAL_RIEMANN_HLLC_H

#include "IO/ParmParse.H"
#include "Solver/Local/Riemann/Riemann.H"
#include "Util/Util.H"

namespace Solver {
namespace Local {
namespace Riemann {

class HLLC : public Riemann
{
public:
    static constexpr const char* name = "hllc";
    HLLC(IO::ParmParse &pp, std::string name) { pp_queryclass(name, *this); }
    HLLC(IO::ParmParse &pp) { pp_queryclass(*this); }
    HLLC() { IO::ParmParse pp; pp_queryclass(*this); }

    static void Parse(HLLC & /*value*/, IO::ParmParse & /*pp*/)
    {
        // pp.query_default("lowmach", value.lowmach, false);
        // pp.query_default("cutoff", value.cutoff, 0.1);
    }

    virtual Flux Solve(State lo, State hi, Set::Scalar p_ref, Set::Scalar small, std::vector<double> &species_mw, std::vector<double> &species_cp) override
    {
        // get densities (and regularize)
        size_t nspecies = lo.rho.size();
        std::vector<double> rho_L = lo.rho;
        std::vector<double> rho_R = hi.rho;
        double rhobar = 0.0;
        double rho_sum_L = 0.0;
        double rho_sum_R = 0.0;
        for (size_t n=0; n<nspecies; ++n) {
            rho_L[n] = std::max(rho_L[n], small);
            rho_R[n] = std::max(rho_R[n], small);
            rho_sum_L += rho_L[n];
            rho_sum_R += rho_R[n];
        }
        // density average
        rhobar = 0.5*(rho_sum_L + rho_sum_R);
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
        std::cout << "gamma: " << gamma_L << " " << gamma_R << std::endl;
            // density ratio
        /*Set::Scalar Rrho = std::sqrt(rho_R / rho_L);*/

        // primitive velocities
        Set::Scalar u_L = lo.M_normal / rho_sum_L;
        Set::Scalar u_R = hi.M_normal / rho_sum_R;
        Set::Scalar v_L = lo.M_tangent / rho_sum_L;
        Set::Scalar v_R = hi.M_tangent / rho_sum_R;
        // roe average velocity
        /*Set::Scalar ubar = (u_L + u_R*Rrho) / (1.0 + Rrho);*/

        // label the energies
        Set::Scalar E_L = lo.E, E_R = hi.E;

        // primitive pressures
        Set::Scalar p_L = (gamma_L-1.0) * (E_L - 0.5*rho_sum_L*(u_L*u_L + v_L*v_L)) + p_ref;
        Set::Scalar p_R = (gamma_R-1.0) * (E_R - 0.5*rho_sum_R*(u_R*u_R + v_R*v_R)) + p_ref;
        
        // sound speeds
        Set::Scalar a_L = std::sqrt(gamma_L*p_L / rho_sum_L);
        Set::Scalar a_R = std::sqrt(gamma_R*p_R / rho_sum_R);
        // sound speed average
        Set::Scalar abar = 0.5*(a_L + a_R);

        Set::Scalar pstar = std::max(0.0, 0.5*(p_L + p_R) - 0.5*(u_R - u_L)*rhobar*abar);
        
        Set::Scalar q_L = pstar <= p_L ?
            1.0 :
            std::sqrt(1.0 + ((gamma_L+1)/(2*gamma_L)) * (pstar/p_L - 1.0));
        Set::Scalar q_R = pstar <= p_R ?
            1.0 :
            std::sqrt(1.0 + ((gamma_R+1)/(2*gamma_R)) * (pstar/p_R - 1.0));
        
        // wave speed estimate
        Set::Scalar S_L = u_L - a_L*q_L;
        Set::Scalar S_R = u_R + a_R*q_R;

        // intermediate speed
        Set::Scalar Sstar =
            (p_R - p_L + rho_sum_L*u_L*(S_L-u_L) - rho_sum_R*u_R*(S_R-u_R)) /
            (      rho_sum_L*(S_L-u_L)  - rho_sum_R*(S_R - u_R)    );

        //if ( Sstar != Sstar ) Sstar = 0.5 * (u_L + u_R);
        
        std::vector<double> rhoYu_L(nspecies);
        std::vector<double> rhoYu_R(nspecies);
        std::vector<double> rhostar_L(nspecies);
        std::vector<double> rhostar_R(nspecies);
        double factor_L = (S_L - u_L) / (S_L - Sstar);
        double factor_R = (S_R - u_R) / (S_R - Sstar);
        for (size_t n=0; n<nspecies; ++n) {
            rhoYu_L[n] = u_L*rho_L[n];
            rhoYu_R[n] = u_R*rho_R[n];
            rhostar_L[n] = rho_L[n] * factor_L;
            rhostar_R[n] = rho_R[n] * factor_R;
        }

        if (0.0 <= S_L) // F_L region
        {
            Flux F_L (  rhoYu_L,
                        rho_sum_L*u_L*u_L + p_L,
                        rho_sum_L*u_L*v_L,
                        u_L * (E_L + p_L) );
            return F_L;
        }
        else if (S_L <= 0.0 && 0.0 <= Sstar)  // Fstar_L region
        {
            Flux F_L (  rhoYu_L, 
                        rho_sum_L*u_L*u_L + p_L,
                        rho_sum_L*u_L*v_L,
                        u_L * (E_L + p_L) );
            State U_L ( rho_L,
                        rho_sum_L*u_L,
                        rho_sum_L*v_L,
                        E_L);

            State Ustar_L ( rhostar_L,
                            rho_sum_L * factor_L * Sstar,
                            rho_sum_L * factor_L * v_L,
                            ((S_L - u_L)*E_L - p_L*u_L + pstar*Sstar + (pstar - p_L)*S_L)/(S_L - Sstar));

            return F_L + S_L*(Ustar_L - U_L);
        }
        else if (Sstar <= 0.0 && 0.0 <= S_R )
        {
            Flux F_R (  rhoYu_R, 
                        rho_sum_R*u_R*u_R + p_R,
                        rho_sum_R*u_R*v_R,
                        u_R * (E_R + p_R) );
            State U_R ( rho_R,
                        rho_sum_R*u_R,
                        rho_sum_R*v_R,
                        E_R);


            State Ustar_R ( rhostar_R,
                            rho_sum_R * factor_R * Sstar,
                            rho_sum_R * factor_R * v_R,
                            ((S_R - u_R)*E_R - p_R*u_R + pstar*Sstar + (pstar - p_R)*S_R)/(S_R - Sstar));

            return F_R + S_R*(Ustar_R - U_R);
        }
        else if (S_R <= 0.0)
        {
            Flux F_R (  rhoYu_R, 
                        rho_sum_R*u_R*u_R + p_R,
                        rho_sum_R*u_R*v_R,
                        u_R * (E_R + p_R) );
            return F_R;
        }
        Util::Message(INFO,lo);
        Util::Message(INFO,hi);
        Util::Message(INFO,S_R);
        Util::Message(INFO,Sstar);
        Util::Message(INFO,S_L);
        Util::Exception(INFO);
        return Flux(std::vector<double>(nspecies,NAN),NAN,NAN,NAN);
    }
};

} // namespace Riemann
} // namespace Local
} // namespace Solver

#endif
