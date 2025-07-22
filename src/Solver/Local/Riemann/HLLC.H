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

class HLLC
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

    virtual Flux Solve(const State &lo, const State &hi,
                       const Set::Scalar gamma,
                       const Set::Scalar p_ref, const Set::Scalar small) const
    {
        // get densities (and regularize)
        Set::Scalar rho_L = std::max(lo.rho, small);
        Set::Scalar rho_R = std::max(hi.rho, small);
        // density average
        Set::Scalar rhobar = 0.5*(rho_L + rho_R);
        // density ratio
        /*Set::Scalar Rrho = std::sqrt(rho_R / rho_L);*/

        // primitive velocities
        Set::Scalar u_L = lo.M_normal / rho_L;
        Set::Scalar u_R = hi.M_normal / rho_R;
        Set::Scalar v_L = lo.M_tangent / rho_L;
        Set::Scalar v_R = hi.M_tangent / rho_R;
        // roe average velocity
        /*Set::Scalar ubar = (u_L + u_R*Rrho) / (1.0 + Rrho);*/

        // label the energies
        Set::Scalar E_L = lo.E, E_R = hi.E;

        // primitive pressures
        Set::Scalar p_L = (gamma-1.0) * (E_L - 0.5*rho_L*(u_L*u_L + v_L*v_L)) + p_ref;
        Set::Scalar p_R = (gamma-1.0) * (E_R - 0.5*rho_R*(u_R*u_R + v_R*v_R)) + p_ref;
        
        // sound speeds
        Set::Scalar a_L = std::sqrt(gamma*p_L / rho_L);
        Set::Scalar a_R = std::sqrt(gamma*p_R / rho_R);
        // sound speed average
        Set::Scalar abar = 0.5*(a_L + a_R);

        Set::Scalar pstar = std::max(0.0, 0.5*(p_L + p_R) - 0.5*(u_R - u_L)*rhobar*abar);
        
        Set::Scalar q_L = pstar <= p_L ?
            1.0 :
            std::sqrt(1.0 + ((gamma+1)/(2*gamma)) * (pstar/p_L - 1.0));
        Set::Scalar q_R = pstar <= p_R ?
            1.0 :
            std::sqrt(1.0 + ((gamma+1)/(2*gamma)) * (pstar/p_R - 1.0));
        
        // wave speed estimate
        Set::Scalar S_L = u_L - a_L*q_L;
        Set::Scalar S_R = u_R + a_R*q_R;

        // intermediate speed
        Set::Scalar Sstar =
            (p_R - p_L + rho_L*u_L*(S_L-u_L) - rho_R*u_R*(S_R-u_R)) /
            (      rho_L*(S_L-u_L)  - rho_R*(S_R - u_R)           );
        

        if (0.0 < S_L) // F_L region
        {
            Flux F_L (  rho_L*u_L, 
                        rho_L*u_L*u_L + p_L,
                        rho_L*u_L*v_L,
                        u_L * (E_L + p_L) );
            return F_L;
        }
        else if (S_L <= 0.0 && 0.0 < Sstar)  // Fstar_L region
        {
            Flux F_L (  rho_L*u_L, 
                        rho_L*u_L*u_L + p_L,
                        rho_L*u_L*v_L,
                        u_L * (E_L + p_L) );
            State U_L ( rho_L,
                        rho_L*u_L,
                        rho_L*v_L,
                        E_L);

            State Ustar_L ( 1.0,
                            Sstar,
                            v_L,
                            (E_L/rho_L) + (Sstar-u_L) * (Sstar + p_L/(rho_L*(S_L-u_L))));
            Ustar_L *= rho_L*(S_L - u_L)/(S_L-Sstar);

            return F_L + S_L*(Ustar_L - U_L);
        }
        else if (Sstar <= 0.0 && 0.0 <= S_R )
        {
            Flux F_R (  rho_R*u_R, 
                        rho_R*u_R*u_R + p_R,
                        rho_R*u_R*v_R,
                        u_R * (E_R + p_R) );
            State U_R ( rho_R,
                        rho_R*u_R,
                        rho_R*v_R,
                        E_R);


            State Ustar_R ( 1.0,
                            Sstar,
                            v_R,
                            (E_R/rho_R) + (Sstar - u_R)*(Sstar + p_R/(rho_R*(S_R-u_R))));
            Ustar_R *= rho_R*(S_R-u_R)/(S_R-Sstar);

            return F_R + S_R*(Ustar_R - U_R);
        }
        else if (S_R < 0.0)
        {
            Flux F_R (  rho_R*u_R, 
                        rho_R*u_R*u_R + p_R,
                        rho_R*u_R*v_R,
                        u_R * (E_R + p_R) );
            return F_R;
        }
        Util::Message(INFO,lo);
        Util::Message(INFO,hi);
        Util::Message(INFO,S_R);
        Util::Message(INFO,Sstar);
        Util::Message(INFO,S_L);
        Util::Abort(INFO);
        return Flux(NAN,NAN,NAN,NAN);
    }
};

} // namespace Riemann
} // namespace Local
} // namespace Solver

#endif
