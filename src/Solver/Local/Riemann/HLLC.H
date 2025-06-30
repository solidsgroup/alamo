#ifndef SOLVER_LOCAL_RIEMANN_HLLC_H
#define SOLVER_LOCAL_RIEMANN_HLLC_H

#include "IO/ParmParse.H"
#include "Solver/Local/Riemann/Riemann.H"

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

    // int lowmach = 0;
    // Set::Scalar cutoff = NAN;

    static void Parse(HLLC & value, IO::ParmParse & pp)
    {
        // pp.query_default("lowmach", value.lowmach, false);
        // pp.query_default("cutoff", value.cutoff, 0.1);
    }

    virtual Flux Solve(State lo, State hi, Set::Scalar gamma, Set::Scalar p_ref, Set::Scalar small) override
    {

        //
        // Compute left/right states and Roe averages
        // 

        // get densities (and regularize)
        Set::Scalar rho_L = std::max(lo.rho, small);
        Set::Scalar rho_R = std::max(hi.rho, small);
        // density ratio
        Set::Scalar Rrho = std::sqrt(rho_R / rho_L);

        // primitive velocities
        Set::Scalar u_L = lo.M_normal / rho_L;
        Set::Scalar u_R = hi.M_normal / rho_R;
        Set::Scalar v_L = lo.M_tangent / rho_L;
        Set::Scalar v_R = hi.M_tangent / rho_R;
        // roe average velocity
        Set::Scalar ubar = (u_L + u_R*Rrho) / (1 + Rrho);

        // label the energies
        Set::Scalar E_L = lo.E, E_R = hi.E;

        // primitive pressures
        Set::Scalar p_L = (gamma-1.0) * (E_L - 0.5*rho_L*u_L*u_L);
        Set::Scalar p_R = (gamma-1.0) * (E_R - 0.5*rho_R*u_R*u_R);
        
        // enthalpiess
        Set::Scalar H_R = (E_R + p_R) / rho_R;
        Set::Scalar H_L = (E_L + p_L) / rho_L;
        // roe average enthalpy
        Set::Scalar Hbar = (H_L + H_R*Rrho) / (1.0 + Rrho);

        // sound speeds
        Set::Scalar c_L = std::sqrt(gamma*p_L / rho_L);
        Set::Scalar c_R = std::sqrt(gamma*p_R / rho_R);
        // roe average sound speed
        Set::Scalar cbar = std::sqrt((gamma-1.0) * (Hbar - 0.5*ubar*ubar));

        

        //
        // Compute wave speeds
        //

        Set::Scalar S_L = std::min(u_L - c_L, ubar - cbar);
        Set::Scalar S_R = std::max(u_R + c_R, ubar + cbar);
        
        Set::Scalar S_M = 
            ( rho_R*u_R*(S_R - u_R) - rho_L*u_L*(S_L - u_L) + p_L - p_R ) / 
            (         rho_R * (S_R - u_R) - rho_L*(S_L - u_L)           ) ;

        
        
        //
        // Determine our state
        //

        if (S_L > 0.0)
        {
            Flux F_L ( rho_L*v_L, 
                       rho_L*v_L*v_L + p_L,
                       0.5 * (rho_L * u_L * v_L + rho_R * u_R * v_R),
                       v_L * (E_L + p_L) );
            return F_L;
        }
        if (S_L < 0.0 && S_M > 0.0)
        {
            Flux F_L ( rho_L*v_L, 
                       rho_L*v_L*v_L + p_L,
                       0.5 * (rho_L * u_L * v_L + rho_R * u_R * v_R),
                       v_L * (E_L + p_L) );

            


        }


        else if (S_R < 0.0)
        {
            Flux F_R ( rho_R*v_R, 
                       rho_R*v_R*v_R + p_L,
                       0.5 * (rho_L * u_L * v_L + rho_R * u_R * v_R),
                       v_R * (E_R + p_R) );
            return F_R;
        }
        
        

    }
};

} // namespace Riemann
} // namespace Local
} // namespace Solver

#endif
