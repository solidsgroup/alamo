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

class HLLE
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

    Flux Solve(State lo, State hi, Model::Gas::Gas& gas, Set::Patch<const Set::Scalar>& X, int i, int j, int k, int side, Set::Scalar small) const
    {
        // Grab and label the input states
        Set::Scalar rho_L = lo.rho       ,  rho_R = hi.rho;
        Set::Scalar Mn_L  = lo.M_normal  ,  Mn_R  = hi.M_normal  ;
        //Set::Scalar Mt_L  = lo.M_tangent ,  Mt_R  = hi.M_tangent ;
        Set::Scalar E_L   = lo.E         ,  E_R   = hi.E         ;

        // Ensure no negative densities
        rho_L = std::max(small,rho_L);
        rho_R = std::max(small,rho_R);


        // Calculate primitives
        Set::Scalar u_L   = Mn_L/(rho_L+small),  u_R   = Mn_R/(rho_R+small);

        Set::Scalar T_L=0, T_R=0, p_L=0, p_R=0, gamma_L=0, gamma_R=0;
        switch (side) {
            case 0: // xlo
                T_L = gas.ComputeT(rho_L, lo.M_normal, lo.M_tangent, E_L, 500.0, X, i-1, j ,k);
                p_L = gas.ComputeP(rho_L, T_L, X, i-1, j, k);
                gamma_L = gas.gamma(T_L, X, i-1, j, k);
                T_R = gas.ComputeT(rho_R, hi.M_normal, hi.M_tangent, E_R, 500.0, X, i, j ,k);
                p_R = gas.ComputeP(rho_R, T_R, X, i, j, k);
                gamma_R = gas.gamma(T_R, X, i, j, k);
                break;
            case 1: // xhi
                T_L = gas.ComputeT(rho_L, lo.M_normal, lo.M_tangent, E_L, 500.0, X, i, j ,k);
                p_L = gas.ComputeP(rho_L, T_L, X, i, j, k);
                gamma_L = gas.gamma(T_L, X, i, j, k);
                T_R = gas.ComputeT(rho_R, hi.M_normal, hi.M_tangent, E_R, 500.0, X, i+1, j ,k);
                p_R = gas.ComputeP(rho_R, T_R, X, i+1, j, k);
                gamma_R = gas.gamma(T_R, X, i+1, j, k);
                break;
            case 2: // ylo
                T_L = gas.ComputeT(rho_L, lo.M_normal, lo.M_tangent, E_L, 500.0, X, i, j-1 ,k);
                p_L = gas.ComputeP(rho_L, T_L, X, i, j-1, k);
                gamma_L = gas.gamma(T_L, X, i, j-1, k);
                T_R = gas.ComputeT(rho_R, hi.M_normal, hi.M_tangent, E_R, 500.0, X, i, j ,k);
                p_R = gas.ComputeP(rho_R, T_R, X, i, j, k);
                gamma_R = gas.gamma(T_R, X, i, j, k);
                break;
            case 3: // yhi
                T_L = gas.ComputeT(rho_L, lo.M_normal, lo.M_tangent, E_L, 500.0, X, i, j ,k);
                p_L = gas.ComputeP(rho_L, T_L, X, i, j, k);
                gamma_L = gas.gamma(T_L, X, i, j, k);
                T_R = gas.ComputeT(rho_R, hi.M_normal, hi.M_tangent, E_R, 500.0, X, i, j+1 ,k);
                p_R = gas.ComputeP(rho_R, T_R, X, i, j+1, k);
                gamma_R = gas.gamma(T_R, X, i, j+1, k);
                break;
            case 4: // zlo
                T_L = gas.ComputeT(rho_L, lo.M_normal, lo.M_tangent, E_L, 500.0, X, i, j ,k-1);
                p_L = gas.ComputeP(rho_L, T_L, X, i, j, k-1);
                gamma_L = gas.gamma(T_L, X, i, j, k-1);
                T_R = gas.ComputeT(rho_R, hi.M_normal, hi.M_tangent, E_R, 500.0, X, i, j ,k);
                p_R = gas.ComputeP(rho_R, T_R, X, i, j, k);
                gamma_R = gas.gamma(T_R, X, i, j, k);
                break;
            case 5: // zhi
                T_L = gas.ComputeT(rho_L, lo.M_normal, lo.M_tangent, E_L, 500.0, X, i, j ,k);
                p_L = gas.ComputeP(rho_L, T_L, X, i, j, k);
                gamma_L = gas.gamma(T_L, X, i, j, k);
                T_R = gas.ComputeT(rho_R, hi.M_normal, hi.M_tangent, E_R, 500.0, X, i, j ,k+1);
                p_R = gas.ComputeP(rho_R, T_R, X, i, j, k+1);
                gamma_R = gas.gamma(T_R, X, i, j, k+1);
                break;
            default:
                Util::Abort(INFO, "Unknown side in Riemann solver, ",side);
        }

        Set::Scalar c_L = std::sqrt(gamma_L * p_L / (rho_L + small));
        Set::Scalar c_R = std::sqrt(gamma_R * p_R / (rho_R + small));

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

        Flux F_L (
            lo.M_normal,
            lo.M_normal * u_L + p_L,
            lo.M_tangent * u_L,
            u_L * (lo.E + p_L)
            );
        Flux F_R (
            hi.M_normal,
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
