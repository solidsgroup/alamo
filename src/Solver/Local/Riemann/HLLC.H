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
#include "Model/Gas/Gas.H"
#include <memory>

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

    virtual Flux Solve(const State& lo, const State& hi, Model::Gas::Gas& gas, Set::Patch<const Set::Scalar>& X, int i, int j, int k, int side, Set::Scalar small) override
    {
        // compute densities (and regularize)
        Set::Scalar rho_sum_L = 0.0;
        Set::Scalar rho_sum_R = 0.0;
        for (int n=0; n<gas.nspecies; ++n)
        {
            rho_sum_L += std::max(lo.rho[n], small);
            rho_sum_R += std::max(hi.rho[n], small);
        }
        // density average
        Set::Scalar rhobar = 0.5*(rho_sum_L + rho_sum_R);

        // primitive velocities
        Set::Scalar u_L = lo.M_normal / rho_sum_L;
        Set::Scalar u_R = hi.M_normal / rho_sum_R;
        Set::Scalar v_L = lo.M_tangent / rho_sum_L;
        Set::Scalar v_R = hi.M_tangent / rho_sum_R;

        // label the energies
        Set::Scalar E_L = lo.E, E_R = hi.E;

        // primitive pressures
        Set::Scalar T_L=0, T_R=0, p_L=0, p_R=0, gamma_L=0, gamma_R=0;
        switch (side) {
            case 0: // xlo
                T_L = std::max(small, gas.ComputeT(rho_sum_L, lo.M_normal, lo.M_tangent, E_L, 500.0, X, i-1, j ,k));
                p_L = std::max(small, gas.ComputeP(rho_sum_L, T_L, X, i-1, j, k));
                gamma_L = gas.gamma(T_L, X, i-1, j, k);
                T_R = std::max(small, gas.ComputeT(rho_sum_R, hi.M_normal, hi.M_tangent, E_R, 500.0, X, i, j ,k));
                p_R = std::max(small, gas.ComputeP(rho_sum_R, T_R, X, i, j, k));
                gamma_R = gas.gamma(T_R, X, i, j, k);
                break;
            case 1: // xhi
                T_L = std::max(small, gas.ComputeT(rho_sum_L, lo.M_normal, lo.M_tangent, E_L, 500.0, X, i, j ,k));
                p_L = std::max(small, gas.ComputeP(rho_sum_L, T_L, X, i, j, k));
                gamma_L = gas.gamma(T_L, X, i, j, k);
                T_R = std::max(small, gas.ComputeT(rho_sum_R, hi.M_normal, hi.M_tangent, E_R, 500.0, X, i+1, j ,k));
                p_R = std::max(small, gas.ComputeP(rho_sum_R, T_R, X, i+1, j, k));
                gamma_R = gas.gamma(T_R, X, i+1, j, k);
                break;
            case 2: // ylo
                T_L = std::max(small, gas.ComputeT(rho_sum_L, lo.M_normal, lo.M_tangent, E_L, 500.0, X, i, j-1 ,k));
                p_L = std::max(small, gas.ComputeP(rho_sum_L, T_L, X, i, j-1, k));
                gamma_L = gas.gamma(T_L, X, i, j-1, k);
                T_R = std::max(small, gas.ComputeT(rho_sum_R, hi.M_normal, hi.M_tangent, E_R, 500.0, X, i, j ,k));
                p_R = std::max(small, gas.ComputeP(rho_sum_R, T_R, X, i, j, k));
                gamma_R = gas.gamma(T_R, X, i, j, k);
                break;
            case 3: // yhi
                T_L = std::max(small, gas.ComputeT(rho_sum_L, lo.M_normal, lo.M_tangent, E_L, 500.0, X, i, j ,k));
                p_L = std::max(small, gas.ComputeP(rho_sum_L, T_L, X, i, j, k));
                gamma_L = gas.gamma(T_L, X, i, j, k);
                T_R = std::max(small, gas.ComputeT(rho_sum_R, hi.M_normal, hi.M_tangent, E_R, 500.0, X, i, j+1 ,k));
                p_R = std::max(small, gas.ComputeP(rho_sum_R, T_R, X, i, j+1, k));
                gamma_R = gas.gamma(T_R, X, i, j+1, k);
                break;
            case 4: // zlo
                T_L = std::max(small, gas.ComputeT(rho_sum_L, lo.M_normal, lo.M_tangent, E_L, 500.0, X, i, j ,k-1));
                p_L = std::max(small, gas.ComputeP(rho_sum_L, T_L, X, i, j, k-1));
                gamma_L = gas.gamma(T_L, X, i, j, k-1);
                T_R = std::max(small, gas.ComputeT(rho_sum_R, hi.M_normal, hi.M_tangent, E_R, 500.0, X, i, j ,k));
                p_R = std::max(small, gas.ComputeP(rho_sum_R, T_R, X, i, j, k));
                gamma_R = gas.gamma(T_R, X, i, j, k);
                break;
            case 5: // zhi
                T_L = std::max(small, gas.ComputeT(rho_sum_L, lo.M_normal, lo.M_tangent, E_L, 500.0, X, i, j ,k));
                p_L = std::max(small, gas.ComputeP(rho_sum_L, T_L, X, i, j, k));
                gamma_L = gas.gamma(T_L, X, i, j, k);
                T_R = std::max(small, gas.ComputeT(rho_sum_R, hi.M_normal, hi.M_tangent, E_R, 500.0, X, i, j ,k+1));
                p_R = std::max(small, gas.ComputeP(rho_sum_R, T_R, X, i, j, k+1));
                gamma_R = gas.gamma(T_R, X, i, j, k+1);
                break;
            default:
                Util::Abort(INFO, "Unknown side in Riemann solver, ",side);
        }

        // sound speeds
        Set::Scalar a_L = std::sqrt(gamma_L*p_L / rho_sum_L);
        Set::Scalar a_R = std::sqrt(gamma_R*p_R / rho_sum_R);
        // sound speed average
        Set::Scalar abar = 0.5*(a_L + a_R);

        Set::Scalar pstar = std::max(0.0, 0.5*(p_L + p_R) - 0.5*(u_R - u_L)*rhobar*abar);
        
        Set::Scalar q_L = pstar <= p_L ?
            1.0 :
            std::sqrt(1.0 + ((gamma_L+1.0)/(2.0*gamma_L)) * (pstar/p_L - 1.0));
        Set::Scalar q_R = pstar <= p_R ?
            1.0 :
            std::sqrt(1.0 + ((gamma_R+1.0)/(2.0*gamma_R)) * (pstar/p_R - 1.0));
        
        // wave speed estimate
        Set::Scalar S_L = u_L - a_L*q_L;
        Set::Scalar S_R = u_R + a_R*q_R;

        // intermediate speed
        Set::Scalar Sstar =
            (p_R - p_L + rho_sum_L*u_L*(S_L-u_L) - rho_sum_R*u_R*(S_R-u_R)) /
            (      rho_sum_L*(S_L-u_L)  - rho_sum_R*(S_R - u_R)           );

        Flux result;
        result.nspecies = gas.nspecies;

        if (0.0 < S_L) // F_L region
        {
            for (int n=0; n<gas.nspecies; ++n) result.mass[n] = std::max(lo.rho[n], small) * u_L;
            result.momentum_normal = rho_sum_L*u_L*u_L + p_L;
            result.momentum_tangent = rho_sum_L*u_L*v_L;
            result.energy = u_L * (E_L + p_L);
            return result;
        }
        else if (S_L <= 0.0 && 0.0 < Sstar)  // Fstar_L region
        {
            double factor_L = (S_L - u_L) / (S_L - Sstar);
            for (int n=0; n<gas.nspecies; ++n) {
                Set::Scalar fL_mass = std::max(lo.rho[n], small) * u_L;
                Set::Scalar Ustar_L_rho = std::max(lo.rho[n], small) * factor_L;
                result.mass[n] = fL_mass + S_L * (Ustar_L_rho - lo.rho[n]);
            }
            Set::Scalar fL_Mn = rho_sum_L*u_L*u_L + p_L;
            Set::Scalar Ustar_L_Mn = Sstar * rho_sum_L * factor_L;
            result.momentum_normal = fL_Mn + S_L * (Ustar_L_Mn - rho_sum_L*u_L);

            Set::Scalar fL_Mt = rho_sum_L*u_L*v_L;
            Set::Scalar Ustar_L_Mt = v_L * rho_sum_L * factor_L;
            result.momentum_tangent = fL_Mt + S_L * (Ustar_L_Mt - rho_sum_L*v_L);

            Set::Scalar fL_E = u_L * (E_L + p_L);
            Set::Scalar Ustar_L_E = ((E_L/rho_sum_L) + (Sstar-u_L) * (Sstar + p_L/(rho_sum_L*(S_L-u_L)))) * rho_sum_L * factor_L;
            result.energy = fL_E + S_L * (Ustar_L_E - E_L);

            return result;
        }
        else if (Sstar <= 0.0 && 0.0 <= S_R )
        {
            double factor_R = (S_R - u_R) / (S_R - Sstar);
            for (int n=0; n<gas.nspecies; ++n) {
                Set::Scalar fR_mass = std::max(hi.rho[n], small) * u_R;
                Set::Scalar Ustar_R_rho = std::max(hi.rho[n], small) * factor_R;
                result.mass[n] = fR_mass + S_R * (Ustar_R_rho - hi.rho[n]);
            }
            Set::Scalar fR_Mn = rho_sum_R*u_R*u_R + p_R;
            Set::Scalar Ustar_R_Mn = Sstar * rho_sum_R * factor_R;
            result.momentum_normal = fR_Mn + S_R * (Ustar_R_Mn - rho_sum_R*u_R);

            Set::Scalar fR_Mt = rho_sum_R*u_R*v_R;
            Set::Scalar Ustar_R_Mt = v_R * rho_sum_R * factor_R;
            result.momentum_tangent = fR_Mt + S_R * (Ustar_R_Mt - rho_sum_R*v_R);

            Set::Scalar fR_E = u_R * (E_R + p_R);
            Set::Scalar Ustar_R_E = ((E_R/rho_sum_R) + (Sstar - u_R)*(Sstar + p_R/(rho_sum_R*(S_R-u_R)))) * rho_sum_R * factor_R;
            result.energy = fR_E + S_R * (Ustar_R_E - E_R);

            return result;
        }
        else if (S_R < 0.0)
        {
            for (int n=0; n<gas.nspecies; ++n) result.mass[n] = std::max(hi.rho[n], small) * u_R;
            result.momentum_normal = rho_sum_R*u_R*u_R + p_R;
            result.momentum_tangent = rho_sum_R*u_R*v_R;
            result.energy = u_R * (E_R + p_R);
            return result;
        }
        Util::Message(INFO,lo);
        Util::Message(INFO,hi);
        Util::Message(INFO,S_R);
        Util::Message(INFO,Sstar);
        Util::Message(INFO,S_L);
        Util::Exception(INFO);
        std::array<double, NSPECIES_MAX> nan_mass; nan_mass.fill(NAN);
        return Flux(nan_mass, gas.nspecies, NAN, NAN, NAN);
    }
};

} // namespace Riemann
} // namespace Local
} // namespace Solver

#endif
