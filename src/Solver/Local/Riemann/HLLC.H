#include <eigen3/Eigen/Core>
#include "Set/Set.H"

/// A bunch of solvers
namespace Solver
{
/// Local solvers
namespace Local
{

namespace Riemann
{

/// HLLC Riemann Solver
class HLLC
{
public:
    struct State {
        Set::Scalar rho;
        //Set::Vector momentum;
        Set::Scalar momentum_normal;
        Set::Scalar momentum_tangent;
        Set::Scalar energy;
        Set::Scalar eta;
        // Construtor for convenience
        State()
        {
            rho = 0.0;
            momentum_normal = 0.0;
            momentum_tangent = 0.0;
            energy = 0.0;
            eta = 0.0;
        }
        State(Set::Scalar a_rho, Set::Scalar a_m_normal, Set::Scalar a_m_tangent, Set::Scalar a_energy, Set::Scalar a_eta)
            : rho(a_rho),
            momentum_normal(a_m_normal), momentum_tangent(a_m_tangent),
            energy(a_energy), eta(a_eta) {}
    };

    struct Flux {
        Set::Scalar mass;
        Set::Scalar energy;
        Set::Scalar momentum_normal;
        Set::Scalar momentum_tangent;
        //Set::Vector momentum;
    };


    // todo: convert to static functor by replacing with static () when c++23 is released
    static Flux Solve(State lo, State hi, Set::Scalar gamma) {

        //lo primitive variables
        Set::Scalar uNorm_lo = lo.momentum_normal / lo.rho;
        Set::Scalar uTang_lo = lo.momentum_tangent / lo.rho;
        Set::Scalar ke_lo = 0.5 * lo.rho * (uNorm_lo * uNorm_lo + uTang_lo * uTang_lo);
        Set::Scalar ue_lo = lo.energy - ke_lo;
        Set::Scalar p_lo = (gamma - 1) * ue_lo;

        //hi primitive variables
        Set::Scalar uNorm_hi = hi.momentum_normal / hi.rho;
        Set::Scalar uTang_hi = hi.momentum_tangent / hi.rho;
        Set::Scalar ke_hi = 0.5 * hi.rho * (uNorm_hi * uNorm_hi + uTang_hi * uTang_hi);
        Set::Scalar ue_hi = hi.energy - ke_hi;
        Set::Scalar p_hi = (gamma - 1) * ue_hi;

        //Find the largest eigenvalues in the normal direction to the interface
        Set::Scalar a_lo = std::sqrt(gamma * p_lo / lo.rho);
        Set::Scalar a_hi = std::sqrt(gamma * p_hi / hi.rho);

        //Compute HLL wave speed
        Set::Scalar S_min = std::min(uNorm_lo, uNorm_hi) - max(a_lo, a_hi);
        Set::Scalar S_max = std::max(uNorm_lo, uNorm_hi) + max(a_lo, a_hi);

        //Compute lagrangian sound speed
        Set::Scalar aLag_min = lo.rho * (uNorm_lo - S_min);
        Set::Scalar aLag_max = hi.rho * (S_max - uNorm_hi);

        //Compute acoustic star state
        Set::Scalar ustar = (aLag_max * uNorm_hi + aLag_min * uNorm_lo + (p_lo - p_hi)) / (aLag_max + aLag_min);
        Set::Scalar pstar = (aLag_max * p_lo + aLag_min * p_hi + aLag_min * aLag_max * (uNorm_lo - uNorm_hi)) / (aLag_max + aLag_min);

        //Left star region variables
        Set::Scalar rstar_lo = lo.rho * (S_min - uNorm_lo) / (S_min - ustar);
        Set::Scalar uestar_lo = ((S_min - uNorm_lo) * ue_lo - p_lo * uNorm_lo + pstar * ustar) / (S_min - ustar);
        Set::Scalar kestar_lo = ((S_min - uNorm_lo) * ke_lo - p_lo * uNorm_lo + pstar * ustar) / (S_min - ustar);
        Set::Scalar etastar_lo = lo.eta * (S_min - uNorm_lo) / (S_min - ustar);

        //Right star region variables
        Set::Scalar rstar_hi = hi.rho * (S_max - uNorm_hi) / (S_max - ustar);
        Set::Scalar uestar_hi = ((S_max - uNorm_hi) * ue_hi - p_hi * uNorm_hi + pstar * ustar) / (S_max - ustar);
        Set::Scalar kestar_hi = ((S_max - uNorm_hi) * ke_hi - p_hi * uNorm_hi + pstar * ustar) / (S_max - ustar);
        Set::Scalar etastar_hi = hi.eta * (S_max - uNorm_hi) / (S_max - ustar);

        //Sample the solution at x/t=0
        Set::Scalar rho_face, u_face, p_face, ke_face, ue_face, eta_face;
        if (S_min > 0.0) {
            rho_face = lo.rho;
            u_face = uNorm_lo;
            p_face = p_lo;
            ue_face = ue_lo;
            ke_face = ke_lo;
            eta_face = lo.eta;
        }
        else if (ustar > 0.0) {
            rho_face = rstar_lo;
            u_face = ustar;
            p_face = pstar;
            ue_face = uestar_lo;
            ke_face = kestar_lo;
            eta_face = etastar_lo;
        }
        else if (S_max > 0.0) {
            rho_face = rstar_hi;
            u_face = ustar;
            p_face = pstar;
            ue_face = uestar_hi;
            ke_face = kestar_hi;
            eta_face = etastar_hi;
        }
        else {
            rho_face = hi.rho;
            u_face = uNorm_hi;
            p_face = p_hi;
            ue_face = ue_hi;
            ke_face = ke_hi;
            eta_face = hi.eta;
        };

        //Compute the Godunov flux
        Flux fl;
        fl.mass = rho_face * u_face;
        fl.momentum_normal = (rho_face * u_face * u_face) * eta_face + p_face;
        if (fl.mass > 0.0) {
            fl.momentum_tangent = (fl.mass * uTang_lo) * u_face * eta_face;
        }
        else {
            fl.momentum_tangent = (fl.mass * uTang_hi) * u_face * eta_face;
        };
        fl.energy = (ke_face + ue_face) * u_face * eta_face;

        return fl;

    }



    static int Test()
    {
        int failed = 0;
        for (int i = 0; i < 10; i++)
        {

            State hi, lo, fixed, pos_vel1, pos_vel2, neg_vel1, neg_vel2;
            hi.rho = Util::Random();
            hi.momentum_normal = Util::Random();
            hi.momentum_tangent = Util::Random();
            hi.energy = Util::Random();
            hi.eta = Util::Random();

            lo.rho = Util::Random();
            lo.momentum_normal = Util::Random();
            lo.momentum_tangent = Util::Random();
            lo.energy = Util::Random();
            lo.eta = Util::Random();

            lo.rho = Util::Random();
            lo.momentum_normal = Util::Random();
            lo.momentum_tangent = Util::Random();
            lo.energy = Util::Random();
            lo.eta = Util::Random();

            pos_vel1.rho = Util::Random();
            pos_vel1.momentum_normal = Util::Random();
            pos_vel1.momentum_tangent = Util::Random();
            pos_vel1.energy = Util::Random();
            pos_vel1.eta = Util::Random();

            neg_vel1.rho = pos_vel1.rho;
            neg_vel1.momentum_normal = -pos_vel1.momentum_normal;
            neg_vel1.momentum_tangent = -pos_vel1.momentum_tangent;
            neg_vel1.energy = pos_vel1.energy;
            neg_vel1.eta = pos_vel1.eta;

            pos_vel2.rho = Util::Random();
            pos_vel2.momentum_normal = Util::Random();
            pos_vel2.momentum_tangent = Util::Random();
            pos_vel2.energy = Util::Random();
            pos_vel2.eta = Util::Random();

            neg_vel2.rho = pos_vel2.rho;
            neg_vel2.momentum_normal = -pos_vel2.momentum_normal;
            neg_vel2.momentum_tangent = -pos_vel2.momentum_tangent;
            neg_vel2.energy = pos_vel2.energy;
            neg_vel2.eta = pos_vel2.eta;

            // 
            // Zero flux test: if the state is the same on both sides,
            //                 the flux should be zero.
            // 
            Flux zero = Solve(fixed, fixed, 1.4);
            if (fabs(zero.mass) > 1E-10) { Util::Warning(INFO, "Nonzero mass flux ", zero.mass); failed++; }
            if (fabs(zero.momentum_normal) > 1E-10) { Util::Warning(INFO, "Nonzero x momentum flux ", zero.momentum_normal); failed++; }
            if (fabs(zero.momentum_tangent) > 1E-10) { Util::Warning(INFO, "Nonzero y momentum flux ", zero.momentum_tangent); failed++; }
            if (fabs(zero.energy) > 1E-10) { Util::Warning(INFO, "Nonzero energy flux ", zero.energy); failed++; }

            // 
            // Antisymmetry test: reversing hi and lo should change the sign
            // of the flux.
            // 

            Flux plus = Solve(hi, lo, 1.4);
            Flux minus = Solve(lo, hi, 1.4);

            if (fabs(plus.mass + minus.mass) > 1E-10)
            {                
                Util::Warning(INFO, "mass flux error: plus=", plus.mass, " minus=", minus.mass); failed++;            
            }
            if (fabs(plus.momentum_normal + minus.momentum_normal) > 1E-10)
            {
                Util::Warning(INFO, "momentum x flux error: plus=", plus.momentum_normal, " minus=", minus.momentum_normal); failed++;
            }
            if (fabs(plus.momentum_tangent + minus.momentum_tangent) > 1E-10)
            {
                Util::Warning(INFO, "momentum y flux error: plus=", plus.momentum_tangent, " minus=", minus.momentum_tangent); failed++;
            }
            if (fabs(plus.energy + minus.energy) > 1E-10)
            {
                Util::Warning(INFO, "energy flux error: plus=", plus.energy, " minus=", minus.energy); failed++;
            }


            // 
            // Wave Propogation Test: changing the sign of velocity should produce negative mass
            //                   and energy fluxes and the same momentum flux
            // 

            Flux pos = Solve(pos_vel1, pos_vel2, 1.4);
            Flux neg = Solve(neg_vel1, neg_vel2, 1.4);

            if (fabs(pos.mass + neg.mass) > 1E-10)
            {
                Util::Warning(INFO, "mass flux error: positive vel=", pos.mass, " negative vel=", neg.mass); failed++;
            }
            if (fabs(pos.momentum_normal - neg.momentum_normal) > 1E-10)
            {
                Util::Warning(INFO, "momentum x flux error: positive vel=", pos.momentum_normal, " negative vel=", neg.momentum_normal); failed++;            
            }
            if (fabs(pos.momentum_tangent - neg.momentum_tangent) > 1E-10)
            {
                Util::Warning(INFO, "momentum y flux error: positive vel=", pos.momentum_tangent, " negative vel=", neg.momentum_tangent); failed++;
            }
            if (fabs(pos.energy + neg.energy) > 1E-10)
            {
                Util::Warning(INFO, "energy flux error: positive vel=", pos.energy, " negative vel=", neg.energy); failed++;
            }
        }

        return failed;
    }

};
}
}
}
