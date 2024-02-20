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

/// Roe Riemann Solver based on Gas Dynamics - Culbert B. Laney
class Roe
{
public:
    struct State {
        Set::Scalar rho;
        Set::Scalar velocity_normal;
        Set::Scalar velocity_tangent;
        Set::Scalar pressure;
        Set::Scalar eta;
        // Construtor for convenience
        State() { rho = 0.0; velocity_normal = 0.0, velocity_tangent = 0.0; pressure = 0.0; eta = 0.0; }
        State(Set::Scalar a_rho, Set::Scalar a_velocity_normal, Set::Scalar a_velocity_tangent, Set::Scalar a_pressure, Set::Scalar a_eta)
            : rho(a_rho), velocity_normal(a_velocity_normal), velocity_tangent(a_velocity_tangent), pressure(a_pressure), eta(a_eta) {}
    };

    struct Flux {
        Set::Scalar mass;
        Set::Scalar energy;
        Set::Scalar momentum_normal;
        Set::Scalar momentum_tangent;
    };

    // todo: convert to static functor by replacing with static () when c++23 is released
    static Flux Solve(State lo, State hi, Set::Scalar gamma, Set::Scalar cell_eta) {
        Set::Scalar ke_lo = 0.5 * lo.rho * lo.velocity_normal * lo.velocity_normal;
        Set::Scalar ue_lo = lo.pressure / (gamma - 1);
        Set::Scalar h_lo = (ke_lo + ue_lo + lo.pressure) / lo.rho;

        Set::Scalar ke_hi = 0.5 * hi.rho * hi.velocity_normal * hi.velocity_normal;
        Set::Scalar ue_hi = hi.pressure / (gamma - 1);
        Set::Scalar h_hi = (ke_hi + ue_hi + hi.pressure) / hi.rho;

        //Util::Message(INFO, "hi.rho ", hi.rho);

        // roe averages
        Set::Scalar rho_RL = std::sqrt(lo.rho) * std::sqrt(hi.rho);
        Set::Scalar u_RL = (std::sqrt(hi.rho) * hi.velocity_normal + std::sqrt(lo.rho) * lo.velocity_normal) / (std::sqrt(hi.rho) + std::sqrt(lo.rho));
        Set::Scalar h_RL = (std::sqrt(hi.rho) * h_hi + std::sqrt(lo.rho) * h_lo) / (std::sqrt(hi.rho) + std::sqrt(lo.rho));
        Set::Scalar a_RL2 = (gamma - 1) * (h_RL - 0.5 * u_RL * u_RL);
        if (std::isnan(a_RL2)) Util::Abort(INFO);
        if (a_RL2 < 0.0) Util::Abort(INFO);
        Set::Scalar a_RL = std::sqrt(a_RL2);

        // flux
        Flux fl;
        fl.mass = 0.5 * (lo.rho * lo.velocity_normal * lo.eta + hi.rho * hi.velocity_normal * hi.eta);
        fl.momentum_normal = 0.5 * (lo.rho * lo.velocity_normal * lo.velocity_normal * lo.eta + hi.rho * hi.velocity_normal * hi.velocity_normal * hi.eta) + 0.5 * (lo.pressure + hi.pressure) * cell_eta;
        fl.momentum_tangent = 0.5 * (lo.rho * lo.velocity_tangent * lo.velocity_normal * lo.eta + hi.rho * hi.velocity_tangent * hi.velocity_normal * hi.eta);
        fl.energy = 0.5 * (lo.velocity_normal * ke_lo * lo.eta + hi.velocity_normal * ke_hi * hi.eta) + 0.5 * (lo.velocity_normal * (ue_lo + lo.pressure) + hi.velocity_normal * (ue_hi + hi.pressure)) * cell_eta;

        //compute eigenvalues
        std::array <double, 3> ev;
        ev[0] = std::abs(u_RL);
        ev[1] = std::abs(u_RL + a_RL);
        ev[2] = std::abs(u_RL - a_RL);

        //compute wave strengths

        Set::Scalar drho, dp, dvel;
        dvel = hi.velocity_normal * hi.eta - lo.velocity_normal * lo.eta;
        dp = (hi.pressure - lo.pressure) * cell_eta;
        drho = hi.rho * hi.eta - lo.rho * lo.eta;

        std::array <double, 3> dv;
        dv[0] = drho - dp / a_RL2;
        dv[1] = dvel + dp / (rho_RL * a_RL);
        dv[2] = dvel - dp / (rho_RL * a_RL);

        //construct right characteristic eigenvectors
        std::array <double, 3> rev_0;
        rev_0[0] = 1;
        rev_0[1] = u_RL;
        rev_0[2] = 0.5 * u_RL * u_RL;

        std::array <double, 3> rev_1;
        rev_1[0] = (rho_RL / 2 / a_RL);
        rev_1[1] = (rho_RL / 2 / a_RL) * (u_RL + a_RL);
        rev_1[2] = (rho_RL / 2 / a_RL) * (h_RL + a_RL * u_RL);

        std::array <double, 3> rev_2;
        rev_2[0] = -(rho_RL / 2 / a_RL);
        rev_2[1] = -(rho_RL / 2 / a_RL) * (u_RL - a_RL);
        rev_2[2] = -(rho_RL / 2 / a_RL) * (h_RL - a_RL * u_RL);

        //upwinding term 
        fl.mass += -0.5 * (rev_0[0] * ev[0] * dv[0] + rev_1[0] * ev[1] * dv[1] + rev_2[0] * ev[2] * dv[2]);
        fl.momentum_normal += -0.5 * (rev_0[1] * ev[0] * dv[0] + rev_1[1] * ev[1] * dv[1] + rev_2[1] * ev[2] * dv[2]);
        fl.momentum_tangent += 0.0;
        fl.energy += -0.5 * (rev_0[2] * ev[0] * dv[0] + rev_1[2] * ev[1] * dv[1] + rev_2[2] * ev[2] * dv[2]);

        return fl;
    }

    static int Test()
    {
        int failed = 0;
        for (int i = 0; i < 10; i++)
        {

            State hi, lo, fixed, pos_vel1, pos_vel2, neg_vel1, neg_vel2;
            hi.rho = Util::Random();
            hi.velocity_normal = Util::Random();
            hi.velocity_tangent = Util::Random();
            hi.pressure = Util::Random();
            hi.eta = Util::Random();

            lo.rho = Util::Random();
            lo.velocity_normal = Util::Random();
            lo.velocity_tangent = Util::Random();
            lo.pressure = Util::Random();
            lo.eta = Util::Random();

            fixed.rho = Util::Random();
            fixed.velocity_normal = Util::Random();
            fixed.velocity_tangent = Util::Random();
            fixed.pressure = Util::Random();
            fixed.eta = Util::Random();

            pos_vel1.rho = Util::Random();
            pos_vel1.velocity_normal = Util::Random();
            pos_vel1.velocity_tangent = Util::Random();
            pos_vel1.pressure = Util::Random();
            pos_vel1.eta = Util::Random();

            neg_vel1.rho = pos_vel1.rho;
            neg_vel1.velocity_normal = -pos_vel1.velocity_normal;
            neg_vel1.velocity_tangent = -pos_vel1.velocity_tangent;
            neg_vel1.pressure = pos_vel1.pressure;
            neg_vel1.eta = pos_vel1.eta;

            pos_vel2.rho = Util::Random();
            pos_vel2.velocity_normal = Util::Random();
            pos_vel2.velocity_tangent = Util::Random();
            pos_vel2.pressure = Util::Random();
            pos_vel2.eta = Util::Random();

            neg_vel2.rho = pos_vel2.rho;
            neg_vel2.velocity_normal = -pos_vel2.velocity_normal;
            neg_vel2.velocity_tangent = -pos_vel2.velocity_tangent;
            neg_vel2.pressure = pos_vel2.pressure;
            neg_vel2.eta = pos_vel2.eta;

            // 
            // Zero flux test: if the state is the same on both sides,
            //                 the flux should be zero.
            // 
            // Flux zero = Solve(fixed, fixed, 1.4, 0.5);

            // if (fabs(zero.mass) > 1E-10) { Util::Warning(INFO, "Nonzero mass flux ", zero.mass); failed++; }
            // if (fabs(zero.velocity_normal) > 1E-10) { Util::Warning(INFO, "Nonzero x.velocity flux ", zero.velocity_normal); failed++; }
            // if (fabs(zero.velocity_tangent) > 1E-10) { Util::Warning(INFO, "Nonzero y.velocity flux ", zero.velocity_tangent); failed++; }
            // if (fabs(zero.pressure) > 1E-10) { Util::Warning(INFO, "Nonzero.pressure flux ", zero.pressure); failed++; }

            // // 
            // // Antisymmetry test: reversing hi and lo should produce the same
            // //                    output 
            // // 

            // Flux plus = Solve(hi, lo, 1.4);
            // Flux minus = Solve(lo, hi, 1.4);

            // if (fabs(plus.mass + minus.mass) > 1E-10)
            // {
            //     Util::Warning(INFO, "mass flux error: plus=", plus.mass, " minus=", minus.mass); failed++;
            // }
            // if (fabs(plus.velocity_normal - minus.velocity_normal) > 1E-10)
            // {
            //     Util::Warning(INFO, .velocity x flux error: plus=", plus.velocity_normal, " minus=", minus.velocity_normal); failed++;
            // }
            // if (fabs(plus.velocity_tangent - minus.velocity_tangent) > 1E-10)
            // {                
            //     Util::Warning(INFO, .velocity y flux error: plus=", plus.velocity_tangent, " minus=", minus.velocity_tangent); failed++;
            // }
            // if (fabs(plus.pressure + minus.pressure) > 1E-10)
            // {
            //     Util::Warning(INFO, .pressure flux error: plus=", plus.pressure, " minus=", minus.pressure); failed++;
            // }


            // // 
            //     // Wave Propogation Test: changing the sign of velocity should produce negative mass
            //     //                   and.pressure fluxes and the same.velocity flux
            //     // 

            // Flux pos = Solve(pos_vel1, pos_vel2, 1.4);
            // Flux neg = Solve(neg_vel1, neg_vel2, 1.4);

            // if (fabs(pos.mass + neg.mass) > 1E-10)
            // {                
            //     Util::Warning(INFO, "mass flux error: positive vel=", pos.mass, " negative vel=", neg.mass); failed++;
            // }
            // if (fabs(pos.velocity_normal - neg.velocity_normal) > 1E-10)
            // {
            //     Util::Warning(INFO, .velocity x flux error: positive vel=", pos.velocity_normal, " negative vel=", neg.velocity_normal); failed++;
            // }
            // if (fabs(pos.velocity_tangent - neg.velocity_tangent) > 1E-10)
            // {
            //     Util::Warning(INFO, .velocity y flux error: positive vel=", pos.velocity_tangent, " negative vel=", neg.velocity_tangent); failed++;
            // }
            // if (fabs(pos.pressure + neg.pressure) > 1E-10)
            // {
            //     Util::Warning(INFO, .pressure flux error: positive vel=", pos.pressure, " negative vel=", neg.pressure); failed++;
            // }


        }

        return failed;
    }
};
}
}
}
