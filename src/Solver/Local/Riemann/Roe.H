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
        Set::Scalar momentum_normal;
        Set::Scalar momentum_tangent;
        Set::Scalar energy;
        Set::Scalar eta;
        // Construtor for convenience
        State() { rho = 0.0; momentum_normal = 0.0, momentum_tangent = 0.0; energy = 0.0; eta = 0.0; }
        State(Set::Scalar a_rho, Set::Scalar a_momentum_normal, Set::Scalar a_momentum_tangent, Set::Scalar a_energy, Set::Scalar a_eta)
            : rho(a_rho), momentum_normal(a_momentum_normal), momentum_tangent(a_momentum_tangent), energy(a_energy), eta(a_eta) {}
    };

    struct Flux {
        Set::Scalar mass;
        Set::Scalar energy;
        Set::Scalar momentum_normal;
        Set::Scalar momentum_tangent;
    };

    // todo: convert to static functor by replacing with static () when c++23 is released
    static Flux Solve(State lo, State hi, Set::Scalar gamma, Set::Scalar cell_eta) {
        Set::Scalar uNorm_lo = lo.momentum_normal / lo.rho;
        Set::Scalar uTang_lo = lo.momentum_tangent / lo.rho;
        Set::Scalar ke_lo = 0.5 * lo.rho * uNorm_lo * uNorm_lo;
        Set::Scalar ue_lo = lo.energy - ke_lo;
        Set::Scalar p_lo = (gamma - 1) * ue_lo;
        Set::Scalar h_lo = (lo.energy + p_lo) / lo.rho;
	Set::Scalar eta_lo = lo.eta;

	//Util::Message(INFO, "lo.rho ", lo.rho);
	//Util::Message(INFO, "hi.rho ", hi.rho);

        Set::Scalar uNorm_hi = hi.momentum_normal / hi.rho;
        Set::Scalar uTang_hi = hi.momentum_tangent / hi.rho;
        Set::Scalar ke_hi = 0.5 * hi.rho * uNorm_hi * uNorm_hi;
        Set::Scalar ue_hi = hi.energy - ke_hi;
        Set::Scalar p_hi = (gamma - 1) * ue_hi;
        Set::Scalar h_hi = (hi.energy + p_hi) / hi.rho;
	Set::Scalar eta_hi = hi.eta;

	//Util::Message(INFO, "hi.rho ", hi.rho);

        // roe averages
        Set::Scalar rho_RL = std::sqrt(lo.rho) * std::sqrt(hi.rho);
        Set::Scalar u_RL = (std::sqrt(hi.rho) * uNorm_hi + std::sqrt(lo.rho) * uNorm_lo) / (std::sqrt(hi.rho) + std::sqrt(lo.rho));
        Set::Scalar h_RL = (std::sqrt(hi.rho) * h_hi + std::sqrt(lo.rho) * h_lo) / (std::sqrt(hi.rho) + std::sqrt(lo.rho));
        Set::Scalar a_RL2 = (gamma - 1) * (h_RL - 0.5 * u_RL * u_RL);
        if (std::isnan(a_RL2)) Util::Abort(INFO);
        if (a_RL2 < 0.0) Util::Abort(INFO);
        Set::Scalar a_RL = std::sqrt(a_RL2);

        // flux
        Flux fl;
        fl.mass = 0.5 * (lo.rho * uNorm_lo * eta_lo + hi.rho * uNorm_hi * eta_hi);
        fl.momentum_normal = 0.5 * (lo.rho * uNorm_lo * uNorm_lo * eta_lo + hi.rho * uNorm_hi * uNorm_hi * eta_hi) + 0.5 * (p_lo + p_hi) * cell_eta;
        fl.momentum_tangent = 0.5 * (lo.rho * uTang_lo * uNorm_lo * eta_lo + hi.rho * uTang_hi * uNorm_hi * eta_hi);
        fl.energy = 0.5 * (uNorm_lo * ke_lo * eta_lo + uNorm_hi * ke_hi * eta_hi) + 0.5 * (uNorm_lo * (ue_lo + p_lo) + uNorm_hi * (ue_hi + p_hi)) * cell_eta;

        //compute eigenvalues
        std::array <double, 3> ev;
        ev[0] = std::abs(u_RL);
        ev[1] = std::abs(u_RL + a_RL);
        ev[2] = std::abs(u_RL - a_RL);

	//compute wave strengths

        Set::Scalar drho, dp, dvel;
        dvel = uNorm_hi * hi.eta - uNorm_lo * lo.eta;
        dp = (p_hi - p_lo);
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
            hi.momentum_normal = Util::Random();
            hi.momentum_tangent = Util::Random();
            hi.energy = Util::Random();
            hi.eta = Util::Random();

            lo.rho = Util::Random();
            lo.momentum_normal = Util::Random();
            lo.momentum_tangent = Util::Random();
            lo.energy = Util::Random();
            lo.eta = Util::Random();

            fixed.rho = Util::Random();
            fixed.momentum_normal = Util::Random();
            fixed.momentum_tangent = Util::Random();
            fixed.energy = Util::Random();
            fixed.eta = Util::Random();

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
            Flux zero = Solve(fixed, fixed, 1.4, 0.5);

            if (fabs(zero.mass) > 1E-10) { Util::Warning(INFO, "Nonzero mass flux ", zero.mass); failed++; }
            if (fabs(zero.momentum_normal) > 1E-10) { Util::Warning(INFO, "Nonzero x momentum flux ", zero.momentum_normal); failed++; }
            if (fabs(zero.momentum_tangent) > 1E-10) { Util::Warning(INFO, "Nonzero y momentum flux ", zero.momentum_tangent); failed++; }
            if (fabs(zero.energy) > 1E-10) { Util::Warning(INFO, "Nonzero energy flux ", zero.energy); failed++; }

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
            // if (fabs(plus.momentum_normal - minus.momentum_normal) > 1E-10)
            // {
            //     Util::Warning(INFO, "momentum x flux error: plus=", plus.momentum_normal, " minus=", minus.momentum_normal); failed++;
            // }
            // if (fabs(plus.momentum_tangent - minus.momentum_tangent) > 1E-10)
            // {                
            //     Util::Warning(INFO, "momentum y flux error: plus=", plus.momentum_tangent, " minus=", minus.momentum_tangent); failed++;
            // }
            // if (fabs(plus.energy + minus.energy) > 1E-10)
            // {
            //     Util::Warning(INFO, "energy flux error: plus=", plus.energy, " minus=", minus.energy); failed++;
            // }


            // // 
            //     // Wave Propogation Test: changing the sign of velocity should produce negative mass
            //     //                   and energy fluxes and the same momentum flux
            //     // 

            // Flux pos = Solve(pos_vel1, pos_vel2, 1.4);
            // Flux neg = Solve(neg_vel1, neg_vel2, 1.4);

            // if (fabs(pos.mass + neg.mass) > 1E-10)
            // {                
            //     Util::Warning(INFO, "mass flux error: positive vel=", pos.mass, " negative vel=", neg.mass); failed++;
            // }
            // if (fabs(pos.momentum_normal - neg.momentum_normal) > 1E-10)
            // {
            //     Util::Warning(INFO, "momentum x flux error: positive vel=", pos.momentum_normal, " negative vel=", neg.momentum_normal); failed++;
            // }
            // if (fabs(pos.momentum_tangent - neg.momentum_tangent) > 1E-10)
            // {
            //     Util::Warning(INFO, "momentum y flux error: positive vel=", pos.momentum_tangent, " negative vel=", neg.momentum_tangent); failed++;
            // }
            // if (fabs(pos.energy + neg.energy) > 1E-10)
            // {
            //     Util::Warning(INFO, "energy flux error: positive vel=", pos.energy, " negative vel=", neg.energy); failed++;
            // }


        }

        return failed;
    }
};
}
}
}
