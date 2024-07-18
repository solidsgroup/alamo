//
// This implements the Riemann Roe solver.
//
// Notation and algorithm follow the presentation in Section 5.3.3
// of *Computational Gasdynamics* by Culbert B. Laney (page 88)
//

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
        // Construtor for convenience
        State() { rho = 0.0; velocity_normal = 0.0, velocity_tangent = 0.0; pressure = 0.0; }
        State(Set::Scalar a_rho, Set::Scalar a_velocity_normal, Set::Scalar a_velocity_tangent, Set::Scalar a_pressure)
            : rho(a_rho), velocity_normal(a_velocity_normal), velocity_tangent(a_velocity_tangent), pressure(a_pressure) {}
    };

    struct Flux {
        Set::Scalar mass;
        Set::Scalar momentum_normal;
        Set::Scalar momentum_tangent;
        Set::Scalar energy;
    };

    // todo: convert to static functor by replacing with static () when c++23 is released
    static Flux Solve(State lo, State hi, Set::Scalar gamma)
    {
        // STEP 1: Compute (or rename) primitives 
        Set::Scalar rho_L = lo.rho,              rho_R = hi.rho;
        Set::Scalar u_L   = lo.velocity_normal,  u_R   = hi.velocity_normal;
        Set::Scalar v_L   = lo.velocity_tangent, v_R   = hi.velocity_tangent;
        Set::Scalar p_L   = lo.pressure,         p_R   = hi.pressure;

        Set::Scalar ke_L = 0.5 * rho_L * u_L * u_L;
        Set::Scalar ue_L = p_L / (gamma - 1);
        Set::Scalar h_TL = (ke_L + ue_L + p_L) / rho_L;

        Set::Scalar ke_R = 0.5 * rho_R * u_R * u_R;
        Set::Scalar ue_R = p_R / (gamma - 1);
        Set::Scalar h_TR = (ke_R + ue_R + p_R) / rho_R;
        
        //
        // STEP 2: Compute Roe-averaged quantities
        // 
        Set::Scalar rho_RL  = std::sqrt(rho_L * rho_R);
        Set::Scalar u_RL    = (std::sqrt(rho_L) * u_L  + std::sqrt(rho_R) * u_R ) / (std::sqrt(rho_L) + std::sqrt(rho_R));
        Set::Scalar h_RL    = (std::sqrt(rho_L) * h_TL + std::sqrt(rho_R) * h_TR) / (std::sqrt(rho_L) + std::sqrt(rho_R));
        Set::Scalar a_RL_sq = (gamma - 1.0) * (h_RL - 0.5 * u_RL * u_RL);

        if ((a_RL_sq<0) || (a_RL_sq!=a_RL_sq))
        {
            Util::Message(INFO,rho_L, " ", rho_R);
            Util::Message(INFO,u_L, " ", u_R);
            Util::Message(INFO,v_L, " ", v_R);
            Util::Message(INFO,p_L, " ", p_R);
        }
        Util::AssertException(INFO,TEST(a_RL_sq==a_RL_sq)," a_RL_sq is nan/inf; (a_RL_sq=", a_RL_sq,")");
        Util::AssertException(INFO,TEST(a_RL_sq>=0),      " a_RL_sq is negative; (a_RL_sq=(",a_RL_sq,")");

        Set::Scalar a_RL = std::sqrt(a_RL_sq);

        //
        // STEP 3: Compute Roe-averaged wave speeds
        //
        Set::Scalar lambda1 = u_RL;          // 5.53a
        Set::Scalar lambda2 = u_RL + a_RL;   // 5.53b
        Set::Scalar lambda3 = u_RL - a_RL;   // 5.53c

        //
        // STEP 4: Compute wave strengths
        //
        Set::Scalar deltarho= rho_R - rho_L;
        Set::Scalar deltap  = p_R - p_L;
        Set::Scalar deltau  = u_R - u_L;

        Set::Scalar deltav1 = deltarho - deltap / a_RL_sq;         // 5.54a
        Set::Scalar deltav2 = deltau   + deltap / (rho_RL * a_RL); // 5.54b
        Set::Scalar deltav3 = deltau   - deltap / (rho_RL * a_RL); // 5.54c

        //
        // STEP 7: Compute fluxes
        //
        Flux fl;
        
        fl.mass =
            rho_L*u_L
            + std::min(0.0,lambda1)*deltav1
            + (0.5*rho_RL/a_RL)*(   std::min(0.0,lambda2)*deltav2 -
                                    std::min(0.0,lambda3)*deltav3);
        
        fl.momentum_normal =
            rho_L*u_L*u_L + p_L
            + u_RL*std::min(0.0,lambda1)*deltav1
            + (0.5*rho_RL/a_RL)*(   lambda2*std::min(0.0,lambda2)*deltav2 -
                                    lambda3*std::min(0.0,lambda3)*deltav3);
        
        fl.energy =
            rho_L*h_TL*u_L
            + 0.5*u_RL*u_RL*std::min(0.0,lambda1)*deltav1
            + (0.5*rho_RL/a_RL)*(   (h_RL + a_RL*u_RL)*std::min(0.0,lambda2)*deltav2 -
                                    (h_RL - a_RL*u_RL)*std::min(0.0,lambda3)*deltav3);

        //
        // (Update the tangential momentum flux)
        //
        fl.momentum_tangent = 0.5 * (rho_L*u_L*v_L + rho_R*u_R*v_R);


        return fl;
    }

    static int Test()
    {
        //
        // !! DO NOT DISABLE THIS TEST !!
        //
        // If you **must** disable this test for some reason, then make sure that it returns
        // a nonzero value - otherwise future users may falsely assume that it is passing tests
        // when it is not.
        //

        int failed = 0;
        for (int i = 0; i < 10; i++)
        {

            State hi, lo, fixed, pos_vel1, pos_vel2, neg_vel1, neg_vel2;
            hi.rho = Util::Random();
            hi.velocity_normal = Util::Random();
            hi.velocity_tangent = Util::Random();
            hi.pressure = Util::Random();

            lo.rho = Util::Random();
            lo.velocity_normal = Util::Random();
            lo.velocity_tangent = Util::Random();
            lo.pressure = Util::Random();

            fixed.rho = Util::Random();
            fixed.velocity_normal = 0.0; // Util::Random();
            fixed.velocity_tangent = 0.0; //Util::Random();
            fixed.pressure = Util::Random();

            pos_vel1.rho = Util::Random();
            pos_vel1.velocity_normal = Util::Random();
            pos_vel1.velocity_tangent = Util::Random();
            pos_vel1.pressure = Util::Random();

            neg_vel1.rho = pos_vel1.rho;
            neg_vel1.velocity_normal = -pos_vel1.velocity_normal;
            neg_vel1.velocity_tangent = -pos_vel1.velocity_tangent;
            neg_vel1.pressure = pos_vel1.pressure;

            pos_vel2.rho = Util::Random();
            pos_vel2.velocity_normal = Util::Random();
            pos_vel2.velocity_tangent = Util::Random();
            pos_vel2.pressure = Util::Random();

            neg_vel2.rho = pos_vel2.rho;
            neg_vel2.velocity_normal = -pos_vel2.velocity_normal;
            neg_vel2.velocity_tangent = -pos_vel2.velocity_tangent;
            neg_vel2.pressure = pos_vel2.pressure;

            // 
            // Zero flux test: if the state is the same on both sides,
            //                 the flux should be zero.
            // 
            Flux zero = Solve(fixed, fixed, 1.4);

            

            if (fabs(zero.mass) > 1E-10) { Util::Warning(INFO, "Nonzero mass flux ", zero.mass); failed++; }
            if (fabs(zero.momentum_normal) > 1E-10) { Util::Warning(INFO, "Nonzero x.momentum flux ", zero.momentum_normal); failed++; }
            if (fabs(zero.momentum_tangent) > 1E-10) { Util::Warning(INFO, "Nonzero y.momentum flux ", zero.momentum_tangent); failed++; }
            if (fabs(zero.energy) > 1E-10) { Util::Warning(INFO, "Nonzero.pressure flux ", zero.energy); failed++; }

            // 
            // Antisymmetry test: reversing hi and lo should produce the same
            //                    output 
            // 

            Flux plus = Solve(hi, lo, 1.4);
            Flux minus = Solve(lo, hi, 1.4);

            if (fabs(plus.mass + minus.mass) > 1E-10)
            {
                Util::Warning(INFO, "mass flux error: plus=", plus.mass, " minus=", minus.mass);
                failed++;
            }
            if (fabs(plus.momentum_normal - minus.momentum_normal) > 1E-10)
            {
                Util::Warning(INFO, ".momentum x flux error: plus=", plus.momentum_normal, " minus=", minus.momentum_normal);
                failed++;
            }
            if (fabs(plus.momentum_tangent - minus.momentum_tangent) > 1E-10)
            {                
                Util::Warning(INFO, ".momentum y flux error: plus=", plus.momentum_tangent, " minus=", minus.momentum_tangent);
                failed++;
            }
            if (fabs(plus.energy + minus.energy) > 1E-10)
            {
                Util::Warning(INFO, ".energy flux error: plus=", plus.energy, " minus=", minus.energy);
                failed++;
            }


            // 
            // Wave Propagation Test: changing the sign of velocity should produce negative mass
            //                   and.pressure fluxes and the same.velocity flux
            // 

            Flux pos = Solve(pos_vel1, pos_vel2, 1.4);
            Flux neg = Solve(neg_vel1, neg_vel2, 1.4);

            if (fabs(pos.mass + neg.mass) > 1E-10)
            {                
                Util::Warning(INFO, "mass flux error: positive vel=", pos.mass, " negative vel=", neg.mass); failed++;
            }
            if (fabs(pos.momentum_normal - neg.momentum_normal) > 1E-10)
            {
                Util::Warning(INFO, ".momentum x flux error: positive momentum=", pos.momentum_normal, " negative momentum=", neg.momentum_normal); failed++;
            }
            if (fabs(pos.momentum_tangent - neg.momentum_tangent) > 1E-10)
            {
                Util::Warning(INFO, ".momentum y flux error: positive momentum=", pos.momentum_tangent, " negative momentum=", neg.momentum_tangent); failed++;
            }
            if (fabs(pos.energy + neg.energy) > 1E-10)
            {
                Util::Warning(INFO, ".energy flux error: positive energy=", pos.energy, " negative vel=", neg.energy); failed++;
            }


        }

        return failed;
    }
};
}
}
}
