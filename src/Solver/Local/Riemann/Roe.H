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
        Set::Scalar etarho;
        Set::Scalar velocity_normal;
        Set::Scalar velocity_tangent;
        Set::Scalar etapressure;
        Set::Scalar eta;
        // Construtor for convenience
        State() { etarho = 0.0; velocity_normal = 0.0, velocity_tangent = 0.0; etapressure = 0.0; eta = 0.0; }
        State(Set::Scalar a_etarho, Set::Scalar a_velocity_normal, Set::Scalar a_velocity_tangent, Set::Scalar a_etapressure, Set::Scalar a_eta)
            : etarho(a_etarho), velocity_normal(a_velocity_normal), velocity_tangent(a_velocity_tangent), etapressure(a_etapressure), eta(a_eta) {}
    };

    struct Flux {
        Set::Scalar etaMass             = 0.0;
        Set::Scalar etaEnergy           = 0.0;
        Set::Scalar etaMomentum_normal  = 0.0;
        Set::Scalar etaMomentum_tangent = 0.0;
    };

    // todo: convert to static functor by replacing with static () when c++23 is released
    static Flux Solve(State lo, State hi, Set::Scalar gamma, Set::Scalar small=1.0e-8) {

        Set::Scalar eta_RL = (std::sqrt(hi.etarho) * hi.eta + std::sqrt(lo.etarho) * lo.eta) / (std::sqrt(hi.etarho) + std::sqrt(lo.etarho) + small);

        if (eta_RL < small) {Flux Zero; return Zero;};

        Set::Scalar etake_lo = 0.5 * lo.etarho * lo.velocity_normal * lo.velocity_normal;
        Set::Scalar etaue_lo = lo.etapressure / (gamma - 1.0);
        Set::Scalar etah_lo = (etake_lo + etaue_lo + lo.etapressure) / (lo.etarho + small);

        Set::Scalar etake_hi = 0.5 * hi.etarho * hi.velocity_normal * hi.velocity_normal;
        Set::Scalar etaue_hi = hi.etapressure / (gamma - 1.0);
        Set::Scalar etah_hi = (etake_hi + etaue_hi + hi.etapressure) / (hi.etarho + small);

        // roe averages
        Set::Scalar etarho_RL =  std::sqrt(lo.etarho) * std::sqrt(hi.etarho);
        Set::Scalar u_RL      = (std::sqrt(hi.etarho) * hi.velocity_normal + std::sqrt(lo.etarho) * lo.velocity_normal) / (std::sqrt(hi.etarho) + std::sqrt(lo.etarho) + small);
        Set::Scalar etah_RL   = (std::sqrt(hi.etarho) * etah_hi            + std::sqrt(lo.etarho) * etah_lo           ) / (std::sqrt(hi.etarho) + std::sqrt(lo.etarho) + small);
        Set::Scalar etap_RL   = (std::sqrt(hi.etarho) * hi.etapressure     + std::sqrt(lo.etarho) * lo.etapressure    ) / (std::sqrt(hi.etarho) + std::sqrt(lo.etarho) + small);
        Set::Scalar a_RL2     = gamma * etap_RL / (etarho_RL + small);

        if (std::isnan(a_RL2)) Util::Abort(INFO, hi.etapressure, " ", lo.etapressure, " ", hi.etarho, " ", lo.etarho);
        if (a_RL2 < 0.0) Util::Abort(INFO, hi.etapressure, " ", lo.etapressure, " ", hi.etarho, " ", lo.etarho);

        Set::Scalar a_RL      = std::sqrt(a_RL2);

        Flux fl;
        fl.etaMass             = 0.5 * (lo.etarho * lo.velocity_normal + hi.etarho * hi.velocity_normal);
        fl.etaMomentum_normal  = 0.5 * (lo.etarho * lo.velocity_normal  * lo.velocity_normal + hi.etarho * hi.velocity_normal  * hi.velocity_normal);// + 0.5 * (lo.etapressure + hi.etapressure); 
        fl.etaMomentum_tangent = 0.5 * (lo.etarho * lo.velocity_tangent * lo.velocity_normal + hi.etarho * hi.velocity_tangent * hi.velocity_normal);
        fl.etaEnergy           = 0.5 * (lo.velocity_normal * etake_lo + hi.velocity_normal * etake_hi)                                              ;// + 0.5 * (lo.velocity_normal * (etaue_lo + lo.etapressure)  + hi.velocity_normal * (etaue_hi + hi.etapressure));

        if (eta_RL > small) {
        fl.etaMomentum_normal += 0.5 * (lo.etapressure + hi.etapressure)/(eta_RL + small);
        fl.etaMomentum_normal += 0.5 * (lo.velocity_normal * (etaue_lo + lo.etapressure)  + hi.velocity_normal * (etaue_hi + hi.etapressure));
        };        

        //compute eigenvalues
        std::array <double, 3> ev;
        ev[0] = std::abs(u_RL);
        ev[1] = std::abs(u_RL + a_RL);
        ev[2] = std::abs(u_RL - a_RL);

        Set::Scalar detarho, detap, detavel;
        detavel = hi.velocity_normal * hi.eta - lo.velocity_normal * lo.eta;
        detap   = hi.etapressure              - lo.etapressure             ;
        if (eta_RL < small) { detap = 0.0; };
        detarho = hi.etarho                   - lo.etarho                  ;

        //compute wave strengths
        std::array <double, 3> detav;
        detav[0] = detarho - detap / (a_RL2 + small);
        detav[1] = detavel + detap * eta_RL / (etarho_RL * a_RL + small);
        detav[2] = detavel - detap * eta_RL / (etarho_RL * a_RL + small);

        //construct right characteristic eigenvectors
        std::array <double, 3> rev_0;
        rev_0[0] = 1;
        rev_0[1] = u_RL;
        rev_0[2] = 0.5 * u_RL * u_RL;

        std::array <double, 3> rev_1;
        rev_1[0] = (etarho_RL / (2.0 * a_RL * eta_RL + small));
        rev_1[1] = rev_1[0] * (u_RL + a_RL);
        rev_1[2] = rev_1[0] * (etah_RL/(eta_RL + small) + a_RL * u_RL);

        std::array <double, 3> rev_2;
        rev_2[0] = -rev_1[0];
        rev_2[1] = rev_2[0] * (u_RL - a_RL);
        rev_2[2] = rev_2[0] * (etah_RL/(eta_RL + small) - a_RL * u_RL);

        //upwinding term 
        fl.etaMass             -= 0.5 * (rev_0[0] * ev[0] * detav[0] + rev_1[0] * ev[1] * detav[1] + rev_2[0] * ev[2] * detav[2]);
        fl.etaMomentum_normal  -= 0.5 * (rev_0[1] * ev[0] * detav[0] + rev_1[1] * ev[1] * detav[1] + rev_2[1] * ev[2] * detav[2]);
        fl.etaMomentum_tangent -= 0.0;
        fl.etaEnergy           -= 0.5 * (rev_0[2] * ev[0] * detav[0] + rev_1[2] * ev[1] * detav[1] + rev_2[2] * ev[2] * detav[2]);

        return fl;
    }

    static int Test()
    {
        int failed = 0;
        for (int i = 0; i < 10; i++)
        {

            // State hi, lo, fixed, pos_vel1, pos_vel2, neg_vel1, neg_vel2;
            // hi.rho = Util::Random();
            // hi.velocity_normal = Util::Random();
            // hi.velocity_tangent = Util::Random();
            // hi.pressure = Util::Random();
            // hi.eta = Util::Random();

            // lo.rho = Util::Random();
            // lo.velocity_normal = Util::Random();
            // lo.velocity_tangent = Util::Random();
            // lo.pressure = Util::Random();
            // lo.eta = Util::Random();

            // fixed.rho = Util::Random();
            // fixed.velocity_normal = Util::Random();
            // fixed.velocity_tangent = Util::Random();
            // fixed.pressure = Util::Random();
            // fixed.eta = Util::Random();

            // pos_vel1.rho = Util::Random();
            // pos_vel1.velocity_normal = Util::Random();
            // pos_vel1.velocity_tangent = Util::Random();
            // pos_vel1.pressure = Util::Random();
            // pos_vel1.eta = Util::Random();

            // neg_vel1.rho = pos_vel1.rho;
            // neg_vel1.velocity_normal = -pos_vel1.velocity_normal;
            // neg_vel1.velocity_tangent = -pos_vel1.velocity_tangent;
            // neg_vel1.pressure = pos_vel1.pressure;
            // neg_vel1.eta = pos_vel1.eta;

            // pos_vel2.rho = Util::Random();
            // pos_vel2.velocity_normal = Util::Random();
            // pos_vel2.velocity_tangent = Util::Random();
            // pos_vel2.pressure = Util::Random();
            // pos_vel2.eta = Util::Random();

            // neg_vel2.rho = pos_vel2.rho;
            // neg_vel2.velocity_normal = -pos_vel2.velocity_normal;
            // neg_vel2.velocity_tangent = -pos_vel2.velocity_tangent;
            // neg_vel2.pressure = pos_vel2.pressure;
            // neg_vel2.eta = pos_vel2.eta;

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
