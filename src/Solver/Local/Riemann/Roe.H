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
        Set::Scalar rho = NAN;
        Set::Scalar M_normal = NAN;
        Set::Scalar M_tangent = NAN;
        Set::Scalar E = NAN;
        Set::Scalar eta = NAN;
        // Construtor for convenience
        State() { rho = 0.0; M_normal = 0.0, M_tangent = 0.0; E = 0.0; eta = 0.0;}
        State(Set::Scalar a_rho, Set::Scalar a_M_normal, Set::Scalar a_M_tangent, Set::Scalar a_E, Set::Scalar a_eta)
            : rho(a_rho), M_normal(a_M_normal), M_tangent(a_M_tangent), E(a_E), eta(a_eta) {}
        friend std::ostream &operator<<(std::ostream &os, const State &state) 
        {
            os << "rho=" << state.rho << ", "; 
            os << "Mn=" << state.M_normal << ", ";
            os << "Mt=" << state.M_tangent << ", ";
            os << "E=" << state.E << ", ";
            os << "eta=" << state.eta;
            // do stuf
            return os;
        }
    };


    struct Flux {
        Set::Scalar mass;
        Set::Scalar momentum_normal;
        Set::Scalar momentum_tangent;
        Set::Scalar energy;
        friend std::ostream &operator<<(std::ostream &os, const Flux &flux) 
        {
            os << "mass=" << flux.mass << ", "; 
            os << "Mn=" << flux.momentum_normal << ", ";
            os << "Mt=" << flux.momentum_tangent << ", ";
            os << "E=" << flux.energy << ", ";
            // do stuff
            return os;
        }
    };

    // todo: convert to static functor by replacing with static () when c++23 is released
    static Flux Solve(State lo, State hi, State solid_lo, State solid_hi, Set::Scalar gamma, Set::Scalar cell_eta, Set::Scalar p_ref, Set::Scalar small)
    {
        // STEP 0: Compute fluid fields 
        Set::Scalar rho_L = (lo.rho       - (1.0 - lo.eta) * solid_lo.rho)      /(lo.eta + small),  rho_R = (hi.rho       - (1.0 - hi.eta) * solid_hi.rho)      /(hi.eta + small) + small;
        Set::Scalar Mn_L  = (lo.M_normal  - (1.0 - lo.eta) * solid_lo.M_normal) /(lo.eta + small),  Mn_R  = (hi.M_normal  - (1.0 - hi.eta) * solid_hi.M_normal) /(hi.eta + small);
        Set::Scalar Mt_L  = (lo.M_tangent - (1.0 - lo.eta) * solid_lo.M_tangent)/(lo.eta + small),  Mt_R  = (hi.M_tangent - (1.0 - hi.eta) * solid_hi.M_tangent)/(hi.eta + small);
        Set::Scalar E_L   = (lo.E         - (1.0 - lo.eta) * solid_lo.E)        /(lo.eta + small),  E_R   = (hi.E         - (1.0 - hi.eta) * solid_hi.E)        /(hi.eta + small);

        // Ensure no negative densities
        rho_L = std::max(0.0,rho_L);
        rho_R = std::max(0.0,rho_R);

        // STEP 1: Compute fluid primitives 
        Set::Scalar ke_L = 0.5 * (Mn_L * Mn_L /*+ Mt_L * Mt_L*/) / (rho_L+small); // KE per unit volume
        Set::Scalar ue_L = E_L - ke_L;                                            // IE per unit volume
        Set::Scalar p_L  = (gamma - 1.0) * ue_L + p_ref;                          // pressure
        Set::Scalar h_TL = (ke_L + ue_L + p_L) / (rho_L+small);                   // specific stagnation enthalpy (per unit mass)

        Set::Scalar ke_R = 0.5 * (Mn_R * Mn_R /*+ Mt_R * Mt_R*/) / (rho_R+small);
        Set::Scalar ue_R = E_R - ke_R;
        Set::Scalar p_R  = (gamma - 1.0) * ue_R + p_ref;
        Set::Scalar h_TR = (ke_R + ue_R + p_R) / (rho_R+small);

        Set::Scalar u_L   = Mn_L/(rho_L+small),  u_R   = Mn_R/(rho_R+small);
        Set::Scalar v_L   = Mt_L/(rho_L+small),  v_R   = Mt_R/(rho_R+small);
        
        //
        // STEP 2: Compute Roe-averaged quantities
        // 
        Set::Scalar rho_RL  = std::sqrt(rho_L * rho_R);
        Set::Scalar u_RL    = (std::sqrt(rho_L) * u_L  + std::sqrt(rho_R) * u_R ) / (std::sqrt(rho_L) + std::sqrt(rho_R) + small);
        Set::Scalar h_RL    = (std::sqrt(rho_L) * h_TL + std::sqrt(rho_R) * h_TR) / (std::sqrt(rho_L) + std::sqrt(rho_R) + small);
        Set::Scalar a_RL_sq = std::max(0.0,(gamma - 1.0) * (h_RL - 0.5 * u_RL * u_RL));

        if ((a_RL_sq<0) || (a_RL_sq!=a_RL_sq))
        {   
            Util::Message(INFO, "sound speed ", a_RL_sq);
            Util::Message(INFO, lo.eta, " ", hi.eta);

            Util::Message(INFO, "solid rho ", solid_lo.rho, " ",       solid_hi.rho);
            Util::Message(INFO, "solid Mn ", solid_lo.M_normal, " ",  solid_hi.M_normal);
            Util::Message(INFO, "solid Mt ", solid_lo.M_tangent, " ", solid_hi.M_tangent);
            Util::Message(INFO, "solid E ", solid_lo.E, " ",         solid_hi.E);

            Util::Message(INFO, "mixed rho ", lo.rho, " ", hi.rho);
            Util::Message(INFO, "mixed Mn ", lo.M_normal, " ", hi.M_normal);
            Util::Message(INFO, "mixed Mt ", lo.M_tangent, " ", hi.M_tangent);
            Util::Message(INFO, "mixed E ", lo.E, " ", hi.E);

            Util::Message(INFO, "fluid rho ", rho_L, " ", rho_R);
            Util::Message(INFO, "fluid Mn ", Mn_L, " ", Mn_R);
            Util::Message(INFO, "fluid Mt ", Mt_L, " ", Mt_R);
            Util::Message(INFO, "fluid E ", E_L, " ", E_R);

            Util::Message(INFO, "fluid rho ", rho_L, " ", rho_R);
            Util::Message(INFO, "fluid u ", u_L, " ", u_R);
            Util::Message(INFO, "fluid v ", v_L, " ", v_R);
            Util::Message(INFO, "fluid p ", p_L, " ", p_R);
        }
        Util::AssertException(INFO,TEST(a_RL_sq==a_RL_sq)," a_RL_sq is nan/inf; (a_RL_sq=", a_RL_sq,")");
        Util::AssertException(INFO,TEST(a_RL_sq>=0),      " a_RL_sq is negative; (a_RL_sq=(",a_RL_sq,")");

        Set::Scalar a_RL = std::sqrt(a_RL_sq) + small;

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

        Set::Scalar deltav1 = deltarho - deltap / (a_RL_sq + small);       // 5.54a
        Set::Scalar deltav2 = deltau   + deltap / (rho_RL * a_RL + small); // 5.54b
        Set::Scalar deltav3 = deltau   - deltap / (rho_RL * a_RL + small); // 5.54c

        //
        // STEP 7: Compute fluxes
        //
        Flux fl;
        
        fl.mass =
            (
                rho_L * u_L 
                + (std::min(0.0,lambda1)*deltav1
                + (0.5*rho_RL/(a_RL))*( std::min(0.0,lambda2)*deltav2 -
                                        std::min(0.0,lambda3)*deltav3)) 
                ) * cell_eta
            ;
        
        if (fl.mass != fl.mass)
        {
            // u_RL    = (std::sqrt(rho_L) * u_L  + std::sqrt(rho_R) * u_R ) / (std::sqrt(rho_L) + std::sqrt(rho_R));
            Util::ParallelMessage(INFO,u_RL); // nan
            Util::ParallelMessage(INFO,rho_L); // 1.15
            Util::ParallelMessage(INFO,rho_R); // -0.067
            Util::ParallelMessage(INFO,u_R); // -0.01


            Util::ParallelMessage(INFO,lo);
            Util::ParallelMessage(INFO,hi);
            Util::ParallelMessage(INFO,solid_lo);
            Util::ParallelMessage(INFO,solid_hi);

            Util::ParallelMessage(INFO,rho_L);
            Util::ParallelMessage(INFO,u_L);
            Util::ParallelMessage(INFO,lambda1); 
            Util::ParallelMessage(INFO,lambda2); 
            Util::ParallelMessage(INFO,lambda3); 
            Util::ParallelMessage(INFO,deltav1); 
            Util::ParallelMessage(INFO,deltav2); // nan
            Util::ParallelMessage(INFO,deltav3); // nan
            Util::ParallelMessage(INFO,rho_RL); // nan
            Util::ParallelMessage(INFO,a_RL);
            Util::ParallelMessage(INFO,cell_eta);
            Util::Exception(INFO);
        }


        fl.momentum_normal =
            (
                rho_L * u_L * u_L + p_L
                + (u_RL*std::min(0.0,lambda1)*deltav1
                + (0.5*rho_RL/(a_RL))*( lambda2*std::min(0.0,lambda2)*deltav2 -
                                        lambda3*std::min(0.0,lambda3)*deltav3)) 
                ) * cell_eta
            ;
        
        fl.energy =
            (
                ke_L * u_L + (p_L + ue_L) * u_L
                + (0.5*u_RL*u_RL*std::min(0.0,lambda1)*deltav1
                + (0.5*rho_RL/(a_RL))*( (h_RL + a_RL*u_RL)*std::min(0.0,lambda2)*deltav2 -
                                        (h_RL - a_RL*u_RL)*std::min(0.0,lambda3)*deltav3))               
                ) * cell_eta
            ;

        //
        // (Update the tangential momentum flux)
        //
        fl.momentum_tangent = 0.5 * (rho_L * u_L * v_L + rho_R * u_R * v_R) * cell_eta;

        return fl;
    }

    static int Test()
    {
        int failed = 0;

        const State solid(1.0, 0.0, 0.0, 0.0, 1.0);
        const State leftsol(1.0, 0.0, 0.0, 0.0, 1.0);
        const State rightsol(1.0, 0.0, 0.0, 0.0, 1.0);
        const Set::Scalar gamma = 1.4;
        const Set::Scalar pref = 10.0;
        const Set::Scalar small = 1E-10;
        const Set::Scalar cell_eta = 1.0;

        // Test 1: Tangential Velocity Difference - No Normal Flux
        try {
            State left  (1.0, 1.0, 0.0, 1.0, 1.0);
            State center(1.0, 1.0, 1.0, 1.0, 1.0);
            State right (1.0, 1.0, 2.0, 1.0, 1.0);
            Flux fluxlo = Solve(center, right, solid, solid, gamma, cell_eta, pref, small);
            Flux fluxhi = Solve(left, center,  solid, solid, gamma, cell_eta, pref, small);

            if (fabs(fluxhi.mass - fluxlo.mass) > 1E-10
                || fabs(fluxhi.momentum_normal - fluxlo.momentum_normal) > 1E-10
                || fabs(fluxhi.energy - fluxlo.energy) > 1E-10) {
                Util::Warning(INFO,   "left:    ",left);
                Util::Warning(INFO,   "center:  ",center);
                Util::Warning(INFO,   "right:   ",right);
                Util::Warning(INFO,   "Fluxlo:  ",fluxlo);
                Util::Warning(INFO,   "Fluxhi:  ",fluxhi);
                Util::Exception(INFO, "Tangential velocity difference incorrectly affecting normal flux.");
            }
            Util::Test::SubMessage("Test 1: Tangential velocity should induce no normal flux",0);
        } catch ( std::_Nested_exception<std::runtime_error> &e)
        {
            failed++;
            Util::Test::SubMessage("Test 1: Tangential velocity should induce no normal flux",1);
        }

        // Test 2: Pure Transverse Velocity Difference
        try {
            State left  (1.0, 0.0, 0.0, 1.0, 1.0);
            State center(1.0, 0.0, 1.0, 1.0, 1.0);
            State right (1.0, 0.0, 2.0, 1.0, 1.0);
            Flux fluxlo = Solve(left, center, leftsol, rightsol, gamma, cell_eta, pref, small);
            Flux fluxhi = Solve(center, right, leftsol, rightsol, gamma, cell_eta, pref, small);

            if (fabs(fluxhi.mass - fluxlo.mass) > 1E-10
                || fabs(fluxhi.momentum_normal - fluxlo.momentum_normal) > 1E-10
                || fabs(fluxhi.energy - fluxlo.energy) > 1E-10) {
                Util::Warning(INFO,   "left:  ",left);
                Util::Warning(INFO,   "center: ",center);
                Util::Warning(INFO,   "right: ",right);
                Util::Warning(INFO,   "Fluxhi:  ",fluxhi);
                Util::Warning(INFO,   "Fluxlo:  ",fluxlo);
                Util::Exception(INFO, "Pure transverse velocity difference affecting normal flux.");
            }
            Util::Test::SubMessage("Test 2: Pure transverse velocity difference",0);
        } catch ( std::_Nested_exception<std::runtime_error> &e)
        {
            failed++;
            Util::Test::SubMessage("Test 2: Pure transverse velocity difference",1);
        }

        // Test 3: Symmetry Test (no flux across identical states)
        try {
            State left(1.0, 0.0, 0.0, 1.0, 1.0);
            State center(1.0, 0.0, 0.0, 1.0, 1.0);
            State right(1.0, 0.0, 0.0, 1.0, 1.0);
            Flux fluxhi = Solve(center, right, solid, solid, gamma, cell_eta, pref, small);
            Flux fluxlo = Solve(left, center, solid, solid, gamma, cell_eta, pref, small);

            if (fabs(fluxhi.mass - fluxlo.mass) > 1E-10 // no change in mass flux
                || fabs(fluxhi.momentum_normal - fluxlo.momentum_normal) > 1E-10 // no change in momentum flux
                || fabs(fluxhi.momentum_tangent) > 1E-10 // zero tangent flux
                || fabs(fluxlo.momentum_tangent) > 1E-10 // zero tangent flux
                || fabs(fluxhi.energy-fluxlo.energy) > 1E-10 // no change in energy flux
                ) {
                Util::Warning(INFO,   "left:  ",left);
                Util::Warning(INFO,   "right: ",right);
                Util::Warning(INFO,   "Fluxhi:  ",fluxhi);
                Util::Warning(INFO,   "Fluxlo:  ",fluxlo);
                Util::Exception(INFO, "Symmetric states should result in zero flux.");
            }
            Util::Test::SubMessage("Test 3: Constant states induces no flux difference",0);
        } catch ( std::_Nested_exception<std::runtime_error> &e)
        {
            failed++;
            Util::Test::SubMessage("Test 3: Constant states induces no flux difference",1);
        }

        // Test 4: Uniform Flow Test (no flux across uniform flow)
        try {
            State left (1.0, 1.0, 0.5, 1.0, 1.0);
            State center(1.0, 1.0, 0.5, 1.0, 1.0);
            State right (1.0, 1.0, 0.5, 1.0, 1.0);
            Flux fluxhi = Solve(center, right, leftsol, rightsol, gamma, cell_eta, pref, small);
            Flux fluxlo = Solve(left, center, leftsol, rightsol, gamma, cell_eta, pref, small);

            if (fabs(fluxhi.mass - fluxlo.mass) > 1E-10 ||
                fabs(fluxhi.momentum_normal - fluxlo.momentum_normal) > 1E-10 ||
                fabs(fluxhi.energy - fluxlo.energy) > 1E-10) {
                Util::Warning(INFO,   "left:   ",left);
                Util::Warning(INFO,   "center: ",center);
                Util::Warning(INFO,   "right:  ",right);
                Util::Warning(INFO,   "Fluxlo: ",fluxlo);
                Util::Warning(INFO,   "Fluxhi: ",fluxhi);
                Util::Exception(INFO, "Uniform flow should result in no flux.");
            }
            Util::Test::SubMessage("Test 4: Uniform flow should maintain constant flux",0);
        } catch ( std::_Nested_exception<std::runtime_error> &e)
        {
            failed++;
            Util::Test::SubMessage("Test 4: Uniform flow should maintain constant flux",1);
        }

        return failed;
    }
};
}
}
}

