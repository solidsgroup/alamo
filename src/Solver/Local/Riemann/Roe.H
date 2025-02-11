#ifndef SOLVER_LOCAL_RIEMANN_ROE_H
#define SOLVER_LOCAL_RIEMANN_ROE_H

//
// This implements the Riemann Roe solver.
//
// Notation and algorithm follow the presentation in Section 5.3.3
// of *Computational Gasdynamics* by Culbert B. Laney (page 88)
//

#include "Set/Set.H"
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

/// Roe Riemann Solver based on Gas Dynamics - Culbert B. Laney
class Roe
{
public:


    static constexpr const char* name = "roe";
    Roe (IO::ParmParse &pp, std::string name) 
    {pp_queryclass(name,*this);}
    Roe (IO::ParmParse &pp) 
    {pp_queryclass(*this);}
    Roe () 
    {
        IO::ParmParse pp;
        pp_queryclass(*this);
    }

    int verbose = 0;

    static void Parse(Roe & value, IO::ParmParse & pp)
    {
        // enable to dump diagnostic data if the roe solver fails
        pp.query_default("verbose",value.verbose,1);
    }

    Flux Solve(State lo, State hi, Set::Scalar gamma, Set::Scalar p_ref, Set::Scalar small, int roefix)
    {
        Set::Scalar rho_L = lo.rho       ,  rho_R = hi.rho;
        Set::Scalar Mn_L  = lo.M_normal  ,  Mn_R  = hi.M_normal  ;
        Set::Scalar Mt_L  = lo.M_tangent ,  Mt_R  = hi.M_tangent ;
        Set::Scalar E_L   = lo.E         ,  E_R   = hi.E         ;

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

        if (verbose && ((a_RL_sq<0) || (a_RL_sq!=a_RL_sq)))
        {   
            Util::Message(INFO, "sound speed ", a_RL_sq);

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
        // STEP 5: Compute the right eigenvectors
        //
        Set::Scalar r11 = 1.0;
        Set::Scalar r12 = u_RL;
        Set::Scalar r13 = 0.5*u_RL*u_RL;
        Set::Scalar r21 = 0.5*rho_RL/a_RL;
        Set::Scalar r22 = 0.5*rho_RL/a_RL * ( u_RL + a_RL );
        Set::Scalar r23 = 0.5*rho_RL/a_RL * ( h_RL + a_RL*u_RL );
        Set::Scalar r31 = -0.5*rho_RL/a_RL;
        Set::Scalar r32 = -0.5*rho_RL/a_RL * ( u_RL - a_RL );
        Set::Scalar r33 = -0.5*rho_RL/a_RL * ( h_RL - a_RL*u_RL );

        //
        // STEP 6: Compute solution - not needed since fluxes will be computed in STEP 7
        //

        //
        // ROE ENTROPY FIX (not included in cited paper)
        //
        if (roefix == 1) {
            Util::ParallelMessage(INFO,"Entropy fix is still experimental and needs to be investigated further");
            //  Eq. 4.3.67 in "Computational Fluid Dynamics for Engineers and Scientists" by Sreenivas Jayanti
            //Set::Scalar a_L = std::sqrt(gamma * p_L / (rho_L + small)); // sound speed
            //Set::Scalar a_R = std::sqrt(gamma * p_R / (rho_R + small));
            //Set::Scalar lambda1_L = u_L;         Set::Scalar lambda1_R = u_R; // eigenvalues
            //Set::Scalar lambda2_L = u_L + a_L;   Set::Scalar lambda2_R = u_R + a_R;
            //Set::Scalar lambda3_L = u_L - a_L;   Set::Scalar lambda3_R = u_R - a_R;
            //Set::Scalar fix1 = std::max(0.0, std::max(lambda1 - lambda1_L, lambda1_R - lambda1));
            //Set::Scalar fix2 = std::max(0.0, std::max(lambda2 - lambda2_L, lambda2_R - lambda2));
            //Set::Scalar fix3 = std::max(0.0, std::max(lambda3 - lambda3_L, lambda3_R - lambda3));
            //if ( lambda1 < fix1 ) lambda1 = fix1;
            //if ( lambda2 < fix2 ) lambda2 = fix2;
            //if ( lambda3 < fix3 ) lambda3 = fix3;

            // chimeracfd method https://chimeracfd.com/programming/gryphon/fluxroe.html
            lambda1 = fabs(lambda1);
            lambda2 = fabs(lambda2);
            lambda3 = fabs(lambda3);
            if ( lambda1 < deltau ) lambda1 = 0.5*(lambda1*lambda1 + deltau*deltau)/deltau;
            if ( lambda2 < deltau ) lambda2 = 0.5*(lambda2*lambda2 + deltau*deltau)/deltau;
            if ( lambda3 < deltau ) lambda3 = 0.5*(lambda3*lambda3 + deltau*deltau)/deltau;
        }

        //
        // STEP 7: Compute fluxes
        //
        Flux fl;
        
        //fl.mass =
        //    (
        //        rho_L * u_L 
        //        + (std::min(0.0,lambda1)*deltav1
        //        + (0.5*rho_RL/(a_RL))*( std::min(0.0,lambda2)*deltav2 -
        //                                std::min(0.0,lambda3)*deltav3)) 
        //        )
        //    ;
        fl.mass = (0.5*(rho_L*u_L + rho_R*u_R) - 0.5*(
                        r11*fabs(lambda1)*deltav1 +
                        r21*fabs(lambda2)*deltav2 +
                        r31*fabs(lambda3)*deltav3)
            );
        
        if (fl.mass != fl.mass)
        {
            if (verbose)
            {
                Util::ParallelMessage(INFO,"hi ", hi);
                Util::ParallelMessage(INFO,"lo ", lo);
                Util::ParallelMessage(INFO,"rho_R ", rho_R);
                Util::ParallelMessage(INFO,"rho_L ", rho_L);
                Util::ParallelMessage(INFO,"rho_RL ", rho_RL); 
                Util::ParallelMessage(INFO,"u_R ", u_R);
                Util::ParallelMessage(INFO,"u_L ", u_L);
                Util::ParallelMessage(INFO,"u_RL ", u_RL); 
                //Util::ParallelMessage(INFO,"a_R ", a_R);
                //Util::ParallelMessage(INFO,"a_L ", a_L);
                Util::ParallelMessage(INFO,"a_RL ", a_RL);
                Util::ParallelMessage(INFO,"lambda1 ", lambda1); 
                Util::ParallelMessage(INFO,"lambda2 ", lambda2); 
                Util::ParallelMessage(INFO,"lambda3 ", lambda3); 
                Util::ParallelMessage(INFO,"deltav1 ", deltav1); 
                Util::ParallelMessage(INFO,"deltav2 ", deltav2);
                Util::ParallelMessage(INFO,"deltav3 ", deltav3);
            }
            Util::Exception(INFO);
        }


        //fl.momentum_normal =
        //    (
        //        rho_L * u_L * u_L + p_L
        //        + (u_RL*std::min(0.0,lambda1)*deltav1
        //        + (0.5*rho_RL/(a_RL))*( lambda2*std::min(0.0,lambda2)*deltav2 -
        //                                lambda3*std::min(0.0,lambda3)*deltav3)) 
        //        ) 
        //    ;
        fl.momentum_normal = ( 0.5*(rho_L*u_L*u_L + p_L + rho_R*u_R*u_R + p_R) - 0.5*(
                                    r12*fabs(lambda1)*deltav1 +
                                    r22*fabs(lambda2)*deltav2 +
                                    r32*fabs(lambda3)*deltav3)
            );

        //fl.energy =
        //    (
        //        ke_L * u_L + (p_L + ue_L) * u_L
        //        + (0.5*u_RL*u_RL*std::min(0.0,lambda1)*deltav1
        //        + (0.5*rho_RL/(a_RL))*( (h_RL + a_RL*u_RL)*std::min(0.0,lambda2)*deltav2 -
        //                                (h_RL - a_RL*u_RL)*std::min(0.0,lambda3)*deltav3))               
        //        ) 
        //    ;
        fl.energy = (   0.5*(u_L*(ke_L + p_L + ue_L) + u_R*(ke_R + p_R + ue_R)) - 0.5*
                        (
                            r13*fabs(lambda1)*deltav1 +
                            r23*fabs(lambda2)*deltav2 +
                            r33*fabs(lambda3)*deltav3)
            );

        //
        // (Update the tangential momentum flux)
        //
        fl.momentum_tangent = 0.5 * (rho_L * u_L * v_L + rho_R * u_R * v_R);

        return fl;
    }


    static int Test()
    {
        Roe solver;
        

        int failed = 0;

        Set::Scalar gamma = 1.4;
        Set::Scalar pref = 10.0;
        Set::Scalar small = 1E-10;
        int roefix = 0;

        // Test 1: Tangential Velocity Difference - No Normal Flux
        try {
            State left  (1.0, 1.0, 0.0, 1.0);
            State center(1.0, 1.0, 1.0, 1.0);
            State right (1.0, 1.0, 2.0, 1.0);
            Flux fluxlo = solver.Solve(center, right, gamma, pref, small, roefix);
            Flux fluxhi = solver.Solve(left, center,  gamma, pref, small, roefix);

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
            State left  (1.0, 0.0, 0.0, 1.0);
            State center(1.0, 0.0, 1.0, 1.0);
            State right (1.0, 0.0, 2.0, 1.0);
            Flux fluxlo = solver.Solve(left, center,  gamma, pref, small, roefix);
            Flux fluxhi = solver.Solve(center, right, gamma, pref, small, roefix);
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
            State left(1.0, 0.0, 0.0, 1.0);
            State center(1.0, 0.0, 0.0, 1.0);
            State right(1.0, 0.0, 0.0, 1.0);
            Flux fluxhi = solver.Solve(center, right, gamma, pref, small, roefix);
            Flux fluxlo = solver.Solve(left, center, gamma, pref, small, roefix);
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
            State left (1.0, 1.0, 0.5, 1.0);
            State center(1.0, 1.0, 0.5, 1.0);
            State right (1.0, 1.0, 0.5, 1.0);
            Flux fluxhi = solver.Solve(center, right, gamma, pref, small, roefix);
            Flux fluxlo = solver.Solve(left, center, gamma, pref, small, roefix);
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

#endif
