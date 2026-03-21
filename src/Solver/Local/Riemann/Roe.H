//
// This implements the Riemann Roe solver.
//
// Notation and algorithm follow the presentation in Section 5.3.3
// of *Computational Gasdynamics* by Culbert B. Laney (page 88)
//
// This solver uses an optional entropy fix
//   Option 1: chimeracfd method https://chimeracfd.com/programming/gryphon/fluxroe.html
//   Option 2: Eq. 4.3.67 in *Computational Fluid Dynamics for Engineers and Scientists* by Sreenivas Jayanti
//

#ifndef SOLVER_LOCAL_RIEMANN_ROE_H
#define SOLVER_LOCAL_RIEMANN_ROE_H

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
class Roe : public Riemann
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
    int entropy_fix = 0;
    int lowmach = 0;
    Set::Scalar phi = NAN;

    static void Parse(Roe & value, IO::ParmParse & pp)
    {
        // enable to dump diagnostic data if the roe solver fails
        pp.query_default("verbose", value.verbose, 1);
        // apply entropy fix if tru
        pp.query_default("entropy_fix", value.entropy_fix, false);

        // Apply the lowmach fix descripte in Rieper 2010
        // "A low-Mach number fix for Roe’s approximate Riemann solver"
        pp.query_default("lowmach",value.lowmach,false);


        if (value.entropy_fix == 1)
            Util::Warning(INFO,"The entropy fix is experimental and should be used with caution");
        else if (value.entropy_fix == 2)
            Util::Warning(INFO,"The entropy fix is experimental and should be used with caution. Has previously caused errors with FlowDrivenCavity regression test");
    }

    virtual Flux Solve(State lo, State hi, Model::Gas::Gas& gas, Set::Patch<const Set::Scalar>& X, int i, int j, int k, int side, Set::Scalar small) override
    {
        std::vector<double> rho_L = lo.rho;
        std::vector<double> rho_R = hi.rho;
        Set::Scalar Mn_L  = lo.M_normal  ,  Mn_R  = hi.M_normal  ;
        Set::Scalar Mt_L  = lo.M_tangent ,  Mt_R  = hi.M_tangent ;
        Set::Scalar E_L   = lo.E         ,  E_R   = hi.E         ;

        // Ensure no negative densities
        Set::Scalar rho_sum_L = 0.0;
        Set::Scalar rho_sum_R = 0.0;
        for (int n=0; n<gas.nspecies; ++n)
        {
            rho_L[n] = std::max(small,rho_L[n]);
            rho_R[n] = std::max(small,rho_R[n]);
            rho_sum_L += rho_L[n];
            rho_sum_R += rho_R[n];
        }

        // STEP 1: Compute fluid primitives

        // This is a modification to Laney allowing for an arbitrary equation of state
        //////////////////////////////////////////////////////////////////////////////////////////////////
        Set::Scalar T_L=0, T_R=0, p_L=0, p_R=0, gamma_L=0, gamma_R=0;
        switch (side) {
            case 0: // xlo
                T_L = gas.ComputeT(rho_sum_L, lo.M_normal, lo.M_tangent, E_L, 500.0, X, i-1, j ,k);
                p_L = gas.ComputeP(rho_sum_L, T_L, X, i-1, j, k);
                gamma_L = gas.gamma(T_L, X, i-1, j, k);
                T_R = gas.ComputeT(rho_sum_R, hi.M_normal, hi.M_tangent, E_R, 500.0, X, i, j ,k);
                p_R = gas.ComputeP(rho_sum_R, T_R, X, i, j, k);
                gamma_R = gas.gamma(T_R, X, i, j, k);
                break;
            case 1: // xhi
                T_L = gas.ComputeT(rho_sum_L, lo.M_normal, lo.M_tangent, E_L, 500.0, X, i, j ,k);
                p_L = gas.ComputeP(rho_sum_L, T_L, X, i, j, k);
                gamma_L = gas.gamma(T_L, X, i, j, k);
                T_R = gas.ComputeT(rho_sum_R, hi.M_normal, hi.M_tangent, E_R, 500.0, X, i+1, j ,k);
                p_R = gas.ComputeP(rho_sum_R, T_R, X, i+1, j, k);
                gamma_R = gas.gamma(T_R, X, i+1, j, k);
                break;
            case 2: // ylo
                T_L = gas.ComputeT(rho_sum_L, lo.M_normal, lo.M_tangent, E_L, 500.0, X, i, j-1 ,k);
                p_L = gas.ComputeP(rho_sum_L, T_L, X, i, j-1, k);
                gamma_L = gas.gamma(T_L, X, i, j-1, k);
                T_R = gas.ComputeT(rho_sum_R, hi.M_normal, hi.M_tangent, E_R, 500.0, X, i, j ,k);
                p_R = gas.ComputeP(rho_sum_R, T_R, X, i, j, k);
                gamma_R = gas.gamma(T_R, X, i, j, k);
                break;
            case 3: // yhi
                T_L = gas.ComputeT(rho_sum_L, lo.M_normal, lo.M_tangent, E_L, 500.0, X, i, j ,k);
                p_L = gas.ComputeP(rho_sum_L, T_L, X, i, j, k);
                gamma_L = gas.gamma(T_L, X, i, j, k);
                T_R = gas.ComputeT(rho_sum_R, hi.M_normal, hi.M_tangent, E_R, 500.0, X, i, j+1 ,k);
                p_R = gas.ComputeP(rho_sum_R, T_R, X, i, j+1, k);
                gamma_R = gas.gamma(T_R, X, i, j+1, k);
                break;
            case 4: // zlo
                T_L = gas.ComputeT(rho_sum_L, lo.M_normal, lo.M_tangent, E_L, 500.0, X, i, j ,k-1);
                p_L = gas.ComputeP(rho_sum_L, T_L, X, i, j, k-1);
                gamma_L = gas.gamma(T_L, X, i, j, k-1);
                T_R = gas.ComputeT(rho_sum_R, hi.M_normal, hi.M_tangent, E_R, 500.0, X, i, j ,k);
                p_R = gas.ComputeP(rho_sum_R, T_R, X, i, j, k);
                gamma_R = gas.gamma(T_R, X, i, j, k);
                break;
            case 5: // zhi
                T_L = gas.ComputeT(rho_sum_L, lo.M_normal, lo.M_tangent, E_L, 500.0, X, i, j ,k);
                p_L = gas.ComputeP(rho_sum_L, T_L, X, i, j, k);
                gamma_L = gas.gamma(T_L, X, i, j, k);
                T_R = gas.ComputeT(rho_sum_R, hi.M_normal, hi.M_tangent, E_R, 500.0, X, i, j ,k+1);
                p_R = gas.ComputeP(rho_sum_R, T_R, X, i, j, k+1);
                gamma_R = gas.gamma(T_R, X, i, j, k+1);
                break;
            default:
                Util::Abort(INFO, "Unknown side in Riemann solver, ",side);
        }
        //////////////////////////////////////////////////////////////////////////////////////////////////

        Set::Scalar ke_L = 0.5 * (Mn_L * Mn_L /*+ Mt_L * Mt_L*/) / (rho_sum_L); // KE per unit volume
        Set::Scalar ue_L = E_L - ke_L;                                            // IE per unit volume
        Set::Scalar h_TL = (ke_L + ue_L + p_L) / (rho_sum_L);                   // specific stagnation enthalpy (per unit mass)

        Set::Scalar ke_R = 0.5 * (Mn_R * Mn_R /*+ Mt_R * Mt_R*/) / (rho_sum_R);
        Set::Scalar ue_R = E_R - ke_R;
        Set::Scalar h_TR = (ke_R + ue_R + p_R) / (rho_sum_R);

        Set::Scalar u_L   = Mn_L/(rho_sum_L),  u_R   = Mn_R/(rho_sum_R);
        Set::Scalar v_L   = Mt_L/(rho_sum_L),  v_R   = Mt_R/(rho_sum_R);
        
        //
        // STEP 2: Compute Roe-averaged quantities
        // 
        Set::Scalar rho_RL  = std::sqrt(rho_sum_L * rho_sum_R);
        Set::Scalar u_RL    = (std::sqrt(rho_sum_L) * u_L  + std::sqrt(rho_sum_R) * u_R ) / (std::sqrt(rho_sum_L) + std::sqrt(rho_sum_R));
        Set::Scalar h_RL    = (std::sqrt(rho_sum_L) * h_TL + std::sqrt(rho_sum_R) * h_TR) / (std::sqrt(rho_sum_L) + std::sqrt(rho_sum_R));
        // This is not thermodynamically consistent to compute gamma_RL - need to compute a Roe-Average composition, temperature, and then gamma for consistency
        // But is should be a semi-close appproximation for now.
        Set::Scalar gamma_RL    = (std::sqrt(rho_sum_L) * gamma_L + std::sqrt(rho_sum_R) * gamma_R) / (std::sqrt(rho_sum_L) + std::sqrt(rho_sum_R)); 
        Set::Scalar a_RL_sq = std::max(0.0,(gamma_RL - 1.0) * (h_RL - 0.5 * u_RL * u_RL));

#ifdef AMREX_DEBUG
        if (verbose && ((a_RL_sq<0) || (a_RL_sq!=a_RL_sq)))
        {   
            Util::Message(INFO, "sound speed ", a_RL_sq);

            //Util::Message(INFO, "mixed rho ", lo.rho, " ", hi.rho);
            Util::Message(INFO, "mixed Mn ", lo.M_normal, " ", hi.M_normal);
            Util::Message(INFO, "mixed Mt ", lo.M_tangent, " ", hi.M_tangent);
            Util::Message(INFO, "mixed E ", lo.E, " ", hi.E);

            //Util::Message(INFO, "fluid rho ", rho_L, " ", rho_R);
            Util::Message(INFO, "fluid Mn ", Mn_L, " ", Mn_R);
            Util::Message(INFO, "fluid Mt ", Mt_L, " ", Mt_R);
            Util::Message(INFO, "fluid E ", E_L, " ", E_R);

            //Util::Message(INFO, "fluid rho ", rho_L, " ", rho_R);
            Util::Message(INFO, "fluid u ", u_L, " ", u_R);
            Util::Message(INFO, "fluid v ", v_L, " ", v_R);
            Util::Message(INFO, "fluid p ", p_L, " ", p_R);
        }
        Util::AssertException(INFO,TEST(a_RL_sq==a_RL_sq)," a_RL_sq is nan/inf; (a_RL_sq=", a_RL_sq,")");
        Util::AssertException(INFO,TEST(a_RL_sq>=0),      " a_RL_sq is negative; (a_RL_sq=(",a_RL_sq,")");
#endif

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
        Set::Scalar deltarho= rho_sum_R - rho_sum_L;
        Set::Scalar deltap  = p_R - p_L;
        Set::Scalar deltau  = u_R - u_L;

        if (lowmach)
        {
            Set::Scalar v_RL    = (std::sqrt(rho_sum_L) * v_L  + std::sqrt(rho_sum_R) * v_R ) / 
                (std::sqrt(rho_sum_L) + std::sqrt(rho_sum_R));
            Set::Scalar Ma = (std::abs(u_RL) + std::abs(v_RL)) / (a_RL + small);
            Ma = std::min(Ma, 1.0);
            deltau *= Ma;
        }

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
        // ROE ENTROPY FIX (Source cited in header comments)
        //
        if (entropy_fix == 1) { // chimeracfd
            lambda1 = fabs(lambda1);
            lambda2 = fabs(lambda2);
            lambda3 = fabs(lambda3);
            if ( lambda1 < deltau ) lambda1 = 0.5*(lambda1*lambda1 + deltau*deltau)/deltau;
            if ( lambda2 < deltau ) lambda2 = 0.5*(lambda2*lambda2 + deltau*deltau)/deltau;
            if ( lambda3 < deltau ) lambda3 = 0.5*(lambda3*lambda3 + deltau*deltau)/deltau;
        }
        else if (entropy_fix == 2) { // Jayanti
            Set::Scalar a_L = std::sqrt(gamma_L * p_L / (rho_sum_L)); // sound speed
            Set::Scalar a_R = std::sqrt(gamma_R * p_R / (rho_sum_R));
            Set::Scalar lambda1_L = u_L;         Set::Scalar lambda1_R = u_R; // eigenvalues
            Set::Scalar lambda2_L = u_L + a_L;   Set::Scalar lambda2_R = u_R + a_R;
            Set::Scalar lambda3_L = u_L - a_L;   Set::Scalar lambda3_R = u_R - a_R;
            Set::Scalar fix1 = std::max(0.0, std::max(lambda1 - lambda1_L, lambda1_R - lambda1));
            Set::Scalar fix2 = std::max(0.0, std::max(lambda2 - lambda2_L, lambda2_R - lambda2));
            Set::Scalar fix3 = std::max(0.0, std::max(lambda3 - lambda3_L, lambda3_R - lambda3));
            if ( lambda1 < fix1 ) lambda1 = fix1;
            if ( lambda2 < fix2 ) lambda2 = fix2;
            if ( lambda3 < fix3 ) lambda3 = fix3;
        }

        //
        // STEP 7: Compute fluxes
        //
        Flux fl;

        Set::Scalar mixed_mass_flux = (0.5*(rho_sum_L*u_L + rho_sum_R*u_R) - 0.5*(
                    r11*fabs(lambda1)*deltav1 +
                    r21*fabs(lambda2)*deltav2 +
                    r31*fabs(lambda3)*deltav3)
        );
        for (int n=0; n<gas.nspecies; ++n)
        {
            Set::Scalar Y_L = rho_L[n] / rho_sum_L;
            Set::Scalar Y_R = rho_R[n] / rho_sum_R;
            Set::Scalar Y_face = (mixed_mass_flux >= 0.0) ? Y_L : Y_R;
            fl.mass[n] = mixed_mass_flux * Y_face;
        }
        
        if (fl.mass != fl.mass)
        {
            if (verbose)
            {
                Util::ParallelMessage(INFO,"hi ", hi);
                Util::ParallelMessage(INFO,"lo ", lo);
                //Util::ParallelMessage(INFO,"rho_R ", rho_R);
                //Util::ParallelMessage(INFO,"rho_L ", rho_L);
                Util::ParallelMessage(INFO,"rho_RL ", rho_RL); 
                Util::ParallelMessage(INFO,"u_R ", u_R);
                Util::ParallelMessage(INFO,"u_L ", u_L);
                Util::ParallelMessage(INFO,"u_RL ", u_RL); 
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


        fl.momentum_normal = ( 0.5*(rho_sum_L*u_L*u_L + p_L + rho_sum_R*u_R*u_R + p_R) - 0.5*(
                                    r12*fabs(lambda1)*deltav1 +
                                    r22*fabs(lambda2)*deltav2 +
                                    r32*fabs(lambda3)*deltav3)
            );


        fl.energy = (   0.5*(u_L*(ke_L + p_L + ue_L) + u_R*(ke_R + p_R + ue_R)) - 0.5*
                        (
                            r13*fabs(lambda1)*deltav1 +
                            r23*fabs(lambda2)*deltav2 +
                            r33*fabs(lambda3)*deltav3)
            );

        //
        // (Update the tangential momentum flux)
        //
        fl.momentum_tangent = 0.5 * (rho_sum_L * u_L * v_L + rho_sum_R * u_R * v_R);

        return fl;
    }


    //static int Test()
    //{
    //    Roe solver;
    //    

    //    int failed = 0;

    //    Set::Scalar gamma = 1.4;
    //    Set::Scalar pref = 10.0;
    //    Set::Scalar small = 1E-10;

    //    // Test 1: Tangential Velocity Difference - No Normal Flux
    //    try {
    //        State left  (1.0, 1.0, 0.0, 1.0);
    //        State center(1.0, 1.0, 1.0, 1.0);
    //        State right (1.0, 1.0, 2.0, 1.0);
    //        Flux fluxlo = solver.Solve(center, right, gamma, pref, small);
    //        Flux fluxhi = solver.Solve(left, center,  gamma, pref, small);

    //        if (fabs(fluxhi.mass - fluxlo.mass) > 1E-10
    //            || fabs(fluxhi.momentum_normal - fluxlo.momentum_normal) > 1E-10
    //            || fabs(fluxhi.energy - fluxlo.energy) > 1E-10) {
    //            Util::Warning(INFO,   "left:    ",left);
    //            Util::Warning(INFO,   "center:  ",center);
    //            Util::Warning(INFO,   "right:   ",right);
    //            Util::Warning(INFO,   "Fluxlo:  ",fluxlo);
    //            Util::Warning(INFO,   "Fluxhi:  ",fluxhi);
    //            Util::Exception(INFO, "Tangential velocity difference incorrectly affecting normal flux.");
    //        }
    //        Util::Test::SubMessage("Test 1: Tangential velocity should induce no normal flux",0);
    //    } catch (const std::runtime_error& e)
    //    {
    //        failed++;
    //        Util::Test::SubMessage("Test 1: Tangential velocity should induce no normal flux",1);
    //    }

    //    // Test 2: Pure Transverse Velocity Difference
    //    try {
    //        State left  (1.0, 0.0, 0.0, 1.0);
    //        State center(1.0, 0.0, 1.0, 1.0);
    //        State right (1.0, 0.0, 2.0, 1.0);
    //        Flux fluxlo = solver.Solve(left, center,  gamma, pref, small);
    //        Flux fluxhi = solver.Solve(center, right, gamma, pref, small);
    //        if (fabs(fluxhi.mass - fluxlo.mass) > 1E-10
    //            || fabs(fluxhi.momentum_normal - fluxlo.momentum_normal) > 1E-10
    //            || fabs(fluxhi.energy - fluxlo.energy) > 1E-10) {
    //            Util::Warning(INFO,   "left:  ",left);
    //            Util::Warning(INFO,   "center: ",center);
    //            Util::Warning(INFO,   "right: ",right);
    //            Util::Warning(INFO,   "Fluxhi:  ",fluxhi);
    //            Util::Warning(INFO,   "Fluxlo:  ",fluxlo);
    //            Util::Exception(INFO, "Pure transverse velocity difference affecting normal flux.");
    //        }
    //        Util::Test::SubMessage("Test 2: Pure transverse velocity difference",0);
    //    } catch (const std::runtime_error& e)
    //    {
    //        failed++;
    //        Util::Test::SubMessage("Test 2: Pure transverse velocity difference",1);
    //    }

    //    // Test 3: Symmetry Test (no flux across identical states)
    //    try {
    //        State left(1.0, 0.0, 0.0, 1.0);
    //        State center(1.0, 0.0, 0.0, 1.0);
    //        State right(1.0, 0.0, 0.0, 1.0);
    //        Flux fluxhi = solver.Solve(center, right, gamma, pref, small);
    //        Flux fluxlo = solver.Solve(left, center, gamma, pref, small);
    //        if (fabs(fluxhi.mass - fluxlo.mass) > 1E-10 // no change in mass flux
    //            || fabs(fluxhi.momentum_normal - fluxlo.momentum_normal) > 1E-10 // no change in momentum flux
    //            || fabs(fluxhi.momentum_tangent) > 1E-10 // zero tangent flux
    //            || fabs(fluxlo.momentum_tangent) > 1E-10 // zero tangent flux
    //            || fabs(fluxhi.energy-fluxlo.energy) > 1E-10 // no change in energy flux
    //            ) {
    //            Util::Warning(INFO,   "left:  ",left);
    //            Util::Warning(INFO,   "right: ",right);
    //            Util::Warning(INFO,   "Fluxhi:  ",fluxhi);
    //            Util::Warning(INFO,   "Fluxlo:  ",fluxlo);
    //            Util::Exception(INFO, "Symmetric states should result in zero flux.");
    //        }
    //        Util::Test::SubMessage("Test 3: Constant states induces no flux difference",0);
    //    } catch (const std::runtime_error& e)
    //    {
    //        failed++;
    //        Util::Test::SubMessage("Test 3: Constant states induces no flux difference",1);
    //    }

    //    // Test 4: Uniform Flow Test (no flux across uniform flow)
    //    try {
    //        State left (1.0, 1.0, 0.5, 1.0);
    //        State center(1.0, 1.0, 0.5, 1.0);
    //        State right (1.0, 1.0, 0.5, 1.0);
    //        Flux fluxhi = solver.Solve(center, right, gamma, pref, small);
    //        Flux fluxlo = solver.Solve(left, center, gamma, pref, small);
    //        if (fabs(fluxhi.mass - fluxlo.mass) > 1E-10 ||
    //            fabs(fluxhi.momentum_normal - fluxlo.momentum_normal) > 1E-10 ||
    //            fabs(fluxhi.energy - fluxlo.energy) > 1E-10) {
    //            Util::Warning(INFO,   "left:   ",left);
    //            Util::Warning(INFO,   "center: ",center);
    //            Util::Warning(INFO,   "right:  ",right);
    //            Util::Warning(INFO,   "Fluxlo: ",fluxlo);
    //            Util::Warning(INFO,   "Fluxhi: ",fluxhi);
    //            Util::Exception(INFO, "Uniform flow should result in no flux.");
    //        }
    //        Util::Test::SubMessage("Test 4: Uniform flow should maintain constant flux",0);
    //    } catch (const std::runtime_error& e)
    //    {
    //        failed++;
    //        Util::Test::SubMessage("Test 4: Uniform flow should maintain constant flux",1);
    //    }

    //    return failed;
    //}
};
}
}
}

#endif
