// This implements a modified Riemann HLLC solver for weakly incompressible flow using the Tammann EOS.
//      HLLC solver uses the fastest left and right wave speeds, along with the middle contact wave.
//      Restores some accuracy when compared to HLLE with a slightly increased computational time, while providing same robustness.
//
// Source: V.A. Titarev, E.F. Toro, Finite-volume WENO schemes for three-dimensional conservation laws, Journal of Computational Physics, Volume 201, Issue 1, 2004, Pages 238-260, ISSN 0021-9991, https://doi.org/10.1016/j.jcp.2004.05.015.
//          URL: https://www.sciencedirect.com/science/article/pii/S0021999104002281
//
// Source: M. J. Del Razo and R. J. LeVeque, "Numerical Methods for Interface Coupling of Compressible and Almost Incompressible Media," SIAM Journal on Scientific Computing, vol. 39, no. 3, pp. B486–B507, Jan. 2017, doi: 10.1137/16m1067834.
//          URL: https://arxiv.org/abs/1604.08544
//          If you wish to use the regular gas HLLC sovler and not the stiffened gas modeled by the Tammann EOS, set P0 = 0

#ifndef SOLVER_LOCAL_RIEMANN_HLLC_H
#define SOLVER_LOCAL_RIEMANN_HLLC_H

#include "IO/ParmParse.H"
#include "Solver/Local/Riemann/Riemann.H"
#include "Util/Util.H"


namespace Solver {
namespace Local {
namespace Riemann {

class HLLC : public Riemann
{
public:
    static constexpr const char *name = "hllc";
    HLLC(IO::ParmParse &pp, std::string name) { pp_queryclass(name, *this); }
    HLLC(IO::ParmParse &pp) { pp_queryclass(*this); }
    HLLC()
    {
        IO::ParmParse pp;
        pp_queryclass(*this);
    }

    static void Parse(HLLC & /*value*/, IO::ParmParse & /*pp*/)
    {
        // pp.query_default("lowmach", value.lowmach, false);
        // pp.query_default("cutoff", value.cutoff, 0.1);
    }

    virtual Flux Solve(State lo, State hi, Set::Scalar gamma, Set::Scalar p_ref, Set::Scalar small) override
    {
        Set::Scalar rho_L = lo.rho, rho_R = hi.rho;
        Set::Scalar Mn_L = lo.M_normal, Mn_R = hi.M_normal;
        Set::Scalar Mt_L = lo.M_tangent, Mt_R = hi.M_tangent;
        Set::Scalar E_L = lo.E, E_R = hi.E;

        // Tamman EOS (NOT INTEGRATED HERE)
        Set::Scalar p0 = 0.0;

        // Ensure no negative densities
        rho_L = std::max(small, rho_L);
        rho_R = std::max(small, rho_R);

        // STEP 1: Compute fluid primitives
        Set::Scalar ke_L, ke_R;
        Set::Scalar ue_L, ue_R;
        Set::Scalar p_L, p_R;
        Set::Scalar h_TL, h_TR;
        
        // Left
        ke_L = 0.5 * (Mn_L * Mn_L + Mt_L * Mt_L) / (rho_L + small); // KE per unit volume
        ue_L = E_L - ke_L;                                          // IE per unit volume
        p_L = (gamma - 1.0) * ue_L - gamma * p0 + p_ref;      // pressure Tammann EOS modification ( e[J/m3] = P[N/m2]/(gamma-1) )
        h_TL = (ke_L + ue_L + p_L);
        // Right
        ke_R = 0.5 * (Mn_R * Mn_R + Mt_R * Mt_R) / (rho_R + small);
        ue_R = E_R - ke_R;
        p_R = (gamma - 1.0) * ue_R - gamma * p0 + p_ref; // Tammann
        h_TR = (ke_R + ue_R + p_R);
      
        Set::Scalar u_L = Mn_L / (rho_L + small), u_R = Mn_R / (rho_R + small);
        Set::Scalar v_L = Mt_L / (rho_L + small), v_R = Mt_R / (rho_R + small);

        // STEP 2: Compute sound speeds
        Set::Scalar a_L = std::sqrt(gamma * (p_L + p0) / (rho_L + small));
        Set::Scalar a_R = std::sqrt(gamma * (p_R + p0) / (rho_R + small));

        // STEP 3: Compute wave speed estimates
        Set::Scalar S_L = std::min(u_L - a_L, u_R - a_R);
        Set::Scalar S_R = std::max(u_L + a_L, u_R + a_R);

        // STEP 4: Compute middle wave speed (contact discontinuity)
        Set::Scalar S_star = (p_R - p_L + rho_L * u_L * (S_L - u_L) - rho_R * u_R * (S_R - u_R)) / (rho_L * (S_L - u_L) - rho_R * (S_R - u_R) + small);

        // STEP 5: Compute physical fluxes for left and right states
        Flux f_L, f_R;

        // STEP 5: Compute left/right fluxes
        // Energies per unit volume [J/m^3]
        f_L.mass = rho_L * u_L;
        f_L.momentum_normal = rho_L * u_L * u_L + p_L;
        f_L.momentum_tangent = rho_L * u_L * v_L;
        f_L.energy = u_L * (E_L + p_L); // [W/m^2]

        f_R.mass = rho_R * u_R;
        f_R.momentum_normal = rho_R * u_R * u_R + p_R;
        f_R.momentum_tangent = rho_R * u_R * v_R;
        f_R.energy = u_R * (E_R + p_R);

        // STEP 6: Compute star state values
        Set::Scalar rho_star_L = rho_L * (S_L - u_L) / (S_L - S_star + small);
        Set::Scalar rho_star_R = rho_R * (S_R - u_R) / (S_R - S_star + small);

        Set::Scalar E_star_L = rho_star_L * (E_L / rho_L + (S_star - u_L) * (S_star + p_L / (rho_L * (S_L - u_L))));
        Set::Scalar E_star_R = rho_star_R * (E_R / rho_R + (S_star - u_R) * (S_star + p_R / (rho_R * (S_R - u_R))));

        // STEP 7: Compute star state fluxes
        Flux f_star_L, f_star_R;

        // Volumetric star fluxes [J/m^3]
        f_star_L.mass = f_L.mass + S_L * (rho_star_L - rho_L);
        f_star_L.momentum_normal = f_L.momentum_normal + S_L * (rho_star_L * S_star - rho_L * u_L);
        f_star_L.momentum_tangent = f_L.momentum_tangent + S_L * (rho_star_L * v_L - rho_L * v_L);
        // f_star_L.energy = f_L.energy + S_L * (rho_star_L * (E_L / rho_L + (S_star - u_L) * (S_star + p_L / (rho_L * (S_L - u_L) + small))) - E_L);
        f_star_L.energy = f_L.energy + S_L * (E_star_L - E_L); // CORRECTED

        f_star_R.mass = f_R.mass + S_R * (rho_star_R - rho_R);
        f_star_R.momentum_normal = f_R.momentum_normal + S_R * (rho_star_R * S_star - rho_R * u_R);
        f_star_R.momentum_tangent = f_R.momentum_tangent + S_R * (rho_star_R * v_R - rho_R * v_R);
        // f_star_R.energy = f_R.energy + S_R * (rho_star_R * (E_R / rho_R + (S_star - u_R) * (S_star + p_R / (rho_R * (S_R - u_R) + small))) - E_R);
        f_star_R.energy = f_R.energy + S_R * (E_star_R - E_R); // CORRECTED

        // STEP 8: Compute HLLC flux based on wave positions
        Flux fl;

        if (0.0 <= S_L)
        {
            fl = f_L;
        }
        else if (S_L <= 0.0 && 0.0 <= S_star)
        {
            fl = f_star_L;
        }
        else if (S_star <= 0.0 && 0.0 <= S_R)
        {
            fl = f_star_R;
        }
        else
        {
            fl = f_R;
        }

        // ERROR CHECKING
        if ((fl.mass != fl.mass)
            or (fl.momentum_normal != fl.momentum_normal)
            or (fl.momentum_tangent != fl.momentum_tangent)
            or (fl.energy != fl.energy)
            or (a_L != a_L)
            or (a_R != a_R))
        {
            bool verbose = false;
            if (verbose)
            {
                Util::ParallelMessage(INFO, "-------------------------------");
                Util::ParallelMessage(INFO, "ERROR IN HLLC SOLVER:");
                Util::ParallelMessage(INFO, "hi:   ", hi);
                Util::ParallelMessage(INFO, "lo:   ", lo);
                Util::ParallelMessage(INFO, "flux: ", fl);
                Util::ParallelMessage(INFO, "rho_L=", rho_L);
                Util::ParallelMessage(INFO, "rho_R=", rho_R);
                Util::ParallelMessage(INFO, "u_L=  ", u_L);
                Util::ParallelMessage(INFO, "u_R=  ", u_R);
                Util::ParallelMessage(INFO, "ke_R= ", ke_R);
                Util::ParallelMessage(INFO, "ke_L= ", ke_L);
                Util::ParallelMessage(INFO, "ue_L= ", ue_L);
                Util::ParallelMessage(INFO, "ue_R= ", ue_R);
                Util::ParallelMessage(INFO, "gamma= ", gamma);
                Util::ParallelMessage(INFO, "gamma= ", gamma);
                Util::ParallelMessage(INFO, "p_L=  ", p_L);
                Util::ParallelMessage(INFO, "p_R=  ", p_R);
                Util::ParallelMessage(INFO, "a_L=  ", a_L);
                Util::ParallelMessage(INFO, "a_R=  ", a_R);
                Util::ParallelMessage(INFO, "S_L=  ", S_L);
                Util::ParallelMessage(INFO, "S_R=  ", S_R);
                Util::ParallelMessage(INFO, "S_star ", S_star);
            }
            Util::Exception(INFO);
        }

        return fl;
    }

    static int
    Test()
    {
        HLLC solver;

        int failed = 0;
        Set::Scalar gamma = 1.4;
        Set::Scalar pref = 10.0;
        Set::Scalar small = 1E-10;

        // Test 1: Tangential Velocity Difference - No Normal Flux
        try
        {
            State left(1.0, 1.0, 0.0, 1.0);
            State center(1.0, 1.0, 1.0, 1.0);
            State right(1.0, 1.0, 2.0, 1.0);
            Flux fluxlo = solver.Solve(center, right, gamma, pref, small);
            Flux fluxhi = solver.Solve(left, center, gamma, pref, small);
            if (fabs(fluxhi.mass - fluxlo.mass) > 1E-10
                || fabs(fluxhi.momentum_normal - fluxlo.momentum_normal) > 1E-10
                || fabs(fluxhi.energy - fluxlo.energy) > 1E-10)
            {
                Util::Warning(INFO, "left:    ", left);
                Util::Warning(INFO, "center:  ", center);
                Util::Warning(INFO, "right:   ", right);
                Util::Warning(INFO, "Fluxlo:  ", fluxlo);
                Util::Warning(INFO, "Fluxhi:  ", fluxhi);
                Util::Exception(INFO, "Tangential velocity difference incorrectly affecting normal flux.");
            }
            Util::Test::SubMessage("Test 1: Tangential velocity should induce no normal flux", 0);
        }
        catch (const std::runtime_error &e)
        {
            failed++;
            Util::Test::SubMessage("Test 1: Tangential velocity should induce no normal flux", 1);
        }

        // Test 2: Pure Transverse Velocity Difference
        try
        {
            State left(1.0, 0.0, 0.0, 1.0);
            State center(1.0, 0.0, 1.0, 1.0);
            State right(1.0, 0.0, 2.0, 1.0);
            Flux fluxlo = solver.Solve(left, center, gamma, pref, small);
            Flux fluxhi = solver.Solve(center, right, gamma, pref, small);
            if (fabs(fluxhi.mass - fluxlo.mass) > 1E-10
                || fabs(fluxhi.momentum_normal - fluxlo.momentum_normal) > 1E-10
                || fabs(fluxhi.energy - fluxlo.energy) > 1E-10)
            {
                Util::Warning(INFO, "left:  ", left);
                Util::Warning(INFO, "center: ", center);
                Util::Warning(INFO, "right: ", right);
                Util::Warning(INFO, "Fluxhi:  ", fluxhi);
                Util::Warning(INFO, "Fluxlo:  ", fluxlo);
                Util::Exception(INFO, "Pure transverse velocity difference affecting normal flux.");
            }
            Util::Test::SubMessage("Test 2: Pure transverse velocity difference", 0);
        }
        catch (const std::runtime_error &e)
        {
            failed++;
            Util::Test::SubMessage("Test 2: Pure transverse velocity difference", 1);
        }

        // Test 3: Symmetry Test (no flux across identical states)
        try
        {
            State left(1.0, 0.0, 0.0, 1.0);
            State center(1.0, 0.0, 0.0, 1.0);
            State right(1.0, 0.0, 0.0, 1.0);
            Flux fluxhi = solver.Solve(center, right, gamma, pref, small);
            Flux fluxlo = solver.Solve(left, center, gamma, pref, small);
            if (fabs(fluxhi.mass - fluxlo.mass) > 1E-10                          // no change in mass flux
                || fabs(fluxhi.momentum_normal - fluxlo.momentum_normal) > 1E-10 // no change in momentum flux
                || fabs(fluxhi.momentum_tangent) > 1E-10                         // zero tangent flux
                || fabs(fluxlo.momentum_tangent) > 1E-10                         // zero tangent flux
                || fabs(fluxhi.energy - fluxlo.energy) > 1E-10                   // no change in energy flux
            )
            {
                Util::Warning(INFO, "left:  ", left);
                Util::Warning(INFO, "right: ", right);
                Util::Warning(INFO, "Fluxhi:  ", fluxhi);
                Util::Warning(INFO, "Fluxlo:  ", fluxlo);
                Util::Exception(INFO, "Symmetric states should result in zero flux.");
            }
            Util::Test::SubMessage("Test 3: Constant states induces no flux difference", 0);
        }
        catch (const std::runtime_error &e)
        {
            failed++;
            Util::Test::SubMessage("Test 3: Constant states induces no flux difference", 1);
        }

        // Test 4: Uniform Flow Test (no flux across uniform flow)
        try
        {
            State left(1.0, 1.0, 0.5, 1.0);
            State center(1.0, 1.0, 0.5, 1.0);
            State right(1.0, 1.0, 0.5, 1.0);
            Flux fluxhi = solver.Solve(center, right, gamma, pref, small);
            Flux fluxlo = solver.Solve(left, center, gamma, pref, small);
            if (fabs(fluxhi.mass - fluxlo.mass) > 1E-10 || fabs(fluxhi.momentum_normal - fluxlo.momentum_normal) > 1E-10 || fabs(fluxhi.energy - fluxlo.energy) > 1E-10)
            {
                Util::Warning(INFO, "left:   ", left);
                Util::Warning(INFO, "center: ", center);
                Util::Warning(INFO, "right:  ", right);
                Util::Warning(INFO, "Fluxlo: ", fluxlo);
                Util::Warning(INFO, "Fluxhi: ", fluxhi);
                Util::Exception(INFO, "Uniform flow should result in no flux.");
            }
            Util::Test::SubMessage("Test 4: Uniform flow should maintain constant flux", 0);
        }
        catch (const std::runtime_error &e)
        {
            failed++;
            Util::Test::SubMessage("Test 4: Uniform flow should maintain constant flux", 1);
        }

        return failed;
    }
};
}
}
}
#endif
