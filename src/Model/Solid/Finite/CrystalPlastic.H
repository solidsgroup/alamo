#ifndef MODEL_SOLID_FINITE_CRYSTALPLASTIC_H_
#define MODEL_SOLID_FINITE_CRYSTALPLASTIC_H_

#include <cmath>
#include <cstdio>
#include <algorithm>
#include </usr/include/eigen3/Eigen/Core>
#include </usr/include/eigen3/Eigen/LU>
#include </usr/include/eigen3/Eigen/SVD>
#include "IO/ParmParse.H"
#include "Model/Solid/Solid.H"
#include "PseudoLinear/Cubic.H"
namespace Model {
namespace Solid {
namespace Finite {

struct DragResult {
    Set::Scalar v_r;
    Set::Scalar dvr_dtau;
};

//
// Relativistic phonon-drag velocity.
//
// Physics:  B(v) = B0 / (1 - v^2/cs^2)   (drag DIVERGES as v -> cs)
//           tau* b = B(v) v
//   =>      tau* b (1 - u^2) = B0 u cs,     u = v/cs
//   =>      eta (1 - u^2) = u,              eta = tau* b/(B0 cs) = v_nr/cs
//
// This is a QUADRATIC in u:   eta u^2 + u - eta = 0
//   =>  u = ( -1 + sqrt(1 + 4 eta^2) ) / (2 eta)
// with the correct limits  u -> eta  (eta -> 0)  and  u -> 1  (eta -> inf).
//
// ### OLD: the previous implementation solved the TRANSPOSED equation
// ###          u (1 - u^2) = eta      (i.e. eta and u swapped)
// ### via the trigonometric depressed-cubic form u = (2/sqrt3) cos(theta/3 - 2pi/3),
// ### theta = acos(-1.5 sqrt3 eta), which is the root of  u^3 - u + eta = 0.
// ### That equation implies B(v) = B0 (1 - u^2), i.e. drag VANISHING at v -> cs,
// ### and it has no real ascending-branch root past eta_max = 2/(3 sqrt3) = 0.385,
// ### so the old code clamped u at u_max = 1/sqrt3 = 0.577.
// ###
// ### Three symptoms of the transposition, all now gone:
// ###   (1) the correct law is not a cubic at all;
// ###   (2) the old small-eta branch read u = eta(1 + eta^2) = eta + eta^3, whereas
// ###       eta(1-u^2) = u gives u = eta(1 - eta^2) = eta - eta^3 (sign flip);
// ###   (3) du/deta = 1/(1 - 3u^2) is singular at u = 1/sqrt3, which is why the old
// ###       code needed the 0.95*u_max derivative clamp.  The correct derivative
// ###       du/deta = (1-u^2)/(1 + 2 eta u) is regular for all u < 1, so no clamp.
// ###
// ### Impact: for Ta (b=2.86e-10 m, B0=2.96e-5 Pa.s, cs=1937 m/s) eta exceeds the
// ### old eta_max at tau* ~ 0.077 GPa, but slip does not activate until tau* ~ 0.57
// ### GPa -- so u was pinned at u_max for EVERY active slip event, giving a constant
// ### v_r = 1118 m/s.  The correct law gives 1641 m/s at tau*=0.6 GPa and 1753 m/s at
// ### tau*=1.0 GPa.  Since t_r >> t_w in that regime (5.2e-10 s vs 1.7e-11 s), drag is
// ### rate-controlling and gdot = b rho_M v_mean ~ b rho_M v_r inherited the full
// ### 32-36% deficit at the shock front.
//
AMREX_GPU_HOST_DEVICE
static inline DragResult
relativistic_drag_velocity(Set::Scalar tau_star_Pa,
                           Set::Scalar b_m,
                           Set::Scalar B0_Pas,
                           Set::Scalar cs_m)
{
    DragResult out = {0.0, 0.0};
    if (tau_star_Pa <= 0.0 || B0_Pas <= 0.0 || cs_m <= 0.0) return out;

    const Set::Scalar v_nr = (b_m / B0_Pas) * tau_star_Pa;
    const Set::Scalar eta  = v_nr / cs_m;

    Set::Scalar u;
    if (eta < 1.0e-6)
        u = eta * (1.0 - eta * eta);          // series root of eta(1-u^2) = u
    else
        u = (std::sqrt(1.0 + 4.0 * eta * eta) - 1.0) / (2.0 * eta);

    if (u < 0.0) u = 0.0;
    if (u > 1.0 - 1.0e-12) u = 1.0 - 1.0e-12; // u -> 1 asymptotically; never reached

    out.v_r = u * cs_m;

    // dv_r/dtau* = cs * (du/deta) * (deta/dtau*) = (b/B0) * (1-u^2)/(1 + 2 eta u).
    // Regular everywhere on u in [0,1): -> b/B0 as eta->0, -> 0 as eta->inf.
    out.dvr_dtau = (b_m / B0_Pas) * (1.0 - u * u) / (1.0 + 2.0 * eta * u);

    return out;
}

class CrystalPlastic : public PseudoLinear::Cubic
{
public:
    static constexpr int NSYS = 24;
    static constexpr Set::Scalar kB = 1.381e-23;

    static constexpr Set::Scalar RHO_MAX = 1.0e12;

    CrystalPlastic() {};
    virtual ~CrystalPlastic() {};

    static Eigen::Matrix3d expand_F(const Set::Matrix &F2d)
    {
#if AMREX_SPACEDIM == 2
        Eigen::Matrix3d F3 = Eigen::Matrix3d::Identity();
        F3(0,0) = F2d(0,0); F3(0,1) = F2d(0,1);
        F3(1,0) = F2d(1,0); F3(1,1) = F2d(1,1);
        return F3;
#else
        return F2d;
#endif
    }

    void build_CC3d(Set::Scalar CC[3][3][3][3]) const
    {
        Eigen::Matrix3d R = q.normalized().toRotationMatrix();
        Set::Scalar Ctmp[3][3][3][3];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    for (int l = 0; l < 3; l++) {
                        if (i==j && j==k && k==l) Ctmp[i][j][k][l] = C11;
                        else if (i==k && j==l) Ctmp[i][j][k][l] = C44;
                        else if (i==l && j==k) Ctmp[i][j][k][l] = C44;
                        else if (i==j && k==l) Ctmp[i][j][k][l] = C12;
                        else Ctmp[i][j][k][l] = 0.0;
                    }
        for (int p = 0; p < 3; p++)
            for (int q_ = 0; q_ < 3; q_++)
                for (int s = 0; s < 3; s++)
                    for (int t = 0; t < 3; t++) {
                        CC[p][q_][s][t] = 0.0;
                        for (int i = 0; i < 3; i++)
                            for (int j = 0; j < 3; j++)
                                for (int k = 0; k < 3; k++)
                                    for (int l = 0; l < 3; l++)
                                        CC[p][q_][s][t] += R(p,i)*R(q_,j)*Ctmp[i][j][k][l]*R(s,k)*R(t,l);
                    }
    }

    Set::Scalar W(const Set::Matrix &F) const
    {
        Eigen::Matrix3d F3 = expand_F(F);
        Eigen::Matrix3d Fpinv = Fp.inverse();
        Eigen::Matrix3d Fe = F3 * Fpinv;
        Eigen::Matrix3d Ce = Fe.transpose() * Fe;
        Eigen::Matrix3d Ee = 0.5 * (Ce - Eigen::Matrix3d::Identity());

        Set::Scalar CC[3][3][3][3];
        build_CC3d(CC);
        Set::Scalar w = 0.0;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                for (int k = 0; k < 3; k++)
                    for (int l = 0; l < 3; l++)
                        w += 0.5 * Ee(i,j) * CC[i][j][k][l] * Ee(k,l);
        return w;
    }

    Set::Matrix DW(const Set::Matrix &F) const
    {
        Eigen::Matrix3d F3 = expand_F(F);
        Eigen::Matrix3d Fpinv = Fp.inverse();
        Eigen::Matrix3d Fe = F3 * Fpinv;
        Eigen::Matrix3d Ce = Fe.transpose() * Fe;

        Set::Scalar CC[3][3][3][3];
        build_CC3d(CC);

        Eigen::Matrix3d Te = Eigen::Matrix3d::Zero();
        for (int I = 0; I < 3; I++)
            for (int J = 0; J < 3; J++)
                for (int K = 0; K < 3; K++)
                    for (int L = 0; L < 3; L++)
                        Te(I,J) += 0.5 * CC[I][J][K][L] * (Ce(K,L) - (K==L ? 1.0 : 0.0));

        Eigen::Matrix3d FpinvT = Fpinv.transpose();
        Eigen::Matrix3d P3 = Fe * Te * FpinvT;

#if AMREX_SPACEDIM == 2
        Set::Matrix P2;
        P2(0,0) = P3(0,0); P2(0,1) = P3(0,1);
        P2(1,0) = P3(1,0); P2(1,1) = P3(1,1);
        return P2;
#else
        return P3;
#endif
    }

    Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major> DDW(const Set::Matrix &F) const
    {
        Eigen::Matrix3d F3    = expand_F(F);
        Eigen::Matrix3d Fpinv = Fp.inverse();
        Eigen::Matrix3d Fe    = F3 * Fpinv;
        Eigen::Matrix3d Ce    = Fe.transpose() * Fe;
        Eigen::Matrix3d I3    = Eigen::Matrix3d::Identity();

        Set::Scalar CC[3][3][3][3];
        build_CC3d(CC);

        Set::Scalar ddw3[3][3][3][3];
        for (int i = 0; i < 3; i++)
        for (int J = 0; J < 3; J++)
        for (int k = 0; k < 3; k++)
        for (int L = 0; L < 3; L++)
        {
            Set::Scalar acc = 0.0;
            for (int P = 0; P < 3; P++)
            for (int R = 0; R < 3; R++)
            {
                Set::Scalar term1 = 0.0;
                for (int A = 0; A < 3; A++)
                    for (int B = 0; B < 3; B++)
                        term1 += CC[P][A][B][R] * 0.5 * (Ce(A,B) - I3(A,B));

                Set::Scalar term2 = 0.0;
                for (int A = 0; A < 3; A++)
                    for (int B = 0; B < 3; B++)
                        term2 += Fe(i,A) * Fe(k,B) * CC[A][P][B][R];

                acc += Fpinv(J,P) * Fpinv(L,R)
                     * ((i == k ? 1.0 : 0.0) * term1 + term2);
            }
            ddw3[i][J][k][L] = acc;
        }

        auto ret = Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Major>::Zero();
        for (int i = 0; i < AMREX_SPACEDIM; i++)
            for (int J = 0; J < AMREX_SPACEDIM; J++)
                for (int k = 0; k < AMREX_SPACEDIM; k++)
                    for (int L = 0; L < AMREX_SPACEDIM; L++)
                        ret(i,J,k,L) = ddw3[i][J][k][L];
        return ret;
    }

    static void build_slip_systems(const Eigen::Quaterniond &q_rot,
                                   Set::Vector3d m0[NSYS],
                                   Set::Vector3d n0[NSYS],
                                   Set::Scalar S0_3d[NSYS][3][3],
                                   int fam[NSYS])
    {
        Eigen::Quaterniond qn = q_rot.normalized();
        Set::Matrix3d Rmat = qn.toRotationMatrix();
        const Set::Vector3d a111(1,1,1), b111(-1,1,1), c111(1,-1,1), d111(-1,-1,1);

        m0[0]=a111; n0[0]=Set::Vector3d( 0,-1, 1);
        m0[1]=a111; n0[1]=Set::Vector3d(-1, 0, 1);
        m0[2]=a111; n0[2]=Set::Vector3d( 1,-1, 0);
        m0[3]=b111; n0[3]=Set::Vector3d( 0,-1, 1);
        m0[4]=b111; n0[4]=Set::Vector3d( 1, 0, 1);
        m0[5]=b111; n0[5]=Set::Vector3d( 1, 1, 0);
        m0[6]=c111; n0[6]=Set::Vector3d( 0, 1, 1);
        m0[7]=c111; n0[7]=Set::Vector3d(-1, 0, 1);
        m0[8]=c111; n0[8]=Set::Vector3d( 1, 1, 0);
        m0[9]=d111; n0[9]=Set::Vector3d( 0, 1, 1);
        m0[10]=d111; n0[10]=Set::Vector3d( 1, 0, 1);
        m0[11]=d111; n0[11]=Set::Vector3d(-1, 1, 0);
        for (int a = 0; a < 12; a++) fam[a] = 0;

        m0[12]=a111; n0[12]=Set::Vector3d( 2,-1,-1);
        m0[13]=a111; n0[13]=Set::Vector3d( 1, 1,-2);
        m0[14]=a111; n0[14]=Set::Vector3d( 1,-2, 1);
        m0[15]=b111; n0[15]=Set::Vector3d( 2, 1, 1);
        m0[16]=b111; n0[16]=Set::Vector3d( 1,-1, 2);
        m0[17]=b111; n0[17]=Set::Vector3d( 1, 2,-1);
        m0[18]=c111; n0[18]=Set::Vector3d( 1, 2, 1);
        m0[19]=c111; n0[19]=Set::Vector3d( 1,-1,-2);
        m0[20]=c111; n0[20]=Set::Vector3d( 2, 1,-1);
        m0[21]=d111; n0[21]=Set::Vector3d( 1, 1, 2);
        m0[22]=d111; n0[22]=Set::Vector3d( 1,-2,-1);
        m0[23]=d111; n0[23]=Set::Vector3d( 2,-1, 1);
        for (int a = 12; a < 24; a++) fam[a] = 1;

        for (int a = 0; a < NSYS; a++) {
            Set::Vector3d mv = Rmat * m0[a].normalized();
            Set::Vector3d nv = Rmat * n0[a].normalized();
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    S0_3d[a][i][j] = mv[i] * nv[j];
        }
    }

    static void build_forest_coplanar(const Set::Vector3d m0[NSYS],
                                      const Set::Vector3d n0[NSYS],
                                      Set::Scalar Afor[NSYS][NSYS],
                                      Set::Scalar Acop[NSYS][NSYS])
    {
        for (int a = 0; a < NSYS; a++) {
            Set::Vector3d na = n0[a].normalized();
            for (int b = 0; b < NSYS; b++) {
                Set::Vector3d sb = m0[b].normalized();
                Set::Vector3d nb = n0[b].normalized();
                Set::Vector3d nb_x_sb = nb.cross(sb);
                Afor[a][b] = 0.5 * std::fabs(na.dot(sb))
                           + 0.5 * std::fabs(na.dot(nb_x_sb));
                Acop[a][b] = 0.5 * na.cross(sb).norm()
                           + 0.5 * na.cross(nb_x_sb).norm();
            }
        }
    }

    void enforce_Fp_incompressibility()
    {
        Set::Scalar dJp = Fp.determinant();
        if (!std::isfinite(dJp) || dJp < 0.1 || dJp > 10.0) {

            fprintf(stderr,
                "[CP] CRITICAL: det(Fp)=%g out of [0.1,10], "
                "resetting Fp to identity (plastic history at this node lost)\n",
                dJp);
            Fp = Set::Matrix3d::Identity();
            return;
        }
        Set::Scalar scale = std::pow(dJp, -1.0/3.0);
        Fp *= scale;
    }

    void floor_and_cap_densities()
    {
        Set::Scalar rhoM_floor = mob_frac0 * rho0_total / NSYS;
        Set::Scalar rhoI_floor = (1.0 - mob_frac0) * rho0_total / NSYS;
        for (int a = 0; a < NSYS; a++) {
            if (!std::isfinite(rho_M[a]) || rho_M[a] < rhoM_floor)
                rho_M[a] = rhoM_floor;
            if (!std::isfinite(rho_I[a]) || rho_I[a] < rhoI_floor)
                rho_I[a] = rhoI_floor;
            if (rho_M[a] > RHO_MAX) rho_M[a] = RHO_MAX;
            if (rho_I[a] > RHO_MAX) rho_I[a] = RHO_MAX;
        }
    }

    //
    // NOTE ON wp_inc / edot_p_eq SEMANTICS:
    //   wp_inc    is an INCREMENT (plastic work over THIS call, GPa).
    //   edot_p_eq is a RATE (equivalent plastic strain rate over THIS call, 1/s).
    // Both are reset at entry, so a caller that substeps Advance() must SUM wp_inc
    // across substeps to obtain the step total.  edot_p_eq, being a rate, is
    // meaningful as-is from the final substep.  See
    // LagrangianShockCP::AdvanceCPModelAllNodes.
    //
    virtual void Advance(Set::Scalar dt, Set::Matrix F_in, Set::Matrix ,
                         Set::Scalar time)
    {
        //
        // ### FIX: `static int adv_count = 0;` was a function-local static
        // ### incremented from inside amrex::ParallelFor -- a non-atomic
        // ### read-modify-write shared across OpenMP threads (data race), and
        // ### unusable on device.  It fed only the debug fprintf below.  Removed.
        //
        wp_inc = 0.0;
        edot_p_eq = 0.0;

        if (time < tstart) return;

        floor_and_cap_densities();

        Set::Vector3d m0[NSYS], n0[NSYS];
        Set::Scalar S0_3d[NSYS][3][3];
        int fam[NSYS];
        build_slip_systems(q, m0, n0, S0_3d, fam);

#ifndef NDEBUG

        {
            Set::Scalar max_mn = 0.0;
            for (int a = 0; a < NSYS; a++) {
                Set::Scalar mn = 0.0;
                for (int i = 0; i < 3; i++) mn += S0_3d[a][i][i];
                max_mn = std::max(max_mn, std::fabs(mn));
            }
            AMREX_ASSERT_WITH_MESSAGE(max_mn < 1e-10,
                "CrystalPlastic: Schmid tensor not traceless — m and n are not "
                "orthogonal.  Check quaternion q and slip system construction.");
        }
#endif

        Eigen::Matrix3d F3 = expand_F(F_in);

        if (use_implicit > 0.5)
            AdvanceNguyenImplicit3D(dt, F3, time, m0, n0, S0_3d, fam, 0);
        else
            AdvanceNguyen3D(dt, F3, time, m0, n0, S0_3d, fam, 0);

        enforce_Fp_incompressibility();

#ifndef NDEBUG
        {
            Eigen::JacobiSVD<Eigen::Matrix3d> svd(Fp);
            Set::Scalar smax = svd.singularValues()(0);
            Set::Scalar smin = svd.singularValues()(2);
            Set::Scalar cond = (smin > 1e-30) ? (smax / smin) : 1.0e30;
            Set::Scalar detFp = Fp.determinant();
            if (!std::isfinite(detFp) || cond > 1.0e6) {
                fprintf(stderr,
                    "[CP-Fp] WARN det(Fp)=%g cond(Fp)=%g smin=%g smax=%g\n",
                    detFp, cond, smin, smax);
            }
            AMREX_ASSERT_WITH_MESSAGE(std::isfinite(detFp) && cond < 1.0e10,
                "CrystalPlastic: Fp ill-conditioned, would produce NaN downstream");
        }
#endif
    }

    void AdvanceNguyen3D(Set::Scalar dt, const Eigen::Matrix3d &F3,
                         Set::Scalar ,
                         const Set::Vector3d m0[NSYS],
                         const Set::Vector3d n0[NSYS],
                         const Set::Scalar S0_3d[NSYS][3][3],
                         const int fam[NSYS],
                         int )
    {
        Set::Scalar mu_star = std::sqrt(C44 * (C11 - C12) / 2.0);

        constexpr Set::Scalar C44_0K = 87.1;
        constexpr Set::Scalar C11_0K = 268.5;
        constexpr Set::Scalar C12_0K = 159.9;
        const Set::Scalar mu_0K = std::sqrt(C44_0K * (C11_0K - C12_0K) / 2.0);
        Set::Scalar mu_ratio = mu_star / mu_0K;

        Set::Scalar Afor_loc[NSYS][NSYS], Acop_loc[NSYS][NSYS];
        build_forest_coplanar(m0, n0, Afor_loc, Acop_loc);

        Eigen::Matrix3d Fpinv = Fp.inverse();
        Eigen::Matrix3d Fe = F3 * Fpinv;
        Eigen::Matrix3d Ce = Fe.transpose() * Fe;
        Eigen::Matrix3d Ee;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Ee(i,j) = 0.5 * (Ce(i,j) - (i==j ? 1.0 : 0.0));

        Set::Scalar CC[3][3][3][3];
        build_CC3d(CC);

        Set::Scalar Te[3][3];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {
                Te[i][j] = 0.0;
                for (int k = 0; k < 3; k++)
                    for (int l = 0; l < 3; l++)
                        Te[i][j] += CC[i][j][k][l] * Ee(k,l);
            }

        Set::Scalar tau[NSYS];
        for (int a = 0; a < NSYS; a++) {
            tau[a] = 0.0;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    tau[a] += Te[i][j] * S0_3d[a][i][j];
        }

        Set::Scalar rho_total[NSYS], rho_for[NSYS], rho_cop[NSYS];
        for (int a = 0; a < NSYS; a++)
            rho_total[a] = rho_M[a] + rho_I[a];
        for (int a = 0; a < NSYS; a++) {
            rho_for[a] = 0.0; rho_cop[a] = 0.0;
            for (int b = 0; b < NSYS; b++) {
                rho_for[a] += Afor_loc[a][b] * rho_total[b];
                rho_cop[a] += Acop_loc[a][b] * rho_total[b];
            }
        }

        Set::Scalar dg[NSYS];
        for (int a = 0; a < NSYS; a++) dg[a] = 0.0;

        for (int a = 0; a < NSYS; a++) {
            Set::Scalar tau_a_val = tau_a0 + Cfor * burg * mu_star
                                  * std::sqrt(std::max(rho_for[a], 0.0));
            Set::Scalar tau_p;
            if (fam[a] == 0)
                tau_p = tau_p_110;
            else
                tau_p = (tau[a] >= 0.0) ? tau_p_112T : tau_p_112AT;

            Set::Scalar tau_r_val;
            if (use_modscaled_taur) {
                tau_r_val = tau_p * mu_ratio;
            } else {
                tau_r_val = tau_p + Ccop * burg * mu_star
                          * std::sqrt(std::max(rho_cop[a], 0.0));
            }

            Set::Scalar tau_star = std::max(std::fabs(tau[a]) - tau_a_val, 0.0);
            if (tau_star <= 0.0 || tau_r_val <= 0.0) continue;
            Set::Scalar sgn = (tau[a] >= 0.0) ? 1.0 : -1.0;
            Set::Scalar Lbar = 1.0 / std::sqrt(std::max(rho_cop[a], 1.0));

            Set::Scalar r = std::min(tau_star / tau_r_val, 1.0 - 1e-12);
            Set::Scalar pw_p = std::pow(r, p_exp);
            Set::Scalar pw_q = std::pow(std::max(1.0 - pw_p, 1e-30), q_exp);
            Set::Scalar DG_eff = DeltaG0_nguyen * pw_q;
            Set::Scalar exp_arg = std::min(DG_eff / (kB * theta), 60.0);
            Set::Scalar t_w = (1.0 / fD) * (std::exp(exp_arg) - 1.0);
            t_w = std::max(t_w, 1e-30);

            const Set::Scalar tau_star_Pa = tau_star * 1.0e9;
            const Set::Scalar b_m         = burg * 1.0e-3;
            const Set::Scalar cs_m        = cs   * 1.0e-3;
            DragResult dr = relativistic_drag_velocity(
                tau_star_Pa, b_m, B0, cs_m);
            Set::Scalar v_r = dr.v_r * 1.0e3;
            Set::Scalar t_r = (v_r > 1e-30) ? Lbar / v_r : 1e30;

            Set::Scalar v_mean = Lbar / (t_w + t_r);
            Set::Scalar rhoM_a = std::max(rho_M[a], 1e-10);
            Set::Scalar gdot_a = burg * rhoM_a * v_mean;
            dg[a] = dt * gdot_a * sgn;
            if (std::fabs(dg[a]) > 0.01) dg[a] = 0.01 * sgn;

            wp_inc += tau[a] * dg[a];
        }

        Set::Scalar Lp[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
        Set::Scalar total_dg = 0.0;
        for (int a = 0; a < NSYS; a++) {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Lp[i][j] += dg[a] * S0_3d[a][i][j];
            total_dg += std::fabs(dg[a]);
        }
        if (!std::isfinite(total_dg) || total_dg > 5.0) return;

        // Equivalent plastic strain rate edot_p_eq = sqrt(2/3 D^p:D^p),
        // with D^p = sym(L^p) 
        {
            Set::Scalar Dp[3][3];
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Dp[i][j] = 0.5 * (Lp[i][j] + Lp[j][i]);

            Set::Scalar DpDp = 0.0;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    DpDp += Dp[i][j] * Dp[i][j];

            Set::Scalar dep_eq_step =
                std::sqrt(std::max((2.0 / 3.0) * DpDp, 0.0));
            edot_p_eq = (dt > 1e-30) ? (dep_eq_step / dt) : 0.0;
            if (!std::isfinite(edot_p_eq) || edot_p_eq < 0.0)
                edot_p_eq = 0.0;
        }

        Eigen::Matrix3d Fp_old = Fp;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {
                Set::Scalar val = 0.0;
                for (int k = 0; k < 3; k++)
                    val += ((i==k ? 1.0 : 0.0) + Lp[i][k]) * Fp_old(k,j);
                Fp(i,j) = val;
            }
        for (int a = 0; a < NSYS; a++)
            gamma[a] += dg[a];

        Set::Scalar CT_val = CT_ratio * CM;
        Set::Scalar CA_val = 4.0 * (CM - CT_val);

        for (int a = 0; a < NSYS; a++) {
            Set::Scalar actual_gdot = std::fabs(dg[a]) / std::max(dt, 1e-30);
            if (actual_gdot < 1e-20) continue;

            Set::Scalar rhoM_a = std::max(rho_M[a], 1e-10);
            Set::Scalar rhoI_a = std::max(rho_I[a], 1e-10);
            Set::Scalar sqrt_rf = std::sqrt(std::max(rho_for[a], 1e-10));
            Set::Scalar va = actual_gdot / (burg * rhoM_a);

            Set::Scalar gdot_ratio = actual_gdot / std::max(gdot0_sat, 1e-30);
            Set::Scalar exponent = kB * theta / std::max(A_sat, 1e-30);
            Set::Scalar rhoI_sat = rho0_sat
                * std::pow(std::max(gdot_ratio, 1e-30), exponent);
            Set::Scalar rhoI_annealed = (1.0 - mob_frac0) * rho0_total / NSYS;
            rhoI_sat = std::max(rhoI_sat, rhoI_annealed);
            Set::Scalar rhoM_sat = varsigma * rhoI_sat;
            Set::Scalar d_A = sqrt_rf / std::max(rhoM_sat, 1e-10);

            Set::Scalar rate_mult = CM * sqrt_rf * rhoM_a * va;
            Set::Scalar rate_ann  = 0.25 * CA_val * d_A * rhoM_a * rhoM_a * va;
            Set::Scalar rate_trap = CT_val * sqrt_rf * rhoM_a * va;

            Set::Scalar rate_nuc = 0.0;
            if (enable_nucleation) {
                Set::Scalar tau_a_val = tau_a0 + Cfor * burg * mu_star
                                      * std::sqrt(std::max(rho_for[a], 0.0));
                Set::Scalar tau_eff_nuc = std::max(std::fabs(tau[a]) - tau_a_val, 0.0);
                if (tau_eff_nuc > 0.0 && tau_n0 > 0.0) {
                    Set::Scalar r_n = std::min(tau_eff_nuc / tau_n0, 1.0 - 1e-12);
                    Set::Scalar pw_pn = std::pow(r_n, pn_exp);
                    Set::Scalar pw_qn = std::pow(std::max(1.0 - pw_pn, 1e-30), qn_exp);
                    Set::Scalar DG_nuc = DeltaGn0 * pw_qn;
                    Set::Scalar nuc_arg = std::min(DG_nuc / (kB * theta), 60.0);
                    rate_nuc = rhodot_n0 * std::exp(-nuc_arg);
                }
            }

            Set::Scalar rhoI_ratio = std::min(rhoI_a / std::max(rhoI_sat, 1e-10), 1.0);
            Set::Scalar f_rec = 1.0 - std::pow(rhoI_ratio, 1.0 / nrec);

            Set::Scalar rhoM_new = rhoM_a + dt * (rate_mult - rate_ann - rate_trap + rate_nuc);
            Set::Scalar rhoI_new = rhoI_a + dt * f_rec * rate_trap;

            Set::Scalar drhoM = rhoM_new - rhoM_a;
            if (std::fabs(drhoM) > 0.05 * rhoM_a)
                drhoM = 0.05 * rhoM_a * (drhoM > 0 ? 1.0 : -1.0);
            Set::Scalar drhoI = rhoI_new - rhoI_a;
            if (std::fabs(drhoI) > 0.05 * rhoI_a)
                drhoI = 0.05 * rhoI_a * (drhoI > 0 ? 1.0 : -1.0);

            rho_M[a] = std::max(rhoM_a + drhoM, 1e-10);
            rho_I[a] = std::max(rhoI_a + drhoI, 1e-10);

            if (rho_M[a] > RHO_MAX) rho_M[a] = RHO_MAX;
            if (rho_I[a] > RHO_MAX) rho_I[a] = RHO_MAX;
            rho_vec[a] = rho_M[a] + rho_I[a];
        }
    }

    void AdvanceNguyenImplicit3D(Set::Scalar dt, const Eigen::Matrix3d &F3,
                                Set::Scalar ,
                                const Set::Vector3d m0[NSYS],
                                const Set::Vector3d n0[NSYS],
                                const Set::Scalar S0_3d[NSYS][3][3],
                                const int fam[NSYS],
                                int )
    {
        Set::Scalar mu_star = std::sqrt(C44 * (C11 - C12) / 2.0);

        constexpr Set::Scalar C44_0K = 87.1;
        constexpr Set::Scalar C11_0K = 268.5;
        constexpr Set::Scalar C12_0K = 159.9;
        const Set::Scalar mu_0K = std::sqrt(C44_0K * (C11_0K - C12_0K) / 2.0);
        Set::Scalar mu_ratio = mu_star / mu_0K;

        Set::Scalar Afor_loc[NSYS][NSYS], Acop_loc[NSYS][NSYS];
        build_forest_coplanar(m0, n0, Afor_loc, Acop_loc);

        Set::Scalar CC[3][3][3][3];
        build_CC3d(CC);

        Eigen::Matrix3d Fpinv = Fp.inverse();
        Eigen::Matrix3d Fe_tr = F3 * Fpinv;
        Eigen::Matrix3d Ce_tr = Fe_tr.transpose() * Fe_tr;
        Set::Scalar Ee[3][3];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Ee[i][j] = 0.5 * (Ce_tr(i,j) - (i==j ? 1.0 : 0.0));

        Set::Scalar Te_tr[3][3];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {
                Te_tr[i][j] = 0.0;
                for (int k = 0; k < 3; k++)
                    for (int l = 0; l < 3; l++)
                        Te_tr[i][j] += CC[i][j][k][l] * Ee[k][l];
            }

        Set::Scalar Ca[NSYS][3][3];
        for (int a = 0; a < NSYS; a++) {
            Set::Scalar CeS[3][3];
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++) {
                    CeS[i][j] = 0.0;
                    for (int k = 0; k < 3; k++)
                        CeS[i][j] += Ce_tr(i,k) * S0_3d[a][k][j];
                }
            Set::Scalar Bs[3][3];
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Bs[i][j] = 0.25 * (CeS[i][j] + CeS[j][i]);
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++) {
                    Ca[a][i][j] = 0.0;
                    for (int k = 0; k < 3; k++)
                        for (int l = 0; l < 3; l++)
                            Ca[a][i][j] += CC[i][j][k][l] * Bs[k][l];
                }
        }

        Set::Scalar rhoM_w[NSYS], rhoI_w[NSYS];
        for (int a = 0; a < NSYS; a++) {
            rhoM_w[a] = rho_M[a];
            rhoI_w[a] = rho_I[a];
        }

        const int MAX_OUTER = 15;
        const int MAX_INNER = 50;
        const Set::Scalar tol_Te = 1.0e-8;
        const Set::Scalar tol_rho = 1.0e-4;
        const Set::Scalar ETA_DAMP = 0.25;

        Set::Scalar Te[3][3], dg[NSYS];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Te[i][j] = Te_tr[i][j];
        for (int a = 0; a < NSYS; a++) dg[a] = 0.0;

        for (int outer = 0; outer < MAX_OUTER; outer++)
        {

            Set::Scalar rho_total_w[NSYS], rho_for_w[NSYS], rho_cop_w[NSYS];
            for (int a = 0; a < NSYS; a++)
                rho_total_w[a] = rhoM_w[a] + rhoI_w[a];
            for (int a = 0; a < NSYS; a++) {
                rho_for_w[a] = 0.0; rho_cop_w[a] = 0.0;
                for (int b = 0; b < NSYS; b++) {
                    rho_for_w[a] += Afor_loc[a][b] * rho_total_w[b];
                    rho_cop_w[a] += Acop_loc[a][b] * rho_total_w[b];
                }
            }

            Set::Scalar tau_a_w[NSYS], tau_r_w[NSYS];
            for (int a = 0; a < NSYS; a++) {
                tau_a_w[a] = tau_a0 + Cfor * burg * mu_star
                           * std::sqrt(std::max(rho_for_w[a], 0.0));
                Set::Scalar tau_p;
                if (fam[a] == 0) tau_p = tau_p_110;
                else tau_p = tau_p_112T;
                if (use_modscaled_taur)
                    tau_r_w[a] = tau_p * mu_ratio;
                else
                    tau_r_w[a] = tau_p + Ccop * burg * mu_star
                               * std::sqrt(std::max(rho_cop_w[a], 0.0));
            }

            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Te[i][j] = Te_tr[i][j];

            for (int inner = 0; inner < MAX_INNER; inner++)
            {
                Set::Scalar tau[NSYS];
                for (int a = 0; a < NSYS; a++) {
                    tau[a] = 0.0;
                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                            tau[a] += Te[i][j] * S0_3d[a][i][j];
                }

                Set::Scalar ddg_dtau[NSYS];
                Set::Scalar max_abs_dg = 0.0;
                for (int a = 0; a < NSYS; a++) {
                    Set::Scalar tau_r_eff = tau_r_w[a];
                    if (fam[a] == 1)
                        tau_r_eff = (tau[a] >= 0.0) ?
                            (use_modscaled_taur ? tau_p_112T * mu_ratio
                                : tau_p_112T + Ccop*burg*mu_star*std::sqrt(std::max(rho_cop_w[a],0.0)))
                          : (use_modscaled_taur ? tau_p_112AT * mu_ratio
                                : tau_p_112AT + Ccop*burg*mu_star*std::sqrt(std::max(rho_cop_w[a],0.0)));

                    SlipResult sr = compute_slip_and_deriv(
                        dt, tau[a], tau_a_w[a], tau_r_eff,
                        rhoM_w[a], rho_cop_w[a], mu_star);
                    dg[a] = sr.dg;
                    ddg_dtau[a] = sr.ddg_dtau;
                    max_abs_dg = std::max(max_abs_dg, std::fabs(dg[a]));
                }

                Set::Scalar G[3][3];
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++) {
                        G[i][j] = Te[i][j] - Te_tr[i][j];
                        for (int a = 0; a < NSYS; a++)
                            G[i][j] += dg[a] * Ca[a][i][j];
                    }

                if (max_abs_dg > 0.5) {
                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                            Te[i][j] -= ETA_DAMP * G[i][j];
                    continue;
                }

                Set::Scalar max_G = 0.0;
                for (int i = 0; i < 3; i++)
                    for (int j = 0; j < 3; j++)
                        max_G = std::max(max_G, std::fabs(G[i][j]));
                if (max_G < tol_Te) break;

                Eigen::Matrix<Set::Scalar,9,9> J_eig = Eigen::Matrix<Set::Scalar,9,9>::Identity();
                Eigen::Matrix<Set::Scalar,9,1> G_vec;

                for (int I = 0; I < 9; I++)
                    G_vec[I] = G[I/3][I%3];

                for (int a = 0; a < NSYS; a++) {
                    if (std::fabs(ddg_dtau[a]) < 1e-30) continue;
                    Set::Scalar Ss[3][3];
                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                            Ss[i][j] = 0.5*(S0_3d[a][i][j] + S0_3d[a][j][i]);
                    for (int ii=0; ii<3; ii++) for (int jj=0; jj<3; jj++)
                    for (int kk=0; kk<3; kk++) for (int ll=0; ll<3; ll++)
                        J_eig(3*ii+jj, 3*kk+ll) += Ca[a][ii][jj] * ddg_dtau[a] * Ss[kk][ll];
                }

                Eigen::Matrix<Set::Scalar,9,1> corr = J_eig.partialPivLu().solve(-G_vec);
                for (int I = 0; I < 9; I++)
                    Te[I/3][I%3] += corr[I];
            }

            Set::Scalar rhoM_old[NSYS];
            for (int a = 0; a < NSYS; a++) {
                rhoM_old[a] = rhoM_w[a];
            }

            Set::Scalar CT_val = CT_ratio * CM;
            Set::Scalar CA_val = 4.0 * (CM - CT_val);

            for (int a = 0; a < NSYS; a++) {
                Set::Scalar actual_gdot = std::fabs(dg[a]) / std::max(dt, 1e-30);
                if (actual_gdot < 1e-20) continue;

                Set::Scalar rhoM_a = std::max(rhoM_w[a], 1e-10);
                Set::Scalar rhoI_a = std::max(rhoI_w[a], 1e-10);
                Set::Scalar sqrt_rf = std::sqrt(std::max(rho_for_w[a], 1e-10));
                Set::Scalar va = actual_gdot / (burg * rhoM_a);

                Set::Scalar gdot_ratio = actual_gdot / std::max(gdot0_sat, 1e-30);
                Set::Scalar exponent = kB * theta / std::max(A_sat, 1e-30);
                Set::Scalar rhoI_sat = rho0_sat
                    * std::pow(std::max(gdot_ratio, 1e-30), exponent);
                rhoI_sat = std::max(rhoI_sat, (1.0-mob_frac0)*rho0_total/NSYS);
                Set::Scalar rhoM_sat = varsigma * rhoI_sat;
                Set::Scalar d_A = sqrt_rf / std::max(rhoM_sat, 1e-10);

                Set::Scalar rate_mult = CM * sqrt_rf * rhoM_a * va;
                Set::Scalar rate_ann  = 0.25 * CA_val * d_A * rhoM_a * rhoM_a * va;
                Set::Scalar rate_trap = CT_val * sqrt_rf * rhoM_a * va;

                Set::Scalar rate_nuc = 0.0;
                if (enable_nucleation) {
                    Set::Scalar tau_resolved = 0.0;
                    for (int i = 0; i < 3; i++)
                        for (int j = 0; j < 3; j++)
                            tau_resolved += Te[i][j] * S0_3d[a][i][j];
                    Set::Scalar tau_eff_nuc = std::max(std::fabs(tau_resolved) - tau_a_w[a], 0.0);
                    if (tau_eff_nuc > 0.0 && tau_n0 > 0.0) {
                        Set::Scalar r_n = std::min(tau_eff_nuc / tau_n0, 1.0 - 1e-12);
                        Set::Scalar pw_pn = std::pow(r_n, pn_exp);
                        Set::Scalar pw_qn = std::pow(std::max(1.0-pw_pn, 1e-30), qn_exp);
                        Set::Scalar nuc_arg = std::min(DeltaGn0 * pw_qn / (kB*theta), 60.0);
                        rate_nuc = rhodot_n0 * std::exp(-nuc_arg);
                    }
                }

                Set::Scalar rhoI_ratio = std::min(rhoI_a / std::max(rhoI_sat, 1e-10), 1.0);
                Set::Scalar f_rec = 1.0 - std::pow(rhoI_ratio, 1.0 / nrec);

                Set::Scalar drhoM = dt * (rate_mult - rate_ann - rate_trap + rate_nuc);
                if (std::fabs(drhoM) > 0.1 * rhoM_a)
                    drhoM = 0.1 * rhoM_a * (drhoM > 0 ? 1.0 : -1.0);
                Set::Scalar drhoI = dt * f_rec * rate_trap;
                if (std::fabs(drhoI) > 0.1 * rhoI_a)
                    drhoI = 0.1 * rhoI_a * (drhoI > 0 ? 1.0 : -1.0);

                rhoM_w[a] = std::max(rho_M[a] + drhoM, 1e-10);
                rhoI_w[a] = std::max(rho_I[a] + drhoI, 1e-10);

                if (rhoM_w[a] > RHO_MAX) rhoM_w[a] = RHO_MAX;
                if (rhoI_w[a] > RHO_MAX) rhoI_w[a] = RHO_MAX;
            }

            Set::Scalar max_drho = 0.0;
            for (int a = 0; a < NSYS; a++) {
                Set::Scalar rel = std::fabs(rhoM_w[a] - rhoM_old[a])
                                / std::max(rhoM_old[a], 1e-10);
                max_drho = std::max(max_drho, rel);
            }
            if (max_drho < tol_rho) break;
        }

        for (int a = 0; a < NSYS; a++) {
            Set::Scalar tau_a_final = 0.0;
            for (int ii = 0; ii < 3; ii++)
                for (int jj = 0; jj < 3; jj++)
                    tau_a_final += Te[ii][jj] * S0_3d[a][ii][jj];
            wp_inc += tau_a_final * dg[a];
        }

        Set::Scalar Lp[3][3] = {{0,0,0},{0,0,0},{0,0,0}};
        Set::Scalar total_dg = 0.0;
        for (int a = 0; a < NSYS; a++) {
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Lp[i][j] += dg[a] * S0_3d[a][i][j];
            total_dg += std::fabs(dg[a]);
        }
        if (!std::isfinite(total_dg) || total_dg > 10.0) return;

        // Equivalent plastic strain rate, mirroring AdvanceNguyen3D.
        // Implicit path's outer Newton iterates rho_M/rho_I and dg
        // jointly; once converged, the final dg[a] is the step's slip
        // increment and L^p*dt is exactly Lp[i][j] below.
        {
            Set::Scalar Dp[3][3];
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    Dp[i][j] = 0.5 * (Lp[i][j] + Lp[j][i]);

            Set::Scalar DpDp = 0.0;
            for (int i = 0; i < 3; i++)
                for (int j = 0; j < 3; j++)
                    DpDp += Dp[i][j] * Dp[i][j];

            Set::Scalar dep_eq_step =
                std::sqrt(std::max((2.0 / 3.0) * DpDp, 0.0));
            edot_p_eq = (dt > 1e-30) ? (dep_eq_step / dt) : 0.0;
            if (!std::isfinite(edot_p_eq) || edot_p_eq < 0.0)
                edot_p_eq = 0.0;
        }

        Eigen::Matrix3d Fp_old = Fp;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++) {
                Set::Scalar val = 0.0;
                for (int k = 0; k < 3; k++)
                    val += ((i==k ? 1.0 : 0.0) + Lp[i][k]) * Fp_old(k,j);
                Fp(i,j) = val;
            }
        for (int a = 0; a < NSYS; a++) {
            gamma[a] += dg[a];
            rho_M[a]   = rhoM_w[a];
            rho_I[a]   = rhoI_w[a];
            rho_vec[a]  = rhoM_w[a] + rhoI_w[a];
        }
    }

    struct SlipResult { Set::Scalar dg; Set::Scalar ddg_dtau; };

    SlipResult compute_slip_and_deriv(
        Set::Scalar dt, Set::Scalar tau_a_resolved,
        Set::Scalar tau_a_val, Set::Scalar tau_r_val,
        Set::Scalar rhoM_a, Set::Scalar rho_cop_a,
        Set::Scalar ) const
    {
        SlipResult res = {0.0, 0.0};
        Set::Scalar abs_tau = std::fabs(tau_a_resolved);
        Set::Scalar tau_star = std::max(abs_tau - tau_a_val, 0.0);
        if (tau_star <= 0.0 || tau_r_val <= 0.0) return res;

        Set::Scalar sgn = (tau_a_resolved >= 0.0) ? 1.0 : -1.0;
        Set::Scalar Lbar = 1.0 / std::sqrt(std::max(rho_cop_a, 1.0));

        Set::Scalar r = std::min(tau_star / tau_r_val, 1.0 - 1e-12);
        Set::Scalar pw_p = std::pow(r, p_exp);
        Set::Scalar one_m_pp = std::max(1.0 - pw_p, 1e-30);
        Set::Scalar pw_q = std::pow(one_m_pp, q_exp);
        Set::Scalar DG_eff = DeltaG0_nguyen * pw_q;
        Set::Scalar exp_arg = std::min(DG_eff / (kB * theta), 60.0);
        Set::Scalar exp_val = std::exp(exp_arg);
        Set::Scalar t_w = (1.0 / fD) * (exp_val - 1.0);
        t_w = std::max(t_w, 1e-30);

        const Set::Scalar tau_star_Pa = tau_star * 1.0e9;
        const Set::Scalar b_m         = burg * 1.0e-3;
        const Set::Scalar cs_m        = cs   * 1.0e-3;
        DragResult dr = relativistic_drag_velocity(
            tau_star_Pa, b_m, B0, cs_m);
        Set::Scalar v_r = dr.v_r * 1.0e3;
        if (v_r < 1e-30) v_r = 1e-30;
        Set::Scalar t_r = Lbar / v_r;

        Set::Scalar t_sum = t_w + t_r;
        Set::Scalar v_mean = Lbar / t_sum;
        rhoM_a = std::max(rhoM_a, 1e-10);
        Set::Scalar gdot = burg * rhoM_a * v_mean;
        res.dg = dt * gdot * sgn;

        Set::Scalar dtw_dtaustar = 0.0;
        if (r > 1e-30 && one_m_pp > 1e-30) {
            dtw_dtaustar = -(1.0 / fD) * exp_val
                * (DeltaG0_nguyen / (kB * theta))
                * q_exp * std::pow(one_m_pp, q_exp - 1.0)
                * p_exp * std::pow(r, p_exp - 1.0) / tau_r_val;
        }

        Set::Scalar dvr_dtaustar = dr.dvr_dtau * 1.0e12;

        Set::Scalar dtr_dtaustar = -Lbar / (v_r * v_r) * dvr_dtaustar;
        Set::Scalar dv_dtaustar = -Lbar / (t_sum * t_sum)
                                 * (dtw_dtaustar + dtr_dtaustar);
        res.ddg_dtau = dt * burg * rhoM_a * dv_dtaustar;

        return res;
    }

    Set::Matrix3d Fp = Set::Matrix3d::Identity();
    Eigen::Matrix<Set::Scalar,NSYS,1> gamma   = Eigen::Matrix<Set::Scalar,NSYS,1>::Zero();
    Eigen::Matrix<Set::Scalar,NSYS,1> rho_vec = Eigen::Matrix<Set::Scalar,NSYS,1>::Zero();
    Eigen::Matrix<Set::Scalar,NSYS,1> s_vec   = Eigen::Matrix<Set::Scalar,NSYS,1>::Zero();
    Eigen::Matrix<Set::Scalar,NSYS,1> rho_M   = Eigen::Matrix<Set::Scalar,NSYS,1>::Zero();
    Eigen::Matrix<Set::Scalar,NSYS,1> rho_I   = Eigen::Matrix<Set::Scalar,NSYS,1>::Zero();
    Set::Scalar theta = 300.0;

    Set::Scalar wp_inc = 0.0;

    // Equivalent plastic strain rate, sqrt(2/3 D^p:D^p),
    Set::Scalar edot_p_eq = 0.0;

    Set::Scalar gammadot0 = 1.0e7;
    Set::Scalar DeltaG    = 2.1e-19;
    Set::Scalar p_exp     = 0.28;
    Set::Scalar q_exp     = 1.34;
    Set::Scalar burg      = 2.86e-7;
    Set::Scalar tstart    = 0.0;
    Set::Scalar rho0_total= 2.4e6;

    Set::Scalar tau_a0    = 0.027;
    Set::Scalar tau_r0    = 0.525;
    Set::Scalar Cfor      = 0.0;
    Set::Scalar Ccop      = 0.5;

    Set::Scalar tau_p_110   = 0.655;
    Set::Scalar tau_p_112T  = 0.575;
    Set::Scalar tau_p_112AT = 1.075;

    Set::Scalar fD        = 1.0e11;
    Set::Scalar g0        = 0.2;

    Set::Scalar B0        = 2.9574e-5;
    Set::Scalar cs        = 1.937e6;
    Set::Scalar DeltaG0_nguyen = 0.0;

    Set::Scalar CM        = 0.35;
    Set::Scalar CT_ratio  = 0.88;
    Set::Scalar rho0_sat  = 2.66e9;
    Set::Scalar nrec      = 4.0;
    Set::Scalar varsigma  = 0.1;

    Set::Scalar mob_frac0 = 0.01;

    Set::Scalar A_sat     = 1.5762e-19;
    Set::Scalar gdot0_sat = 2.5143e9;

    Set::Scalar use_modscaled_taur = 0.0;

    Set::Scalar enable_nucleation = 0.0;
    Set::Scalar rhodot_n0  = 1.0e15;
    Set::Scalar DeltaGn0   = 5.0e-19;
    Set::Scalar tau_n0     = 3.0;
    Set::Scalar pn_exp     = 0.3;
    Set::Scalar qn_exp     = 1.5;

    Set::Scalar use_implicit = 0.0;

public:
    //
    // ### FIX: Zero() must zero EVERY member named in OP_VARS.
    // ###
    // ### OP_VARS lists 48 members; the old Zero() set only 15.  The other 33
    // ### (tau_p_110, tau_p_112T, tau_p_112AT, burg, B0, cs, DeltaG0_nguyen, CM, ...)
    // ### retained their member-initializer defaults.  The generated operators in
    // ### InClassOperators.H act on the full OP_VARS list, so ANY zero-then-accumulate
    // ### path -- `acc = Zero(); acc += A; acc += B; acc *= 0.5;` -- returns those 33
    // ### multiplied by the number of accumulations.  That is the yield-threshold
    // ### doubling: tau_p_110 = 0.655 -> (0.655 + 0.655 + 0.655)*0.5 = 0.98.
    // ###
    // ### The current LagrangianShockCP does no MODEL arithmetic (ComputeCellStressField
    // ### copies a node then overwrites Fp), so this is latent today.  It fires the
    // ### moment AMR is enabled, because the interpolators take weighted averages of
    // ### MODEL.  Fix it before the static-AMR work, not during.
    //
    static CrystalPlastic Zero()
    {
        CrystalPlastic ret;
        // --- elastic constants + orientation ---
        ret.C11 = 0; ret.C12 = 0; ret.C44 = 0;
        ret.q = Set::Quaternion(0,0,0,0);
        // --- state ---
        ret.Fp = Set::Matrix3d::Zero();
        ret.gamma   = Eigen::Matrix<Set::Scalar,NSYS,1>::Zero();
        ret.rho_vec = Eigen::Matrix<Set::Scalar,NSYS,1>::Zero();
        ret.s_vec   = Eigen::Matrix<Set::Scalar,NSYS,1>::Zero();
        ret.rho_M   = Eigen::Matrix<Set::Scalar,NSYS,1>::Zero();
        ret.rho_I   = Eigen::Matrix<Set::Scalar,NSYS,1>::Zero();
        ret.theta = 0;
        ret.wp_inc = 0;
        ret.edot_p_eq = 0;
        // --- kinetics ---
        ret.gammadot0 = 0; ret.DeltaG = 0; ret.p_exp = 0; ret.q_exp = 0;
        ret.burg = 0; ret.tstart = 0; ret.rho0_total = 0;
        // --- slip resistance ---
        ret.tau_a0 = 0; ret.tau_r0 = 0; ret.Cfor = 0; ret.Ccop = 0;
        ret.tau_p_110 = 0; ret.tau_p_112T = 0; ret.tau_p_112AT = 0;
        // --- thermal activation / drag ---
        ret.fD = 0; ret.g0 = 0; ret.B0 = 0; ret.cs = 0;
        ret.DeltaG0_nguyen = 0;
        // --- density evolution ---
        ret.CM = 0; ret.CT_ratio = 0; ret.rho0_sat = 0; ret.nrec = 0;
        ret.varsigma = 0; ret.mob_frac0 = 0;
        ret.A_sat = 0; ret.gdot0_sat = 0;
        // --- switches ---
        ret.use_modscaled_taur = 0;
        ret.enable_nucleation = 0;
        ret.rhodot_n0 = 0; ret.DeltaGn0 = 0; ret.tau_n0 = 0;
        ret.pn_exp = 0; ret.qn_exp = 0;
        ret.use_implicit = 0;
        return ret;
    }
    static CrystalPlastic Random()
    { return Random(Util::Random(), Util::Random(), Util::Random()); }
    static CrystalPlastic Random(Set::Scalar c11, Set::Scalar c12, Set::Scalar c44)
    {
        CrystalPlastic ret = CrystalPlastic::Zero();
        ret.C11 = c11; ret.C12 = c12; ret.C44 = c44;
        ret.q = Set::Quaternion::UnitRandom();
        return ret;
    }

    static void Parse(CrystalPlastic &value, IO::ParmParse &pp)
    {
        PseudoLinear::Cubic::Parse(value, pp);

        pp_query_default("gammadot0", value.gammadot0, 1.0e7);
        pp_query_default("DeltaG",    value.DeltaG,    2.1e-19);
        pp_query_default("p_exp",     value.p_exp,     0.28);
        pp_query_default("q_exp",     value.q_exp,     1.34);
        pp_query_default("burg",      value.burg,      2.86e-7);
        pp_query_default("theta",     value.theta,     300.0);
        pp_query_default("tstart",    value.tstart,    0.0);
        pp_query_default("rho0_total",value.rho0_total, 6.0e6);

        pp_query_default("tau_a0",    value.tau_a0,    0.027);
        pp_query_default("tau_r0",    value.tau_r0,    0.525);
        pp_query_default("Cfor",      value.Cfor,      0.0);
        pp_query_default("Ccop",      value.Ccop,      0.5);

        pp_query_default("tau_p_110",   value.tau_p_110,   0.655);
        pp_query_default("tau_p_112T",  value.tau_p_112T,  0.575);
        pp_query_default("tau_p_112AT", value.tau_p_112AT, 1.075);

        pp_query_default("fD",        value.fD,        1.0e11);
        pp_query_default("g0",        value.g0,        0.2);

        pp_query_default("B0",        value.B0,        2.9574e-5);
        pp_query_default("cs",        value.cs,        1.937e6);

        pp_query_default("CM",        value.CM,        0.35);
        pp_query_default("CT_ratio",  value.CT_ratio,  0.88);
        pp_query_default("rho0_sat",  value.rho0_sat,  2.66e9);
        pp_query_default("nrec",      value.nrec,      4.0);
        pp_query_default("varsigma",  value.varsigma,  0.1);
        pp_query_default("mob_frac0", value.mob_frac0, 0.01);

        pp_query_default("A_sat",     value.A_sat,     1.5762e-19);
        pp_query_default("gdot0_sat", value.gdot0_sat, 2.5143e9);

        int modscaled_int = 0;
        pp_query_default("use_modscaled_taur", modscaled_int, 0);
        value.use_modscaled_taur = (Set::Scalar)modscaled_int;

        int nuc_int = 0;
        pp_query_default("enable_nucleation", nuc_int, 0);
        value.enable_nucleation = (Set::Scalar)nuc_int;
        pp_query_default("rhodot_n0", value.rhodot_n0, 1.0e15);
        pp_query_default("DeltaGn0",  value.DeltaGn0,  5.0e-19);
        pp_query_default("tau_n0",    value.tau_n0,    3.0);
        pp_query_default("pn_exp",    value.pn_exp,    0.3);
        pp_query_default("qn_exp",    value.qn_exp,    1.5);

        int implicit_int = 0;
        pp_query_default("use_implicit", implicit_int, 0);
        value.use_implicit = (Set::Scalar)implicit_int;

        Set::Scalar mu_now = std::sqrt(value.C44 * (value.C11 - value.C12) / 2.0);

        constexpr Set::Scalar C44_0K = 87.1;
        constexpr Set::Scalar C11_0K = 268.5;
        constexpr Set::Scalar C12_0K = 159.9;
        const Set::Scalar mu_0K = std::sqrt(C44_0K * (C11_0K - C12_0K) / 2.0);

        Set::Scalar b_m = value.burg * 1e-3;
        Set::Scalar mu_Pa = mu_now * 1e9;
        value.DeltaG0_nguyen = value.g0 * mu_Pa * b_m * b_m * b_m;

        Set::Scalar rho_per_sys = value.rho0_total / NSYS;
        for (int a = 0; a < NSYS; a++) {
            value.rho_M[a] = value.mob_frac0 * rho_per_sys;
            value.rho_I[a] = (1.0 - value.mob_frac0) * rho_per_sys;
            value.rho_vec[a] = rho_per_sys;
        }
        value.wp_inc = 0.0;
        value.edot_p_eq = 0.0;

        Set::Vector3d m0[NSYS], n0[NSYS];
        Set::Scalar S0_dummy[NSYS][3][3];
        int fam_dummy[NSYS];
        build_slip_systems(value.q, m0, n0, S0_dummy, fam_dummy);
        Set::Scalar Afor_init[NSYS][NSYS], Acop_init[NSYS][NSYS];
        build_forest_coplanar(m0, n0, Afor_init, Acop_init);

        Set::Scalar rho_tot_arr[NSYS];
        for (int a = 0; a < NSYS; a++) rho_tot_arr[a] = value.rho_vec[a];
        for (int a = 0; a < NSYS; a++) {
            Set::Scalar rho_for = 0.0, rho_cop = 0.0;
            for (int b = 0; b < NSYS; b++) {
                rho_for += Afor_init[a][b] * rho_tot_arr[b];
                rho_cop += Acop_init[a][b] * rho_tot_arr[b];
            }
            Set::Scalar tau_a = value.tau_a0 + value.Cfor * value.burg * mu_now
                * std::sqrt(std::max(rho_for, 0.0));
            Set::Scalar tau_p = (fam_dummy[a] == 0) ? value.tau_p_110 : value.tau_p_112T;
            Set::Scalar tau_r;
            if (modscaled_int) {
                tau_r = tau_p * (mu_now / mu_0K);
            } else {
                tau_r = tau_p + value.Ccop * value.burg * mu_now
                    * std::sqrt(std::max(rho_cop, 0.0));
            }
            value.s_vec[a] = tau_a + tau_r;
        }

        fprintf(stderr, "[CP-Nguyen] Parse: tau_a0=%g tau_p_110=%g tau_p_112T=%g "
                "tau_p_112AT=%g DG0=%g rho_M/sys=%g rho_I/sys=%g B0=%g Pa.s "
                "nucleation=%s modscaled_taur=%s implicit=%s\n",
            value.tau_a0, value.tau_p_110, value.tau_p_112T, value.tau_p_112AT,
            value.DeltaG0_nguyen, value.rho_M[0], value.rho_I[0], value.B0,
            nuc_int ? "ON" : "OFF",
            modscaled_int ? "ON" : "OFF",
            implicit_int ? "ON" : "OFF");
    }

#define OP_CLASS CrystalPlastic
#define OP_VARS X(C11) X(C12) X(C44) X(q) X(Fp) X(gamma) X(rho_vec) X(s_vec) \
    X(rho_M) X(rho_I) \
    X(theta) X(wp_inc) X(edot_p_eq) X(gammadot0) X(DeltaG) X(p_exp) X(q_exp) X(burg) \
    X(rho0_total) X(tstart) \
    X(tau_a0) X(tau_r0) X(Cfor) X(Ccop) X(fD) X(g0) X(B0) X(cs) \
    X(CM) X(CT_ratio) X(rho0_sat) X(nrec) X(varsigma) X(mob_frac0) X(DeltaG0_nguyen) \
    X(A_sat) X(gdot0_sat) \
    X(tau_p_110) X(tau_p_112T) X(tau_p_112AT) \
    X(use_modscaled_taur) \
    X(enable_nucleation) X(rhodot_n0) X(DeltaGn0) X(tau_n0) X(pn_exp) X(qn_exp) \
    X(use_implicit)
#include "Model/Solid/InClassOperators.H"
};
#include "Model/Solid/ExtClassOperators.H"
}
}
}

template<>
inline int Set::Field<Model::Solid::Finite::CrystalPlastic>::NComp() const { return 120; }

template<>
inline std::string Set::Field<Model::Solid::Finite::CrystalPlastic>::Name(int i) const
{
    if (i < 24)  return name + "gamma" + std::to_string(i+1);
    if (i < 48)  return name + "rho"   + std::to_string(i-24+1);
    if (i < 72)  return name + "s"     + std::to_string(i-48+1);
    if (i < 96)  return name + "rhoM"  + std::to_string(i-72+1);
    if (i < 120) return name + "rhoI"  + std::to_string(i-96+1);
    return name;
}

template<>
inline void Set::Field<Model::Solid::Finite::CrystalPlastic>::Copy(
    int a_lev, amrex::MultiFab &a_dst, int a_dstcomp, int a_nghost) const
{
    for (amrex::MFIter mfi(a_dst, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.growntilebox(amrex::IntVect(a_nghost));
        if (bx.ok()) {
            amrex::Array4<const Model::Solid::Finite::CrystalPlastic> const& src
                = ((*this)[a_lev])->array(mfi);
            amrex::Array4<Set::Scalar> const& dst = a_dst.array(mfi);
            for (int n = 0; n < 24; n++)
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    dst(i,j,k,a_dstcomp+n) = src(i,j,k).gamma(n);
                });
            for (int n = 0; n < 24; n++)
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    dst(i,j,k,a_dstcomp+24+n) = src(i,j,k).rho_vec(n);
                });
            for (int n = 0; n < 24; n++)
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    dst(i,j,k,a_dstcomp+48+n) = src(i,j,k).s_vec(n);
                });
            for (int n = 0; n < 24; n++)
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    dst(i,j,k,a_dstcomp+72+n) = src(i,j,k).rho_M(n);
                });
            for (int n = 0; n < 24; n++)
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                    dst(i,j,k,a_dstcomp+96+n) = src(i,j,k).rho_I(n);
                });
        }
    }
}
#endif