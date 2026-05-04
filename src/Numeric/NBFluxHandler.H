// ========================= NBFluxHandler.H =========================
#ifndef NUMERIC_NBFLUXHANDLER_H_
#define NUMERIC_NBFLUXHANDLER_H_

#include <AMReX.H>
#include <AMReX_Array4.H>
#include <AMReX_REAL.H>
#include <AMReX_MultiFab.H>

#include "Numeric/NBStencil.H"
#include "Util/NB_Util.H"
#include "Set/Set.H"

namespace NB {

struct FluxHandler
{
    static constexpr int NL = 3;
    static constexpr int NR = 3;

    // ============================================================
    // HJ-WENO5 derivative (cell-centered, Jiang-Peng style)
    // ============================================================
    template<int DIR, typename PhiArr, typename MaskArr, typename VelArr>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static amrex::Real weno5_derivative(
        int i, int j, int k,
        PhiArr const& phi,
        MaskArr const& mask,
        VelArr const& vel,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx)
    {
        using DirT = NBStencil::Dir<DIR>;

        // Velocity at cell center
        amrex::Real u = vel(i,j,k,DIR);

        // Upwind direction: dir=-1 means left-biased (u>0), dir=+1 means right-biased (u<0)
        int dir = (u > 0.0) ? -1 : 1;

        // Load 7-point stencil centered at (i,j,k)
        // NB_MASK aware: stencil clamps at narrowband boundary
        amrex::Real s[7];
        NBStencil::load_stencil<DirT, NL, NR>(i, j, k, phi, mask, s, 0);

        // Map stencil s[0:6] → Fortran f(-3:3) via f(m) = s[m+3]
        auto f = [&](int m) -> amrex::Real { return s[m + 3]; };

        amrex::Real dxi = 1.0 / dx[DIR];

        // Upwind-biased index offsets (Fortran ord(:,dir) logic)
        // dir=-1 (left-biased):  o1=-2, o2=-1, o3=0, o4=1, o5=2
        // dir=+1 (right-biased): o1= 3, o2= 2, o3=1, o4=0, o5=-1
        int o1 = (dir == -1) ? -2 :  3;
        int o2 = (dir == -1) ? -1 :  2;
        int o3 = (dir == -1) ?  0 :  1;
        int o4 = (dir == -1) ?  1 :  0;
        int o5 = (dir == -1) ?  2 : -1;

        // One-sided finite difference slopes
        amrex::Real v1 = (f(o1) - f(o1-1)) * dxi;
        amrex::Real v2 = (f(o2) - f(o2-1)) * dxi;
        amrex::Real v3 = (f(o3) - f(o3-1)) * dxi;
        amrex::Real v4 = (f(o4) - f(o4-1)) * dxi;
        amrex::Real v5 = (f(o5) - f(o5-1)) * dxi;

        // Candidate 3rd-order derivatives from each sub-stencil
        amrex::Real phi1 =  (1.0/3.0)*v1 - (7.0/6.0)*v2 + (11.0/6.0)*v3;
        amrex::Real phi2 = -(1.0/6.0)*v2 + (5.0/6.0)*v3 +  (1.0/3.0)*v4;
        amrex::Real phi3 =  (1.0/3.0)*v3 + (5.0/6.0)*v4 -  (1.0/6.0)*v5;

        // Jiang-Shu smoothness indicators
        amrex::Real S1a = v1 - 2.0*v2 + v3;
        amrex::Real S1b = v1 - 4.0*v2 + 3.0*v3;
        amrex::Real S1  = (13.0/12.0)*S1a*S1a + 0.25*S1b*S1b;

        amrex::Real S2a = v2 - 2.0*v3 + v4;
        amrex::Real S2b = v2 - v4;
        amrex::Real S2  = (13.0/12.0)*S2a*S2a + 0.25*S2b*S2b;

        amrex::Real S3a = v3 - 2.0*v4 + v5;
        amrex::Real S3b = 3.0*v3 - 4.0*v4 + v5;
        amrex::Real S3  = (13.0/12.0)*S3a*S3a + 0.25*S3b*S3b;

        // WENO weights (ideal: 0.1, 0.6, 0.3)
        constexpr amrex::Real eps = 1.0e-6;
        amrex::Real a1 = 0.1 / ((S1 + eps)*(S1 + eps));
        amrex::Real a2 = 0.6 / ((S2 + eps)*(S2 + eps));
        amrex::Real a3 = 0.3 / ((S3 + eps)*(S3 + eps));
        amrex::Real asum = 1.0 / (a1 + a2 + a3);

        amrex::Real w1 = a1 * asum;
        amrex::Real w2 = a2 * asum;
        amrex::Real w3 = 1.0 - w1 - w2;

        return w1*phi1 + w2*phi2 + w3*phi3;
    }

    // ============================================================
    // Per-cell WENO5 RHS:  L(phi) = -u . grad(phi)
    //
    // Returns 0 for cells outside the narrowband (NB_MASK == 0).
    // Made static so it can be called inside ParallelFor lambdas
    // on GPU without capturing `this`.
    // ============================================================
    template<typename PhiArr, typename MaskArr, typename VelArr>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static Set::Scalar computeWENO5RHS(
        int i, int j, int k,
        PhiArr const& phi,
        MaskArr const& mask,
        VelArr const& vel,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx)
    {
        const int bandval =
            static_cast<int>(mask(i,j,k,NarrowBandData::TubeIdx::NB_MASK));

        if (bandval == 0) return Set::Scalar(0.0);

        Set::Scalar rhs = 0.0;

        rhs -= vel(i,j,k,0) *
            weno5_derivative<NBStencil::XDir::value>(i,j,k,phi,mask,vel,dx);

    #if (AMREX_SPACEDIM >= 2)
        rhs -= vel(i,j,k,1) *
            weno5_derivative<NBStencil::YDir::value>(i,j,k,phi,mask,vel,dx);
    #endif

    #if (AMREX_SPACEDIM == 3)
        rhs -= vel(i,j,k,2) *
            weno5_derivative<NBStencil::ZDir::value>(i,j,k,phi,mask,vel,dx);
    #endif

        return rhs;
    }

    // ============================================================
    // Per-cell SSP-RK3 stage update:
    //
    //   phi_out(i,j,k) = a * phi_old(i,j,k)
    //                  + b * (phi_stage(i,j,k) + dt * L(phi_stage))
    //
    // Stage coefficients (a, b):
    //   Stage 1:  a=0,    b=1       ->  phi* = phi^n + dt*L(phi^n)
    //   Stage 2:  a=3/4,  b=1/4     ->  phi** = 3/4*phi^n + 1/4*(phi* + dt*L(phi*))
    //   Stage 3:  a=1/3,  b=2/3     ->  phi^n+1 = 1/3*phi^n + 2/3*(phi** + dt*L(phi**))
    //
    // Non-narrowband cells are always assigned phi_old unchanged.
    // Made static so it can be called inside ParallelFor lambdas on GPU.
    // ============================================================
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static void rk_stage_adv(
        int i, int j, int k,
        amrex::Array4<const amrex::Real> const& phi_old,    // phi^n (read-only throughout)
        amrex::Array4<const amrex::Real> const& phi_stage,  // stage input (read for RHS + update)
        amrex::Array4<      amrex::Real> const& phi_out,    // stage output (must differ from phi_stage)
        amrex::Array4<const Set::Scalar> const& vel,
        amrex::Array4<const amrex::Real> const& mask,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx,
        amrex::Real dt,
        amrex::Real a,
        amrex::Real b)
    {
        const int m =
            static_cast<int>(mask(i,j,k,NarrowBandData::TubeIdx::NB_MASK));

        // Outside narrowband: preserve phi^n exactly, no RHS needed
        if (!m) {
            phi_out(i,j,k) = phi_old(i,j,k);
            return;
        }

        amrex::Real rhs = computeWENO5RHS(i,j,k,phi_stage,mask,vel,dx);

        phi_out(i,j,k) = a * phi_old(i,j,k)
                       + b * (phi_stage(i,j,k) + dt * rhs);
    }

    // ============================================================
    // Field-level SSP-RK3 advection (single scratch fab)
    //
    // Memory layout (safe: read and write fabs are always distinct):
    //
    //   Stage 1:  phi_old(LSold) + phi_stage(LSold)  -> phi_out(scratch)   [phi*]
    //   Stage 2:  phi_old(LSold) + phi_stage(scratch) -> phi_out(LS)        [phi**]
    //   Stage 3:  phi_old(LSold) + phi_stage(LS)      -> phi_out(scratch)   [phi^n+1]
    //   Swap:     std::swap(LS[lev], scratch)
    //
    // Ghost fill strategy between stages:
    //   - FillBoundary()    : halo exchange between MPI ranks / tiles (interior ghosts)
    //   - bc.FillBoundary() : physical BCs applied at domain boundaries
    //   Both are required. bc.define(geom) must have been called before Advance().
    //
    // The caller is responsible for:
    //   - Registering scratch as a Set::Field with the same BC/nghost as LS.
    //   - Copying LS -> LSold before calling this (i.e. LSold holds phi^n).
    //   - Calling update_Tube() after this returns if the tube needs refreshing.
    // ============================================================
    static void advanceRK3(
        int lev,
        Set::Field<Set::Scalar>      & LS,       // current level set (phi^n+1 on exit)
        Set::Field<Set::Scalar> const& LSold,     // phi^n (read-only, never written)
        Set::Field<Set::Scalar> const& TUBE,      // NB_MASK and tube type flags
        Set::Field<Set::Scalar> const& vel,       // advection velocity
        amrex::Vector<char>          & tiles,
        BC::BC<Set::Scalar>          & bc,        // physical BC for domain boundary ghosts
        const amrex::Geometry        & geom,
        amrex::Real                    time,      // current sim time (needed by bc.FillBoundary)
        amrex::Real                    dt)
    {
        BL_PROFILE("NB::FluxHandler::advanceRK3");

        auto const dx = geom.CellSizeArray();

        // SSP-RK3 Shu-Osher coefficients
        // stage: (a,   b  )
        //   1:   (0,   1  )
        //   2:   (3/4, 1/4)
        //   3:   (1/3, 2/3)
        constexpr amrex::Real a1 = 0.0,        b1 = 1.0;
        constexpr amrex::Real a2 = 3.0/4.0,    b2 = 1.0/4.0;
        constexpr amrex::Real a3 = 1.0/3.0,    b3 = 2.0/3.0;

        // Define scratch multifab for RHS initialized to LS
        amrex::MultiFab scratch;
        scratch.define(
            LS[lev]->boxArray(),
            LS[lev]->DistributionMap(),
            1,
            0
        );

        amrex::MultiFab::Copy(
            scratch,
            *LS[lev],
            0,
            0,
            1,
            0
        );

        // --------------------------------------------------------
        // Stage 1:  phi* = phi^n + dt*L(phi^n)
        //   phi_old   = LSold   (phi^n)
        //   phi_stage = LSold   (phi^n, same pointer — safe because
        //                        phi_out goes to scratch, not LSold)
        //   phi_out   = scratch
        // --------------------------------------------------------
        int tile_idx = 0;
        for (amrex::MFIter mfi(*LS[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi, ++tile_idx)
        {
            if (!tiles[tile_idx]) continue;

            const auto& tilebox = mfi.tilebox();

            Set::Patch<const Set::Scalar> phi_old_arr = LSold.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> tube_arr    = TUBE.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> vel_arr     = vel.Patch(lev, mfi);
            auto       scratch_arr   = scratch.array(mfi);

            amrex::ParallelFor(tilebox, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                rk_stage_adv(i, j, k,
                    phi_old_arr,    // phi^n
                    phi_old_arr,    // stage input = phi^n (same fab, different ptr role)
                    scratch_arr,    // phi* written here
                    vel_arr, tube_arr, dx, dt,
                    a1, b1);
            });
        }

        // Fill interior ghosts (MPI halo exchange), then domain boundary ghosts
        scratch.FillBoundary();
        bc.FillBoundary(scratch, 0, 0, time, 0);

        // --------------------------------------------------------
        // Stage 2:  phi** = 3/4*phi^n + 1/4*(phi* + dt*L(phi*))
        //   phi_old   = LSold   (phi^n)
        //   phi_stage = scratch (phi*)
        //   phi_out   = LS[lev] (phi** written here — distinct from scratch)
        // --------------------------------------------------------
        tile_idx = 0;
        for (amrex::MFIter mfi(*LS[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi, ++tile_idx)
        {
            if (!tiles[tile_idx]) continue;

            const auto& tilebox = mfi.tilebox();

            Set::Patch<const Set::Scalar> phi_old_arr  = LSold.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> phi_star_arr = scratch.const_array(mfi);
            Set::Patch<const Set::Scalar> tube_arr     = TUBE.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> vel_arr      = vel.Patch(lev, mfi);
            Set::Patch<      Set::Scalar> ls_arr       = LS.Patch(lev, mfi);

            amrex::ParallelFor(tilebox, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                rk_stage_adv(i, j, k,
                    phi_old_arr,    // phi^n
                    phi_star_arr,   // stage input = phi*
                    ls_arr,         // phi** written into LS (distinct from scratch)
                    vel_arr, tube_arr, dx, dt,
                    a2, b2);
            });
        }

        // Fill interior ghosts (MPI halo exchange), then domain boundary ghosts
        LS[lev]->FillBoundary();
        bc.FillBoundary(*LS[lev], 0, 0, time, 0);

        // --------------------------------------------------------
        // Stage 3:  phi^n+1 = 1/3*phi^n + 2/3*(phi** + dt*L(phi**))
        //   phi_old   = LSold   (phi^n)
        //   phi_stage = LS[lev] (phi**)
        //   phi_out   = scratch (phi^n+1 written here — distinct from LS)
        // --------------------------------------------------------
        tile_idx = 0;
        for (amrex::MFIter mfi(*LS[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi, ++tile_idx)
        {
            if (!tiles[tile_idx]) continue;

            const auto& tilebox = mfi.tilebox();

            Set::Patch<const Set::Scalar> phi_old_arr   = LSold.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> phi_dstar_arr = LS.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> tube_arr      = TUBE.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> vel_arr       = vel.Patch(lev, mfi);
            Set::Patch<      Set::Scalar> scratch_arr   = scratch.array(mfi);

            amrex::ParallelFor(tilebox, [=] AMREX_GPU_DEVICE (int i, int j, int k)
            {
                rk_stage_adv(i, j, k,
                    phi_old_arr,    // phi^n
                    phi_dstar_arr,  // stage input = phi**
                    scratch_arr,    // phi^n+1 written into scratch
                    vel_arr, tube_arr, dx, dt,
                    a3, b3);
            });
        }

        // phi^n+1 now lives in scratch. Swap so LS[lev] holds the final result.
        // std::swap exchanges the unique_ptr internals — O(1), no data copy.
        std::swap(*LS[lev], scratch);

        LS[lev]->FillBoundary();
        bc.FillBoundary(*LS[lev], 0, 0, time, 0);
    }

    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    static Set::Scalar minmod(Set::Scalar a, Set::Scalar b)
    {
        const Set::Scalar prod  = a * b;
        const Set::Scalar abs_a = amrex::Math::abs(a);
        const Set::Scalar abs_b = amrex::Math::abs(b);

        return (prod <= 0.0) ? Set::Scalar(0.0)
                            : ((abs_a < abs_b) ? a : b);
    }

    template<typename DIR>
    AMREX_GPU_DEVICE AMREX_FORCE_INLINE
    static amrex::Real eno2_flux(
        int i, int j, int k,
        const Set::Scalar sgn,
        amrex::Array4<const amrex::Real> const& reinitold,
        amrex::Array4<const amrex::Real> const& mask,
        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx)
    {
        using Set::Scalar;

        constexpr int dim = DIR::value;

        const int comp = 0;

        // --- stencil ---
        Scalar s[5];
        NBStencil::load_stencil<DIR,2,2>(i,j,k,reinitold,mask,s,comp);

        // --- spacing ---
        const Scalar d = dx[dim];

        // --- first derivatives ---
        Scalar dm = (s[2] - s[1]) / d;
        Scalar dp = (s[3] - s[2]) / d;

        // --- second derivatives ---
        Scalar d2c = (s[1] - 2.0*s[2] + s[3]) / (2.0*d);
        Scalar d2l = (s[0] - 2.0*s[1] + s[2]) / (2.0*d);
        Scalar d2r = (s[2] - 2.0*s[3] + s[4]) / (2.0*d);

        // --- ENO correction ---
        dm += 0.5 * d * minmod(d2c, d2l);
        dp -= 0.5 * d * minmod(d2c, d2r);

        // --- upwinding ---
        const Scalar dmS = dm * sgn;
        const Scalar dpS = dp * sgn;

        return (dpS < 0.0 && dmS < -dpS) ? dp :
            (dmS >= 0.0 && dpS >= -dmS) ? dm : 0.0;
    }
};

} // namespace NB

#endif