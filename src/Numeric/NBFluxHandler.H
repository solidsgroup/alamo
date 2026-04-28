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
        Set::Field<Set::Scalar>      & scratch,   // one scratch field, same layout as LS
        BC::BC<Set::Scalar>          & bc,        // physical BC for domain boundary ghosts
        const amrex::Box             & bbox,      // narrowband bounding box for this object
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

        // --------------------------------------------------------
        // Stage 1:  phi* = phi^n + dt*L(phi^n)
        //   phi_old   = LSold   (phi^n)
        //   phi_stage = LSold   (phi^n, same pointer — safe because
        //                        phi_out goes to scratch, not LSold)
        //   phi_out   = scratch
        // --------------------------------------------------------
        for (amrex::MFIter mfi(*LS[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box workbox = mfi.tilebox() & bbox;
            if (!workbox.ok()) continue;

            auto const phi_old_arr   = LSold[lev]->const_array(mfi);
            auto const tube_arr      = TUBE [lev]->const_array(mfi);
            auto const vel_arr       = vel  [lev]->const_array(mfi);
            auto       scratch_arr   = scratch[lev]->array(mfi);

            amrex::ParallelFor(workbox, [=] AMREX_GPU_DEVICE (int i, int j, int k)
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
        scratch[lev]->FillBoundary(geom.periodicity());
        bc.FillBoundary(*scratch[lev], 0, 0, time, 0);

        // --------------------------------------------------------
        // Stage 2:  phi** = 3/4*phi^n + 1/4*(phi* + dt*L(phi*))
        //   phi_old   = LSold   (phi^n)
        //   phi_stage = scratch (phi*)
        //   phi_out   = LS[lev] (phi** written here — distinct from scratch)
        // --------------------------------------------------------
        for (amrex::MFIter mfi(*LS[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box workbox = mfi.tilebox() & bbox;
            if (!workbox.ok()) continue;

            auto const phi_old_arr   = LSold  [lev]->const_array(mfi);
            auto const phi_star_arr  = scratch[lev]->const_array(mfi);
            auto const tube_arr      = TUBE   [lev]->const_array(mfi);
            auto const vel_arr       = vel    [lev]->const_array(mfi);
            auto       ls_arr        = LS     [lev]->array(mfi);

            amrex::ParallelFor(workbox, [=] AMREX_GPU_DEVICE (int i, int j, int k)
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
        LS[lev]->FillBoundary(geom.periodicity());
        bc.FillBoundary(*LS[lev], 0, 0, time, 0);

        // --------------------------------------------------------
        // Stage 3:  phi^n+1 = 1/3*phi^n + 2/3*(phi** + dt*L(phi**))
        //   phi_old   = LSold   (phi^n)
        //   phi_stage = LS[lev] (phi**)
        //   phi_out   = scratch (phi^n+1 written here — distinct from LS)
        // --------------------------------------------------------
        for (amrex::MFIter mfi(*LS[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box workbox = mfi.tilebox() & bbox;
            if (!workbox.ok()) continue;

            auto const phi_old_arr    = LSold  [lev]->const_array(mfi);
            auto const phi_dstar_arr  = LS     [lev]->const_array(mfi);
            auto const tube_arr       = TUBE   [lev]->const_array(mfi);
            auto const vel_arr        = vel    [lev]->const_array(mfi);
            auto       scratch_arr    = scratch[lev]->array(mfi);

            amrex::ParallelFor(workbox, [=] AMREX_GPU_DEVICE (int i, int j, int k)
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
        std::swap(LS[lev], scratch[lev]);

        // Caller should call FillBoundary + bc.FillBoundary on LS[lev],
        // then update_Tube() and compute_geometry() as needed.
    }
};

} // namespace NB

#endif