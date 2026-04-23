#ifndef Numeric_NBStencil_H_
#define Numeric_NBStencil_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_REAL.H>
#include <AMReX_Array4.H>
#include <AMReX_Extension.H>

namespace NBStencil
{
    //==================================================
    // Direction tags
    //==================================================
    template<int DIR>
    struct Dir { static constexpr int value = DIR; };

    using XDir = Dir<0>;
    using YDir = Dir<1>;
    using ZDir = Dir<2>;

    //==================================================
    // Direction mapping (compile-time)
    //==================================================
    template<int DIR> struct DirMap;

    template<> struct DirMap<0> { static constexpr int dx = 1, dy = 0, dz = 0; };
    template<> struct DirMap<1> { static constexpr int dx = 0, dy = 1, dz = 0; };
    template<> struct DirMap<2> { static constexpr int dx = 0, dy = 0, dz = 1; };

    //==================================================
    // Mask helper (handles Real mask cleanly)
    //==================================================
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    int is_nb(amrex::Real m)
    {
        return m != amrex::Real(0.0);
    }

    //==================================================
    // Generic stencil loader (directional)
    //==================================================
    template<typename DIR, int NL, int NR>
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    void load_stencil (int i, int j, int k,
                       amrex::Array4<const amrex::Real> const& phi,
                       amrex::Array4<const amrex::Real> const& mask,
                       amrex::Real (&s)[NL+NR+1])
    {
        constexpr int dx = DirMap<DIR::value>::dx;
        constexpr int dy = DirMap<DIR::value>::dy;
        constexpr int dz = DirMap<DIR::value>::dz;

        // Center index
        constexpr int C = NL;

        // Center value
        amrex::Real center = phi(i,j,k);
        s[C] = center;

        // ---- negative direction ----
        {
            int ii = i;
            int jj = j;
            int kk = k;

            amrex::Real last_val = center;

            #pragma unroll
            for (int n = 1; n <= NL; ++n)
            {
                ii -= dx;
                jj -= dy;
                kk -= dz;

                int m = is_nb(mask(ii,jj,kk));

                // branchless update (GPU-friendly)
                amrex::Real val = phi(ii,jj,kk);
                last_val = m * val + (1-m) * last_val;

                s[C - n] = last_val;
            }
        }

        // ---- positive direction ----
        {
            int ii = i;
            int jj = j;
            int kk = k;

            amrex::Real last_val = center;

            #pragma unroll
            for (int n = 1; n <= NR; ++n)
            {
                ii += dx;
                jj += dy;
                kk += dz;

                int m = is_nb(mask(ii,jj,kk));

                amrex::Real val = phi(ii,jj,kk);
                last_val = m * val + (1-m) * last_val;

                s[C + n] = last_val;
            }
        }
    }


    //==================================================
    // Fused gradient (x/y/z in one pass)
    //==================================================
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    Set::Vector norm_gradient (int i, int j, int k,
                        amrex::Array4<const amrex::Real> const& phi,
                        amrex::Array4<const amrex::Real> const& mask,
                        amrex::GpuArray<amrex::Real, AMREX_SPACEDIM> const& dx)
    {
        Set::Vector g;
        g.setZero();

        // -------------------------------------------------
        // Load all stencils ONCE
        // -------------------------------------------------
        amrex::Real sx[5];
        load_stencil<XDir,2,2>(i,j,k,phi,mask,sx);

    #if AMREX_SPACEDIM >= 2
        amrex::Real sy[5];
        load_stencil<YDir,2,2>(i,j,k,phi,mask,sy);
    #endif

    #if AMREX_SPACEDIM == 3
        amrex::Real sz[5];
        load_stencil<ZDir,2,2>(i,j,k,phi,mask,sz);
    #endif

        // -------------------------------------------------
        // Load ALL mask values once (reuse across logic)
        // -------------------------------------------------
        // X-direction
        int xm2 = is_nb(mask(i-2,j,k));
        int xm1 = is_nb(mask(i-1,j,k));
        int xp1 = is_nb(mask(i+1,j,k));
        int xp2 = is_nb(mask(i+2,j,k));

    #if AMREX_SPACEDIM >= 2
        int ym2 = is_nb(mask(i,j-2,k));
        int ym1 = is_nb(mask(i,j-1,k));
        int yp1 = is_nb(mask(i,j+1,k));
        int yp2 = is_nb(mask(i,j+2,k));
    #endif

    #if AMREX_SPACEDIM == 3
        int zm2 = is_nb(mask(i,j,k-2));
        int zm1 = is_nb(mask(i,j,k-1));
        int zp1 = is_nb(mask(i,j,k+1));
        int zp2 = is_nb(mask(i,j,k+2));
    #endif

        // -------------------------------------------------
        // Helper lambda (inlined) for derivative logic
        // -------------------------------------------------
        auto compute_dir = [] AMREX_GPU_DEVICE (
            const amrex::Real* s,
            int m2m, int m1m, int m1p, int m2p,
            amrex::Real invdx)
        {
            int use4 = m2m & m1m & m1p & m2p;
            int use2 = m1m & m1p;
            int useF = m1p & m2p;
            int useB = m1m & m2m;

            // priority: 4th > 2nd > one-sided
            use2 = use2 & (1 - use4);
            useF = useF & (1 - use4) & (1 - use2);
            useB = useB & (1 - use4) & (1 - use2) & (1 - useF);

            // 4th-order
            amrex::Real d4 =
                (-s[4] + 8.0*s[3] - 8.0*s[1] + s[0]) * (invdx / 12.0);

            // 2nd central
            amrex::Real d2 =
                (s[3] - s[1]) * (invdx / 2.0);

            // forward
            amrex::Real dF =
                (s[3] - s[2]) * invdx;

            // backward
            amrex::Real dB =
                (s[2] - s[1]) * invdx;

            amrex::Real d =
                use4 * d4
                + use2 * d2
                + useF * dF
                + useB * dB;

            int any = use4 | use2 | useF | useB;
            return any * d;
        };

        // -------------------------------------------------
        // Compute derivatives
        // -------------------------------------------------
        g[0] = compute_dir(sx, xm2, xm1, xp1, xp2, 1.0/dx[0]);

    #if AMREX_SPACEDIM >= 2
        g[1] = compute_dir(sy, ym2, ym1, yp1, yp2, 1.0/dx[1]);
    #endif

    #if AMREX_SPACEDIM == 3
        g[2] = compute_dir(sz, zm2, zm1, zp1, zp2, 1.0/dx[2]);
    #endif

        return g;
    }
}

#endif