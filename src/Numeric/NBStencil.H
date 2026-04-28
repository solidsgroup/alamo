#ifndef NUMERIC_NBSTENCIL_H_
#define NUMERIC_NBSTENCIL_H_

#include <AMReX_GpuQualifiers.H>
#include <AMReX_Array4.H>
#include <AMReX_REAL.H>
#include <AMReX_Math.H>

#include "Set/Base.H"

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
    // Direction mapping
    //==================================================
    template<int DIR> struct DirMap;

    template<> struct DirMap<0> { static constexpr int dx=1, dy=0, dz=0; };
    template<> struct DirMap<1> { static constexpr int dx=0, dy=1, dz=0; };
    template<> struct DirMap<2> { static constexpr int dx=0, dy=0, dz=1; };

    //==================================================
    // Mask helper
    //==================================================
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    int is_nb(amrex::Real m)
    {
        return m != amrex::Real(0.0);
    }

    struct Mask4
    {
        int m2m, m1m, m1p, m2p;
    };

    template<typename DIR>
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    Mask4 get_mask4(int i, int j, int k,
                    amrex::Array4<const amrex::Real> const& mask)
    {
        constexpr int dx = DirMap<DIR::value>::dx;
        constexpr int dy = DirMap<DIR::value>::dy;
        constexpr int dz = DirMap<DIR::value>::dz;

        Mask4 m;

        m.m2m = is_nb(mask(i - 2*dx, j - 2*dy, k - 2*dz));
        m.m1m = is_nb(mask(i - 1*dx, j - 1*dy, k - 1*dz));
        m.m1p = is_nb(mask(i + 1*dx, j + 1*dy, k + 1*dz));
        m.m2p = is_nb(mask(i + 2*dx, j + 2*dy, k + 2*dz));

        return m;
    }

    //==================================================
    // Generic stencil loader
    //==================================================
    template<typename DIR, int NL, int NR>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void load_stencil(int i,int j,int k,
        amrex::Array4<const amrex::Real> const& phi,
        amrex::Array4<const amrex::Real> const& mask,
        amrex::Real (&s)[NL+NR+1])
    {
        constexpr int dx = DirMap<DIR::value>::dx;
        constexpr int dy = DirMap<DIR::value>::dy;
        constexpr int dz = DirMap<DIR::value>::dz;

        constexpr int C = NL;

        amrex::Real center = phi(i,j,k);
        s[C] = center;

        // negative
        {
            int ii=i,jj=j,kk=k;
            amrex::Real last=center;

            #pragma unroll
            for(int n=1;n<=NL;++n){
                ii-=dx; jj-=dy; kk-=dz;
                int m = is_nb(mask(ii,jj,kk));
                amrex::Real val = phi(ii,jj,kk);
                last = m*val + (1-m)*last;
                s[C-n] = last;
            }
        }

        // positive
        {
            int ii=i,jj=j,kk=k;
            amrex::Real last=center;

            #pragma unroll
            for(int n=1;n<=NR;++n){
                ii+=dx; jj+=dy; kk+=dz;
                int m = is_nb(mask(ii,jj,kk));
                amrex::Real val = phi(ii,jj,kk);
                last = m*val + (1-m)*last;
                s[C+n] = last;
            }
        }
    }

    template <class DIR, int NL, int NR, int COMP>
    AMREX_GPU_HOST_DEVICE 
    AMREX_FORCE_INLINE
    void load_stencil_vec(
        int i, int j, int k,
        Set::Patch<const Set::Vector> const& v,
        amrex::Real (&s)[(NL + NR) + 1])
    {
        constexpr int N = NL + NR + 1;

        for (int n = 0; n < N; ++n) {
            int ii = i;
            int jj = j;
            int kk = k;

            if constexpr (std::is_same_v<DIR, XDir>) ii += n - NL;
    #if AMREX_SPACEDIM >= 2
            if constexpr (std::is_same_v<DIR, YDir>) jj += n - NL;
    #endif
    #if AMREX_SPACEDIM == 3
            if constexpr (std::is_same_v<DIR, ZDir>) kk += n - NL;
    #endif

            s[n] = v(ii,jj,kk)(COMP);
        }
    }

    //==================================================
    // Core derivative (reusable everywhere)
    //==================================================
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real compute_derivative(
        const amrex::Real* s,
        int m2m,int m1m,int m1p,int m2p,
        amrex::Real invdx)
    {
        int use4 = m2m & m1m & m1p & m2p;
        int use2 = m1m & m1p;
        int useF = m1p & m2p;
        int useB = m1m & m2m;

        use2 &= (1-use4);
        useF &= (1-use4)&(1-use2);
        useB &= (1-use4)&(1-use2)&(1-useF);

        amrex::Real d4 = (-s[4]+8*s[3]-8*s[1]+s[0])*(invdx/12.0);
        amrex::Real d2 = (s[3]-s[1])*(invdx/2.0);
        amrex::Real dF = (s[3]-s[2])*invdx;
        amrex::Real dB = (s[2]-s[1])*invdx;

        amrex::Real d =
            use4*d4 + use2*d2 + useF*dF + useB*dB;

        int any = use4|use2|useF|useB;
        return any*d;
    }

    //==================================================
    // Gradient
    //==================================================
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Vector gradient(int i,int j,int k,
        amrex::Array4<const amrex::Real> const& phi,
        amrex::Array4<const amrex::Real> const& mask,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx)
    {
        Set::Vector g; g.setZero();

        amrex::Real s[5];
        load_stencil<XDir,2,2>(i,j,k,phi,mask,s);

        auto m = get_mask4<XDir>(i,j,k,mask);

        g[0] = compute_derivative(s,m.m2m,m.m1m,m.m1p,m.m2p,1.0/dx[0]);

#if AMREX_SPACEDIM>=2
        load_stencil<YDir,2,2>(i,j,k,phi,mask,s);

        m = get_mask4<YDir>(i,j,k,mask);

        g[1] = compute_derivative(s,m.m2m,m.m1m,m.m1p,m.m2p,1.0/dx[1]);
#endif

#if AMREX_SPACEDIM==3
        load_stencil<ZDir,2,2>(i,j,k,phi,mask,s);

        m = get_mask4<ZDir>(i,j,k,mask);

        g[2] = compute_derivative(s,m.m2m,m.m1m,m.m1p,m.m2p,1.0/dx[2]);
#endif

        return g;
    }

    //==================================================
    // Normal = normalized gradient
    //==================================================
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Vector normal(int i,int j,int k,
        amrex::Array4<const amrex::Real> const& phi,
        amrex::Array4<const amrex::Real> const& mask,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx)
    {
        Set::Vector g = gradient(i,j,k,phi,mask,dx);
        amrex::Real mag = amrex::max(g.norm(), amrex::Real(1e-12));
        return g / mag;
    }

    //==================================================
    // Curvature = div(n)
    //==================================================
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::Real curvature(int i,int j,int k,
        Set::Patch<const Set::Vector> const& n,
        amrex::Array4<const amrex::Real> const& mask,
        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dx)
    {
        amrex::Real dnx_dx = 0.0;
        amrex::Real dny_dy = 0.0;
        amrex::Real dnz_dz = 0.0;

        amrex::Real s[5];

        // X
        load_stencil_vec<XDir,2,2,0>(i,j,k,n,s);


        auto m = get_mask4<XDir>(i,j,k,mask);
        dnx_dx = compute_derivative(s,m.m2m,m.m1m,m.m1p,m.m2p,1.0/dx[0]);

#if AMREX_SPACEDIM>=2
        load_stencil_vec<YDir,2,2,1>(i,j,k,n,s);

        m = get_mask4<YDir>(i,j,k,mask);
        dny_dy = compute_derivative(s,m.m2m,m.m1m,m.m1p,m.m2p,1.0/dx[1]);
#endif

#if AMREX_SPACEDIM==3
        load_stencil_vec<ZDir,2,2,2>(i,j,k,n,s);

        m = get_mask4<ZDir>(i,j,k,mask);
        dnz_dz = compute_derivative(s,m.m2m,m.m1m,m.m1p,m.m2p,1.0/dx[2]);
#endif

        return dnx_dx + dny_dy + dnz_dz;
    }

} // namespace NBStencil

#endif