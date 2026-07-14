//
// This header defines :code:`Numeric::Advect::MUSCL`, a second-order finite-volume
// advection operator using piecewise-linear MUSCL reconstruction with a selectable
// slope limiter.
//
// The operator is written as a GPU-safe functor:
//
// .. code-block:: c++
//
//    rhs = muscl(phi, velocity, i, j, k, comp, DX, options, stencil);
//
// It reconstructs left and right states of :math:`\phi` at cell faces, selects
// the upwind state using face-centered velocities obtained by averaging neighboring
// cell-centered velocity values, and returns either
// :math:`-\nabla\cdot(\phi\mathbf{u})` or
// :math:`-\nabla\cdot(\phi\mathbf{u}) + \phi\nabla\cdot\mathbf{u}` depending on
// :code:`Numeric::Advect::Options::form`.
//
// This implementation assumes collocated cell-centered data:
// :code:`phi_location = Set::HC::Cell` and
// :code:`velocity_location = Set::HC::Cell`.  It requires two ghost cells in each
// coordinate direction because the limited reconstruction reaches to
// :math:`i\pm2`, :math:`j\pm2`, and :math:`k\pm2`.
//

#ifndef NUMERIC_ADVECT_MUSCL_H
#define NUMERIC_ADVECT_MUSCL_H

#include "Advect.H"
#include "Limiter/Limiter.H"
#include "Limiter/Koren.H"
#include "Limiter/MC.H"
#include "Limiter/MinMod.H"
#include "Limiter/Superbee.H"
#include "Limiter/UMIST.H"
#include "Limiter/VanAlbada.H"
#include "Limiter/VanLeer.H"

namespace Numeric
{
namespace Advect
{
class MUSCL
{
public:
    static constexpr const char* name = "muscl";
    static constexpr Set::HC phi_location = Set::HC::Cell;
    static constexpr Set::HC velocity_location = Set::HC::Cell;

    Limiter::Limiter<Limiter::MC,
                     Limiter::MinMod,
                     Limiter::Superbee,
                     Limiter::VanLeer,
                     Limiter::VanAlbada,
                     Limiter::Koren,
                     Limiter::UMIST> limiter;

    static void Parse(MUSCL& value, IO::ParmParse& pp)
    {
        pp.select<Limiter::MC,
                  Limiter::MinMod,
                  Limiter::Superbee,
                  Limiter::VanLeer,
                  Limiter::VanAlbada,
                  Limiter::Koren,
                  Limiter::UMIST>("limiter",value.limiter);
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    int NGhost() const
    {
        return 2;
    }

    template<class Phi>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar Slope(Phi const& phi,
                      int ilo, int jlo, int klo,
                      int ic, int jc, int kc,
                      int ihi, int jhi, int khi,
                      int comp) const
    {
        return limiter(phi(ilo,jlo,klo,comp),
                       phi(ic,jc,kc,comp),
                       phi(ihi,jhi,khi,comp));
    }

    template<class Phi, class Velocity>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar operator()(Phi const& phi,
                           Velocity const& velocity,
                           int i, int j, int k,
                           int comp,
                           const Set::Scalar* DX,
                           Options options = {},
                           std::array<Numeric::StencilType, AMREX_SPACEDIM> stencil = Numeric::DefaultType) const
    {
        (void)stencil;

        Set::Scalar div_phi_u = 0.0;
        Set::Scalar div_u = 0.0;

        Set::Scalar phi_c = phi(i,j,k,comp);

        Set::Scalar phi_xm = phi(i-1,j,k,comp);
        Set::Scalar phi_xp = phi(i+1,j,k,comp);
        Set::Scalar slope_xm = Slope(phi, i-2,j,k, i-1,j,k, i,j,k, comp);
        Set::Scalar slope_x  = Slope(phi, i-1,j,k, i,j,k, i+1,j,k, comp);
        Set::Scalar slope_xp = Slope(phi, i,j,k, i+1,j,k, i+2,j,k, comp);

        Set::Scalar ux_xlo = 0.5 * (velocity(i-1,j,k,0) + velocity(i,j,k,0));
        Set::Scalar ux_xhi = 0.5 * (velocity(i,j,k,0) + velocity(i+1,j,k,0));

        Set::Scalar phi_xlo_l = Util::Clamp(phi_xm + 0.5 * slope_xm, phi_xm, phi_c);
        Set::Scalar phi_xlo_r = Util::Clamp(phi_c  - 0.5 * slope_x,  phi_xm, phi_c);
        Set::Scalar phi_xhi_l = Util::Clamp(phi_c  + 0.5 * slope_x,  phi_c, phi_xp);
        Set::Scalar phi_xhi_r = Util::Clamp(phi_xp - 0.5 * slope_xp, phi_c, phi_xp);

        Set::Scalar flux_xlo = ux_xlo * (ux_xlo >= 0.0 ? phi_xlo_l : phi_xlo_r);
        Set::Scalar flux_xhi = ux_xhi * (ux_xhi >= 0.0 ? phi_xhi_l : phi_xhi_r);
        div_phi_u += (flux_xhi - flux_xlo) / DX[0];
        div_u += (ux_xhi - ux_xlo) / DX[0];

#if AMREX_SPACEDIM >= 2
        Set::Scalar phi_ym = phi(i,j-1,k,comp);
        Set::Scalar phi_yp = phi(i,j+1,k,comp);
        Set::Scalar slope_ym = Slope(phi, i,j-2,k, i,j-1,k, i,j,k, comp);
        Set::Scalar slope_y  = Slope(phi, i,j-1,k, i,j,k, i,j+1,k, comp);
        Set::Scalar slope_yp = Slope(phi, i,j,k, i,j+1,k, i,j+2,k, comp);

        Set::Scalar uy_ylo = 0.5 * (velocity(i,j-1,k,1) + velocity(i,j,k,1));
        Set::Scalar uy_yhi = 0.5 * (velocity(i,j,k,1) + velocity(i,j+1,k,1));

        Set::Scalar phi_ylo_l = Util::Clamp(phi_ym + 0.5 * slope_ym, phi_ym, phi_c);
        Set::Scalar phi_ylo_r = Util::Clamp(phi_c  - 0.5 * slope_y,  phi_ym, phi_c);
        Set::Scalar phi_yhi_l = Util::Clamp(phi_c  + 0.5 * slope_y,  phi_c, phi_yp);
        Set::Scalar phi_yhi_r = Util::Clamp(phi_yp - 0.5 * slope_yp, phi_c, phi_yp);

        Set::Scalar flux_ylo = uy_ylo * (uy_ylo >= 0.0 ? phi_ylo_l : phi_ylo_r);
        Set::Scalar flux_yhi = uy_yhi * (uy_yhi >= 0.0 ? phi_yhi_l : phi_yhi_r);
        div_phi_u += (flux_yhi - flux_ylo) / DX[1];
        div_u += (uy_yhi - uy_ylo) / DX[1];
#endif

#if AMREX_SPACEDIM == 3
        Set::Scalar phi_zm = phi(i,j,k-1,comp);
        Set::Scalar phi_zp = phi(i,j,k+1,comp);
        Set::Scalar slope_zm = Slope(phi, i,j,k-2, i,j,k-1, i,j,k, comp);
        Set::Scalar slope_z  = Slope(phi, i,j,k-1, i,j,k, i,j,k+1, comp);
        Set::Scalar slope_zp = Slope(phi, i,j,k, i,j,k+1, i,j,k+2, comp);

        Set::Scalar uz_zlo = 0.5 * (velocity(i,j,k-1,2) + velocity(i,j,k,2));
        Set::Scalar uz_zhi = 0.5 * (velocity(i,j,k,2) + velocity(i,j,k+1,2));

        Set::Scalar phi_zlo_l = Util::Clamp(phi_zm + 0.5 * slope_zm, phi_zm, phi_c);
        Set::Scalar phi_zlo_r = Util::Clamp(phi_c  - 0.5 * slope_z,  phi_zm, phi_c);
        Set::Scalar phi_zhi_l = Util::Clamp(phi_c  + 0.5 * slope_z,  phi_c, phi_zp);
        Set::Scalar phi_zhi_r = Util::Clamp(phi_zp - 0.5 * slope_zp, phi_c, phi_zp);

        Set::Scalar flux_zlo = uz_zlo * (uz_zlo >= 0.0 ? phi_zlo_l : phi_zlo_r);
        Set::Scalar flux_zhi = uz_zhi * (uz_zhi >= 0.0 ? phi_zhi_l : phi_zhi_r);
        div_phi_u += (flux_zhi - flux_zlo) / DX[2];
        div_u += (uz_zhi - uz_zlo) / DX[2];
#endif

        if (options.form == Form::Conservative) return -div_phi_u;
        return -div_phi_u + phi_c * div_u;
    }

    template<class Phi, class Velocity>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar Scalar(Phi const& phi,
                       Velocity const& velocity,
                       int i, int j, int k,
                       int comp,
                       const Set::Scalar* DX,
                       Options options = {},
                       std::array<Numeric::StencilType, AMREX_SPACEDIM> stencil = Numeric::DefaultType) const
    {
        return (*this)(phi, velocity, i, j, k, comp, DX, options, stencil);
    }
};
}
}

#endif
