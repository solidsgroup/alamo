//
// This header defines :code:`Numeric::Advect::WENO5`, a fifth-order finite-volume
// advection operator using Jiang-Shu weighted essentially non-oscillatory
// reconstruction.
//
// The operator is written as a GPU-safe functor:
//
// .. code-block:: c++
//
//    rhs = weno5(phi, velocity, i, j, k, comp, DX, options, stencil);
//
// It computes face velocities by averaging neighboring cell-centered velocity
// values, reconstructs left and right states of :math:`\phi` at each face using
// WENO5, selects the upwind state, and returns either
// :math:`-\nabla\cdot(\phi\mathbf{u})` or
// :math:`-\nabla\cdot(\phi\mathbf{u}) + \phi\nabla\cdot\mathbf{u}` depending on
// :code:`Numeric::Advect::Options::form`.
//
// This implementation assumes collocated cell-centered data:
// :code:`phi_location = Set::HC::Cell` and
// :code:`velocity_location = Set::HC::Cell`.  It requires three ghost cells in each
// coordinate direction.
//

#ifndef NUMERIC_ADVECT_WENO5_H
#define NUMERIC_ADVECT_WENO5_H

#include "Advect.H"

namespace Numeric
{
namespace Advect
{
class WENO5
{
public:
    static constexpr const char* name = "weno5";
    static constexpr Set::HC phi_location = Set::HC::Cell;
    static constexpr Set::HC velocity_location = Set::HC::Cell;

    Set::Scalar eps = 1.0e-12;

    static void Parse(WENO5& value, IO::ParmParse& pp)
    {
        pp.query_default("eps", value.eps, 1.0e-12);
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    int NGhost() const
    {
        return 3;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ReconstructLeft(Set::Scalar fm2,
                                Set::Scalar fm1,
                                Set::Scalar f0,
                                Set::Scalar fp1,
                                Set::Scalar fp2) const
    {
        Set::Scalar p0 = (2.0 * fm2 - 7.0 * fm1 + 11.0 * f0) / 6.0;
        Set::Scalar p1 = (-fm1 + 5.0 * f0 + 2.0 * fp1) / 6.0;
        Set::Scalar p2 = (2.0 * f0 + 5.0 * fp1 - fp2) / 6.0;

        Set::Scalar b0 = (13.0 / 12.0) * (fm2 - 2.0 * fm1 + f0) * (fm2 - 2.0 * fm1 + f0) +
                         0.25 * (fm2 - 4.0 * fm1 + 3.0 * f0) * (fm2 - 4.0 * fm1 + 3.0 * f0);
        Set::Scalar b1 = (13.0 / 12.0) * (fm1 - 2.0 * f0 + fp1) * (fm1 - 2.0 * f0 + fp1) +
                         0.25 * (fm1 - fp1) * (fm1 - fp1);
        Set::Scalar b2 = (13.0 / 12.0) * (f0 - 2.0 * fp1 + fp2) * (f0 - 2.0 * fp1 + fp2) +
                         0.25 * (3.0 * f0 - 4.0 * fp1 + fp2) * (3.0 * f0 - 4.0 * fp1 + fp2);

        Set::Scalar a0 = 0.1 / ((eps + b0) * (eps + b0));
        Set::Scalar a1 = 0.6 / ((eps + b1) * (eps + b1));
        Set::Scalar a2 = 0.3 / ((eps + b2) * (eps + b2));
        Set::Scalar asum = a0 + a1 + a2;
        return (a0 * p0 + a1 * p1 + a2 * p2) / asum;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ReconstructRight(Set::Scalar fm1,
                                 Set::Scalar f0,
                                 Set::Scalar fp1,
                                 Set::Scalar fp2,
                                 Set::Scalar fp3) const
    {
        Set::Scalar p0 = (2.0 * fp3 - 7.0 * fp2 + 11.0 * fp1) / 6.0;
        Set::Scalar p1 = (-fp2 + 5.0 * fp1 + 2.0 * f0) / 6.0;
        Set::Scalar p2 = (2.0 * fp1 + 5.0 * f0 - fm1) / 6.0;

        Set::Scalar b0 = (13.0 / 12.0) * (fp3 - 2.0 * fp2 + fp1) * (fp3 - 2.0 * fp2 + fp1) +
                         0.25 * (fp3 - 4.0 * fp2 + 3.0 * fp1) * (fp3 - 4.0 * fp2 + 3.0 * fp1);
        Set::Scalar b1 = (13.0 / 12.0) * (fp2 - 2.0 * fp1 + f0) * (fp2 - 2.0 * fp1 + f0) +
                         0.25 * (fp2 - f0) * (fp2 - f0);
        Set::Scalar b2 = (13.0 / 12.0) * (fp1 - 2.0 * f0 + fm1) * (fp1 - 2.0 * f0 + fm1) +
                         0.25 * (3.0 * fp1 - 4.0 * f0 + fm1) * (3.0 * fp1 - 4.0 * f0 + fm1);

        Set::Scalar a0 = 0.1 / ((eps + b0) * (eps + b0));
        Set::Scalar a1 = 0.6 / ((eps + b1) * (eps + b1));
        Set::Scalar a2 = 0.3 / ((eps + b2) * (eps + b2));
        Set::Scalar asum = a0 + a1 + a2;
        return (a0 * p0 + a1 * p1 + a2 * p2) / asum;
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

        const Set::Scalar phi_c = phi(i,j,k,comp);
        Set::Scalar div_phi_u = 0.0;
        Set::Scalar div_u = 0.0;

        Set::Scalar ux_xlo = 0.5 * (velocity(i-1,j,k,0) + velocity(i,j,k,0));
        Set::Scalar ux_xhi = 0.5 * (velocity(i,j,k,0) + velocity(i+1,j,k,0));
        Set::Scalar phi_xlo_l = ReconstructLeft(phi(i-3,j,k,comp), phi(i-2,j,k,comp), phi(i-1,j,k,comp), phi_c, phi(i+1,j,k,comp));
        Set::Scalar phi_xlo_r = ReconstructRight(phi(i-2,j,k,comp), phi(i-1,j,k,comp), phi_c, phi(i+1,j,k,comp), phi(i+2,j,k,comp));
        Set::Scalar phi_xhi_l = ReconstructLeft(phi(i-2,j,k,comp), phi(i-1,j,k,comp), phi_c, phi(i+1,j,k,comp), phi(i+2,j,k,comp));
        Set::Scalar phi_xhi_r = ReconstructRight(phi(i-1,j,k,comp), phi_c, phi(i+1,j,k,comp), phi(i+2,j,k,comp), phi(i+3,j,k,comp));
        Set::Scalar flux_xlo = ux_xlo * (ux_xlo >= 0.0 ? phi_xlo_l : phi_xlo_r);
        Set::Scalar flux_xhi = ux_xhi * (ux_xhi >= 0.0 ? phi_xhi_l : phi_xhi_r);
        div_phi_u += (flux_xhi - flux_xlo) / DX[0];
        div_u += (ux_xhi - ux_xlo) / DX[0];

#if AMREX_SPACEDIM >= 2
        Set::Scalar uy_ylo = 0.5 * (velocity(i,j-1,k,1) + velocity(i,j,k,1));
        Set::Scalar uy_yhi = 0.5 * (velocity(i,j,k,1) + velocity(i,j+1,k,1));
        Set::Scalar phi_ylo_l = ReconstructLeft(phi(i,j-3,k,comp), phi(i,j-2,k,comp), phi(i,j-1,k,comp), phi_c, phi(i,j+1,k,comp));
        Set::Scalar phi_ylo_r = ReconstructRight(phi(i,j-2,k,comp), phi(i,j-1,k,comp), phi_c, phi(i,j+1,k,comp), phi(i,j+2,k,comp));
        Set::Scalar phi_yhi_l = ReconstructLeft(phi(i,j-2,k,comp), phi(i,j-1,k,comp), phi_c, phi(i,j+1,k,comp), phi(i,j+2,k,comp));
        Set::Scalar phi_yhi_r = ReconstructRight(phi(i,j-1,k,comp), phi_c, phi(i,j+1,k,comp), phi(i,j+2,k,comp), phi(i,j+3,k,comp));
        Set::Scalar flux_ylo = uy_ylo * (uy_ylo >= 0.0 ? phi_ylo_l : phi_ylo_r);
        Set::Scalar flux_yhi = uy_yhi * (uy_yhi >= 0.0 ? phi_yhi_l : phi_yhi_r);
        div_phi_u += (flux_yhi - flux_ylo) / DX[1];
        div_u += (uy_yhi - uy_ylo) / DX[1];
#endif

#if AMREX_SPACEDIM == 3
        Set::Scalar uz_zlo = 0.5 * (velocity(i,j,k-1,2) + velocity(i,j,k,2));
        Set::Scalar uz_zhi = 0.5 * (velocity(i,j,k,2) + velocity(i,j,k+1,2));
        Set::Scalar phi_zlo_l = ReconstructLeft(phi(i,j,k-3,comp), phi(i,j,k-2,comp), phi(i,j,k-1,comp), phi_c, phi(i,j,k+1,comp));
        Set::Scalar phi_zlo_r = ReconstructRight(phi(i,j,k-2,comp), phi(i,j,k-1,comp), phi_c, phi(i,j,k+1,comp), phi(i,j,k+2,comp));
        Set::Scalar phi_zhi_l = ReconstructLeft(phi(i,j,k-2,comp), phi(i,j,k-1,comp), phi_c, phi(i,j,k+1,comp), phi(i,j,k+2,comp));
        Set::Scalar phi_zhi_r = ReconstructRight(phi(i,j,k-1,comp), phi_c, phi(i,j,k+1,comp), phi(i,j,k+2,comp), phi(i,j,k+3,comp));
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
