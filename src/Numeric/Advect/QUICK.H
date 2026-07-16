//
// This header defines :code:`Numeric::Advect::QUICK`, a third-order upwind-biased
// advection operator for collocated cell-centered fields.
//
// The operator is written as a GPU-safe functor:
//
// .. code-block:: c++
//
//    rhs = quick(phi, velocity, i, j, k, comp, DX, options, stencil);
//
// It computes face velocities by averaging neighboring cell-centered velocity
// values, reconstructs face values of :math:`\phi` with the QUICK stencil, and
// returns either :math:`-\nabla\cdot(\phi\mathbf{u})` or
// :math:`-\nabla\cdot(\phi\mathbf{u}) + \phi\nabla\cdot\mathbf{u}` depending on
// :code:`Numeric::Advect::Options::form`.
//
// This implementation assumes collocated cell-centered data:
// :code:`phi_location = Set::HC::Cell` and
// :code:`velocity_location = Set::HC::Cell`.  It requires two ghost cells in each
// coordinate direction.
//

#ifndef NUMERIC_ADVECT_QUICK_H
#define NUMERIC_ADVECT_QUICK_H

#include "Advect.H"

namespace Numeric
{
namespace Advect
{
class QUICK
{
public:
    static constexpr const char* name = "quick";
    static constexpr Set::HC phi_location = Set::HC::Cell;
    static constexpr Set::HC velocity_location = Set::HC::Cell;

    static void Parse(QUICK&, IO::ParmParse&) {}

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    int NGhost() const
    {
        return 2;
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar Reconstruct(Set::Scalar up2,
                            Set::Scalar up1,
                            Set::Scalar down1) const
    {
        return (-up2 + 6.0 * up1 + 3.0 * down1) / 8.0;
    }

    template<class Phi, class Velocity>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar operator()( Phi const& phi,
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
        Set::Scalar flux_xlo = ux_xlo * (ux_xlo >= 0.0 ?
            Reconstruct(phi(i-2,j,k,comp), phi(i-1,j,k,comp), phi_c) :
            Reconstruct(phi(i+1,j,k,comp), phi_c, phi(i-1,j,k,comp)));
        Set::Scalar flux_xhi = ux_xhi * (ux_xhi >= 0.0 ?
            Reconstruct(phi(i-1,j,k,comp), phi_c, phi(i+1,j,k,comp)) :
            Reconstruct(phi(i+2,j,k,comp), phi(i+1,j,k,comp), phi_c));
        div_phi_u += (flux_xhi - flux_xlo) / DX[0];
        div_u += (ux_xhi - ux_xlo) / DX[0];

#if AMREX_SPACEDIM >= 2
        Set::Scalar uy_ylo = 0.5 * (velocity(i,j-1,k,1) + velocity(i,j,k,1));
        Set::Scalar uy_yhi = 0.5 * (velocity(i,j,k,1) + velocity(i,j+1,k,1));
        Set::Scalar flux_ylo = uy_ylo * (uy_ylo >= 0.0 ?
            Reconstruct(phi(i,j-2,k,comp), phi(i,j-1,k,comp), phi_c) :
            Reconstruct(phi(i,j+1,k,comp), phi_c, phi(i,j-1,k,comp)));
        Set::Scalar flux_yhi = uy_yhi * (uy_yhi >= 0.0 ?
            Reconstruct(phi(i,j-1,k,comp), phi_c, phi(i,j+1,k,comp)) :
            Reconstruct(phi(i,j+2,k,comp), phi(i,j+1,k,comp), phi_c));
        div_phi_u += (flux_yhi - flux_ylo) / DX[1];
        div_u += (uy_yhi - uy_ylo) / DX[1];
#endif

#if AMREX_SPACEDIM == 3
        Set::Scalar uz_zlo = 0.5 * (velocity(i,j,k-1,2) + velocity(i,j,k,2));
        Set::Scalar uz_zhi = 0.5 * (velocity(i,j,k,2) + velocity(i,j,k+1,2));
        Set::Scalar flux_zlo = uz_zlo * (uz_zlo >= 0.0 ?
            Reconstruct(phi(i,j,k-2,comp), phi(i,j,k-1,comp), phi_c) :
            Reconstruct(phi(i,j,k+1,comp), phi_c, phi(i,j,k-1,comp)));
        Set::Scalar flux_zhi = uz_zhi * (uz_zhi >= 0.0 ?
            Reconstruct(phi(i,j,k-1,comp), phi_c, phi(i,j,k+1,comp)) :
            Reconstruct(phi(i,j,k+2,comp), phi(i,j,k+1,comp), phi_c));
        div_phi_u += (flux_zhi - flux_zlo) / DX[2];
        div_u += (uz_zhi - uz_zlo) / DX[2];
#endif

        if (options.form == Form::Conservative) return -div_phi_u;
        return -div_phi_u + phi_c * div_u;
    }

    template<class Phi, class Velocity>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar Scalar( Phi const& phi,
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
