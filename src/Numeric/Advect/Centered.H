//
// This header defines :code:`Numeric::Advect::Centered`, a second-order centered
// finite-volume advection operator for collocated cell-centered fields.
//
// The operator is written as a GPU-safe functor:
//
// .. code-block:: c++
//
//    rhs = centered(phi, velocity, i, j, k, comp, DX, options, stencil);
//
// It computes face velocities by averaging neighboring cell-centered velocity
// values, reconstructs face values of :math:`\phi` by arithmetic averaging, and
// returns either :math:`-\nabla\cdot(\phi\mathbf{u})` or
// :math:`-\nabla\cdot(\phi\mathbf{u}) + \phi\nabla\cdot\mathbf{u}` depending on
// :code:`Numeric::Advect::Options::form`.
//
// This implementation assumes collocated cell-centered data:
// :code:`phi_location = Set::HC::Cell` and
// :code:`velocity_location = Set::HC::Cell`.  It requires one ghost cell in each
// coordinate direction.
//

#ifndef NUMERIC_ADVECT_CENTERED_H
#define NUMERIC_ADVECT_CENTERED_H

#include "Advect.H"

namespace Numeric
{
namespace Advect
{
class Centered
{
public:
    static constexpr const char* name = "centered";
    static constexpr Set::HC phi_location = Set::HC::Cell;
    static constexpr Set::HC velocity_location = Set::HC::Cell;

    static void Parse(Centered&, IO::ParmParse&) {}

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    int NGhost() const
    {
        return 1;
    }

    template<class Phi>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    FaceStates ReconstructFace(Phi const& phi,
                               int i, int j, int k,
                               int direction,
                               int comp) const
    {
        const int di = direction == 0 ? 1 : 0;
        const int dj = direction == 1 ? 1 : 0;
        const int dk = direction == 2 ? 1 : 0;
        return {phi(i,j,k,comp), phi(i+di,j+dj,k+dk,comp)};
    }

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar FaceValue(FaceStates const& state,
                          Set::Scalar /*face_velocity*/) const
    {
        return 0.5 * (state.left + state.right);
    }

    template<class Phi>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar ReconstructTransportedFaceValue(
        Phi const& phi,
        int i, int j, int k,
        int direction,
        int comp,
        Set::Scalar /*face_velocity*/) const
    {
        const int di = direction == 0 ? 1 : 0;
        const int dj = direction == 1 ? 1 : 0;
        const int dk = direction == 2 ? 1 : 0;
        return 0.5 *
            (phi(i,j,k,comp) + phi(i+di,j+dj,k+dk,comp));
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
        Set::Scalar flux_xlo = ux_xlo * 0.5 * (phi(i-1,j,k,comp) + phi_c);
        Set::Scalar flux_xhi = ux_xhi * 0.5 * (phi_c + phi(i+1,j,k,comp));
        div_phi_u += (flux_xhi - flux_xlo) / DX[0];
        div_u += (ux_xhi - ux_xlo) / DX[0];

#if AMREX_SPACEDIM >= 2
        Set::Scalar uy_ylo = 0.5 * (velocity(i,j-1,k,1) + velocity(i,j,k,1));
        Set::Scalar uy_yhi = 0.5 * (velocity(i,j,k,1) + velocity(i,j+1,k,1));
        Set::Scalar flux_ylo = uy_ylo * 0.5 * (phi(i,j-1,k,comp) + phi_c);
        Set::Scalar flux_yhi = uy_yhi * 0.5 * (phi_c + phi(i,j+1,k,comp));
        div_phi_u += (flux_yhi - flux_ylo) / DX[1];
        div_u += (uy_yhi - uy_ylo) / DX[1];
#endif

#if AMREX_SPACEDIM == 3
        Set::Scalar uz_zlo = 0.5 * (velocity(i,j,k-1,2) + velocity(i,j,k,2));
        Set::Scalar uz_zhi = 0.5 * (velocity(i,j,k,2) + velocity(i,j,k+1,2));
        Set::Scalar flux_zlo = uz_zlo * 0.5 * (phi(i,j,k-1,comp) + phi_c);
        Set::Scalar flux_zhi = uz_zhi * 0.5 * (phi_c + phi(i,j,k+1,comp));
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
