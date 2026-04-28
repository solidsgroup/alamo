#ifndef BC_LEVELSETNEUMANN_H_
#define BC_LEVELSETNEUMANN_H_

#include "BC/BC.H"
#include "Set/Set.H"
#include "AMReX_LO_BCTYPES.H"

namespace BC {

class LevelSetNeumann final : public BC<Set::Scalar>
{
public:
    // IMPORTANT: no geometry passed here
    LevelSetNeumann() = default;

    // Homogeneous Neumann: copy nearest interior value
    void FillBoundary(amrex::BaseFab<Set::Scalar>& in,
                      const amrex::Box& box,
                      int ngrow,
                      int dcomp,
                      int ncomp,
                      amrex::Real /*time*/,
                      Orientation /*face*/ = Orientation::All,
                      const amrex::Mask* /*mask*/ = nullptr) override
    {
        const amrex::Box& domain = m_geom.Domain();

        // Tile + ghost region
        amrex::Box gbox = amrex::grow(box, ngrow);

        auto arr = in.array();

        const amrex::Dim3 lo = amrex::lbound(domain);
        const amrex::Dim3 hi = amrex::ubound(domain);

        for (int n = 0; n < ncomp; ++n)
        {
            amrex::ParallelFor(gbox,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                // Clamp to valid interior cell
                int ii = amrex::max(lo.x, amrex::min(i, hi.x));
                int jj = amrex::max(lo.y, amrex::min(j, hi.y));
                int kk = amrex::max(lo.z, amrex::min(k, hi.z));

                // If inside domain, do nothing
                if (i == ii && j == jj && k == kk) return;

                // Homogeneous Neumann (zero normal gradient)
                arr(i,j,k,dcomp+n) = arr(ii,jj,kk,dcomp+n);
            });
        }
    }

    // Tell AMReX what type of BC this is
    amrex::BCRec GetBCRec() override
    {
        int lo[AMREX_SPACEDIM];
        int hi[AMREX_SPACEDIM];

        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
            // For ghost filling (not MLMG), this is the correct type
            lo[d] = amrex::BCType::foextrap;
            hi[d] = amrex::BCType::foextrap;
        }

        return amrex::BCRec(lo, hi);
    }
};

} // namespace BC

#endif