//
// Provide an interpolator function that can work with templated fields.
//
// :bdg-danger-line:`Warning: This is a low-level class used to add templating functionality to unerlying AMReX classes. Edit at your own risk!`
//    
//


#ifndef NUMERIC_INTERPOLATOR_NODEBILINEAR_H
#define NUMERIC_INTERPOLATOR_NODEBILINEAR_H

#include <AMReX_Box.H>
#include <AMReX_BCRec.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuControl.H>
#include <AMReX_Interp_C.H>
#include <AMReX_Interpolater.H>

namespace Numeric
{
namespace Interpolator
{

template <class T>
class NodeBilinear: public amrex::NodeBilinear
{
public:
    void interp(const amrex::BaseFab<T>& crse,
        int crse_comp,
        amrex::BaseFab<T>& fine,
        int fine_comp,
        int ncomp,
        const amrex::Box& fine_region,
        const amrex::IntVect& ratio,
        const amrex::Geometry& /*crse_geom */,
        const amrex::Geometry& /*fine_geom */,
        amrex::Vector<amrex::BCRec> const& /*bcr*/,
        int /*actual_comp*/,
        int /*actual_state*/,
        amrex::RunOn runon)
    {
        // bool run_on_gpu = (runon == amrex::RunOn::Gpu && amrex::Gpu::inLaunchRegion());

        int num_slope = ncomp * (AMREX_D_TERM(2, *2, *2) - 1);
        const amrex::Box cslope_bx = amrex::enclosedCells(CoarseBox(fine_region, ratio));
        amrex::BaseFab<T> slopefab(cslope_bx, num_slope);
        // amrex::Elixir slopeeli;
        // if (run_on_gpu) slopeeli = slopefab.elixir();

        amrex::Array4<T const> const& crsearr = crse.const_array();
        amrex::Array4<T> const& finearr = fine.array();
        amrex::Array4<T> const& slopearr = slopefab.array();

        // AMReX's launch macro uses RunOn and Gpu without namespace qualifiers.
        using namespace amrex;

        AMREX_LAUNCH_HOST_DEVICE_LAMBDA_FLAG(runon, cslope_bx, tbx,
            {
                amrex::nodebilin_slopes<T>(tbx, slopearr, crsearr, crse_comp, ncomp, ratio);
            });

        AMREX_LAUNCH_HOST_DEVICE_LAMBDA_FLAG(runon, fine_region, tbx,
            {
                amrex::nodebilin_interp<T>(tbx, finearr, fine_comp, ncomp, slopearr, crsearr, crse_comp, ratio);
            });
    }
    virtual void interp_face(const amrex::BaseFab<T>& crse,
        const int        crse_comp,
        amrex::BaseFab<T>&       fine,
        const int        fine_comp,
        const int        ncomp,
        const amrex::Box&       fine_region,
        const amrex::IntVect&   ratio,
        const amrex::IArrayBox& /*solve_mask*/,
        const amrex::Geometry&  crse_geom,
        const amrex::Geometry&  fine_geom,
        amrex::Vector<amrex::BCRec> const& bcr,
        const int        bccomp,
        amrex::RunOn            gpu_or_cpu)
    {
#if AMREX_SPACEDIM == 1
        interp(crse, crse_comp, fine, fine_comp, ncomp, fine_region, ratio,
            crse_geom, fine_geom, bcr, bccomp, 0, gpu_or_cpu);
#else
        (void)crse;
        (void)crse_comp;
        (void)fine;
        (void)fine_comp;
        (void)ncomp;
        (void)fine_region;
        (void)ratio;
        (void)crse_geom;
        (void)fine_geom;
        (void)bcr;
        (void)bccomp;
        (void)gpu_or_cpu;
        Util::Abort("The version of this Interpolater for face-based data is not implemented or does not apply. Call 'interp' instead.");
#endif
    }

    virtual void interp_face(const amrex::FArrayBox& crse,
        const int        crse_comp,
        amrex::FArrayBox&       fine,
        const int        fine_comp,
        const int        ncomp,
        const amrex::Box&       fine_region,
        const amrex::IntVect&   ratio,
        const amrex::IArrayBox& /*solve_mask*/,
        const amrex::Geometry&  crse_geom,
        const amrex::Geometry&  fine_geom,
        amrex::Vector<amrex::BCRec> const& bcr,
        const int        bccomp,
        amrex::RunOn            gpu_or_cpu) override
    {
#if AMREX_SPACEDIM == 1
        amrex::NodeBilinear::interp(crse, crse_comp, fine, fine_comp, ncomp,
            fine_region, ratio, crse_geom, fine_geom, bcr, bccomp, 0,
            gpu_or_cpu);
#else
        Util::Abort("The version of this Interpolater for face-based data is not implemented or does not apply. Call 'interp' instead.");
#endif
        interp(crse, crse_comp, fine, fine_comp, ncomp, fine_region, ratio,
            crse_geom, fine_geom, bcr, bccomp, 0, gpu_or_cpu);
    }
private:
    using amrex::NodeBilinear::interp;
    using amrex::NodeBilinear::interp_face;

};
}
}

#endif
