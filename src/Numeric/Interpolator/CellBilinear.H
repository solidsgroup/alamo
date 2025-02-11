//
// Provide an interpolator function that can work with templated fields.
//
// :bdg-danger-line:`Warning: This is a low-level class used to add templating functionality to unerlying AMReX classes. Edit at your own risk!`
//    
//


#ifndef NUMERIC_INTERPOLATOR_CELLBILINEAR_H
#define NUMERIC_INTERPOLATOR_CELLBILINEAR_H

#include <AMReX_Box.H>
#include <AMReX_BCRec.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuControl.H>
#include <AMReX_MFInterp_C.H>
#include <AMReX_Interpolater.H>

namespace Numeric
{
namespace Interpolator
{

template <class T>
class CellBilinear: public amrex::CellBilinear
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
        amrex::Array4<T const> const& crsearr = crse.const_array();
        amrex::Array4<T> const& finearr = fine.array();

        AMREX_HOST_DEVICE_PARALLEL_FOR_4D_FLAG(runon, fine_region, ncomp, i,j,k,n,
        {
            mf_cell_bilin_interp(i,j,k,n,finearr,fine_comp,crsearr,crse_comp,ratio);
        });
    }
    virtual void interp_face(const amrex::BaseFab<T>& /*crse*/,
        const int        /*crse_comp*/,
        amrex::BaseFab<T>&       /*fine*/,
        const int        /*fine_comp*/,
        const int        /*ncomp*/,
        const amrex::Box&       /*fine_region*/,
        const amrex::IntVect&   /*ratio*/,
        const amrex::IArrayBox& /*solve_mask*/,
        const amrex::Geometry&  /*crse_geom*/,
        const amrex::Geometry&  /*fine_geom*/,
        amrex::Vector<amrex::BCRec> const& /*bcr*/,
        const int        /*bccomp*/,
        amrex::RunOn            /*gpu_or_cpu*/)
    {
        Util::Abort("The version of this Interpolater for face-based data is not implemented or does not apply. Call 'interp' instead.");
    }
private:
    using amrex::CellBilinear::interp;
    using amrex::CellBilinear::interp_face;

};
}
}

#endif
