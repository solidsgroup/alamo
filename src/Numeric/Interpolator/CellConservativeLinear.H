#ifndef NUMERIC_INTERPOLATOR_CELLCONSERVATIVELINEAR_H
#define NUMERIC_INTERPOLATOR_CELLCONSERVATIVELINEAR_H

#include <AMReX_Box.H>
#include <AMReX_BCRec.H>
#include <AMReX_REAL.H>
#include <AMReX_GpuControl.H>
#include <AMReX_Interp_C.H>
#include <AMReX_MFInterp_C.H>
#include <AMReX_Interpolater.H>

namespace Numeric
{
namespace Interpolator
{

template <class T>
class CellConservativeLinear: public amrex::CellConservativeLinear
{
public:
    void interp(const amrex::BaseFab<T>& crse,
                int crse_comp,
                amrex::BaseFab<T>& fine,
                int fine_comp,
                int ncomp,
                const amrex::Box& fine_region,
                const amrex::IntVect& ratio,
                const amrex::Geometry& crse_geom ,
                const amrex::Geometry& fine_geom ,
                amrex::Vector<amrex::BCRec> const& bcr,
                int /*actual_comp*/,
                int /*actual_state*/,
                amrex::RunOn runon)
    {


        Util::Abort(INFO);

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
                             amrex::RunOn            /*runon*/)
    {
        Util::Abort(INFO);
        

    }
private:
    using amrex::CellConservativeLinear::interp;
    using amrex::CellConservativeLinear::interp_face;

};
}
}

#endif
