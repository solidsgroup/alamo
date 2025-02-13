//
// This is the mechanism for impositing boundary conditions on :code:`Set::Field` objects of 
// scalar type.
// Typical convention is for :code:`[prefix]` to be given by the field name. For instance,
// 
// .. code-block:: make
// 
//                 bc.temp.type.xhi = dirichlet dirichlet dirichlet
//                 bc.temp.val.xhi  = 0.0 1.0 0.0
//
// corresponds to the boundary condition for temperature. See the specific integrator
// for details.
//

#ifndef BC_BC_H_
#define BC_BC_H_

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BCRec.H>
#include <AMReX_PhysBCFunct.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_Mask.H>
#include <AMReX_Periodicity.H>

#include "Util/Util.H"
#include "Set/Set.H"

/// \brief Collection of boundary condition (BC) objects
namespace BC
{
enum Orientation {
    All = -1,
    AMREX_D_DECL(xlo = 0,
    ylo = 1,
    zlo = 2),
    AMREX_D_DECL(xhi = 0 + AMREX_SPACEDIM,
    yhi = 1 + AMREX_SPACEDIM,
    zhi = 2 + AMREX_SPACEDIM)
};

template<class T>
class BC
{

public:
    virtual ~BC() {};

    void define(const amrex::Geometry& a_geom) { m_geom = a_geom; };

    virtual void FillBoundary(amrex::BaseFab<T>& in,
        const amrex::Box& box,
        int ngrow, int dcomp, int ncomp,
        amrex::Real time,
        Orientation face = Orientation::All,
        const amrex::Mask* mask = nullptr) = 0;

    virtual void FillBoundary(amrex::FabArray<amrex::BaseFab<T>>& mf,
        int dcomp, int ncomp,
        amrex::Real time, int /*bccomp*/)
    {
        mf.FillBoundary(m_geom.periodicity());
        for (amrex::MFIter mfi(mf, true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& box = mfi.tilebox();
            amrex::BaseFab<T>& in = mf[mfi];
            FillBoundary(in, box, mf.nGrow(), dcomp, ncomp, time);
        }
    }

    virtual void FillBoundary(amrex::FabArray<amrex::BaseFab<T>>& mf,
        int dcomp, int ncomp, amrex::IntVect const& /*nghost*/,
        amrex::Real time, int bccomp) //override
    {
        FillBoundary(mf, dcomp, ncomp, time, bccomp);
    }

    void operator () (amrex::FabArray<amrex::BaseFab<T>>& mf,
        int dcomp, int ncomp, amrex::IntVect const& /*nghost*/,
        amrex::Real time, int bccomp)
    {
        FillBoundary(mf, dcomp, ncomp, time, bccomp);
    }



    // template <class Q = T>
    // typename std::enable_if<std::is_same<Q, amrex::Real>::value>::type
    void FillBoundary(amrex::MultiFab& mf,
            int dcomp, int ncomp,
            amrex::Real time, int /*bccomp*/)
    {
        mf.FillBoundary(m_geom.periodicity());
        for (amrex::MFIter mfi(mf, true); mfi.isValid(); ++mfi)
        {
            const amrex::Box& box = mfi.tilebox();
            amrex::BaseFab<T>& in = mf[mfi];
            FillBoundary(in, box, mf.nGrow(), dcomp, ncomp, time);
        }
    }

    // template <class Q = T>
    // typename std::enable_if<std::is_same<Q, amrex::Real>::value>::type
    void operator () (amrex::MultiFab& mf,
            int dcomp, int ncomp, amrex::IntVect const& /*nghost*/,
            amrex::Real time, int bccomp)
    {
        FillBoundary(mf, dcomp, ncomp, time, bccomp);
    }


    virtual amrex::BCRec GetBCRec() = 0;

    virtual amrex::Array<int, AMREX_SPACEDIM> IsPeriodic()
    {
        return { {AMREX_D_DECL(m_geom.isPeriodic(0),m_geom.isPeriodic(1),m_geom.isPeriodic(2))} };
    }

    virtual amrex::Periodicity Periodicity() const {
        return m_geom.periodicity();
    }

    virtual amrex::Periodicity Periodicity(const amrex::Box& b) {
        return m_geom.periodicity(b);
    }


protected:
    amrex::Geometry m_geom;
    //std::vector<amrex::Geometry> &geom;
    //int lev=0;
};

namespace BCUtil
{
int ReadString(std::string bcstring);
bool IsPeriodic(int bctype);
bool IsNeumann(int bctype);
bool IsReflectEven(int bctype);
bool IsReflectOdd(int bctype);
bool IsDirichlet(int bctype);
}

}
#endif
