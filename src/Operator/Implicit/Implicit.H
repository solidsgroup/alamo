#ifndef OPERATOR_IMPLICIT_IMPLICIT_H_
#define OPERATOR_IMPLICIT_IMPLICIT_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

#include "Operator/Operator.H"


using namespace amrex;

namespace Operator
{
namespace Implicit
{
class Implicit : public Operator<Grid::Cell>
{
public:
    Implicit () {}
    Implicit (const Vector<Geometry>& a_geom,
        const Vector<BoxArray>& a_grids,
        const Vector<DistributionMapping>& a_dmap,
        BC::BC<Set::Scalar>& a_bc,
        const LPInfo& a_info);
    virtual ~Implicit () {};
    Implicit (const Implicit&) = delete;
    Implicit (Implicit&&) = delete;
    Implicit& operator= (const Implicit&) = delete;
    Implicit& operator= (Implicit&&) = delete;

protected:

    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
    virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh, int redblack) const final;
    virtual void FFlux (int amrlev, const MFIter& mfi,
                const Array<FArrayBox*,AMREX_SPACEDIM>& flux,
                const FArrayBox& sol, Location loc, const int face_only=0) const;
    virtual int getNComp() const {return AMREX_SPACEDIM;};
    virtual bool isCrossStencil () const { return false; }
private:
};
}
}
#endif
