#ifndef AMREX_OPERATOR_DIAGONAL_H_
#define AMREX_OPERATOR_DIAGONAL_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

#include "Operator/Operator.H"


using namespace amrex;

namespace Operator
{
namespace Diagonal
{
class Diagonal : public Operator
{
public:
	Diagonal () {}
	Diagonal (const Vector<Geometry>& a_geom,
		  const Vector<BoxArray>& a_grids,
		  const Vector<DistributionMapping>& a_dmap,
		  const LPInfo& a_info);
	virtual ~Diagonal ();
	Diagonal (const Diagonal&) = delete;
	Diagonal (Diagonal&&) = delete;
	Diagonal& operator= (const Diagonal&) = delete;
	Diagonal& operator= (Diagonal&&) = delete;
protected:

	virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
	virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh) const final;
	virtual int getNComp() const {return AMREX_SPACEDIM;};
};
}
}
#endif
