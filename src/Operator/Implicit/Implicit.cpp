#include <AMReX_MultiFabUtil.H>
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>
#include "Set/Set.H"

#include <AMReX_ArrayLim.H>

#include "Util/Util.H"
#include "Set/Set.H"
#include "Implicit.H"

namespace Operator
{
namespace Implicit
{
Implicit::Implicit (const Vector<Geometry>& a_geom,
		  const Vector<BoxArray>& a_grids,
		  const Vector<DistributionMapping>& a_dmap,
		  BC::BC& a_bc,
		  const LPInfo& a_info)
{
	define(a_geom, a_grids, a_dmap, a_bc, a_info);
}

void
Implicit::Fapply (int /*amrlev*/, ///<[in] AMR Level
		  int /*mglev*/,  ///<[in]
		  MultiFab& /*f*/,///<[out] The force vector
		  const MultiFab& /*u*/ ///<[in] The displacements vector
		 ) const
{
	Util::Message(INFO);
}


void
Implicit::Fsmooth (int /*amrlev*/,          ///<[in] AMR level
		   int /*mglev*/,           ///<[in]
		   MultiFab& /*u*/,       ///<[inout] Solution (displacement field)
		  const MultiFab& /*rhs*/, ///<[in] Body force vectors (rhs=right hand side)
		   int /*redblack*/         ///<[in] Smooth even vs. odd modes
		  ) const
{
	Util::Message(INFO);
}

void Implicit::FFlux (int /*amrlev*/, const MFIter& /*mfi*/,
		     const Array<FArrayBox*,AMREX_SPACEDIM>& sigmafab,
		     const FArrayBox& /*sol*/, Location /*loc*/, const int /*face_only*/) const
{
	amrex::BaseFab<amrex::Real> AMREX_D_DECL( &fxfab = *sigmafab[0],
	 					  &fyfab = *sigmafab[1],
	 					  &fzfab = *sigmafab[2] ) ;
	AMREX_D_TERM(fxfab.setVal(0.0);,
	 	     fyfab.setVal(0.0);,
	 	     fzfab.setVal(0.0););

}

}
}
