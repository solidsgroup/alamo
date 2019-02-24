#ifndef OPERATOR_ELASTIC_ELASTIC_H_
#define OPERATOR_ELASTIC_ELASTIC_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

#include "Operator/Operator.H"


using namespace amrex;

namespace Operator
{
namespace CellElastic
{
class CellElastic : public Operator<Grid::Cell>
{
public:
	CellElastic () {}
	CellElastic (const Vector<Geometry>& a_geom,
		 const Vector<BoxArray>& a_grids,
		 const Vector<DistributionMapping>& a_dmap,
		 BC::BC& a_bc,
		 const LPInfo& a_info);
	virtual ~CellElastic ();
	CellElastic (const CellElastic&) = delete;
	CellElastic (CellElastic&&) = delete;
	CellElastic& operator= (const CellElastic&) = delete;
	CellElastic& operator= (CellElastic&&) = delete;

	void SetEigenstrain(amrex::Vector<std::unique_ptr<amrex::MultiFab> > &);
	void SetEigenstrain(amrex::Vector<amrex::MultiFab> &);

	void AddEigenstrainToRHS (FArrayBox& esrhs,
									int amrlev, const MFIter& mfi) const;
	void AddEigenstrainToRHS (amrex::Vector<amrex::MultiFab>& rhsfab) const;

	void Stress (FArrayBox& sigmafab,
		     const FArrayBox& ufab,
		     int amrlev, const MFIter& mfi,
		     bool voigt = false) const;
	void Energy (FArrayBox& energyfab,
				 const FArrayBox& ufab,
				 int amrlev, const MFIter& mfi) const;
protected:

	virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
	virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh, int redblack) const final;
	virtual void FFlux (int amrlev, const MFIter& mfi,
			    const Array<FArrayBox*,AMREX_SPACEDIM>& flux,
			    const FArrayBox& sol, Location loc, const int face_only=0) const;

	virtual int getNComp() const {return AMREX_SPACEDIM;};
	virtual bool isCrossStencil () const { return false; }

private:
	virtual amrex::Real C(const int i, const int j, const int k, const int l, const amrex::IntVect loc,
						  const int amrlev, const int mglev, const MFIter &mfi) const = 0;
	bool usingEigenstrain = false;
};
}
}
#endif
