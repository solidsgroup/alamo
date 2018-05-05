#ifndef AMREX_OPERATOR_FEM_H_
#define AMREX_OPERATOR_FEM_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

#include "Operator/Operator.H"

#include "Model/Solid/Solid.H"

#include "Element/Q4.H"

using namespace amrex;

namespace Operator
{
namespace FEM
{

class FEM : public Operator
{
public:
  FEM (Model::Solid::Solid &);
  // FEM (const Vector<Geometry>& a_geom,
  // 	   const Vector<BoxArray>& a_grids,
  // 	   const Vector<DistributionMapping>& a_dmap,
  // 	   const LPInfo& a_info);
  virtual ~FEM ();
  FEM (const FEM&) = delete;
  FEM (FEM&&) = delete;
  FEM& operator= (const FEM&) = delete;
  FEM& operator= (FEM&&) = delete;

  void Stress (FArrayBox& sigmafab,
	       const FArrayBox& ufab,
	       int amrlev, const MFIter& mfi) const;

  template<int n>
  void Energy (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const;

  virtual void smooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh, int redblack) const;

protected:
  virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
  virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh, int redblack) const final;
  virtual void FFlux (int amrlev, const MFIter& mfi,
		      const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
		      const FArrayBox& sol, const int face_only=0) const final;
  virtual int getNComp() const {return AMREX_SPACEDIM;};
  virtual bool isCrossStencil () const { return false; }

private:
  virtual amrex::Real C(const int i, const int j, const int k, const int l,
			const amrex::IntVect loc, const int amrlev, const int mglev, const MFIter &mfi)
    const;

  Model::Solid::Solid &model;
  //mutable amrex::Vector<amrex::Vector<amrex::MultiFab> > rho,aa;

};
}
}
#endif
