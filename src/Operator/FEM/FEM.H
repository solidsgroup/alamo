#ifndef AMREX_OPERATOR_FEM_H_
#define AMREX_OPERATOR_FEM_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

#include "Operator/Operator.H"


using namespace amrex;

namespace Operator
{
namespace FEM
{
class FEM : public Operator
{
public:
  FEM () {}
  FEM (const Vector<Geometry>& a_geom,
	   const Vector<BoxArray>& a_grids,
	   const Vector<DistributionMapping>& a_dmap,
	   const LPInfo& a_info);
  virtual ~FEM ();
  FEM (const FEM&) = delete;
  FEM (FEM&&) = delete;
  FEM& operator= (const FEM&) = delete;
  FEM& operator= (FEM&&) = delete;

  void Stress (FArrayBox& sigmafab,
	       const FArrayBox& ufab,
	       int amrlev, const MFIter& mfi) const;
  void Energy (FArrayBox& energyfab,
	       const FArrayBox& ufab,
	       int amrlev, const MFIter& mfi) const;

protected:
  virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
  virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh, int redblack) const final;
  virtual void FFlux (int amrlev, const MFIter& mfi,
		      const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
		      const FArrayBox& sol, const int face_only=0) const final;
  virtual int getNComp() const {return AMREX_SPACEDIM;};
  virtual bool isCrossStencil () const { return false; }

private:
  virtual amrex::Real C(const int i, const int j, const int k, const int l, const amrex::IntVect loc,
			const int amrlev, const int mglev, const MFIter &mfi) const = 0;

  constexpr static amrex::Real GradPhiQ[4][2][4] = // GradPhiQ[phi_n][d/dx_i][quadrature point]
    {{{0.2113248655, -0.2113248655, 0.7886751345, -0.7886751345}, {0.2113248655, 0.7886751345, -0.2113248655, -0.7886751345}},
     {{0.2113248655, -0.2113248655, 0.7886751345, -0.7886751345}, {0.7886751345, 0.2113248655, -0.7886751345, -0.2113248655}},
     {{0.7886751345, -0.7886751345, 0.2113248655, -0.2113248655}, {0.2113248655, 0.7886751345, -0.2113248655, -0.7886751345}},
     {{0.7886751345, -0.7886751345, 0.2113248655, -0.2113248655}, {0.7886751345, 0.2113248655, -0.7886751345, -0.2113248655}}};

  constexpr static amrex::Real WeightQ[4] = // GradPhiQ[phi_n][d/dx_i][quadrature point]
    {0.25, 0.25, 0.25, 0.25};


};
}
}
#endif
