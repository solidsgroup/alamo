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
enum ShapeFunctionType
  {
    Q4
  };

template<ShapeFunctionType sftype> struct Element{};


template <> struct Element<Q4>
{
public:
  Element<Q4>(amrex::Real _dx, amrex::Real _dy): dx(_dx), dy(_dy) {} ;

  template<int n>
  inline amrex::Real Phi(std::array<amrex::Real,2>);

  template<int n>
  inline std::array<amrex::Real,2> DPhi(std::array<amrex::Real,2>);

  template<int n>
  inline std::array<amrex::Real,2> QPoint();

  template<int n>
  inline amrex::Real QWeight();


private:
  amrex::Real dx, dy;
};

template <>
inline amrex::Real Element<Q4>::Phi<1>(std::array<amrex::Real,2> x)
{return (x[0]-0.0)*(x[1]-0.0)/dx/dy;} 
 
template <> inline amrex::Real Element<Q4>::Phi<2>(std::array<amrex::Real,2> x)
{return (dx-x[0])*(x[1]-0.0)/dx/dy;} 

template <> inline amrex::Real Element<Q4>::Phi<3>(std::array<amrex::Real,2> x)
{return (x[0]-0.0)*(dy-x[1])/dx/dy;} 

template <> inline amrex::Real Element<Q4>::Phi<4>(std::array<amrex::Real,2> x)
{return (dx-x[0])*(dy-x[1])/dx/dy;}

template <> inline std::array<amrex::Real,2> Element<Q4>::DPhi<1>(std::array<amrex::Real,2> x)
{return {x[1]/dx/dy, x[0]/dx/dy};} 

template <> inline std::array<amrex::Real,2> Element<Q4>::DPhi<2>(std::array<amrex::Real,2> x)
{return {-x[1]/dx/dy,(dx-x[0])/dx/dy};} 

template <> inline std::array<amrex::Real,2> Element<Q4>::DPhi<3>(std::array<amrex::Real,2> x)
{return {(dy-x[1])/dx/dy,-x[0]/dx/dy};} 

template <> inline std::array<amrex::Real,2> Element<Q4>::DPhi<4>(std::array<amrex::Real,2> x) 
{return {-(dy-x[1])/dx/dy,-(dx-x[0])/dx/dy};}

template <> amrex::Real Element<Q4>::QWeight<1>() {return 0.25*dx*dy; }
template <> amrex::Real Element<Q4>::QWeight<2>() {return 0.25*dx*dy; }
template <> amrex::Real Element<Q4>::QWeight<3>() {return 0.25*dx*dy; }
template <> amrex::Real Element<Q4>::QWeight<4>() {return 0.25*dx*dy; }


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
};
}
}
#endif
