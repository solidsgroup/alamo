#ifndef AMREX_ML_ABECLAPLACIAN_H_
#define AMREX_ML_ABECLAPLACIAN_H_

#include <AMReX_MLLinOp.H>
#include <AMReX_Array.H>
#include <limits>

namespace amrex {

// (alpha * a - beta * (del dot b grad)) phi

class MLStiffnessMatrix
  : public MLLinOp
{
public:

  MLStiffnessMatrix () {}
  MLStiffnessMatrix (const Vector<Geometry>& a_geom,
                     const Vector<BoxArray>& a_grids,
                     const Vector<DistributionMapping>& a_dmap,
                     const LPInfo& a_info = LPInfo());
  virtual ~MLStiffnessMatrix ();

  MLStiffnessMatrix (const MLStiffnessMatrix&) = delete;
  MLStiffnessMatrix (MLStiffnessMatrix&&) = delete;
  MLStiffnessMatrix& operator= (const MLStiffnessMatrix&) = delete;
  MLStiffnessMatrix& operator= (MLStiffnessMatrix&&) = delete;

  void define (const Vector<Geometry>& a_geom,
	       const Vector<BoxArray>& a_grids,
	       const Vector<DistributionMapping>& a_dmap,
	       const LPInfo& a_info = LPInfo());

  void setScalars (Real a, Real b);
  void setACoeffs (int amrlev, const MultiFab& alpha);
  void setBCoeffs (int amrlev, const std::array<MultiFab const*,AMREX_SPACEDIM>& beta);

  void setLevelBC (int amrlev, const MultiFab* a_levelbcdata);

protected:

  virtual void prepareForSolve () final;
  virtual bool isSingular (int amrlev) const final { return m_is_singular[amrlev]; }
  virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
  virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh, int redblack) const final;
  virtual void FFlux (int amrlev, const MFIter& mfi,
		      const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
		      const FArrayBox& sol, const int face_only=0) const final;
  virtual Real Anorm (int amrlev, int mglev) const final;
  virtual Real getBScalar () const final { return m_b_scalar; }

private:

  Real m_a_scalar = std::numeric_limits<Real>::quiet_NaN();
  Real m_b_scalar = std::numeric_limits<Real>::quiet_NaN();
  Vector<Vector<MultiFab> > m_a_coeffs;
  Vector<Vector<std::array<MultiFab,AMREX_SPACEDIM> > > m_b_coeffs;

  mutable Vector<Vector<Real> > m_Anorm;

  Vector<int> m_is_singular;

  //
  // functions
  //

  void averageDownCoeffsSameAmrLevel (Vector<MultiFab>& a,
				      Vector<std::array<MultiFab,AMREX_SPACEDIM> >& b);
  void averageDownCoeffs ();
  void averageDownCoeffsToCoarseAmrLevel (int flev);
  void applyMetricTermsCoeffs ();

  amrex::Real K(const amrex::IntVect m, const amrex::IntVect n, const int i, const int k) const;

  // Elastic constants for isotropic
  amrex::Real mu, lambda;
  amrex::Real dx1dx2=1.0, dx2dx1=1.0; 

};

}

#endif
