#ifndef AMREX_MLSTIFFNESSMATRIX_H_
#define AMREX_MLSTIFFNESSMATRIX_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

namespace amrex {

// (alpha * a - beta * (del dot b grad)) phi

class MLStiffnessMatrix
  : public MLCellLinOp
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

protected:

  virtual void prepareForSolve () final;
  virtual bool isSingular (int amrlev) const final { return m_is_singular[amrlev]; }
  virtual bool isBottomSingular () const final { return m_is_singular[0]; }
  virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
  virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh, int redblack) const final;
  virtual void FFlux (int amrlev, const MFIter& mfi,
		      const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
		      const FArrayBox& sol, const int face_only=0) const final;
  virtual int getNComp() const;
  virtual Real getAScalar () const final { return m_a_scalar; }
  //virtual Real Anorm (int amrlev, int mglev) const final;
  virtual Real getBScalar () const final { return m_b_scalar; }
  virtual MultiFab const* getACoeffs (int amrlev, int mglev) const final
  { return &(m_a_coeffs[amrlev][mglev]); }
  virtual std::array<MultiFab const*,AMREX_SPACEDIM> getBCoeffs (int amrlev, int mglev) const final
  { return amrex::GetArrOfConstPtrs(m_b_coeffs[amrlev][mglev]); }

  virtual std::unique_ptr<MLLinOp> makeNLinOp (int /*grid_size*/) const final {
    amrex::Abort("MLABecLaplacian::makeNLinOp: Not implmented");
    return std::unique_ptr<MLLinOp>{};
  }

private:

  Real m_a_scalar = std::numeric_limits<Real>::quiet_NaN();
  Real m_b_scalar = std::numeric_limits<Real>::quiet_NaN();
  Vector<Vector<MultiFab> > m_a_coeffs;
  Vector<Vector<std::array<MultiFab,AMREX_SPACEDIM> > > m_b_coeffs;

  //mutable Vector<Vector<Real> > m_Anorm;

  Vector<int> m_is_singular;

  //
  // functions
  //

  void averageDownCoeffsSameAmrLevel (Vector<MultiFab>& a,
				      Vector<std::array<MultiFab,AMREX_SPACEDIM> >& b);
  void averageDownCoeffs ();
  void averageDownCoeffsToCoarseAmrLevel (int flev);
  void applyMetricTermsCoeffs ();

  amrex::Real K(amrex::IntVect m, amrex::IntVect n, int i, int k);

  // Elastic constants for isotropic
  amrex::Real mu, lambda;
  amrex::Real dx1dx2=1.0, dx2dx1=1.0; 

};

}

#endif
