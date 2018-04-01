#ifndef AMREX_MLPFSTIFFNESSMATRIX_H_
#define AMREX_MLPFSTIFFNESSMATRIX_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

using namespace amrex;

/// \class MLPFStiffnessMatrix
/// \brief Linear elasticity solver
///
/// This class will eventually be a pure abstract class that wraps the amrex MLCellLinOp
/// in a wrapper that translates from the language of elasticity // to amrex language.
// /// Currently, it is implemented as an instantiable class that implements 2D, isotropic
// /// elasticity.
// ///
// /// \bug Currently, Neumann boundary conditions **do not work!**
// ///
/// \bug This is **not** yet generalized to 3D
///
class MLPFStiffnessMatrix
  : public MLCellLinOp
{
public:

  /// \brief Base constructor
  //MLPFStiffnessMatrix () {}
  /// \brief Typical constructor
  MLPFStiffnessMatrix (const Vector<Geometry>& a_geom,
		       const Vector<BoxArray>& a_grids,
		       const Vector<DistributionMapping>& a_dmap,
		       const LPInfo& a_info);
  /// \brief Destructor (does nothing)
  virtual ~MLPFStiffnessMatrix ();

  MLPFStiffnessMatrix (const MLPFStiffnessMatrix&) = delete;
  MLPFStiffnessMatrix (MLPFStiffnessMatrix&&) = delete;
  MLPFStiffnessMatrix& operator= (const MLPFStiffnessMatrix&) = delete;
  MLPFStiffnessMatrix& operator= (MLPFStiffnessMatrix&&) = delete;

  /// \brief Relay to call define for parent class
  void define (const Vector<Geometry>& a_geom,
	       const Vector<BoxArray>& a_grids,
	       const Vector<DistributionMapping>& a_dmap,
	       const LPInfo& a_info = LPInfo());
  
  void setACoeffs (int amrlev, const MultiFab& alpha);

  void test(int amrlev, int mglev, MultiFab& out, const MultiFab& in)
  {Fapply(amrlev, mglev, out, in);}

protected:

  void averageDownCoeffs ();
  void averageDownCoeffsSameAmrLevel (Vector<MultiFab>& a);
  void applyMetricTermsCoeffs ();


  /// \brief Relay to call prepareForSolve for parent class
  virtual void prepareForSolve () final;
  /// \brief Required pure virtual function, but not used here
  virtual bool isSingular (int /*amrlev*/) const final { return false; }
  /// \brief Required pure virtual function, but not used here
  virtual bool isBottomSingular () const final { return false; }
  /// \brief Numerically implement \f$A\mathbf{x} = \mathbf{b}\f$
  virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final;
  /// \brief Compute one half Gauss-Seidel iteration corresponding to the operation described in MLPFStiffnessMatrix::Fapply
  virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh, int redblack) const final;
  /// \brief Compute the flux corresponding to the operator described in MLPFStiffnessMatrix::Fapply. Not used here.
  virtual void FFlux (int amrlev, const MFIter& mfi,
		      const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
		      const FArrayBox& sol, const int face_only=0) const final;
  /// \brief Tell how many fab components this operator is designed to work with
  virtual int getNComp() const;
  /// \brief Required pure virtual function, not used here
  virtual Real getAScalar () const final { return 0.0; }
  /// \brief Required pure virtual function, not used here
  virtual Real getBScalar () const final { return 0.0; }
  /// \brief Required pure virtual function, not used here
  virtual MultiFab const* getACoeffs (int /*amrlev*/, int /*mglev*/) const final { return nullptr;}
  /// \brief Required pure virtual function, not used here
  virtual std::array<MultiFab const*,AMREX_SPACEDIM> getBCoeffs (int /*amrlev*/, int /*mglev*/) const final {
    std::array<MultiFab const*,AMREX_SPACEDIM> ret;
    AMREX_D_TERM(ret[0] = nullptr;, ret[1] = nullptr;,ret[2] = nullptr;);
    return ret;}
  /// \brief Pointer cast utility
  virtual std::unique_ptr<MLLinOp> makeNLinOp (int /*grid_size*/) const final {
    amrex::Abort("MLABecLaplacian::makeNLinOp: Not implmented");
    return std::unique_ptr<MLLinOp>{};
  }
  
  virtual bool isCrossStencil () const { return false; }

private:

  amrex::Real C(const int i, const int j, const int k, const int l, amrex::IntVect loc, int amrlev, int mglev, amrex::MFIter &mfi) const;

  amrex::Real mu1, mu2;     ///< First Lame parameter
  amrex::Real lambda1, lambda2; ///< Second Lame parameter

  Vector<Vector<MultiFab> > m_a_coeffs;
};

#endif
