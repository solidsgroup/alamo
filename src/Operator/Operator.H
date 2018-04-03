#ifndef OPERATOR_H_
#define OPERATOR_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

using namespace amrex;

namespace Operator
{
class Operator
  : public amrex::MLCellLinOp
{
public:

  /// \brief Base constructor
  Operator () {};
  /// \brief Typical constructor
  Operator (const Vector<Geometry>& a_geom,
		       const Vector<BoxArray>& a_grids,
		       const Vector<DistributionMapping>& a_dmap,
		       const LPInfo& a_info);
  /// \brief Destructor (does nothing)
  virtual ~Operator ();

  Operator (const Operator&) = delete;
  Operator (Operator&&) = delete;
  Operator& operator= (const Operator&) = delete;
  Operator& operator= (Operator&&) = delete;

  /// \brief Relay to call define for parent class
  void define (const Vector<Geometry>& a_geom,
	       const Vector<BoxArray>& a_grids,
	       const Vector<DistributionMapping>& a_dmap,
	       const LPInfo& a_info = LPInfo());
  
protected:

  void RegisterNewFab(amrex::Array<std::unique_ptr<amrex::MultiFab> > &input);
  const amrex::FArrayBox & GetFab(const int num, const int amrlev, const int mglev, const amrex::MFIter &mfi) const;


  void averageDownCoeffs ();
  void averageDownCoeffsSameAmrLevel (Vector<MultiFab>& a);
  void applyMetricTermsCoeffs ();


  virtual void prepareForSolve () final;
  virtual bool isSingular (int /*amrlev*/) const final { return false; }
  virtual bool isBottomSingular () const final { return false; }
  virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const = 0;
  virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rsh, int redblack) const  = 0;
  virtual void FFlux (int amrlev, const MFIter& mfi,
		      const std::array<FArrayBox*,AMREX_SPACEDIM>& flux,
		      const FArrayBox& sol, const int face_only=0) const =0;
  virtual Real getAScalar () const final { return 0.0; }
  virtual Real getBScalar () const final { return 0.0; }
  virtual MultiFab const* getACoeffs (int /*amrlev*/, int /*mglev*/) const final { return nullptr;}
  virtual std::array<MultiFab const*,AMREX_SPACEDIM> getBCoeffs (int /*amrlev*/, int /*mglev*/) const final {
    std::array<MultiFab const*,AMREX_SPACEDIM> ret;
    AMREX_D_TERM(ret[0] = nullptr;, ret[1] = nullptr;,ret[2] = nullptr;);
    return ret;}
  /// \brief Pointer cast utility
  virtual std::unique_ptr<MLLinOp> makeNLinOp (int /*grid_size*/) const final {
    amrex::Abort("MLABecLaplacian::makeNLinOp: Not implmented");
    return std::unique_ptr<MLLinOp>{};
  }

private:
  int m_num_a_fabs = 0;
  Vector<Vector<Vector<MultiFab> > > m_a_coeffs; // m_a_coeffs[index of fab][amrlevel][mglevel]

  int current_amrlev = 0;
  int current_mglev = 0;
  amrex::MFIter *current_mfiter;
};
}
#endif
