#ifndef OPERATOR_ELASTIC_ISOTROPIC_H_
#define OPERATOR_ELASTIC_ISOTROPIC_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

#include "Operator/Elastic/Elastic.H"

using namespace amrex;

namespace Operator
{
namespace Elastic
{
class Isotropic : public Elastic
{
public:

  /// \brief Base constructor
  Isotropic ();
  /// \brief Typical constructor
  Isotropic (const Vector<Geometry>& a_geom,
		       const Vector<BoxArray>& a_grids,
		       const Vector<DistributionMapping>& a_dmap,
		       const LPInfo& a_info);
  /// \brief Destructor (does nothing)
  virtual ~Isotropic () {};

  void SetEta(amrex::Array<std::unique_ptr<amrex::MultiFab> > &eta,
	      GeneralAMRIntegratorPhysBC &eta_bc);

protected:

  amrex::Real C(const int i, const int j, const int k, const int l, const amrex::IntVect loc,
		const int amrlev, const int mglev, const MFIter &mfi) const;

  amrex::Real mu1, mu2;     ///< First Lame parameter
  amrex::Real lambda1, lambda2; ///< Second Lame parameter



};
}
}
#endif
