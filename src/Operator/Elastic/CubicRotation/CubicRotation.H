#ifndef OPERATOR_ELASTIC_CUBIC_ROTATION_H_
#define OPERATOR_ELASTIC_CUBIC_ROTATION_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

#include "eigen3/Eigen/Core"

#include "Operator/Elastic/Elastic.H"

using namespace amrex;

namespace Operator
{
namespace Elastic
{
class CubicRotation : public Elastic
{
public:

  /// \brief Base constructor
  CubicRotation ();
  /// \brief Typical constructor
  CubicRotation (const Vector<Geometry>& a_geom,
		       const Vector<BoxArray>& a_grids,
		       const Vector<DistributionMapping>& a_dmap,
		       const LPInfo& a_info);
  /// \brief Destructor (does nothing)
  virtual ~CubicRotation () {};

  void SetEta(amrex::Array<std::unique_ptr<amrex::MultiFab> > &eta, GeneralAMRIntegratorPhysBC &eta_bc);

protected:

  amrex::Real C(const int i, const int j, const int k, const int l, const amrex::IntVect loc,
		const int amrlev, const int mglev, const MFIter &mfi) const;

  amrex::Real E1, E2;       ///< Young's modulus
  amrex::Real nu1, nu2;     ///< Poisson's ratio
  amrex::Real mu1, mu2;     ///< Shear modulus
  eigen::Matrix<amrex::Real, AMREX_SPACEDIM, AMREX_SPACEDIM> R;
  amrex::Real C11, C12, C44;
};
}
}
#endif
