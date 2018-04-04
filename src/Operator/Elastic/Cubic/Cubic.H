#ifndef OPERATOR_ELASTIC_CUBIC_H_
#define OPERATOR_ELASTIC_CUBIC_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

#include "Operator/Elastic/Elastic.H"

using namespace amrex;

namespace Operator
{
namespace Elastic
{
class Cubic : public Elastic
{
public:

  /// \brief Base constructor
  Cubic ();
  /// \brief Typical constructor
  Cubic (const Vector<Geometry>& a_geom,
		       const Vector<BoxArray>& a_grids,
		       const Vector<DistributionMapping>& a_dmap,
		       const LPInfo& a_info);
  /// \brief Destructor (does nothing)
  virtual ~Cubic () {};

  void SetEta(amrex::Array<std::unique_ptr<amrex::MultiFab> > &eta);

protected:

  amrex::Real C(const int i, const int j, const int k, const int l, const amrex::IntVect loc,
		const int amrlev, const int mglev, const MFIter &mfi) const;

  amrex::Real E1, E2;       ///< Young's modulus
  amrex::Real nu1, nu2;     ///< Poisson's ratio
  amrex::Real mu1, mu2;     ///< Shear modulus


};
}
}
#endif
