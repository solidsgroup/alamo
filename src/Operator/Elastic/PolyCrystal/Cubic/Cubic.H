#ifndef OPERATOR_ELASTIC_POLYCRYSTAL_CUBIC_H_
#define OPERATOR_ELASTIC_POLYCRYSTAL_CUBIC_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

#include "Operator/Elastic/PolyCrystal/PolyCrystal.H"

using namespace amrex;

namespace Operator
{
namespace Elastic
{
namespace PolyCrystal
{
class Cubic : public PolyCrystal
{
public:

  /// \brief Base constructor
  Cubic ();

protected:

  amrex::Real C(const int i, const int j, const int k, const int l, const amrex::IntVect loc,
		const int amrlev, const int mglev, const MFIter &mfi) const;

  amrex::Real E1, E2;       ///< Young's modulus
  amrex::Real nu1, nu2;     ///< Poisson's ratio
  amrex::Real mu1, mu2;     ///< Shear modulus


};
}
}
}
#endif
