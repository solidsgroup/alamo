#ifndef OPERATOR_ELASTIC_POLYCRYSTAL_ISOTROPIC_H_
#define OPERATOR_ELASTIC_POLYCRYSTAL_ISOTROPIC_H_

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
class Isotropic : public PolyCrystal
{
public:

  /// \brief Base constructor
  Isotropic ();
  /// \brief Destructor (does nothing)
  virtual ~Isotropic () {};

protected:

  amrex::Vector<amrex::Real> C(const int i, const int j, const int k, const int l) const;

  amrex::Real mu1, mu2;     ///< First Lame parameter
  amrex::Real lambda1, lambda2; ///< Second Lame parameter



};
}
}
}
#endif
