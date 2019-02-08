#ifndef OPERATORCELL_ELASTIC_DEGRADATION_ISOTROPIC_H_
#define OPERATORCELL_ELASTIC_DEGRADATION_ISOTROPIC_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

#include "Operator/Elastic/Degradation/Degradation.H"

using namespace amrex;

namespace Operator
{
namespace Elastic
{
namespace Degradation
{
class Isotropic : public PristineMaterialModel
{
public:

  /// \brief Base constructor
  Isotropic (amrex::Real _lambda1, amrex::Real _mu1);
  /// \brief Destructor (does nothing)
  virtual ~Isotropic () {};

protected:

  amrex::Real C(const int i, const int j, const int k, const int l) const;

  amrex::Real lambda1 = 1.0, mu1 = 1.0;     ///< First Lame parameter



};
}
}
}
#endif
