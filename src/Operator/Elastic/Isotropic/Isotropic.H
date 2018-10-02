#ifndef OPERATOR_ELASTIC_ISOTROPIC_H_
#define OPERATOR_ELASTIC_ISOTROPIC_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

#include "Model/Solid/Elastic/Isotropic/Isotropic.H"
#include "Operator/Elastic/Elastic.H"

using namespace amrex;

namespace Operator
{
namespace Elastic
{
class Isotropic : public Elastic<Model::Solid::Elastic::Isotropic::Isotropic>
{
public:

  Isotropic ();
  /// \brief Destructor (does nothing)
  virtual ~Isotropic () {};

protected:

  amrex::Real C(const int i, const int j, const int k, const int l, const amrex::IntVect loc,
		const int amrlev, const int mglev, const MFIter &mfi) const;

  amrex::Real mu;     ///< First Lame parameter
  amrex::Real lambda; ///< Second Lame parameter

};
}
}
#endif
