#ifndef OPERATORCELL_ELASTIC_ISOTROPIC_H_
#define OPERATORCELL_ELASTIC_ISOTROPIC_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

#include "OperatorCell/Elastic/Elastic.H"

using namespace amrex;

namespace OperatorCell
{
namespace Elastic
{
class Isotropic : public Elastic
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
