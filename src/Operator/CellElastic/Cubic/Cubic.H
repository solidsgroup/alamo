#ifndef OPERATOR_CELLELASTIC_CUBIC_H_
#define OPERATOR_CELLELASTIC_CUBIC_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

#include "Operator/CellElastic/CellElastic.H"

using namespace amrex;

namespace Operator
{
namespace CellElastic
{
class Cubic : public CellElastic
{
public:

  Cubic ();
  virtual ~Cubic () {};

protected:

  amrex::Real C(const int i, const int j, const int k, const int l, const amrex::IntVect loc,
		const int amrlev, const int mglev, const MFIter &mfi) const;

  amrex::Real E,nu,mu;
  amrex::Real C11, C12, C44;

};
}
}
#endif
