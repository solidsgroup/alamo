#ifndef OPERATOR_ELASTIC_POLYCRYSTAL_H_
#define OPERATOR_ELASTIC_POLYCRYSTAL_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

#include "Operator/Elastic/Elastic.H"

using namespace amrex;

namespace Operator
{
namespace Elastic
{
namespace PolyCrystal
{
class PolyCrystal : public Elastic
{
public:

  /// \brief Base constructor
  PolyCrystal () {};
  /// \brief Destructor (does nothing)
  virtual ~PolyCrystal () {};

  void SetEta(amrex::Array<std::unique_ptr<amrex::MultiFab> > &eta, GeneralAMRIntegratorPhysBC &eta_bc);

protected:

  virtual amrex::Real C(const int i, const int j, const int k, const int l, const amrex::IntVect loc,
			const int amrlev, const int mglev, const MFIter &mfi) const = 0;

};
}
}
}
#endif
