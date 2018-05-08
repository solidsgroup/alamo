#ifndef OPERATOR_ELASTIC_POLYCRYSTAL_H_
#define OPERATOR_ELASTIC_POLYCRYSTAL_H_

#include <AMReX_MLCellLinOp.H>
#include <AMReX_Array.H>
#include <limits>

#include "BC/BC.H"
#include "Operator/Elastic/Elastic.H"

using namespace amrex;

namespace Operator
{
namespace Elastic
{
namespace PolyCrystal
{

class PolyCrystalModel
{
  friend class PolyCrystal;
  virtual amrex::Real C(const int i, const int j, const int k, const int l) const = 0;
};

class PolyCrystal : public Elastic
{
public:

  PolyCrystal () {};
  virtual ~PolyCrystal () {};
  void SetEta(amrex::Vector<std::unique_ptr<amrex::MultiFab> > &eta, BC::BC &eta_bc, std::vector<PolyCrystalModel *> &models);
  void Energies (FArrayBox& energyfab,
		 const FArrayBox& ufab,
		 int amrlev, const MFIter& mfi);

protected:

  virtual amrex::Real C(const int i, const int j, const int k, const int l, const amrex::IntVect loc,
			const int amrlev, const int mglev, const MFIter &mfi) const;

  int num_eta = 0;

  std::vector<std::array<std::array<std::array<std::array<amrex::Real,AMREX_SPACEDIM>,AMREX_SPACEDIM>,AMREX_SPACEDIM>,AMREX_SPACEDIM> > Cs;
};
}
}
}
#endif
