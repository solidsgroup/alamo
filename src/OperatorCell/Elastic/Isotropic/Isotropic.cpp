#include <AMReX_MultiFabUtil.H>
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_ArrayLim.H>

#include "Isotropic.H"


OperatorCell::Elastic::Isotropic::Isotropic ()
{
  mu = 26.0;
  lambda = 60.0;
}

amrex::Real
OperatorCell::Elastic::Isotropic::C(const int i, const int j, const int k, const int l,
				const amrex::IntVect /*loc*/,
				const int /*amrlev*/, const int /*mglev*/, const MFIter & /*mfi*/) const
{
  amrex::Real ret = 0.0;
  if (i==k && j==l) ret += mu;
  if (i==l && j==k) ret += mu;
  if (i==j && k==l) ret += lambda;

  return ret;
}
