#include <AMReX_MultiFabUtil.H>
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_ArrayLim.H>

#include "Cubic.H"


Operator::Elastic::Cubic::Cubic()
{
  E = 1.0; nu = 0.25; mu = 2.0;

  C11 = E*(1-nu)/(1-nu-2.0*nu*nu);
  C12 = E*nu/(1-nu-2.0*nu*nu);
  C44 = mu;
}

amrex::Real
Operator::Elastic::Cubic::C(const int i, const int j, const int k, const int l,
			    const amrex::IntVect /*loc*/,
			    const int /*amrlev*/, const int /*mglev*/, const MFIter &/*mfi*/) const
{
  if(i == j && j == k && k == l) return C11;
  if (i==k && j==l) return C44;
  if (i==j && k==l) return C12;

  return 0.0;
}
