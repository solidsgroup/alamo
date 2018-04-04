#include <AMReX_MultiFabUtil.H>
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_ArrayLim.H>

#include "Operator/Elastic/Isotropic/Isotropic.H"


Operator::Elastic::Isotropic::Isotropic()
{
  lambda1 = 1.0;
  lambda2 = 1.0;
  mu1 = 1.0;
  mu2 = 2.0;
}

void
Operator::Elastic::Isotropic::SetEta(amrex::Array<std::unique_ptr<amrex::MultiFab> > &eta)
{
  RegisterNewFab(eta);
}

amrex::Real
Operator::Elastic::Isotropic::C(const int i, const int j, const int k, const int l,
				const amrex::IntVect loc,
				const int amrlev, const int mglev, const MFIter &mfi) const
{
  amrex::Real mu, lambda;

  const FArrayBox &etafab = GetFab(0,amrlev,mglev,mfi);

  mu = (mu1*etafab(loc,0) + mu2*etafab(loc,1)) / (etafab(loc,0) + etafab(loc,1));
  lambda = (lambda1*etafab(loc,0) + lambda2*etafab(loc,1)) / (etafab(loc,0) + etafab(loc,1));

  // TODO: This is a hack, needs to be fixed.

  if (mu != mu) mu = 0.5*(mu1+mu2);
  if (lambda != lambda) lambda = 0.5*(lambda1+lambda2);

  amrex::Real ret = 0.0;
  if (i==k && j==l) ret += mu;
  if (i==l && j==k) ret += mu;
  if (i==j && k==l) ret += lambda;

  //  std::cout << "(eta1="<<etafab(loc,0)<<",eta2="<<etafab(loc,1)<<")";
  return ret;
}
