#include <AMReX_MultiFabUtil.H>
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_ArrayLim.H>

#include "Operator/Elastic/Cubic/Cubic.H"


Operator::Elastic::Cubic::Cubic()
{
  // Default values hard coded for now
  E1 = 1.0; nu1 = 0.25; mu1 = 2.0;
  E2 = 1.0; nu2 = 0.25; mu2 = 2.0;
}

void
Operator::Elastic::Cubic::SetEta(amrex::Array<std::unique_ptr<amrex::MultiFab> > &eta, GeneralAMRIntegratorPhysBC &eta_bc)
{
  RegisterNewFab(eta,eta_bc);
}

amrex::Real
Operator::Elastic::Cubic::C(const int i, const int j, const int k, const int l,
				const amrex::IntVect loc,
				const int amrlev, const int mglev, const MFIter &mfi) const
{
  amrex::Real C11, C12, C44;
  amrex::Real C11a = E1*(1-nu1)/(1-nu1-2.0*nu1*nu1);
  amrex::Real C11b = E2*(1-nu2)/(1-nu2-2.0*nu2*nu2);
  amrex::Real C12a = E1*nu1/(1-nu1-2.0*nu1*nu1);
  amrex::Real C12b = E2*nu2/(1-nu2-2.0*nu2*nu2);

  const FArrayBox &etafab = GetFab(0,amrlev,mglev,mfi);

  C11 = (C11a*etafab(loc,0) + C11b*etafab(loc,1)) / (etafab(loc,0) + etafab(loc,1));
  C12 = (C12a*etafab(loc,0) + C12b*etafab(loc,1)) / (etafab(loc,0) + etafab(loc,1));
  C44 = (mu1*etafab(loc,0) + mu2*etafab(loc,1)) / (etafab(loc,0) + etafab(loc,1));

  // TODO: This is a hack, needs to be fixed.

  if (C11 != C11) C11 = 0.5*(C11a + C11b);
  if (C12 != C12) C12 = 0.5*(C12a + C12b);
  if (C44 != C44) C44 = 0.5*(mu1 + mu2);

  if(i == j && j == k && k == l) return C11;
  if (i==k && j==l) return C44;
  if (i==j && k==l) return C12;

  return 0.0;
}
