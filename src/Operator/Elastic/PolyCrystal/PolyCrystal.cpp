#include <AMReX_MultiFabUtil.H>
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_ArrayLim.H>

#include "PolyCrystal.H"


void
Operator::Elastic::PolyCrystal::PolyCrystal::SetEta(amrex::Array<std::unique_ptr<amrex::MultiFab> > &eta, GeneralAMRIntegratorPhysBC &eta_bc)
{
  RegisterNewFab(eta,eta_bc);
  num_eta = eta[0].get()->nComp();
  Cs.resize(num_eta);
  
  for (int i = 0; i < AMREX_SPACEDIM; i++)
    for (int j = 0; j < AMREX_SPACEDIM; j++)
      for (int k = 0; k < AMREX_SPACEDIM; k++)
	for (int l = 0; l < AMREX_SPACEDIM; l++)
	  for (int n = 0; n<num_eta; n++)
	    {
	      amrex::Vector<amrex::Real> C_tmp = C(i,j,k,l);
	      Cs[n][i][j][k][l] = C_tmp[n];
	    }
}

amrex::Real
Operator::Elastic::PolyCrystal::PolyCrystal::C(const int i, const int j, const int k, const int l, const amrex::IntVect loc,
					       const int amrlev, const int mglev, const MFIter &mfi) const
{
  const FArrayBox &etafab = GetFab(0,amrlev,mglev,mfi);

  amrex::Real etasum = 0.0, etaCsum = 0.0;

  for (int n = 0; n < num_eta; n++)
    {
      etaCsum += Cs[n][i][j][k][l] * etafab(loc,n);
      etasum += etafab(loc,n);
    }

  amrex::Real C = etaCsum / etasum;

  if (C != C)
    {
      amrex::Real Cavg = 0.0;
      for (int n = 0; n < etafab.nComp(); n++) Cavg += Cs[n][i][j][k][l];
      Cavg /= (amrex::Real)num_eta;
      return Cavg;
    }

  return C;
}
