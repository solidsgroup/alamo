#include <AMReX_MultiFabUtil.H>
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_ArrayLim.H>

#include "Set/Set.H"
#include "PolyCrystal.H"


void
Operator::Elastic::PolyCrystal::PolyCrystal::SetEta(amrex::Vector<std::unique_ptr<amrex::MultiFab> > &eta,
						    BC::BC &eta_bc,
						    std::vector<PolyCrystalModel *> &models)
{
  RegisterNewFab(eta,eta_bc);
  num_eta = eta[0]->nComp();
  Cs.resize(num_eta);
  for (int i = 0; i < AMREX_SPACEDIM; i++)
    for (int j = 0; j < AMREX_SPACEDIM; j++)
      for (int k = 0; k < AMREX_SPACEDIM; k++)
	for (int l = 0; l < AMREX_SPACEDIM; l++)
	  for (int n = 0; n<num_eta; n++)
	    Cs[n][i][j][k][l] = models[n]->C(i,j,k,l);
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


void
Operator::Elastic::PolyCrystal::PolyCrystal::Energies (FArrayBox& energyfab,
						       const FArrayBox& ufab,
													   int amrlev, const MFIter& mfi)
{
	const amrex::Real* DX = m_geom[amrlev][0].CellSize();
  
#if AMREX_SPACEDIM == 2
	amrex::IntVect dx(1,0), dy(0,1);
#elif AMREX_SPACEDIM == 3
	amrex::IntVect dx(1,0,0), dy(0,1,0), dz(0,0,1);
#endif 

	const Box& bx = mfi.tilebox();
	for (int n = 0; n < num_eta; n++)
		for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++)
			for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++)
#if AMREX_SPACEDIM == 3
				for (int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++)
#endif
					{
						amrex::IntVect m(AMREX_D_DECL(m1,m2,m3));

						Set::Matrix gradu;
						gradu(0,0) = (ufab(m+dx,0) - ufab(m-dx,0))/(2.0*DX[0]);
						gradu(0,1) = (ufab(m+dy,0) - ufab(m-dy,0))/(2.0*DX[1]);
						gradu(1,0) = (ufab(m+dx,1) - ufab(m-dx,1))/(2.0*DX[0]);
						gradu(1,1) = (ufab(m+dy,1) - ufab(m-dy,1))/(2.0*DX[1]);
#if AMREX_SPACEDIM == 3
						gradu(0,2) = (ufab(m+dz,0) - ufab(m-dz,0))/(2.0*DX[2]);
						gradu(1,2) = (ufab(m+dz,1) - ufab(m-dz,1))/(2.0*DX[2]);
						gradu(2,0) = (ufab(m+dx,2) - ufab(m-dx,2))/(2.0*DX[0]);
						gradu(2,1) = (ufab(m+dy,2) - ufab(m-dy,2))/(2.0*DX[1]);
						gradu(2,2) = (ufab(m+dz,2) - ufab(m-dz,2))/(2.0*DX[2]);
#endif

						energyfab(m,n) = 0.0;

						for (int i=0; i<AMREX_SPACEDIM; i++)
							for (int j=0; j<AMREX_SPACEDIM; j++)
								for (int k=0; k<AMREX_SPACEDIM; k++)
									for (int l=0; l<AMREX_SPACEDIM; l++)
										energyfab(m,n) += gradu(i,j) * Cs[n][i][j][k][l] * gradu(k,l);
						energyfab(m,n) *= 0.5;
					}
}
