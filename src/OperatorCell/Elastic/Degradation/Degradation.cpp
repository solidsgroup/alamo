#include <AMReX_MultiFabUtil.H>
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_ArrayLim.H>

#include "Degradation.H"

OperatorCell::Elastic::Degradation::Degradation::Degradation(const amrex::Real WInput)
{
	W = WInput;
}

void
OperatorCell::Elastic::Degradation::Degradation::SetEta(amrex::Vector<std::unique_ptr<amrex::MultiFab> > &eta,
						    std::vector<PristineMaterialModel *> &models)
{
	RegisterNewFab(eta);
	num_eta = eta[0].get()->nComp();
	std::cout << "number of eta = " << num_eta << std::endl;

	Cs.resize(1);
	//Cs.resize(num_phase);
	for (int i = 0; i < AMREX_SPACEDIM; i++)
		for (int j = 0; j < AMREX_SPACEDIM; j++)
			for (int k = 0; k < AMREX_SPACEDIM; k++)
				for (int l = 0; l < AMREX_SPACEDIM; l++)
					Cs[0][i][j][k][l] = models[0]->C(i,j,k,l);
					//for (int n = 0; n<num_phase; n++)
						//Cs[n][i][j][k][l] = models[n]->C(i,j,k,l);
}

amrex::Real
OperatorCell::Elastic::Degradation::Degradation::C(const int i, const int j, const int k, const int l, const amrex::IntVect loc,
					       const int amrlev, const int mglev, const MFIter &mfi) const
{

	amrex:: Real C = 0.0;

	for (int p = 0; p < AMREX_SPACEDIM; p ++)
		for (int q = 0; q < AMREX_SPACEDIM; q++)
			for (int r = 0; r < AMREX_SPACEDIM; r++)
				for (int s = 0; s < AMREX_SPACEDIM; s++)
					C += Cs[0][i][j][k][l]*M(i,j,p,q,loc,amrlev,mglev,mfi)*M(k,l,r,s,loc,amrlev,mglev,mfi);
					//for (int n = 0; n < num_phase; n++)
						//C += Cs[n][p][q][r][s]*M(i,j,p,q,loc,amrlev,mglev,mfi)*M(k,l,r,s,loc,amrlev,mglev,mfi);

	if (C != C)
	{
		std::cout << "Nans in C detected at ("<<i<< ","<<j<<","<<k<<","<<l<<")" << std::endl;
		//amrex::Real Cavg = 0.0;
		//for (int n = 0; n < etafab.nComp(); n++) Cavg += Cs[n][i][j][k][l];
		//Cavg /= (amrex::Real)num_eta;
		//return Cavg;
	}

	return C;
}

amrex::Real
OperatorCell::Elastic::Degradation::Degradation::M(const int i, const int j, const int k, const int l, const amrex::IntVect loc,
					       const int amrlev, const int mglev, const MFIter &mfi) const
{
	const FArrayBox &etafab = GetFab(0,amrlev,mglev,mfi);
	int m = i*(i==j ? 1:0) + (1- (i==j ? 1:0))*(9 - i - j);
	int n = k*(k==l ? 1:0) + (1- (k==l ? 1:0))*(9 - k - l);

	if (m != n) return 0;
	if (m == 1) return etafab(loc,0);
	if (m == 2) return etafab(loc,1);
	if (m == 3) return etafab(loc,2);
	if (m == 4) return sqrt(etafab(loc,1)*etafab(loc,2));
	if (m == 5) return sqrt(etafab(loc,0)*etafab(loc,2));
	if (m == 6) return sqrt(etafab(loc,0)*etafab(loc,1));
	return 0;
}

std::vector<amrex::Real>
OperatorCell::Elastic::Degradation::Degradation::dMdEta(const int i, const int j, const int k, const int l,
									const amrex::IntVect loc, const int amrlev, const int mglev, const MFIter &mfi) const
{
	const FArrayBox &etafab = GetFab(0,amrlev,mglev,mfi);
	int m = i*(i==j ? 1:0) + (1- (i==j ? 1:0))*(9 - i - j);
	int n = k*(k==l ? 1:0) + (1- (k==l ? 1:0))*(9 - k - l);

	std::vector<amrex::Real> dMdE = {0.0,0.0,0.0}; //assuming num_eta = 3 here.

	if (m != n) return dMdE;
	if (m == 1) return {1.0,0.0,0.0};
	if (m == 2) return {0.0,1.0,0.0};
	if (m == 3) return {0.0,0.0,1.0};
	if (m == 4)
	{
		dMdE[1] = 0.5*sqrt(etafab(loc,2))*pow(etafab(loc,1),-1.5);
		if(dMdE[1] != dMdE[1]) dMdE[1] = 0.0;
		dMdE[2] = 0.5*sqrt(etafab(loc,1))*pow(etafab(loc,2),-1.5);
		if(dMdE[2] != dMdE[2]) dMdE[2] = 0.0;
	}
	if (m == 5)
	{
		dMdE[0] = 0.5*sqrt(etafab(loc,2))*pow(etafab(loc,0),-1.5);
		if(dMdE[0] != dMdE[0]) dMdE[0] = 0.0;
		dMdE[2] = 0.5*sqrt(etafab(loc,0))*pow(etafab(loc,2),-1.5);
		if(dMdE[2] != dMdE[2]) dMdE[2] = 0.0;
	}
	if (m == 6)
	{
		dMdE[0] = 0.5*sqrt(etafab(loc,1))*pow(etafab(loc,0),-1.5);
		if(dMdE[0] != dMdE[0]) dMdE[0] = 0.0;
		dMdE[1] = 0.5*sqrt(etafab(loc,0))*pow(etafab(loc,1),-1.5);
		if(dMdE[1] != dMdE[1]) dMdE[1] = 0.0;
	}
	return dMdE;
}

amrex::Real
OperatorCell::Elastic::Degradation::Degradation::N(const int i, const int j, const amrex::IntVect loc,
					       const int amrlev, const int mglev, const MFIter &mfi) const
{
	const FArrayBox &etafab = GetFab(0,amrlev,mglev,mfi);
	int m = i*(i==j ? 1:0) + (1- (i==j ? 1:0))*(9 - i - j);

	if (m == 1) return (1-etafab(loc,1))*(1+etafab(loc,1));
	if (m == 2) return (1-etafab(loc,2))*(1+etafab(loc,2));
	if (m == 3) return (1-etafab(loc,3))*(1+etafab(loc,3));
	if (m == 4) return sqrt((1-etafab(loc,2))*(1+etafab(loc,2))*(1-etafab(loc,3))*(1+etafab(loc,3)));
	if (m == 5) return sqrt((1-etafab(loc,1))*(1+etafab(loc,1))*(1-etafab(loc,3))*(1+etafab(loc,3)));
	if (m == 6) return sqrt((1-etafab(loc,1))*(1+etafab(loc,1))*(1-etafab(loc,2))*(1+etafab(loc,2)));
	return 0;
}

std::vector<amrex::Real>
OperatorCell::Elastic::Degradation::Degradation::dNdEta(const int i, const int j, const amrex::IntVect loc,
								const int amrlev, const int mglev, const MFIter &mfi) const
{
	const FArrayBox &etafab = GetFab(0,amrlev,mglev,mfi);
	int m = i*(i==j ? 1:0) + (1- (i==j ? 1:0))*(9 - i - j);

	std::vector<amrex::Real> dNdE = {0.0,0.0,0.0}; //assuming num_eta = 3 here.
	if(m==1) dNdE[0] = -2.0*etafab(loc,0);
	if(m==2) dNdE[1] = -2.0*etafab(loc,1);
	if(m==3) dNdE[2] = -2.0*etafab(loc,2);
	if(m==4)
	{
		dNdE[1] = -etafab(loc,1)*sqrt((1-etafab(loc,2)*etafab(loc,2))/(1-etafab(loc,1)*etafab(loc,1)));
		dNdE[2] = -etafab(loc,2)*sqrt((1-etafab(loc,1)*etafab(loc,1))/(1-etafab(loc,2)*etafab(loc,2)));
	}
	if(m==5)
	{
		dNdE[0] = -etafab(loc,0)*sqrt((1-etafab(loc,2)*etafab(loc,2))/(1-etafab(loc,0)*etafab(loc,0)));
		dNdE[2] = -etafab(loc,2)*sqrt((1-etafab(loc,0)*etafab(loc,0))/(1-etafab(loc,2)*etafab(loc,2)));
	}
	if(m==6)
	{
		dNdE[0] = -etafab(loc,0)*sqrt((1-etafab(loc,1)*etafab(loc,1))/(1-etafab(loc,0)*etafab(loc,0)));
		dNdE[1] = -etafab(loc,1)*sqrt((1-etafab(loc,0)*etafab(loc,0))/(1-etafab(loc,1)*etafab(loc,1)));
	}
	return dNdE;
}

void
OperatorCell::Elastic::Degradation::Degradation::Energies (FArrayBox& energyfab,
						       const FArrayBox& ufab,
						       int amrlev, int mglev, const MFIter& mfi)
{
	const amrex::Real* DX = m_geom[amrlev][0].CellSize();

#if AMREX_SPACEDIM == 1
	amrex::IntVect dx(1);
#elif AMREX_SPACEDIM == 2
	amrex::IntVect dx(1,0);
	amrex::IntVect dy(0,1);
#elif AMREX_SPACEDIM > 2
	amrex::IntVect dx(1,0,0);
	amrex::IntVect dy(0,1,0);
	amrex::IntVect dz(0,0,1);
#endif

	const Box& bx = mfi.tilebox();
	for (int n = 0; n < num_eta; n++)
		for (int m1 = bx.loVect()[0]; m1<=bx.hiVect()[0]; m1++)
			for (int m2 = bx.loVect()[1]; m2<=bx.hiVect()[1]; m2++)
#if AMREX_SPACEDIM > 2
				for (int m3 = bx.loVect()[2]; m3<=bx.hiVect()[2]; m3++)
#endif
			{
#if AMREX_SPACEDIM == 1
				amrex::IntVect m(m1);
				amrex::Real gradu[1][1] = {ufab(m+dx,0) - ufab(m-dx,0)/(2.0*DX[0])};
#elif AMREX_SPACEDIM == 2
				amrex::IntVect m(m1,m2);
				amrex::Real gradu[2][2] = {{(ufab(m+dx,0) - ufab(m-dx,0))/(2.0*DX[0]),
					(ufab(m+dy,0) - ufab(m-dy,0))/(2.0*DX[1])},
					{(ufab(m+dx,1) - ufab(m-dx,1))/(2.0*DX[0]),
					(ufab(m+dy,1) - ufab(m-dy,1))/(2.0*DX[1])}};
#elif AMREX_SPACEDIM == 3
				amrex::IntVect m(m1,m2,m3);
				amrex::Real gradu[3][3] = {{(ufab(m+dx,0) - ufab(m-dx,0))/(2.0*DX[0]),
					(ufab(m+dy,0) - ufab(m-dy,0))/(2.0*DX[1]),
					(ufab(m+dz,0) - ufab(m-dz,0))/(2.0*DX[2])},
					{(ufab(m+dx,1) - ufab(m-dx,1))/(2.0*DX[0]),
					(ufab(m+dy,1) - ufab(m-dy,1))/(2.0*DX[1]),
					(ufab(m+dz,1) - ufab(m-dz,1))/(2.0*DX[2])},
					{(ufab(m+dx,2) - ufab(m-dx,2))/(2.0*DX[0]),
					(ufab(m+dy,2) - ufab(m-dy,2))/(2.0*DX[1]),
					(ufab(m+dz,2) - ufab(m-dz,2))/(2.0*DX[2])}};
#endif

				energyfab(m,n) = 0.0;

				for (int i=0; i<AMREX_SPACEDIM; i++)
					for (int j=0; j<AMREX_SPACEDIM; j++)
					{
						for (int k=0; k<AMREX_SPACEDIM; k++)
							for (int l=0; l<AMREX_SPACEDIM; l++)
								for (int p=0; p<AMREX_SPACEDIM; p++)
									for (int q=0; q<AMREX_SPACEDIM; q++)
										for (int r=0; r<AMREX_SPACEDIM; r++)
											for (int s=0; s<AMREX_SPACEDIM; s++)
											{
												std::vector<amrex::Real> dMdE1 = dMdEta(p,q,i,j,m,amrlev,mglev,mfi);
												std::vector<amrex::Real> dMdE2 = dMdEta(r,s,k,l,m,amrlev,mglev,mfi);
												energyfab(m,n) += gradu[i][j]*(dMdE1[n]*Cs[0][p][q][r][s]*M(r,s,k,l,m,amrlev,mglev,mfi) +
													M(p,q,i,j,m,amrlev,mglev,mfi)*Cs[0][p][q][r][s]*dMdE2[n])*gradu[k][l];
												energyfab(m,n) *= 0.5;
											}
						energyfab(m,n) += 2.0*W*dNdEta(i,j,m,amrlev,mglev,mfi)[n]*N(i,j,m,amrlev,mglev,mfi);
					}
			}

}
