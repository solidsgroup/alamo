#include <AMReX_MultiFabUtil.H>
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_ArrayLim.H>

#include "eigen3/Eigen/Core" //not sure if this import is necessary

#include "Operator/Elastic/CubicRotation/CubicRotation.H"
//#include "Operator/Elastic/Cubic/Cubic.H"  //allows use of Cijkl internally


Operator::Elastic::CubicRotation::CubicRotation(Eigen::Matrix<amrex::Real, AMREX_SPACEDIM, AMREX_SPACEDIM> R,
						amrex::Real C11in, amrex::Real C12in, amrex::Real C44in)
{
  // Default values hard coded for now
  E1 = 1.0; nu1 = 0.25; mu1 = 2.0;
  E2 = 1.0; nu2 = 0.25; mu2 = 2.0;
  C11 = C11in;
  C12 = C12in;
  C44 = C44in;
}

void
Operator::Elastic::CubicRotation::SetEta(amrex::Array<std::unique_ptr<amrex::MultiFab> > &eta, GeneralAMRIntegratorPhysBC &eta_bc)
{
  RegisterNewFab(eta,eta_bc);
}

amrex::Real
Operator::Elastic::CubicRotation::C(const int i, const int j, const int k, const int l,
				    const amrex::IntVect loc,
				    const int amrlev, const int mglev, const MFIter &mfi) //const
{
  //begin "move to constructor once figured out"
  amrex::Real C11a = E1*(1-nu1)/(1-nu1-2.0*nu1*nu1);  //taken from Cubic
  amrex::Real C11b = E2*(1-nu2)/(1-nu2-2.0*nu2*nu2);
  amrex::Real C12a = E1*nu1/(1-nu1-2.0*nu1*nu1);
  amrex::Real C12b = E2*nu2/(1-nu2-2.0*nu2*nu2);

  const FArrayBox &etafab = GetFab(0,amrlev,mglev,mfi);

  C11 = (C11a*etafab(loc,0) + C11b*etafab(loc,1)) / (etafab(loc,0) + etafab(loc,1));
  C12 = (C12a*etafab(loc,0) + C12b*etafab(loc,1)) / (etafab(loc,0) + etafab(loc,1));
  C44 = (mu1*etafab(loc,0) + mu2*etafab(loc,1)) / (etafab(loc,0) + etafab(loc,1));

  // TODO: This is a hack, needs to be fixed.
  if (C11 != C11) C11 = 0.5*(C11a + C11b);  //presumably, this checks if NaN and fixes the problem
  if (C12 != C12) C12 = 0.5*(C12a + C12b);
  if (C44 != C44) C44 = 0.5*(mu1 + mu2);

  amrex::Real C[AMREX_SPACEDIM][AMREX_SPACEDIM][AMREX_SPACEDIM][AMREX_SPACEDIM];  //may create naming conflict w/ C()

  for(int i = 0; i < AMREX_SPACEDIM; i++) {
    for(int j = 0; j < AMREX_SPACEDIM; j++) {
      for(int k = 0; k < AMREX_SPACEDIM; k++) {
	for(int l = 0; l < AMREX_SPACEDIM; l++) {
	  if(i == j && j == k && k == l) {
	    C[i][j][k][l] = C11;
	  } else if (i==k && j==l) {
	    C[i][j][k][l] = C44;
	  } else if (i==j && k==l) {
	    C[i][j][k][l] = C12;
	  } else {
	    C[i][j][k][l] = 0.0;
	  }
	  for(int p = 0; p < AMREX_SPACEDIM; p++) {
	    for(int q = 0; q < AMREX_SPACEDIM; q++) {
	      for(int s = 0; s < AMREX_SPACEDIM; s++) {
		for(int t = 0; t < AMREX_SPACEDIM; t++) {
		  C[p][q][s][t] = R(p,i)*R(s,k)*C[i][j][k][l]*R(q,i)*R(t,l);
		}
	      }
	    }
	  }
	}
      }
    }
  }

  //end of "move to constructor once figured out"

  return C[i][j][k][l];
}
