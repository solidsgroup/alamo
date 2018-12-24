#include <AMReX_MultiFabUtil.H>
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_ArrayLim.H>

#include "eigen3/Eigen/Core" //not sure if this import is necessary

#include "OperatorCell/Elastic/CubicRotation/CubicRotation.H"
//#include "OperatorCell/Elastic/Cubic/Cubic.H"  //allows use of Cijkl internally
amrex::Real Cpqst[AMREX_SPACEDIM][AMREX_SPACEDIM][AMREX_SPACEDIM][AMREX_SPACEDIM];

OperatorCell::Elastic::CubicRotation::CubicRotation(Eigen::Matrix<amrex::Real, AMREX_SPACEDIM, AMREX_SPACEDIM> R,
						amrex::Real C11in, amrex::Real C12in, amrex::Real C44in)
{
  // Default values hard coded for now
  E1 = 1.0; nu1 = 0.25; mu1 = 2.0;
  E2 = 1.0; nu2 = 0.25; mu2 = 2.0;
  C11 = C11in;
  C12 = C12in;
  C44 = C44in;
  
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
        }
      }
    }
  }
  
  //amrex::Real Cpqst[AMREX_SPACEDIM][AMREX_SPACEDIM][AMREX_SPACEDIM][AMREX_SPACEDIM];
  amrex::Real SumStore;
  
  for(int p = 0; p < AMREX_SPACEDIM; p++) {
    for(int q = 0; q < AMREX_SPACEDIM; q++) {
      for(int s = 0; s < AMREX_SPACEDIM; s++) {
        for(int t = 0; t < AMREX_SPACEDIM; t++) {
	  SumStore = 0.0;
	  for(int i = 0; i < AMREX_SPACEDIM; i++) {
	    for(int j = 0; j < AMREX_SPACEDIM; j++) {
	      for(int k = 0; k < AMREX_SPACEDIM; k++) {
		for(int l = 0; l < AMREX_SPACEDIM; l++) {
		  SumStore += R(p,i)*R(s,k)*C[i][j][k][l]*R(q,i)*R(t,l);
		}
	      }
	    }
	  }
	  Cpqst[p][q][s][t] = SumStore;
	}
      }
    }
  }

}

// void
// OperatorCell::Elastic::CubicRotation::SetEta(amrex::Vector<std::unique_ptr<amrex::MultiFab> > &eta, BC::BC &eta_bc)
// {
//   RegisterNewFab(eta,eta_bc);
// }

amrex::Real
OperatorCell::Elastic::CubicRotation::C(const int i, const int j, const int k, const int l,
				    const amrex::IntVect /*loc*/,
				    const int /*amrlev*/, const int /*mglev*/, const MFIter &/*mfi*/) const
{
  return Cpqst[i][j][k][l];
}
