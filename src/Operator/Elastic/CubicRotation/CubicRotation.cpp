#include <AMReX_MultiFabUtil.H>
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_ArrayLim.H>
//#include <iostream>

#include "eigen3/Eigen/Core" //not sure if this import is necessary

#include "Operator/Elastic/CubicRotation/CubicRotation.H"

Operator::Elastic::CubicRotation::CubicRotation(Eigen::Matrix<amrex::Real, AMREX_SPACEDIM, AMREX_SPACEDIM> R)
{
  E = 1.0; nu = 0.25; mu = 2.0;

  C11 = E*(1-nu)/(1-nu-2.0*nu*nu);
  C12 = E*nu/(1-nu-2.0*nu*nu);
  C44 = mu;
  //std::cout << "this is C11: " << C11 << "  this is C12: " << C12 << "  this is C44: " << C44 << "\n";
  //std::cout << "this is R: " << R(0,0) << " " << R(0,1) << " " << R(1,0) << " " << R(1,1) << "\n";
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
  
  amrex::Real SumStore;
  
  for(int p = 0; p < AMREX_SPACEDIM; p++) {
    for(int q = 0; q < AMREX_SPACEDIM; q++) {
      for(int s = 0; s < AMREX_SPACEDIM; s++) {
        for(int t = 0; t < AMREX_SPACEDIM; t++) {
	  Ctmp[p][q][s][t] = 0.0;
	  for(int i = 0; i < AMREX_SPACEDIM; i++) {
	    for(int j = 0; j < AMREX_SPACEDIM; j++) {
	      for(int k = 0; k < AMREX_SPACEDIM; k++) {
		for(int l = 0; l < AMREX_SPACEDIM; l++) {
		  Ctmp[p][q][s][t] += R(p,i)*R(s,k)*C[i][j][k][l]*R(q,j)*R(t,l);
		  //std::cout << "this is Ctmp: " << Ctmp[p][q][s][t] << " with " << p << q << s << t;
		}
	      }
	    }
	  }
	  //std::cout << "this is Ctmp: " << Ctmp[p][q][s][t] << " with " << p << q << s << t << "\n";
	  //std::cout << "this is C: " << C[p][q][s][t] << " with " << p << q << s << t << "\n";
	}
      }
    }
  }

}


  

void
Operator::Elastic::CubicRotation::SetEta(amrex::Array<std::unique_ptr<amrex::MultiFab> > &eta, GeneralAMRIntegratorPhysBC &eta_bc)
{
  RegisterNewFab(eta,eta_bc);
}

amrex::Real
Operator::Elastic::CubicRotation::C(const int i, const int j, const int k, const int l,
				    const amrex::IntVect loc,
				    const int amrlev, const int mglev, const MFIter &mfi) const
{
  return Ctmp[i][j][k][l];
}
