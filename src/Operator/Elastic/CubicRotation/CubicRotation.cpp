//add proper Doxygen comments
#include <AMReX_MultiFabUtil.H>
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_ArrayLim.H>

#include "eigen3/Eigen/Core" //not sure if this import is necessary

#include "Operator/Elastic/CubicRotation/CubicRotation.H"

Operator::Elastic::CubicRotation::CubicRotation(Eigen::Matrix<amrex::Real, 3, 3> R,
         amrex::Real C11in, amrex::Real C12in, amrex::Real C44in)
{
  C11 = C11in;
  C12 = C12in;
  C44 = C44in;
  
  amrex::Real C[3][3][3][3];  //may create naming conflict w/ C()
  
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      for(int k = 0; k < 3; k++) {
        for(int l = 0; l < 3; l++) {
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
  
  
  for(int p = 0; p < 3; p++) {
    for(int q = 0; q < 3; q++) {
      for(int s = 0; s < 3; s++) {
        for(int t = 0; t < 3; t++) {
	  Ctmp[p][q][s][t] = 0.0;
	  for(int i = 0; i < 3; i++) {
	    for(int j = 0; j < 3; j++) {
	      for(int k = 0; k < 3; k++) {
		for(int l = 0; l < 3; l++) {
		  Ctmp[p][q][s][t] += R(p,i)*R(s,k)*C[i][j][k][l]*R(q,j)*R(t,l);
		}
	      }
	    }
	  }
	  
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
