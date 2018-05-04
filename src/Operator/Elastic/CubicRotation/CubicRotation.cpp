///
/// \file CubicRotatation.cpp
/// \brief Apply a rotation to the elastic modulus tensor \f$\mathbb{C}_{ijkl}\f$
///

#include <AMReX_MultiFabUtil.H>
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_ArrayLim.H>

#include "eigen3/Eigen/Core" //not sure if this import is necessary

#include "Operator/Elastic/CubicRotation/CubicRotation.H"

amrex::Real C11, C12, C44;
Eigen::Matrix<amrex::Real, 3, 3> R;
amrex::Real Cc[3][3][3][3];

//beginning of Rotation Matrix constructor
Operator::Elastic::CubicRotation::CubicRotation(Eigen::Matrix<amrex::Real, 3, 3> Rin,
         amrex::Real C11, amrex::Real C12, amrex::Real C44)
{
  R = Rin;
  
  amrex::Real Ctmp[3][3][3][3];
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      for(int k = 0; k < 3; k++) {
        for(int l = 0; l < 3; l++) {
          if(i == j && j == k && k == l) {
            Ctmp[i][j][k][l] = C11;
          } else if (i==k && j==l) {
            Ctmp[i][j][k][l] = C44;
          } else if (i==j && k==l) {
            Ctmp[i][j][k][l] = C12;
          } else {
            Ctmp[i][j][k][l] = 0.0;
          }
        }
      }
    }
  }
  
  for(int p = 0; p < 3; p++) {
    for(int q = 0; q < 3; q++) {
      for(int s = 0; s < 3; s++) {
        for(int t = 0; t < 3; t++) {
	  Cc[p][q][s][t] = 0.0;
	  for(int i = 0; i < 3; i++) {
	    for(int j = 0; j < 3; j++) {
	      for(int k = 0; k < 3; k++) {
		for(int l = 0; l < 3; l++) {
		  Cc[p][q][s][t] += R(p,i)*R(s,k)*Ctmp[i][j][k][l]*R(q,j)*R(t,l);
		}
	      }
	    }
	  }
	  
	}
      }
    }
  }

}

//beginning of Bunge Euler Angle constructor (ZXZ, radians)
Operator::Elastic::CubicRotation::CubicRotation(amrex::Real phi1, amrex::Real theta, amrex::Real phi2,
		   amrex::Real C11, amrex::Real C12, amrex::Real C44)
{
  R << std::cos(phi1)*std::cos(phi2) - std::cos(theta)*std::sin(phi1)*std::sin(phi2),
       -std::cos(phi1)*std::sin(phi2) - std::cos(theta)*std::cos(phi2)*std::sin(phi1),
       std::sin(phi1)*std::sin(theta),
       std::cos(phi2)*std::sin(phi1) + std::cos(phi1)*std::cos(theta)*std::sin(phi2),
       std::cos(phi1)*std::cos(theta)*std::cos(phi2) - std::sin(phi1)*std::sin(phi2),
       -std::cos(phi1)*std::sin(theta),
       std::sin(theta)*std::sin(phi2),
       std::cos(phi2)*std::sin(theta),
       std::cos(theta);
	   
  amrex::Real Ctmp[3][3][3][3];
  for(int i = 0; i < 3; i++) {
    for(int j = 0; j < 3; j++) {
      for(int k = 0; k < 3; k++) {
        for(int l = 0; l < 3; l++) {
          if(i == j && j == k && k == l) {
            Ctmp[i][j][k][l] = C11;
          } else if (i==k && j==l) {
            Ctmp[i][j][k][l] = C44;
          } else if (i==j && k==l) {
            Ctmp[i][j][k][l] = C12;
          } else {
            Ctmp[i][j][k][l] = 0.0;
          }
        }
      }
    }
  }
  
  for(int p = 0; p < 3; p++) {
    for(int q = 0; q < 3; q++) {
      for(int s = 0; s < 3; s++) {
        for(int t = 0; t < 3; t++) {
	  Cc[p][q][s][t] = 0.0;
	  for(int i = 0; i < 3; i++) {
	    for(int j = 0; j < 3; j++) {
	      for(int k = 0; k < 3; k++) {
		for(int l = 0; l < 3; l++) {
		  Cc[p][q][s][t] += R(p,i)*R(s,k)*Ctmp[i][j][k][l]*R(q,j)*R(t,l);
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
  return Cc[i][j][k][l];
}
