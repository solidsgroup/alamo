#include <AMReX_MultiFabUtil.H>
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_ArrayLim.H>

#include "Isotropic.H"


OperatorCell::Elastic::Degradation::Isotropic::Isotropic(amrex::Real _lambda1, amrex::Real _mu1):
  lambda1(_lambda1), mu1(_mu1)
{

}

amrex::Real
OperatorCell::Elastic::Degradation::Isotropic::C(const int i, const int j, const int k, const int l) const
{
  amrex::Real C = 0.0;

  if (i==k && j==l) {C += mu1;}
  if (i==l && j==k) {C += mu1;}
  if (i==j && k==l) {C += lambda1;}

  return C;
}
