#include <AMReX_MultiFabUtil.H>
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_ArrayLim.H>

#include "Isotropic.H"


Operator::Elastic::PolyCrystal::Isotropic::Isotropic()
{
  lambda1 = 1000.0;
  lambda2 = 1000.0;
  mu1 = 1000.0;
  mu2 = 2000.0;
}

amrex::Vector<amrex::Real>
Operator::Elastic::PolyCrystal::Isotropic::C(const int i, const int j, const int k, const int l) const
{
  amrex::Vector<amrex::Real> Cs;

  Cs.push_back(0.0); Cs.push_back(0.0);
  if (i==k && j==l) {Cs[0] += mu1; Cs[1] += mu2;} 
  if (i==l && j==k) {Cs[0] += mu1; Cs[1] += mu2;}
  if (i==j && k==l) {Cs[0] += lambda1; Cs[1] += lambda2;}

  return Cs;
}
