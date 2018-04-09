#include <AMReX_MultiFabUtil.H>
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_ArrayLim.H>

#include "Cubic.H"


Operator::Elastic::PolyCrystal::Cubic::Cubic()
{
  // Default values hard coded for now
  E1 = 1.0; nu1 = 0.25; mu1 = 2.0;
  E2 = 1.0; nu2 = 0.25; mu2 = 2.0;
}

amrex::Vector<amrex::Real>
Operator::Elastic::PolyCrystal::Cubic::C(const int i, const int j, const int k, const int l) const
{
  amrex::Real C11, C12, C44;
  amrex::Real C11a = E1*(1-nu1)/(1-nu1-2.0*nu1*nu1);
  amrex::Real C11b = E2*(1-nu2)/(1-nu2-2.0*nu2*nu2);
  amrex::Real C12a = E1*nu1/(1-nu1-2.0*nu1*nu1);
  amrex::Real C12b = E2*nu2/(1-nu2-2.0*nu2*nu2);
  amrex::Real C44a = mu1;
  amrex::Real C44b = mu2;

  amrex::Vector<amrex::Real> Cs;
  if(i == j && j == k && k == l) {Cs.push_back(C11a); Cs.push_back(C11b); }
  else if (i==k && j==l) {Cs.push_back(C44a); Cs.push_back(C44b);}
  else if (i==j && k==l) {Cs.push_back(C12a); Cs.push_back(C12b);} 

  return Cs;
}
