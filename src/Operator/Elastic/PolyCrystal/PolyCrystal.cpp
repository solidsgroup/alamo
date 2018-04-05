#include <AMReX_MultiFabUtil.H>
#include <AMReX_REAL.H>
#include <AMReX_MLCGSolver.H>
#include <AMReX_ArrayLim.H>

#include "PolyCrystal.H"


void
Operator::Elastic::PolyCrystal::PolyCrystal::SetEta(amrex::Array<std::unique_ptr<amrex::MultiFab> > &eta, GeneralAMRIntegratorPhysBC &eta_bc)
{
  RegisterNewFab(eta,eta_bc);
}

