#if BL_SPACEDIM == 2
#ifndef PHASEFIELDMICROSTRUCTURE_H

#define PHASEFIELDMICROSTRUCTURE_H
#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParallelDescriptor.H"

#include "GeneralAMRIntegrator.H"
#include "GeneralAMRIntegratorBC.H"

#include "PFBoundary/PFBoundary.H"
#include "PFBoundary/PFBoundarySin.H"
#include "PFBoundary/PFBoundaryAbsSin.H"
#include "PFBoundary/PFBoundaryRead.H"

class PhaseFieldMicrostructure : public GeneralAMRIntegrator
{
public:
  PhaseFieldMicrostructure();

protected:

  void Advance (int lev, Real time, Real dt);

  void Initialize (int lev);

  void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow);

private:

  int number_of_grains = 2;
  int number_of_ghost_cells = 2;

  amrex::Array<std::unique_ptr<amrex::MultiFab> > eta_new;
  amrex::Array<std::unique_ptr<amrex::MultiFab> > eta_old;

  GeneralAMRIntegratorPhysBC mybc;

  amrex::Real M, mu, gamma, sigma0, l_gb, beta;
  int anisotropy = 0;

  amrex::Real anisotropy_tstart;

  std::string ic_type, gb_type, filename;

  // Voronoi
  amrex::Array<amrex::Real> voronoi_x;
  amrex::Array<amrex::Real> voronoi_y;
#if BL_SPACEDIM==3
  amrex::Array<amrex::Real> voronoi_z;
#endif

  amrex::Real damp;

  PFBoundary *boundary;// = PFBoundarySin();

};

#endif
#endif
