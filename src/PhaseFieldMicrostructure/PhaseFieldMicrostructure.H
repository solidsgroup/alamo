///
/// \file PhaseFieldMicrostructure.H
///
#ifndef PHASEFIELDMICROSTRUCTURE_H
#define PHASEFIELDMICROSTRUCTURE_H

#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "AMReX_ParmParse.H"
#include "AMReX_ParallelDescriptor.H"

#include "GeneralAMRIntegrator/GeneralAMRIntegrator.H"
#include "GeneralAMRIntegrator/GeneralAMRIntegratorBC.H"

#include "InitialConditions/PerturbedInterface.H"

#include "PFBoundary/PFBoundary.H"
#include "PFBoundary/PFBoundarySin.H"
#include "PFBoundary/PFBoundaryAbsSin.H"
#include "PFBoundary/PFBoundaryRead.H"

///
/// \class PhaseFieldMicrostructure
/// \brief Microstructure evolution with grain boundary anisotropy
///
/// Solve the Allen-Cahn evolution equation for microstructure with parameters \f$\eta_1\ldots\eta_n\f$,
/// where n corresponds to the number of grains.
///
class PhaseFieldMicrostructure : public GeneralAMRIntegrator
{
public:
  PhaseFieldMicrostructure();

protected:

  /// \fn    Advance
  /// \brief Evolve phase field in time
  void Advance (int lev, Real time, Real dt);

  void Initialize (int lev);

  void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow);

private:

  int number_of_grains = 2;
  int number_of_ghost_cells = 2;

  amrex::Array<std::unique_ptr<amrex::MultiFab> > eta_new; ///< Multicomponent field variable storing \t$\eta_i\t$ for the __current__ timestep
  amrex::Array<std::unique_ptr<amrex::MultiFab> > eta_old; ///< Multicomponent field variable storing \t$\eta_i\t$ for the __previous__ timestep
  //  amrex::Array<std::unique_ptr<amrex::MultiFab> > Curvature; ///< Multicomponent field variable storing \t$\eta_i\t$ for the __previous__ timestep
  // amrex::Array<std::unique_ptr<amrex::MultiFab> > Boundary_energy; ///< Multicomponent field variable storing \t$\eta_i\t$ for the __previous__ timestep

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

  PFBoundary *boundary;

  GeneralAMRIntegratorIC *ic;

};

#endif
