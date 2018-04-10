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
#include <AMReX_MLMG.H>

#include "Integrator/Integrator.H"

#include "BC/BC.H"
#include "IC/PerturbedInterface/PerturbedInterface.H"

#include "PFBoundary/PFBoundary.H"
#include "PFBoundary/PFBoundarySin.H"
#include "PFBoundary/PFBoundaryAbsSin.H"
#include "PFBoundary/PFBoundaryRead.H"

#include "Operator/Elastic/PolyCrystal/PolyCrystal.H"
#include "Operator/Elastic/Isotropic/Isotropic.H"
#include "Operator/Elastic/PolyCrystal/Isotropic/Isotropic.H"
#include "Operator/Elastic/Cubic/Cubic.H"

///
/// \class PhaseFieldMicrostructure
/// \brief Microstructure evolution with grain boundary anisotropy
///
/// Solve the Allen-Cahn evolution equation for microstructure with parameters \f$\eta_1\ldots\eta_n\f$,
/// where n corresponds to the number of grains.
///
class PhaseFieldMicrostructure : public Integrator::Integrator
{
public:
  PhaseFieldMicrostructure();

protected:

  /// \fn    Advance
  /// \brief Evolve phase field in time
  void Advance (int lev, Real time, Real dt);

  void Initialize (int lev);

  void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow);

  void TimeStepComplete(amrex::Real time, int iter);

private:

  int number_of_grains = 2;
  int number_of_ghost_cells = 2;

  amrex::Array<std::unique_ptr<amrex::MultiFab> > eta_new; ///< Multicomponent field variable storing \t$\eta_i\t$ for the __current__ timestep
  amrex::Array<std::unique_ptr<amrex::MultiFab> > eta_old; ///< Multicomponent field variable storing \t$\eta_i\t$ for the __previous__ timestep
  amrex::Array<std::unique_ptr<amrex::MultiFab> > displacement; ///< Multicomponent field variable storing pointwise displacements
  amrex::Array<std::unique_ptr<amrex::MultiFab> > strain; ///< Multicomponent field variable storing pointwise strains (6 components)
  amrex::Array<std::unique_ptr<amrex::MultiFab> > stress; ///< Multicomponent field variable storing pointwise strains (6 components)
  amrex::Array<std::unique_ptr<amrex::MultiFab> > body_force; ///< Multicomponent field variable storing pointwise displacements
  amrex::Array<std::unique_ptr<amrex::MultiFab> > energy; ///< Multicomponent field variable storing pointwise displacements

  BC::BC mybc;

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

  IC::IC *ic;

  bool        elastic_on = false;
  std::string elastic_type;
  int         elastic_max_iter = 200;
  int         elastic_max_fmg_iter = 0;
  int         elastic_verbose = 3;
  int         elastic_cgverbose = 3;
  amrex::Real elastic_tol_rel = 0.0;
  amrex::Real elastic_tol_abs = 1.0E-10;
  Operator::Elastic::PolyCrystal::PolyCrystal *elastic_operator;

};

#endif
