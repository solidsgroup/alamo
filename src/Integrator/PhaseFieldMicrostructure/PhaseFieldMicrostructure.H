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
#include "IC/Voronoi/Voronoi.H"

#include "Model/Interface/GrainBoundary/GrainBoundary.H"
#include "Model/Interface/GrainBoundary/Sin.H"
#include "Model/Interface/GrainBoundary/AbsSin.H"
#include "Model/Interface/GrainBoundary/Read.H"

#include "Operator/Elastic/PolyCrystal/PolyCrystal.H"
#include "Operator/Elastic/Isotropic/Isotropic.H"
#include "Operator/Elastic/PolyCrystal/Isotropic/Isotropic.H"
#include "Operator/Elastic/PolyCrystal/Cubic/Cubic.H"
#include "Operator/Elastic/Cubic/Cubic.H"

namespace Integrator
{
///
/// \class PhaseFieldMicrostructure
/// \brief Microstructure evolution with grain boundary anisotropy
///
/// Solve the Allen-Cahn evolution equation for microstructure with parameters \f$\eta_1\ldots\eta_n\f$,
/// where n corresponds to the number of grains.
///
class PhaseFieldMicrostructure : public Integrator
{
public:
  PhaseFieldMicrostructure();

protected:

  /// \fn    Advance
  /// \brief Evolve phase field in time
  void Advance (int lev, Real time, Real dt);

  void Initialize (int lev);

  void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow);

  void TimeStepBegin(amrex::Real time, int iter);

private:

  int number_of_grains = 2;
  int number_of_ghost_cells = 2;

  amrex::Vector<std::unique_ptr<amrex::MultiFab> > eta_new; ///< Multicomponent field variable storing \t$\eta_i\t$ for the __current__ timestep
  amrex::Vector<std::unique_ptr<amrex::MultiFab> > eta_old; ///< Multicomponent field variable storing \t$\eta_i\t$ for the __previous__ timestep
  amrex::Vector<std::unique_ptr<amrex::MultiFab> > displacement; ///< Multicomponent field variable storing pointwise displacements
  amrex::Vector<std::unique_ptr<amrex::MultiFab> > strain; ///< Multicomponent field variable storing pointwise strains (6 components)
  amrex::Vector<std::unique_ptr<amrex::MultiFab> > stress; ///< Multicomponent field variable storing pointwise strains (6 components)
  amrex::Vector<std::unique_ptr<amrex::MultiFab> > stress_vm; ///< Multicomponent field variable storing pointwise strains (6 components)
  amrex::Vector<std::unique_ptr<amrex::MultiFab> > body_force; ///< Multicomponent field variable storing pointwise displacements
  amrex::Vector<std::unique_ptr<amrex::MultiFab> > energy; ///< Single component variable storing strain energy
  amrex::Vector<std::unique_ptr<amrex::MultiFab> > energies; ///< Multicomponent field variable storing individual strain energies

  BC::BC *mybc;

  amrex::Real M, mu, gamma, sigma0, l_gb, beta;
  int anisotropy = 0;

  amrex::Real anisotropy_tstart;

  std::string ic_type, gb_type, filename;

  amrex::Real damp;

  Model::Interface::GrainBoundary::GrainBoundary *boundary;

  IC::IC *ic;

  bool        elastic_on = false;
  int         elastic_int = 1;
  std::string elastic_type;
  int         elastic_max_iter = 200;
  int         elastic_max_fmg_iter = 0;
  int         elastic_verbose = 3;
  int         elastic_cgverbose = 3;
  amrex::Real elastic_tol_rel = 0.0;
  amrex::Real elastic_tol_abs = 1.0E-10;
  amrex::Real elastic_tstart = 0.0;
  amrex::Vector<amrex::Real> elastic_load_t;
  amrex::Vector<amrex::Real> elastic_load_disp;
  Operator::Elastic::PolyCrystal::PolyCrystal *elastic_operator;

};
}
#endif
