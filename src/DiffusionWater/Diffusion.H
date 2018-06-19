#if BL_SPACEDIM == 2
#ifndef DIFFUSION_H
#define DIFFUSION_H
#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"

#include "GeneralAMRIntegrator.H"
#include "GeneralAMRIntegratorBC.H"

class Diffusion : public GeneralAMRIntegrator
{
public:
  Diffusion();

protected:

  void Advance (int lev, Real /*time*/, Real dt);

  void Initialize (int lev);

  void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/);

private:
  int number_of_components = 1;
  int number_of_ghost_cells = 2;

  amrex::Array<std::unique_ptr<amrex::MultiFab> > conc;
  amrex::Array<std::unique_ptr<amrex::MultiFab> > conc_old;
  GeneralAMRIntegratorPhysBC mybc;

  amrex::Real mu;	// diffusivity to be read from input file
};

#endif
#endif
