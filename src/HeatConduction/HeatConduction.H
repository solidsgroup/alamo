#if BL_SPACEDIM == 2
#ifndef HEATCONDUCTION_H
#define HEATCONDUCTION_H
#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"

#include "GeneralAMRIntegrator.H"
#include "GeneralAMRIntegratorBC.H"

class HeatConduction : public GeneralAMRIntegrator
{
public:
  HeatConduction();

protected:

  void Advance (int lev, Real /*time*/, Real dt);

  void Initialize (int lev);

  void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/);

private:
  int number_of_components = 1;
  int number_of_ghost_cells = 2;

  amrex::Array<std::unique_ptr<amrex::MultiFab> > Temp;
  amrex::Array<std::unique_ptr<amrex::MultiFab> > Temp_old;
  GeneralAMRIntegratorPhysBC mybc;
};

#endif
#endif
