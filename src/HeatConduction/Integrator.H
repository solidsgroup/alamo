#ifndef HEATCONDUCTION_INTEGRATOR_H
#define HEATCONDUCTION_INTEGRATOR_H
#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"

#include "GeneralAMRIntegrator/GeneralAMRIntegrator.H"
#include "GeneralAMRIntegrator/GeneralAMRIntegratorBC.H"
#include "GeneralAMRIntegrator/GeneralAMRIntegratorIC.H"

#include "HeatConduction/InitialCondition/Cylinder.H"
#include "HeatConduction/InitialCondition/Constant.H"

#define TEMP_OLD(i,j,k) Temp_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k)))
#define TEMP(i,j,k) Temp_box(amrex::IntVect(AMREX_D_DECL(i,j,k)))


namespace HeatConduction
{
/// \class HeatConduction::Integrator
/// \brief Solve the heat diffusion equation in 2D or 3D.
///
/// This is a nontrivial example class that demonstrates how to write an explicit AMR
/// code using the #GeneralAMRIntegrator virtual class that abstracts the AmrBase class
/// from Amrex.
///
/// For more details:
///    - See documentation on #Initialize for input parameters
///    - See documentation on #Advance for equations and discretization
///    - See documentation on #TagCellsForRefinement for cell refinement criteria
/// For boundary conditions:
///    - See #GeneralAMRIntegratorBC
/// For initial conditions:
///    - See #Cylinder or #Constant
class Integrator : public GeneralAMRIntegrator
{
public:

  /// \brief Read in parameters and register field variables
  Integrator();

protected:

  /// \brief Use the #ic object to initialize #Temp
  void Initialize (int lev);

  /// \brief Integrate the heat equation
  void Advance (int lev, Real /*time*/, Real dt);

  /// \brief Tag cells for mesh refinement based on temperature gradient
  void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/);

private:
  int number_of_components = 1;					///< Number of components
  int number_of_ghost_cells = 1;				///< Number of ghost cells

  amrex::Array<std::unique_ptr<amrex::MultiFab> > Temp;		///< Temperature field variable (current timestep)
  amrex::Array<std::unique_ptr<amrex::MultiFab> > Temp_old;	///< Temperature field variable (previous timestep)

  amrex::Real alpha = 1.0;					///< Thermal diffusivity
  amrex::Real refinement_threshold = 0.01;			///< Criterion for cell refinement

  std::string ic_type;						///< String to determine what kind of initial condition to use
  GeneralAMRIntegratorIC *ic;					///< Pointer to abstract IC object
  GeneralAMRIntegratorPhysBC mybc;				///< Stock generic boundary condition object
};
}
#endif
