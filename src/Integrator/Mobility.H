#ifndef INTEGRATOR_MOBILITY_H
#define INTEGRATOR_MOBILITY_H
#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_ParmParse.H"

#include "Integrator/Integrator.H"

#include "IC/IC.H"
#include "BC/BC.H"

#define TEMP_OLD(i,j,k) Temp_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k)))
#define TEMP(i,j,k) Temp_box(amrex::IntVect(AMREX_D_DECL(i,j,k)))


/// \class Mobility::Integrator
/// \brief Solve the heat diffusion equation in 2D or 3D.
///
/// This is a nontrivial example class that demonstrates how to write an explicit AMR
/// code using the #Integrator virtual class that abstracts the AmrBase class
/// from Amrex.
///
/// For more details:
///    - See documentation on #Initialize for input parameters
///    - See documentation on #Advance for equations and discretization
///    - See documentation on #TagCellsForRefinement for cell refinement criteria
/// For boundary conditions:
///    - See #BC
/// For initial conditions:
///    - See #Cylinder or #Constant
namespace Integrator
{
class Mobility : public Integrator
{
public:

  /// \brief Read in parameters and register field variables
  Mobility();
  ~Mobility();

protected:

  /// \brief Use the #ic object to initialize #Temp
  void Initialize (int lev) override;

  void TimeStepBegin(amrex::Real /*time*/, int /*iter*/) override;

  /// \brief Integrate the heat equation
  void Advance (int lev, amrex::Real /*time*/, amrex::Real dt) override;

  /// \brief Tag cells for mesh refinement based on temperature gradient
  void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/) override;

private:
  int number_of_components = 1;					///< Number of components
  int number_of_ghost_cells = 2;				///< Number of ghost cells

  amrex::Vector<std::unique_ptr<amrex::MultiFab> > disp;  	///< Displacement
  amrex::Vector<std::unique_ptr<amrex::MultiFab> > rhs;  	///< Right Hand Side
  amrex::Vector<std::unique_ptr<amrex::MultiFab> > res;  	///< Residual
  amrex::Vector<std::unique_ptr<amrex::MultiFab> > sigma;  	///< Right Hand Side
  amrex::Vector<std::unique_ptr<amrex::MultiFab> > energy_mf;  	///< Right Hand Side
  amrex::Vector<std::unique_ptr<amrex::MultiFab> > gammagb_mf;  	///< Right Hand Side
  amrex::Vector<std::unique_ptr<amrex::MultiFab> > gammagbold_mf;  	///< Right Hand Side

  struct {
    Set::Scalar ref_criterion = 0.01;    
  } amr;

  struct {
    int interval = 0;
    int fixed_iter = -1;
    int verbose = -1;
    int bottom_max_iter = -1;
    int max_coarsening_level = -1;
    std::string bottomsolver = "bicgstab";
  } solver;

  struct {
    Set::Scalar gamma = 1.0;
    Set::Scalar elastic_mult = 1.0;
    Set::Scalar kappa = 1.0;
    Set::Scalar L = 1.0;
    Set::Scalar gammagb0 = 0.1;
  } physics;

  struct {
    struct {
      Set::Scalar yhi_x = 0.0;
    } disp;
  } bc;

  amrex::Real alpha = 1.0;					///< Thermal diffusivity
  amrex::Real refinement_threshold = 0.01;			///< Criterion for cell refinement

  IC::IC *ic;            					///< Pointer to abstract IC object
  BC::BC *mybc;                   				///< Stock generic boundary condition object
};
}
#endif
