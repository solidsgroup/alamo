#ifndef GENERALAMRINTEGRATOR_H
#define GENERALAMRINTEGRATOR_H

#include <string>
#include <limits>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_Utility.H>
#include <AMReX_PlotFileUtil.H>

#include "GeneralAMRIntegrator/GeneralAMRIntegratorBC.H"


///
/// \class GeneralAMRIntegrator
/// \brief Pure abstract class for managing data structures, time integration (with substepping),
///        mesh refinement, and I/O. 
///
class GeneralAMRIntegrator
  : public amrex::AmrCore
{
public:

  /// \function GeneralAMRIntegrator
  /// \brief Constructor
  ///
  /// Does the following things:
  ///    - Read in simulation TIME(STEP) information
  ///    - Read in simulation output and AMR information
  ///    - Initalize timestep substep information
  ///    - Create a clean directory
  /// For derived classes this **must** be called for the derived constructor. For instance: `code`
  /// ```cpp
  /// class MyDerivedClass : GeneralAMRIntegrator
  /// {
  ///    MyDerivedClass() : GeneralAMRIntegrator() { ... }
  ///    ...
  /// }
  /// ```
  GeneralAMRIntegrator ();

  virtual ~GeneralAMRIntegrator();

  /// \fn    FrontData
  /// \brief Front-end method to initialize simulation
  void InitData ();

  /// \fn    Evolve
  /// \brief Front-end method to start simulation
  void Evolve ();

protected:

  /// \fn    Initialize
  /// \brief Apply initial conditions
  ///
  /// You **must** override this function to inherit this class.
  /// This function is called before the simulation begins, and is where
  /// initial conditions should be applied.
  virtual void Initialize (int lev ///<[in] AMR Level
			   ) = 0;

  /// \fn    Advance
  /// \brief Perform computation
  ///
  /// You **must** override this function to inherit this class.
  /// Advance is called every time(sub)step, and implements the evolution of 
  /// the system in time.
  /// 
  virtual void Advance               (int lev,          ///<[in] AMR Level
				      amrex::Real time, ///< [in] System time
				      amrex::Real dt    ///< [in] Timestep for this level
				      )=0;

  /// \fn    TagCellsForRefinement
  /// \brief Tag cells where mesh refinement is needed
  ///
  /// You **must** override this function to inherit this class.
  /// Advance is called every time(sub)step, and implements the evolution of 
  /// the system in time.
  /// 
  virtual void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real time,
				      int ngrow) = 0;

  /// \fn    RegisterNewFab
  /// \brief Register a field variable for AMR with this class 
  ///
  /// You **must** call this function, typically in the constructor, for every field variable 
  /// that you use in any inherited class.
  /// 
  void RegisterNewFab (amrex::Array<std::unique_ptr<amrex::MultiFab> > &new_fab,
		       ///<[in] Reference to field variable that you want to use.
		       ///      **Note** `new_fab` must not fall out of scope! Should be a
		       ///      private member variable.
		       GeneralAMRIntegratorPhysBC &new_bc,
		       ///<[in] Boundary condition object. (Same BC object can be re-used 
		       ///      if applicable
		       int ncomp,
		       ///<[in] Number of components that the fab should have.
		       ///      (e.g. velocity in 2D would have two components.)
		       int nghost,
		       ///<[in] Number of ghost cells. Depends on the order of the
		       ///      maximum numerical derivative. (e.g. a Laplacian would require one.)
		       std::string name
		       ///<[in] The name of the field to be used when dumping output
		       );

private:
  virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
					const amrex::DistributionMapping& dm) override;
  virtual void MakeNewLevelFromCoarse  (int lev, amrex::Real time, const amrex::BoxArray& ba,
					const amrex::DistributionMapping& dm) override;
  virtual void RemakeLevel             (int lev, amrex::Real time, const amrex::BoxArray& ba,
			                const amrex::DistributionMapping& dm) override;
  virtual void ClearLevel              (int lev) override;
  virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

  void FillPatch (int lev, amrex::Real time,
		  amrex::Array<std::unique_ptr<amrex::MultiFab> > &source_mf,
		  amrex::MultiFab &destination_multifab, GeneralAMRIntegratorPhysBC &physbc,
		  int icomp);
  void InitFromCheckpoint ();
  long CountCells (int lev);
  void TimeStep (int lev, amrex::Real time, int iteration);
  void FillCoarsePatch (int lev, amrex::Real time, amrex::Vector<std::unique_ptr<amrex::MultiFab> >& mf, GeneralAMRIntegratorPhysBC &physbc, int icomp, int ncomp);
  void GetData (const int lev, const amrex::Real time, amrex::Array<amrex::MultiFab*>& data, amrex::Array<amrex::Real>& datatime);

  std::vector<std::string> PlotFileName (int lev) const;
  void WritePlotFile () const;
  void WriteMetaData () const;

  //
  // MEMBER VARIABLES
  //

  // TIME (STEP) KEEPING
  amrex::Real timestep;							///< Timestep for the base level of refinement
  amrex::Array<amrex::Real> dt;						///< Timesteps for each level of refinement
  amrex::Array<int> nsubsteps;						///< how many substeps on each level?
  amrex::Array<int> istep;						///< Keep track of where each level is
  amrex::Array<amrex::Real> t_new;					///< Keep track of current old simulation time on each level
  amrex::Array<amrex::Real> t_old;					///< Keep track of current old simulation time on each level
  int max_step = std::numeric_limits<int>::max();			///< Maximum allowable timestep
  amrex::Real tstart = 0;						///< Default start time (default: 0)
  amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();	///< Default stop time (default: max)

  // KEEPING TRACK OF ALL EVOLVING VARIABLES
  int number_of_fabs = 0;									///< Number of evolving variables (updated automatically)
  amrex::Array<amrex::Array<std::unique_ptr<amrex::MultiFab> > *> fab_array;	                ///< List of ptrs to all updated fabs
  amrex::Array<int> ncomp_array;								///< Array containing number of components for each fab array
  amrex::Array<int> nghost_array;								///< Array containing number of ghost cells in each fab array
  amrex::Array<std::string> name_array;								///< Array containing names for each variable in the fab array
  amrex::Array<std::unique_ptr<GeneralAMRIntegratorPhysBC> > physbc_array;			///< Array of pointers to BC objects

  // REGRIDDING
  int regrid_int = 2; ///< Determine how often to regrid (default: 2)

  // RESTART CHECKFILE -- NOT CURRENTLY USED!
  //std::string restart_chkfile;     ///< Checkfile name __Not currently used__
  //std::string check_file {"chk"};  ///< Checkfiles __Not currently used__
  //int check_int = -1;              ///< How frequently to dump checkfile __Not currently used__

  // PLOT FILES
  std::string plot_file {"plt"};   ///< Plotfile name
  int plot_int = -1;               ///< How frequently to dump plot file (default: never)
};

#endif
