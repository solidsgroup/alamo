#ifndef INTEGRATOR_INTEGRATOR_H
#define INTEGRATOR_INTEGRATOR_H

#include <chrono>
#include <ctime>
#include <string>
#include <limits>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_BC_TYPES.H>
#include <AMReX_AmrCore.H>
#include <AMReX_FluxRegister.H>
#include <AMReX_Utility.H>
#include <AMReX_PlotFileUtil.H>

#include "Set/Set.H"
#include "BC/BC.H"
#include "BC/Nothing.H"
#include "IO/WriteMetaData.H"

/// \brief Collection of numerical integrator objects
namespace Integrator
{

///
/// \class Integrator
/// \brief Pure abstract class for managing data structures, time integration (with substepping),
///        mesh refinement, and I/O. 
///
/// ## Input file parameters ##
/// ### Native input file parameters ###
///
///     max_step  = [maximum number of timesteps]
///     stop_time = [maximum simulation time]
///     timestep  = [time step for coarsest level]
///
///     amr.regrid_int = [number of timesteps between regridding]
///     amr.plot_int   = [number of timesteps between dumping output]
///     amr.plot_file  = [base name of output directory]
///     
///     amr.nsubsteps  = [number of temporal substeps at each level. This can be
///                       either a single int (which is then applied to every refinement
///                       level) or an array of ints (equal to amr.max_level) 
///                       corresponding to the refinement for each level.]
///
/// ### Inherited input file parameters (from amrex AmrMesh class) ###
///
///     amr.v                  = [verbosity level]
///     amr.max_level          = [maximum level of refinement]
///     amr.n_proper           = 
///     amr.grid_eff           = 
///     amr.n_error_buff       = 
///     amr.ref_ratio_vect     = [refinement ratios in each direction]
///     amr.ref_ratio          = [refinement ratio in all directions (cannot be used with ref_ratio_vect)]
///     amr.max_grid_x         = 
///     amr.blocking_factor    =
///     amr.n_cell             = [number of cells on coarsest level]
///     amr.refine_grid_layout = 
///     amr.check_input        = 
///
class Integrator
	: public amrex::AmrCore
{
public:

	/// \function Integrator
	/// \brief Constructor
	///
	/// Does the following things:
	///    - Read in simulation TIME(STEP) information
	///    - Read in simulation output and AMR information
	///    - Initalize timestep substep information
	///    - Create a clean directory
	/// For derived classes this **must** be called for the derived constructor. For instance: `code`
	/// ```cpp
	/// class MyDerivedClass : Integrator
	/// {
	///    MyDerivedClass() : Integrator() { ... }
	///    ...
	/// }
	/// ```
	Integrator ();

	~Integrator();

	/// \fn    FrontData
	/// \brief Front-end method to initialize simulation
	void InitData ();

	void Restart (std::string restartfile);

	/// \fn    Evolve
	/// \brief Front-end method to start simulation
	void Evolve ();

	void SetFilename(std::string _plot_file) {plot_file = _plot_file;};
	std::string GetFilename() {return plot_file;};

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

	/// \fn    TimeStepComplete
	/// \brief Run another system calculation (e.g. implicit solve) before integration step
	///
	/// This function is called at the beginning of every timestep. This function can be used
	/// to complete additional global solves, e.g. a MLMG implicit solve.
	///
	/// Overriding is optional; the default is to do nothing.
	///
	virtual void TimeStepBegin(amrex::Real /*time*/, int /*iter*/) {};

	/// \fn    TimeStepComplete
	/// \brief Run another system calculation (e.g. implicit solve) after integration step
	///
	/// This function is called at the end of every timestep. This function can be used
	/// to complete additional global solves, e.g. a MLMG implicit solve.
	///
	/// Overriding is optional; the default is to do nothing.
	///
	virtual void TimeStepComplete(amrex::Real /*time*/, int /*iter*/) {};

	/// \fn    Integrate
	/// \brief Perform an integration to compute integrated quantities
	///
	/// This is a function that is called by `Integrator` to update the variables registered in
	/// RegisterIntegratedVariable.
	/// The following variables are used:
	///   -  amrlev: current amr level
	///   -  time: current simulation time
	///   -  iter: current simulation iteration
	///   -  mfi:  current MFIter object (used to get FArrayBox from MultiFab)
	///   -  box:  Use this box (not mfi.tilebox). This box covers only cells on this level that are
	///            not also on a finer level.
	virtual void Integrate(int /*amrlev*/, Set::Scalar /*time*/, int /*iter*/,
			       const amrex::MFIter &/*mfi*/, const amrex::Box &/*box*/)
	{
		if (thermo.number > 0)
			Util::Warning(INFO,"integrated variables registered, but no integration implemented!"); 
	}
		
	virtual void Regrid(int /* amrlev */, Set::Scalar /* time */)
	{}


	/// \fn    RegisterNewFab
	/// \brief Register a field variable for AMR with this class 
	///
	/// You **must** call this function, typically in the constructor, for every field variable 
	/// that you use in any inherited class.
	/// 
	void RegisterNewFab (amrex::Vector<std::unique_ptr<amrex::MultiFab> > &new_fab,
			     ///<[in] Reference to field variable that you want to use.
			     ///      **Note** `new_fab` must not fall out of scope! Should be a
			     ///      private member variable.
			     BC::BC *new_bc,
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

	void RegisterNewFab (amrex::Vector<std::unique_ptr<amrex::MultiFab> > &new_fab,
			     int ncomp,
			     std::string name
			     );
	void RegisterNodalFab (amrex::Vector<std::unique_ptr<amrex::MultiFab> > &new_fab,
			       int ncomp,
			       int nghost,
			       std::string name
			       );

	void RegisterIntegratedVariable(Set::Scalar *integrated_variable, std::string name);

	void SetTimestep(Set::Scalar _timestep);
	void SetPlotInt(int plot_int);
	void SetThermoInt(int a_thermo_int) {thermo.interval = a_thermo_int;}
	void SetThermoPlotInt(int a_thermo_plot_int) {thermo.plot_int = a_thermo_plot_int;}

	amrex::Vector<amrex::Real> t_new;         ///< Keep track of current old simulation time on each level
	amrex::Vector<int> istep;           ///< Keep track of where each level is
	// PLOT FILES
	std::string plot_file {"plt"};   ///< Plotfile name

private:
	virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
					      const amrex::DistributionMapping& dm) override;
	virtual void MakeNewLevelFromCoarse  (int lev, amrex::Real time, const amrex::BoxArray& ba,
					      const amrex::DistributionMapping& dm) override;
	virtual void RemakeLevel             (int lev, amrex::Real time, const amrex::BoxArray& ba,
					      const amrex::DistributionMapping& dm) override;
	virtual void ClearLevel              (int lev) override;
	virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

	void IntegrateVariables(Set::Scalar cur_time, int step);

	void FillPatch (int lev, amrex::Real time,
			amrex::Vector<std::unique_ptr<amrex::MultiFab> > &source_mf,
			amrex::MultiFab &destination_multifab, BC::BC &physbc,
			int icomp);
	long CountCells (int lev);
	void TimeStep (int lev, amrex::Real time, int iteration);
	void FillCoarsePatch (int lev, amrex::Real time, amrex::Vector<std::unique_ptr<amrex::MultiFab> >& mf, BC::BC &physbc, int icomp, int ncomp);
	void GetData (const int lev, const amrex::Real time, amrex::Vector<amrex::MultiFab*>& data, amrex::Vector<amrex::Real>& datatime);

	std::vector<std::string> PlotFileName (int lev) const;
protected:
	void WritePlotFile (bool initial = false) const;

	//
	// MEMBER VARIABLES
	//

	// TIME (STEP) KEEPING
protected:
	amrex::Real timestep;	///< Timestep for the base level of refinement
private:
	amrex::Vector<amrex::Real> dt;	///< Timesteps for each level of refinement
	amrex::Vector<int> nsubsteps;	///< how many substeps on each level?
  
  
	amrex::Vector<amrex::Real> t_old;///< Keep track of current old simulation time on each level
	int max_step = std::numeric_limits<int>::max();	///< Maximum allowable timestep
	amrex::Real tstart = 0;	///< Default start time (default: 0)
	amrex::Real stop_time = std::numeric_limits<amrex::Real>::max(); ///< Default stop time 

protected:
	struct {
		int number_of_fabs = 0;
		amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab> > *> fab_array;
		amrex::Vector<int> ncomp_array;
		amrex::Vector<int> nghost_array;
		amrex::Vector<std::string> name_array;
		amrex::Vector<BC::BC *> physbc_array;
	} node;

	struct {
		int number_of_fabs = 0;
		amrex::Vector<amrex::Vector<std::unique_ptr<amrex::MultiFab> > *> fab_array;
		amrex::Vector<int> ncomp_array;
		amrex::Vector<int> nghost_array;
		amrex::Vector<std::string> name_array;
		amrex::Vector<BC::BC *> physbc_array;
	} cell;

	BC::Nothing bcnothing;

	// KEEP TRACK OF ALL INTEGRATED VARIABLES
	struct{
		int interval = -1;
		Set::Scalar dt = -1.0;
		int plot_int = -1;
		Set::Scalar plot_dt = -1.0;
		int number = 0;
		std::vector<Set::Scalar *> vars;
		std::vector<std::string> names;
	} thermo;

	// REGRIDDING
	int regrid_int = 2; ///< Determine how often to regrid (default: 2)

	std::string restart_file = "";

protected:
	int plot_int = -1;               ///< How frequently to dump plot file (default: never)
	Set::Scalar plot_dt = -1.0;
};
}
#endif
