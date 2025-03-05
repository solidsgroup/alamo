//
// Pure abstract class for managing data structures, time integration (with substepping),
// mesh refinement, and I/O. 
//
// Native input file parameters:
//
// .. code-block:: make
//
//     max_step  = [maximum number of timesteps]
//     stop_time = [maximum simulation time]
//     timestep  = [time step for coarsest level]
//
//     amr.regrid_int = [number of timesteps between regridding]
//     amr.plot_int   = [number of timesteps between dumping output]
//     amr.plot_file  = [base name of output directory]
//     
//     amr.nsubsteps  = [number of temporal substeps at each level. This can be
//                       either a single int (which is then applied to every refinement
//                       level) or an array of ints (equal to amr.max_level) 
//                       corresponding to the refinement for each level.]
//
// Inherited input file parameters (from amrex AmrMesh class):
//
// .. code-block:: make
//
//     amr.v                  = [verbosity level]
//     amr.max_level          = [maximum level of refinement]
//     amr.n_proper           = 
//     amr.grid_eff           = 
//     amr.n_error_buff       = 
//     amr.ref_ratio_vect     = [refinement ratios in each direction]
//     amr.ref_ratio          = [refinement ratio in all directions (cannot be used with ref_ratio_vect)]
//     amr.max_grid_x         = 
//     amr.blocking_factor    =
//     amr.n_cell             = [number of cells on coarsest level]
//     amr.refine_grid_layout = 
//     amr.check_input        = 
//

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
#include "BaseField.H"

/// \brief Collection of numerical integrator objects
namespace Integrator
{

class Integrator
    : public amrex::AmrCore
{
public:

    /// This is the constructor for the intetgrator class, which reads timestep information,
    /// simulation output and AMR, initialized time substep, and creates a new directory.
    Integrator();

    /// Virtual destructure; make sure delete any pointers that you create here.
    virtual ~Integrator();

    /// Front-end method to initialize simulation on all levels
    void InitData();

    /// Read in output from previous simulation and start simulation at that point -
    /// Not currently tested
    void Restart(std::string restartfile, bool a_node = false);

    /// Front-end method to start simulation
    void Evolve();

    /// Simple setter to set filename
    void SetFilename(std::string _plot_file) { plot_file = _plot_file; };

    /// Simple getter to get filename
    std::string GetFilename() { return plot_file; };

    /// This overrides an AMReX method just to allow for explicit meshing when
    /// desired.
    void regrid(int lbase, Set::Scalar time, bool initial = false) override
    {
        if (!explicitmesh.on)
            AmrCore::regrid(lbase, time, initial);
    }

    /// This creates a new levels that have not previously been
    /// used.
    void InitFromScratch(Set::Scalar time)
    {
        if (!explicitmesh.on) AmrCore::InitFromScratch(time);
        else
        {
            // Generate the coarse level mesh.
            {
                finest_level = 0;
                const amrex::BoxArray& ba = MakeBaseGrids();
                amrex::DistributionMapping dm(ba);
                const auto old_num_setdm = num_setdm;
                const auto old_num_setba = num_setba;
                MakeNewLevelFromScratch(0, time, ba, dm);
                if (old_num_setba == num_setba) SetBoxArray(0, ba);
                if (old_num_setdm == num_setdm) SetDistributionMap(0, dm);
            }
            // Generate subsequent level meshes based on user input
            for (int ilev = 0; ilev < maxLevel(); ++ilev)
            {
                finest_level = ilev + 1;
                amrex::BoxArray grids(explicitmesh.box[ilev]);
                ChopGrids(ilev + 1, grids, amrex::ParallelDescriptor::NProcs());
                amrex::DistributionMapping dmap(grids);
                SetBoxArray(ilev + 1, grids);
                SetDistributionMap(ilev + 1, dmap);
                MakeNewLevelFromScratch(ilev + 1, time, grids, dmap);
            }

        }
    }

protected:

    /// You must override this function to inherit this class;
    /// this function is called before the simulation begins, and is where
    /// initial conditions should be applied.
    virtual void Initialize(int lev ///<[in] AMR Level
    ) = 0;

    /// You must override this function to inherit this class;
    /// Advance is called every time(sub)step, and implements the evolution of 
    /// the system in time.
    virtual void Advance(int lev,          ///<[in] AMR Level
        amrex::Real time, ///< [in] System time
        amrex::Real dt    ///< [in] Timestep for this level
    ) = 0;

    /// You must override this function to inherit this class;
    /// Advance is called every time(sub)step, and implements the evolution of 
    /// the system in time.
    virtual void TagCellsForRefinement(int lev, amrex::TagBoxArray& tags, amrex::Real time,
        int ngrow) = 0;

    /// This optional function is  called at the beginning of every timestep, and can be used
    /// to complete additional global solves, e.g. a MLMG implicit solve.
    virtual void TimeStepBegin(Set::Scalar /*time*/, int /*iter*/) {};

    /// This optional function is called at the end of every timestep, and can be used
    /// to complete additional global solves, e.g. a MLMG implicit solve.
    virtual void TimeStepComplete(Set::Scalar /*time*/, int /*iter*/) {};

    /// This is a function that is called by `Integrator` to update the variables registered in
    /// RegisterIntegratedVariable; you can override this to do your integration.
    virtual void Integrate( int /*amrlev*/,
                            Set::Scalar /*time*/,
                            int /*iter*/,
                            const amrex::MFIter&/*mfi*/,
                            const amrex::Box&/*box*/
        )
    {
        if (thermo.number > 0)
            Util::Warning(INFO, "integrated variables registered, but no integration implemented!");
    }

    /// An optionally overridable method to trigger behavior whenver a regrid occurs.
    virtual void Regrid(int /* amrlev */, Set::Scalar /* time */)
    {}


    /// Add a new cell-based scalar field.
    void RegisterNewFab(Set::Field<Set::Scalar>& new_fab, BC::BC<Set::Scalar>* new_bc, int ncomp, int nghost, std::string name, bool writeout,std::vector<std::string> suffix = {});
    /// Add a new cell-based scalar field (with additional arguments).
    void RegisterNewFab(Set::Field<Set::Scalar>& new_fab, int ncomp, std::string name, bool writeout,std::vector<std::string> suffix = {});
    /// Add a new node-based scalar field
    void RegisterNodalFab(Set::Field<Set::Scalar>& new_fab, int ncomp, int nghost, std::string name, bool writeout,std::vector<std::string> suffix = {});
    /// Add a new node-based scalar field (wtih additional arguments)
    void RegisterNodalFab(Set::Field<Set::Scalar>& new_fab, BC::BC<Set::Scalar>* new_bc, int ncomp, int nghost, std::string name, bool writeout,std::vector<std::string> suffix = {});
    template<class T>
    /// Add a templated nodal field
    void RegisterGeneralFab(Set::Field<T>& new_fab, int ncomp, int nghost, bool evolving = true);
    template<class T>
    /// Add a templated nodal field (additional arguments)
    void RegisterGeneralFab(Set::Field<T>& new_fab, int ncomp, int nghost, std::string a_name, bool evolving = true);
    template<class T>
    /// Add a templated nodal field (additional arguments)
    void RegisterGeneralFab(Set::Field<T>& new_fab, int ncomp, int nghost, bool writeout, std::string a_name, bool evolving = true);

    /// Add a field with arbitrary type (templated with T) and grid location (templated with d).
    template<class T, int d>
    void AddField(  Set::Field<T>& new_field, BC::BC<T>* new_bc, int ncomp, int nghost,
                    std::string, bool writeout, bool evolving, std::vector<std::string> suffix = {});

    template<class T, int d>
    void AddField(Set::Field<T>& new_field, std::string name, bool writeout, bool evolving );

    /// Utility to ensure that all fields know what the finest level is
    void SetFinestLevel(const int a_finestlevel)
    {
        for (unsigned int i = 0; i < cell.fab_array.size(); i++)
            cell.fab_array[i]->finest_level = a_finestlevel;
        for (unsigned int i = 0; i < node.fab_array.size(); i++)
            node.fab_array[i]->finest_level = a_finestlevel;
        for (unsigned int i = 0; i < m_basefields_cell.size(); i++)
            m_basefields_cell[i]->SetFinestLevel(finest_level);
        for (unsigned int i = 0; i < m_basefields.size(); i++)
            m_basefields[i]->SetFinestLevel(finest_level);
    }

    /// Register a variable to be integrated over the spatial domain using
    /// the Integrate function.
    void RegisterIntegratedVariable(Set::Scalar* integrated_variable, std::string name, bool extensive=true);

    /// Utility to set the coarse-grid timestep
    void SetTimestep(Set::Scalar _timestep);
    
    /// Utility to set the frequency (in timesteps) of plotfile dumping
    void SetPlotInt(int plot_int);
    
    /// Utility to set the frequency (in timesteps) of thermo data calculation
    void SetThermoInt(int a_thermo_int) { thermo.interval = a_thermo_int; }
    /// Utility to set the frequency (in timesteps) of thermo data writing to file
    void SetThermoPlotInt(int a_thermo_plot_int) { thermo.plot_int = a_thermo_plot_int; }
    /// Utility to set the global stop time.
    void SetStopTime(Set::Scalar a_stop_time) { stop_time = a_stop_time; }


    // Dynamic timestep adjustment

    struct {
        // user params
        bool on = false;
        int verbose = -1;
        int nprevious = -1;
        Set::Scalar cfl = NAN;
        Set::Scalar min = NAN;
        Set::Scalar max = NAN;
        // internal variables
        std::vector<Set::Scalar> dt_limit_min;
        std::vector<Set::Scalar> previous_timesteps;
    } dynamictimestep; /// Params for the dynamic timestp
    void DynamicTimestep_SyncTimeStep(int lev, Set::Scalar dt_min)
    {
        if (!dynamictimestep.on) return;
        if (!dynamictimestep.dt_limit_min.size())
        {
            dynamictimestep.dt_limit_min.resize(max_level+1,
                                                std::numeric_limits<Set::Scalar>::max());
        }
        dynamictimestep.dt_limit_min[lev] = std::min(dt_min,
                                                    dynamictimestep.dt_limit_min[lev]);
        
        amrex::ParallelDescriptor::ReduceRealMin(dynamictimestep.dt_limit_min[lev]);
    }
    void DynamicTimestep_Reset()
    {
        if (!dynamictimestep.on) return;
        dynamictimestep.previous_timesteps.clear();
    }
    void DynamicTimestep_Update()
    {
        if (!dynamictimestep.on) return;
        if (!dynamictimestep.dt_limit_min.size()) return;

        Set::Scalar final_timestep = NAN;

        if (amrex::ParallelDescriptor::IOProcessor())
        {
            Set::Scalar timestep_average = this->dt[0];
            if (dynamictimestep.previous_timesteps.size() > 0)
            {
                timestep_average = 0.0;
                for (unsigned int d = 0; d < dynamictimestep.previous_timesteps.size(); d++)
                    timestep_average += dynamictimestep.previous_timesteps[d];
                timestep_average /= dynamictimestep.previous_timesteps.size();
            }

            Set::Scalar new_timestep = std::numeric_limits<Set::Scalar>::max();
            for (int lev = 0; lev <= this->max_level; lev++)
            {
                //const Set::Scalar* DX = this->geom[lev].CellSize();
                Set::Scalar dt_lev = dynamictimestep.dt_limit_min[lev];

                Util::Message(INFO,"lev=",lev," ",dt_lev, " (",this->nsubsteps[lev],")");

                for (int ilev = lev; ilev > 0; ilev--) dt_lev *= (Set::Scalar)(this->nsubsteps[ilev]);

                Util::Message(INFO,"lev=",lev,"       -->    ",dt_lev);

                new_timestep = std::min(new_timestep,dt_lev);
            }

            if (new_timestep < timestep_average)
            {
                dynamictimestep.previous_timesteps.clear();

                final_timestep = new_timestep;
                final_timestep = std::max(final_timestep,dynamictimestep.min);
                final_timestep = std::min(final_timestep,dynamictimestep.max);

                dynamictimestep.previous_timesteps.push_back(new_timestep);
            }
            else
            {
                final_timestep = timestep_average;
                final_timestep = std::max(final_timestep,dynamictimestep.min);
                final_timestep = std::min(final_timestep,dynamictimestep.max);

                if ((int)(dynamictimestep.previous_timesteps.size()) > dynamictimestep.nprevious)
                    dynamictimestep.previous_timesteps.erase(dynamictimestep.previous_timesteps.begin()); // pop first
                dynamictimestep.previous_timesteps.push_back(new_timestep); // push back new timestep
            }

        }
        amrex::ParallelDescriptor::Bcast(&final_timestep,1);
        this->SetTimestep(final_timestep);
        dynamictimestep.dt_limit_min.clear();
    }



    amrex::Vector<amrex::Real> t_new;         ///< Keep track of current old simulation time on each level
    amrex::Vector<int> istep;           ///< Keep track of where each level is
    // PLOT FILES
    std::string plot_file{ "plt" };   ///< Plotfile name

private:
    virtual void MakeNewLevelFromScratch(int lev, amrex::Real time, const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override;
    virtual void MakeNewLevelFromCoarse(int lev, amrex::Real time, const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override;
    virtual void RemakeLevel(int lev, amrex::Real time, const amrex::BoxArray& ba,
        const amrex::DistributionMapping& dm) override;
    virtual void ClearLevel(int lev) override;
    virtual void ErrorEst(int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;


    /// This is the function that is responsible for updating patch data.
    void FillPatch(int lev, amrex::Real time,
        amrex::Vector<std::unique_ptr<amrex::MultiFab>>& source_mf,
        amrex::MultiFab& destination_multifab,
        BC::BC<Set::Scalar>& physbc,
        int icomp);
    /// Simple utility to count cells 
    long CountCells(int lev);
    /// Timestep marching
    void TimeStep(int lev, amrex::Real time, int iteration);
    void FillCoarsePatch(int lev, amrex::Real time, Set::Field<Set::Scalar>& mf, BC::BC<Set::Scalar>& physbc, int icomp, int ncomp);
    void GetData(const int lev, const amrex::Real time, amrex::Vector<amrex::MultiFab*>& data, amrex::Vector<amrex::Real>& datatime);

    std::vector<std::string> PlotFileName(int lev, std::string prefix = "") const;
protected:
    void IntegrateVariables(Set::Scalar cur_time, int step);
    void WritePlotFile(bool initial = false) const;
    void WritePlotFile(std::string prefix, Set::Scalar time, int step) const;
    void WritePlotFile(Set::Scalar time, amrex::Vector<int> iter, bool initial = false, std::string prefix = "") const;

    //
    // MEMBER VARIABLES
    //

    // TIME (STEP) KEEPINGamrex::Vector<std::unique_ptr<amrex::MultiFab> >
protected:
    amrex::Real timestep = NAN;   ///< Timestep for the base level of refinement
    amrex::Vector<amrex::Real> dt;  ///< Timesteps for each level of refinement
    amrex::Vector<int> nsubsteps;   ///< how many substeps on each level?
private:
    int max_plot_level = -1;

    amrex::Vector<amrex::Real> t_old;///< Keep track of current old simulation time on each level
    int max_step = std::numeric_limits<int>::max(); ///< Maximum allowable timestep
    amrex::Real tstart = 0; ///< Default start time (default: 0)
    amrex::Real stop_time = NAN; ///< Default stop time 

protected:
    bool integrate_variables_before_advance = true;
    bool integrate_variables_after_advance = false;

protected:
    struct {
        int number_of_fabs = 0;
        std::vector<Set::Field<Set::Scalar>*> fab_array;
        std::vector<int> ncomp_array;
        std::vector<int> nghost_array;
        std::vector<std::vector<std::string>> name_array;
        std::vector<BC::BC<Set::Scalar>*> physbc_array;
        std::vector<bool> writeout_array;
        bool any = true;
        bool all = false;
    } node;

    struct {
        int number_of_fabs = 0;
        std::vector<Set::Field<Set::Scalar>*> fab_array;
        std::vector<int> ncomp_array;
        std::vector<int> nghost_array;
        std::vector<std::vector<std::string>> name_array;
        std::vector<BC::BC<Set::Scalar>*> physbc_array;
        std::vector<bool> writeout_array;
        bool any = true;
        bool all = false;
    } cell;

    std::vector<BaseField*> m_basefields;
    std::vector<BaseField*> m_basefields_cell;

    BC::Nothing<Set::Scalar> bcnothing;

    // KEEP TRACK OF ALL INTEGRATED VARIABLES
    struct {
        int interval = -1;
        Set::Scalar dt = NAN;
        int plot_int = -1;
        Set::Scalar plot_dt = NAN;
        int number = 0;
        std::vector<Set::Scalar*> vars;
        std::vector<std::string> names;
        std::vector<bool> extensives;
    } thermo;

    // REGRIDDING
    int regrid_int = -1;       ///< Determine how often to regrid (default: 2)
    int base_regrid_int = -1; ///< Determine how often to regrid based on coarse level only (default: 0)

    std::string restart_file_cell = "";
    std::string restart_file_node = "";

    struct {
        int on = 0;
        std::vector<amrex::Box> box;
    } explicitmesh;

protected:
    int plot_int = -1;               ///< How frequently to dump plot file (default: never)
    Set::Scalar plot_dt = -1.0;

    int abort_on_nan = true;
};


template<>
ALAMO_SINGLE_DEFINITION
void Integrator::AddField<Set::Scalar, Set::Hypercube::Cell>
(   Set::Field<Set::Scalar>& new_field,
    BC::BC<Set::Scalar>* new_bc,
    int ncomp,
    int nghost,
    std::string name,
    bool writeout,
    bool /*evolving*/,
    std::vector<std::string> suffix)
{
    int nlevs_max = maxLevel() + 1;
    new_field.resize(nlevs_max);
    cell.fab_array.push_back(&new_field);
    if (new_bc != nullptr) cell.physbc_array.push_back(new_bc);
    else                   cell.physbc_array.push_back(&bcnothing);
    cell.ncomp_array.push_back(ncomp);
    cell.nghost_array.push_back(nghost);
    //cell.name_array.push_back(name);
    cell.writeout_array.push_back(writeout);
    cell.number_of_fabs++;

    Util::Assert(INFO,TEST((int)suffix.size() == 0 || (int)suffix.size() == ncomp));
    std::vector<std::string> names;
    if (ncomp == 1)
        names.push_back(name);
    else
    {
        if (suffix.size() == 0)
            for (int j = 0; j < ncomp; j++)
                names.push_back(amrex::Concatenate(name,j+1,3));
        else
            for (int j = 0; j < ncomp; j++)
                names.push_back(name + suffix[j]);
    }
    cell.name_array.push_back(names);
}

template<>
ALAMO_SINGLE_DEFINITION
void Integrator::AddField<Set::Scalar, Set::Hypercube::Node>
(   Set::Field<Set::Scalar>& new_field,
    BC::BC<Set::Scalar>* new_bc,
    int ncomp,
    int nghost,
    std::string name,
    bool writeout,
    bool /*evolving*/,
    std::vector<std::string> suffix)
{
    BL_PROFILE("Integrator::RegisterNodalFab");
    Util::Assert(INFO, TEST(new_bc == nullptr));
    int nlevs_max = maxLevel() + 1;
    new_field.resize(nlevs_max);
    node.fab_array.push_back(&new_field);
    node.physbc_array.push_back(&bcnothing);
    node.ncomp_array.push_back(ncomp);
    node.nghost_array.push_back(nghost);
    //node.name_array.push_back(name);
    node.writeout_array.push_back(writeout);
    node.number_of_fabs++;

    Util::Assert(INFO,TEST((int)suffix.size() == 0 || (int)suffix.size() == ncomp));
    std::vector<std::string> names;
    if (ncomp == 1)
        names.push_back(name);
    else
    {
        if (suffix.size() == 0)
            for (int j = 0; j < ncomp; j++)
                names.push_back(amrex::Concatenate(name,j+1,3));
        else
            for (int j = 0; j < ncomp; j++)
                names.push_back(name + suffix[j]);
    }
    node.name_array.push_back(names);
}

template<class T, int d>
ALAMO_SINGLE_DEFINITION
void Integrator::AddField
(   Set::Field<T>& new_field, BC::BC<T>* new_bc, int ncomp,
    int nghost, std::string name, bool writeout, bool evolving,
    std::vector<std::string> /*suffix*/)
{
    if (d == Set::Hypercube::Node)
    {
        Util::Assert(INFO, TEST(new_bc == nullptr));
        int nlevs_max = maxLevel() + 1;
        new_field.resize(nlevs_max);
        m_basefields.push_back(new Field<T>(new_field, geom, refRatio(), ncomp, nghost));
        m_basefields.back()->evolving = evolving;
        m_basefields.back()->writeout = writeout;
        m_basefields.back()->setName(name);
        m_basefields.back()->evolving = evolving;
        m_basefields.back()->m_gridtype = Set::Hypercube::Node;
    }
    else if (d == Set::Hypercube::Cell)
    {
        int nlevs_max = maxLevel() + 1;
        new_field.resize(nlevs_max);
        m_basefields_cell.push_back(new Field<T>(new_field, geom, refRatio(), ncomp, nghost));
        m_basefields_cell.back()->evolving = evolving;
        m_basefields_cell.back()->writeout = writeout;
        m_basefields_cell.back()->setName(name);
        m_basefields_cell.back()->evolving = evolving;
        if (new_bc) m_basefields_cell.back()->setBC(new_bc);
        else m_basefields_cell.back()->setBC(&bcnothing);
        m_basefields_cell.back()->m_gridtype = Set::Hypercube::Cell;
    }
    else
    {
        Util::Abort(INFO, "Only node and cell based fields can be added at this time");
    }
}

/// Convenience function for fields that have a single component,
/// no ghost nodes/cells, and (therefore) no boundary condition.
template<class T, int d>
ALAMO_SINGLE_DEFINITION
void Integrator::AddField(
    Set::Field<T>& new_field, std::string name, bool writeout, bool evolving )
{
    AddField<T,d>(new_field,new BC::Nothing<T>, 1, 0, name, writeout, evolving);
}


template<class T>
ALAMO_SINGLE_DEFINITION
void Integrator::RegisterGeneralFab(Set::Field<T>& new_fab, int ncomp, int nghost, bool evolving) 
{
    //Util::Warning(INFO, "RegisterGeneralFab is depricated. Please replace with AddField");
    AddField<T, Set::Hypercube::Node>(new_fab, nullptr, ncomp, nghost, "", true, evolving);
}
template<class T>
ALAMO_SINGLE_DEFINITION
void Integrator::RegisterGeneralFab(Set::Field<T>& new_fab, int ncomp, int nghost, std::string a_name, bool evolving)
{
    //Util::Warning(INFO, "RegisterGeneralFab is depricated. Please replace with AddField");
    AddField<T, Set::Hypercube::Node>(new_fab, nullptr, ncomp, nghost, a_name, true, evolving);
}
template<class T>
AMREX_ATTRIBUTE_WEAK
void Integrator::RegisterGeneralFab(Set::Field<T>& new_fab, int ncomp, int nghost, bool writeout, std::string a_name, bool evolving)
{
    //Util::Warning(INFO, "RegisterGeneralFab is depricated. Please replace with AddField");
    AddField<T, Set::Hypercube::Node>(new_fab, nullptr, ncomp, nghost, a_name, writeout, evolving);
}

}
#endif
