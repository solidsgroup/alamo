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

#include "GeneralAMRIntegratorBC.H"

class GeneralAMRIntegrator
  : public amrex::AmrCore
{
public:
  GeneralAMRIntegrator ();
  virtual ~GeneralAMRIntegrator();

  void InitData ();
  void Evolve ();

  virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
					const amrex::DistributionMapping& dm) override;
  virtual void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
				       const amrex::DistributionMapping& dm) override;
  virtual void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
			    const amrex::DistributionMapping& dm) override;
  virtual void ClearLevel (int lev) override;
  virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) override;

  void RegisterNewFab(amrex::Array<std::unique_ptr<amrex::MultiFab> > &new_fab, GeneralAMRIntegratorPhysBC &new_bc, int ncomp, int nghost, std::string name);

  virtual void Initialize (int lev) = 0;

private:

  void InitFromCheckpoint ();

  long CountCells (int lev);

  void TimeStep (int lev, amrex::Real time, int iteration);
protected:
  virtual void Advance (int lev, amrex::Real time, amrex::Real dt) = 0;

  virtual void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow) = 0;
  

  void FillPatch (int lev, amrex::Real time, amrex::Array<std::unique_ptr<amrex::MultiFab> > &source_mf, amrex::MultiFab &destination_multifab, GeneralAMRIntegratorPhysBC &physbc, int icomp);
private:
  void FillCoarsePatch (int lev, amrex::Real time, amrex::MultiFab& mf, GeneralAMRIntegratorPhysBC &physbc, int icomp, int ncomp);
  void GetData (const int lev, const amrex::Real time, amrex::Array<amrex::MultiFab*>& data, amrex::Array<amrex::Real>& datatime);

  void CreateCleanDirectory () const;
  std::vector<std::string> PlotFileName (int lev) const;
  void WritePlotFile () const;

  amrex::Array<int> istep;      // which step?
  amrex::Array<int> nsubsteps;  // how many substeps on each level?

protected:
  amrex::Array<amrex::Real> t_new;
  amrex::Array<amrex::Real> t_old;
  amrex::Array<amrex::Real> dt;

  amrex::Array<std::unique_ptr<amrex::Array<std::unique_ptr<amrex::MultiFab> > > > fab_array;
  amrex::Array<int> ncomp_array;  // Array containing number of components for each fab array
  amrex::Array<int> nghost_array; // Array containing number of ghost cells in each fab array

  amrex::Array<std::string> name_array;
private:

  int max_step   = std::numeric_limits<int>::max();
  amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();

  amrex::Real tstart = 0;
  int regrid_int = 2;

  std::string restart_chkfile;
  std::string check_file {"chk"};
  int check_int = -1;

  std::string plot_file {"plt"};
  int plot_int = -1;

  amrex::Real timestep;
  int number_of_fabs = 0;

  amrex::Array<std::unique_ptr<GeneralAMRIntegratorPhysBC> > physbc_array;
  //GeneralAMRIntegratorPhysBC physbc;
};

#endif
