#ifndef AMR_Adv_C_V2_Adv_H_
#define AMR_Adv_C_V2_Adv_H_

#include <string>
#include <limits>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_FluxRegister.H>

class PFAmr
  : public amrex::AmrCore
{
public:
  PFAmr ();
  virtual ~PFAmr();

  void InitData ();
  void Evolve ();

  virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
					const amrex::DistributionMapping& dm) override;
  virtual void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
				       const amrex::DistributionMapping& dm) override;
  virtual void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
			    const amrex::DistributionMapping& dm) override;
  virtual void ClearLevel (int lev) override;
  virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int /*ngrow*/) override;

private:

  void InitFromCheckpoint ();

  long CountCells (int lev);

  void TimeStep (int lev, amrex::Real time, int iteration);
  void Advance (int lev, amrex::Real time, amrex::Real dt);

  void FillPatch (int lev, amrex::Real time, amrex::Array<std::unique_ptr<amrex::MultiFab> >& mf, int icomp);
  void FillCoarsePatch (int lev, amrex::Real time, amrex::MultiFab& mf, int icomp, int ncomp);
  void GetData (int lev, amrex::Real time, amrex::Array<amrex::Array<amrex::MultiFab*> >& data, amrex::Array<amrex::Real>& datatime);

  std::string PlotFileName (int lev) const;
  amrex::Array<const amrex::MultiFab*> PlotFileMF (int grain) const;
  amrex::Array<std::string> PlotFileVarNames () const;
  void WritePlotFile () const;

  amrex::Array<int> istep;      // which step?
  amrex::Array<int> nsubsteps;  // how many substeps on each level?

  amrex::Array<amrex::Real> t_new;
  amrex::Array<amrex::Real> t_old;
  amrex::Array<amrex::Real> dt;

  amrex::Array<amrex::Array<std::unique_ptr<amrex::MultiFab> > > phi_new;
  amrex::Array<amrex::Array<std::unique_ptr<amrex::MultiFab> > > phi_old;

  // runtime parameters

  int max_step   = std::numeric_limits<int>::max();
  amrex::Real stop_time = std::numeric_limits<amrex::Real>::max();

  int regrid_int = 2;

  std::string restart_chkfile;

  std::string check_file {"chk"};
  std::string plot_file {"plt"};
  int check_int = -1;
  int plot_int = -1;

  // Phase Field Parameters
  amrex::Real timestep;
  int number_of_grains = 2;
  amrex::Real L, mu, gamma, kappa;

};

#endif
