#ifndef PFFLAME_H
#define PFFLAME_H

#include <string>
#include <limits>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_AmrCore.H>
#include <AMReX_FluxRegister.H>

#include "PFAmr.H"


#ifndef DEBUG
#define DEBUG std::cout << __FILE__ << ":" << __LINE__ << std::endl;
#endif

class PFFlame : public PFAmr
{
public:
  PFFlame ();
  // virtual ~PFFlame();

  void InitData ();
  void Evolve ();

  virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
					const amrex::DistributionMapping& dm);
  // virtual void MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& ba,
  // 				       const amrex::DistributionMapping& dm);
  // virtual void RemakeLevel (int lev, amrex::Real time, const amrex::BoxArray& ba,
  // 			    const amrex::DistributionMapping& dm);
  // virtual void ClearLevel (int lev);
  virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int /*ngrow*/);

private:

  void InitFromCheckpoint ();

  // long CountCells (int lev);

  void TimeStep (int lev, amrex::Real time, int iteration);
  void Advance (int lev, amrex::Real time, amrex::Real dt);

  // void FillPatch (int lev, amrex::Real time, amrex::Array<std::unique_ptr<amrex::MultiFab> >& mf, int icomp);
  // void FillCoarsePatch (int lev, amrex::Real time, amrex::MultiFab& mf, int icomp, int ncomp);
  // void GetData (int lev, amrex::Real time, amrex::Array<amrex::Array<amrex::MultiFab*> >& data, amrex::Array<amrex::Real>& datatime);

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
  amrex::Real timestep;

  int regrid_int = 2;
  std::string restart_chkfile;
  std::string check_file {"chk"};
  std::string plot_file {"plt"};
  int check_int = -1;
  int plot_int = -1;

  // Framework Parameters
  int number_of_fabs = 1;
  int number_of_components = 2;

  // Physics parametrs
  amrex::Real M;                  // Burn speed parameter
  amrex::Real kappa;              // Burn width parameter
  amrex::Real w1, w12, w0;        // Energy (before, during, after)
  amrex::Real rho1, rho0;         // Density (before and after)
  amrex::Real k1, k0;             // Thermal conductivity (before and after)
  amrex::Real cp1, cp0;             // Specific heat (before and after)
  amrex::Real qdotburn;           // Heat lost in flame

  //amrex::Real L, mu, gamma, kappa;

  // More phase field parameters


  // Voronoi
  amrex::Array<amrex::Real> voronoi_x;
  amrex::Array<amrex::Real> voronoi_y;
};

#endif
