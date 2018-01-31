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

  virtual void MakeNewLevelFromScratch (int lev, amrex::Real time, const amrex::BoxArray& ba,
   					const amrex::DistributionMapping& dm);
  virtual void ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int /*ngrow*/);

protected:

  virtual void Advance (int lev, amrex::Real time, amrex::Real dt);

private:


  int number_of_fabs = 1;
  int number_of_components = 2;

  amrex::Real M;                  // Burn speed parameter
  amrex::Real kappa;              // Burn width parameter
  amrex::Real w1, w12, w0;        // Energy (before, during, after)
  amrex::Real rho1, rho0;         // Density (before and after)
  amrex::Real k1, k0;             // Thermal conductivity (before and after)
  amrex::Real cp1, cp0;           // Specific heat (before and after)
  amrex::Real qdotburn;           // Heat lost in flame

};

#endif
