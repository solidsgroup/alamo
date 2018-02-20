#ifndef PFFLAME_H
#define PFFLAME_H

#include <string>
#include <limits>
#include <memory>


#include "GeneralAMRIntegrator/GeneralAMRIntegratorBC.H"
#include "GeneralAMRIntegrator/GeneralAMRIntegrator.H"

#include "PFFlame/PFFlameBC.H"

class PFFlame : public GeneralAMRIntegrator
{
public:
  PFFlame ();

protected:
  void Initialize (int lev) ;
  void Advance (int lev, Real time, Real dt);
  void TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/);

private:

  amrex::Array<std::unique_ptr<amrex::MultiFab> > Temp;
  amrex::Array<std::unique_ptr<amrex::MultiFab> > Temp_old;
  amrex::Array<std::unique_ptr<amrex::MultiFab> > Eta;
  amrex::Array<std::unique_ptr<amrex::MultiFab> > Eta_old;
  PFFlamePhysBC TempBC;
  PFFlamePhysBC EtaBC;

  // int number_of_fabs = 1;
  // int number_of_components = 2;

  amrex::Real M;                  // Burn speed parameter
  amrex::Real kappa;              // Burn width parameter
  amrex::Real w1, w12, w0;        // Energy (before, during, after)
  amrex::Real rho1, rho0;         // Density (before and after)
  amrex::Real k1, k0;             // Thermal conductivity (before and after)
  amrex::Real cp1, cp0;           // Specific heat (before and after)
  amrex::Real qdotburn;           // Heat lost in flame

};

#endif
