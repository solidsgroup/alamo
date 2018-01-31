#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>

#include <AMReX_PhysBCFunct.H>

#include <PFFlame.H>
#include <PFFlameBC.H>

using namespace amrex;

PFFlame::PFFlame () :
  PFAmr()
{
  //
  // READ INPUT PARAMETERS
  //

  {
    ParmParse pp("physics"); // Phase-field model parameters
    pp.query("M",M);
    pp.query("kappa",kappa);
    pp.query("w1",w1);
    pp.query("w12",w12);
    pp.query("w0",w0);
    pp.query("rho1",rho1);
    pp.query("rho0",rho0);
    pp.query("k1",k1);
    pp.query("k0",k0);
    pp.query("cp1",cp1);
    pp.query("cp0",cp0);
    pp.query("qdotburn",qdotburn);
  }

  // int nlevs_max = maxLevel() + 1;

  // istep.resize(nlevs_max, 0);
  // nsubsteps.resize(nlevs_max, 1);
  // for (int lev = 1; lev <= maxLevel(); ++lev) {
  //   nsubsteps[lev] = MaxRefRatio(lev-1);
  // }

  // phi_new.resize(number_of_fabs + 1); // <-- resize for number of order parameters
  // phi_old.resize(number_of_fabs + 1); //
  // for (int n = 0; n < number_of_fabs; n++)
  //   {
  //     phi_new[n].resize(nlevs_max);
  //     phi_old[n].resize(nlevs_max);
  //   }

}

