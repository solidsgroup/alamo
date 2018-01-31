#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_BC_TYPES.H>

#include "GeneralAMRIntegrator.H"
#include "GeneralAMRIntegratorBC.H"


using namespace amrex;

GeneralAMRIntegrator::GeneralAMRIntegrator ()
  : physbc(geom)
{
  { 
    ParmParse pp;   // Basic run parameters
    pp.query("max_step", max_step);
    pp.query("stop_time", stop_time);
    pp.query("timestep",timestep);
  }

  {
    ParmParse pp("amr"); // AMR specific parameters
    pp.query("regrid_int", regrid_int);
    pp.query("check_file", check_file);
    pp.query("check_int", check_int);
    pp.query("plot_file", plot_file);
    pp.query("plot_int", plot_int);
    pp.query("restart", restart_chkfile);
  }

  int nlevs_max = maxLevel() + 1;

  istep.resize(nlevs_max, 0);
  nsubsteps.resize(nlevs_max, 1);
  for (int lev = 1; lev <= maxLevel(); ++lev) {
    nsubsteps[lev] = MaxRefRatio(lev-1);
  }

  t_new.resize(nlevs_max, 0.0);
  t_old.resize(nlevs_max, -1.e100);
  dt.resize(nlevs_max, 1.e100);
  for (int i = 0; i < nlevs_max; i++)
    dt[i] = timestep;

  CreateCleanDirectory();
}

void
GeneralAMRIntegrator::RegisterNewFab(amrex::Array<std::unique_ptr<amrex::MultiFab> > &new_fab,
				     int ncomp,
				     int nghost)
{
  int nlevs_max = maxLevel() + 1;

  new_fab.resize(nlevs_max);
  fab_array.push_back((std::unique_ptr<amrex::Array<std::unique_ptr<amrex::MultiFab> > >)&new_fab);
  ncomp_array.push_back(ncomp);
  nghost_array.push_back(nghost);
  number_of_fabs++;
}


GeneralAMRIntegrator::~GeneralAMRIntegrator ()
{
  // TODO
}

// Part of REGRIDDING
void
GeneralAMRIntegrator::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
				const DistributionMapping& dm)
{
  for (int n = 0; n < number_of_fabs; n++)
    {
      const int ncomp = (*fab_array[n])[lev-1]->nComp();
      const int nghost = (*fab_array[n])[lev-1]->nGrow();

      (*fab_array[n])[lev].reset(new MultiFab(ba, dm, ncomp, nghost));

      t_new[lev] = time;
      t_old[lev] = time - 1.e200;

      FillCoarsePatch(lev, time, *(*fab_array[n])[lev], 0, ncomp);
    }
}


// Part of REGRIDDING
void
GeneralAMRIntegrator::RemakeLevel (int lev, Real time, const BoxArray& ba,
		     const DistributionMapping& dm)
{
  amrex::Array<std::unique_ptr<MultiFab> > new_state(number_of_fabs); 
  amrex::Array<std::unique_ptr<MultiFab> > old_state(number_of_fabs);
  for (int n=0; n < number_of_fabs; n++)
    {
      const int ncomp = (*fab_array[n])[lev]->nComp();
      const int nghost = (*fab_array[n])[lev]->nGrow();

      new_state[n].reset(new MultiFab(ba, dm, ncomp, nghost));
    }

  FillPatch(lev, time, new_state, 0);

  for (int n=0; n < number_of_fabs; n++)
    {
      std::swap(new_state[n], (*fab_array[n])[lev]);

      t_new[lev] = time;
      t_old[lev] = time - 1.e200;
    }
}

void
GeneralAMRIntegrator::ClearLevel (int lev)
{
  for (int n = 0; n < number_of_fabs; n++)
    {
      (*fab_array[n])[lev].reset(nullptr);
    }
}

long
GeneralAMRIntegrator::CountCells (int lev)
{
  const int N = grids[lev].size();

  long cnt = 0;

  for (int i = 0; i < N; ++i)
    {
      cnt += grids[lev][i].numPts();
    }

  return cnt;
}

void
GeneralAMRIntegrator::FillPatch (int lev, Real time, Array<std::unique_ptr<MultiFab> >& mf, int icomp)
{
  if (lev == 0)
    {
      Array<Array<MultiFab*> > smf;
      smf.resize(number_of_fabs);
      Array<Real> stime;
      GetData(0, time, smf, stime);

      physbc.SetLevel(lev);
      for (int n = 0; n<number_of_fabs; n++)
	{
	  amrex::FillPatchSingleLevel(*mf[n], time, smf[n], stime, 0, icomp, mf[n]->nComp(),
				      geom[lev], physbc);
	}
    }
  else
    {
      Array<Array<MultiFab*> > cmf, fmf;
      cmf.resize(number_of_fabs);
      fmf.resize(number_of_fabs);
      Array<Real> ctime, ftime;
      GetData(lev-1, time, cmf, ctime);
      GetData(lev  , time, fmf, ftime);

      physbc.SetLevel(lev);
      Interpolater* mapper = &cell_cons_interp;

      for (int n = 0; n<number_of_fabs; n++)
	{
	  Array<BCRec> bcs(mf[n]->nComp(), physbc.GetBCRec()); // todo
	  amrex::FillPatchTwoLevels(*mf[n], time, cmf[n], ctime, fmf[n], ftime,
				    0, icomp, mf[n]->nComp(), geom[lev-1], geom[lev],
				    physbc, physbc,
				    refRatio(lev-1),
				    mapper, bcs);
	}
    }
}

void
GeneralAMRIntegrator::FillCoarsePatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
  BL_ASSERT(lev > 0);

  Array<Array<MultiFab*> > cmf(number_of_fabs);
  Array<Real> ctime;
  GetData(lev-1, time, cmf, ctime);
    
  for (int n = 0; n < number_of_fabs; n++)
    if (cmf.size() != 1) 
      amrex::Abort("FillCoarsePatch: how did this happen?");
  physbc.SetLevel(lev);
  Interpolater* mapper = &cell_cons_interp;
    
  for (int n = 0; n < number_of_fabs; n++)
    {
      Array<BCRec> bcs(ncomp, physbc.GetBCRec());
      amrex::InterpFromCoarseLevel(mf, time, *cmf[n][0], 0, icomp, ncomp, geom[lev-1], geom[lev],
				   physbc, physbc,
				   refRatio(lev-1),
				   mapper, bcs);
    }
}
 
void
GeneralAMRIntegrator::GetData (const int lev, const Real time,
			       Array<Array<MultiFab*> >& data, Array<Real>& datatime)
{
  for (int n = 0; n < number_of_fabs; n++) data[n].clear();
  datatime.clear();

  const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

  // TODO - Fix this
  // if (time > t_new[lev] - teps && time < t_new[lev] + teps)
  //   {
  //     for (int n = 0; n < number_of_fabs; n++) data[n].push_back(phi_new[n][lev].get());
  //     datatime.push_back(t_new[lev]);
  //   }
  // else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
  //   {
  //     for (int n = 0; n < number_of_fabs; n++) data[n].push_back(phi_old[n][lev].get());
  //     datatime.push_back(t_old[lev]);
  //   }
  // else
    {
      for (int n = 0; n < number_of_fabs; n++)
	{
	  //data[n].push_back(phi_old[n][lev].get());
	  data[n].push_back((*fab_array[n])[lev].get());
	}
      //datatime.push_back(t_old[lev]);
      datatime.push_back(t_new[lev]);
    }
}

