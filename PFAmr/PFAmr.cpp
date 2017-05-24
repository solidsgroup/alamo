
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>

#include <PFAmr.H>
#include <PFAmrBC.H>

using namespace amrex;

PFAmr::PFAmr ()
{
  ReadParameters();

  // Geometry on all levels has been defined already.

  // No valid BoxArray and DistributionMapping have been defined.
  // But the arrays for them have been resized.

  int nlevs_max = maxLevel() + 1;

  istep.resize(nlevs_max, 0);
  nsubsteps.resize(nlevs_max, 1);
  for (int lev = 1; lev <= maxLevel(); ++lev) {
    nsubsteps[lev] = MaxRefRatio(lev-1);
  }

  t_new.resize(nlevs_max, 0.0);
  t_old.resize(nlevs_max, -1.e100);
  dt.resize(nlevs_max, 1.e100);

  phi_new.resize(nlevs_max);
  phi_old.resize(nlevs_max);
}

PFAmr::~PFAmr ()
{
  
}

void
PFAmr::ReadParameters ()
{
  {
    ParmParse pp;  // Traditionally, max_step and stop_time do not have prefix.
    pp.query("max_step", max_step);
    pp.query("stop_time", stop_time);
  }

  {
    ParmParse pp("amr"); // Traditionally, these have prefix, amr.

    pp.query("regrid_int", regrid_int);

    pp.query("check_file", check_file);
    pp.query("check_int", check_int);

    pp.query("plot_file", plot_file);
    pp.query("plot_int", plot_int);

    pp.query("restart", restart_chkfile);
  }
}

void
PFAmr::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
				const DistributionMapping& dm)
{
  const int ncomp = phi_new[lev-1]->nComp();
  const int nghost = phi_new[lev-1]->nGrow();
    
  phi_new[lev].reset(new MultiFab(ba, dm, ncomp, nghost));
  phi_old[lev].reset(new MultiFab(ba, dm, ncomp, nghost));

  t_new[lev] = time;
  t_old[lev] = time - 1.e200;

  FillCoarsePatch(lev, time, *phi_new[lev], 0, ncomp);
}


void
PFAmr::RemakeLevel (int lev, Real time, const BoxArray& ba,
		     const DistributionMapping& dm)
{
  const int ncomp = phi_new[lev]->nComp();
  const int nghost = phi_new[lev]->nGrow();

#if __cplusplus >= 201402L
  auto new_state = std::make_unique<MultiFab>(ba, dm, ncomp, nghost);
  auto old_state = std::make_unique<MultiFab>(ba, dm, ncomp, nghost);
#else
  std::unique_ptr<MultiFab> new_state(new MultiFab(ba, dm, ncomp, nghost));
  std::unique_ptr<MultiFab> old_state(new MultiFab(ba, dm, ncomp, nghost));
#endif

  FillPatch(lev, time, *new_state, 0, ncomp);

  std::swap(new_state, phi_new[lev]);
  std::swap(old_state, phi_old[lev]);

  t_new[lev] = time;
  t_old[lev] = time - 1.e200;

}

void
PFAmr::ClearLevel (int lev)
{
  phi_new[lev].reset(nullptr);
  phi_old[lev].reset(nullptr);
}

long
PFAmr::CountCells (int lev)
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
PFAmr::FillPatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
  if (lev == 0)
    {
      Array<MultiFab*> smf;
      Array<Real> stime;
      GetData(0, time, smf, stime);

      PFAmrPhysBC physbc;
      amrex::FillPatchSingleLevel(mf, time, smf, stime, 0, icomp, ncomp,
				  geom[lev], physbc);
    }
  else
    {
      Array<MultiFab*> cmf, fmf;
      Array<Real> ctime, ftime;
      GetData(lev-1, time, cmf, ctime);
      GetData(lev  , time, fmf, ftime);

      PFAmrPhysBC cphysbc, fphysbc;
      Interpolater* mapper = &cell_cons_interp;

      int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR}; // periodic boundaryies
      int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
      Array<BCRec> bcs(1, BCRec(lo_bc, hi_bc));

      amrex::FillPatchTwoLevels(mf, time, cmf, ctime, fmf, ftime,
				0, icomp, ncomp, geom[lev-1], geom[lev],
				cphysbc, fphysbc, refRatio(lev-1),
				mapper, bcs);
    }
}

void
PFAmr::FillCoarsePatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp)
{
  BL_ASSERT(lev > 0);

  Array<MultiFab*> cmf;
  Array<Real> ctime;
  GetData(lev-1, time, cmf, ctime);
    
  if (cmf.size() != 1) {
    amrex::Abort("FillCoarsePatch: how did this happen?");
  }

  PFAmrPhysBC cphysbc, fphysbc;
  Interpolater* mapper = &cell_cons_interp;
    
  int lo_bc[] = {INT_DIR, INT_DIR, INT_DIR}; // periodic boundaryies
  int hi_bc[] = {INT_DIR, INT_DIR, INT_DIR};
  Array<BCRec> bcs(1, BCRec(lo_bc, hi_bc));

  amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
			       cphysbc, fphysbc, refRatio(lev-1),
			       mapper, bcs);
}

void
PFAmr::GetData (int lev, Real time, Array<MultiFab*>& data, Array<Real>& datatime)
{
  data.clear();
  datatime.clear();

  const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

  if (time > t_new[lev] - teps && time < t_new[lev] + teps)
    {
      data.push_back(phi_new[lev].get());
      datatime.push_back(t_new[lev]);
    }
  else if (time > t_old[lev] - teps && time < t_old[lev] + teps)
    {
      data.push_back(phi_old[lev].get());
      datatime.push_back(t_old[lev]);
    }
  else
    {
      data.push_back(phi_old[lev].get());
      data.push_back(phi_new[lev].get());
      datatime.push_back(t_old[lev]);
      datatime.push_back(t_new[lev]);
    }
}

