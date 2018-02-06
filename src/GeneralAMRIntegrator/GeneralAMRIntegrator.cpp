///
/// \file GeneralAMRIntegrator.cpp
/// \brief Compute the volume integral of two multiplied Fourier series
///

#include <chrono>
#include <ctime>

#include <AMReX_ParallelDescriptor.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_FillPatchUtil.H>
#include <AMReX_BC_TYPES.H>

#include "GeneralAMRIntegrator.H"
#include "GeneralAMRIntegratorBC.H"


using namespace amrex;

GeneralAMRIntegrator::GeneralAMRIntegrator ()
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
  dt[0] = timestep;
  for (int i = 1; i < nlevs_max; i++)
    dt[i] = dt[i-1] / (amrex::Real)nsubsteps[i];

  amrex::UtilCreateCleanDirectory(plot_file, false);
  WriteMetaData();
}

///
/// \func  ~GeneralAMRIntegrator
/// \brief Does nothing -- check here first if there are memory leaks
///
GeneralAMRIntegrator::~GeneralAMRIntegrator ()
{
}

/// 
/// \func  MakeNewLevelFromCoarse
/// \brief Wrapper to call FillCoarsePatch
///
/// Note: **THIS OVERRIDES A PURE VIRTUAL METHOD - DO NOT CHANGE**
///
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

      FillCoarsePatch(lev, time, *(*fab_array[n])[lev], *physbc_array[n], 0, ncomp);
    }
}


///
/// RESETS ALL MULTIFABS AT A GIVEN LEVEL
///
/// (OVERRIDES PURE VIRTUAL METHOD - DO NOT CHANGE)
///
void
GeneralAMRIntegrator::RemakeLevel (int lev,       ///<[in] AMR Level
				   Real time,     ///<[in] Simulation time
				   const BoxArray& ba, 
				   const DistributionMapping& dm)
{
  for (int n=0; n < number_of_fabs; n++)
    {

      const int ncomp = (*fab_array[n])[lev]->nComp();
      const int nghost = (*fab_array[n])[lev]->nGrow();

      MultiFab new_state(ba, dm, ncomp, nghost); 

      FillPatch(lev, time, *fab_array[n], new_state, *physbc_array[n], 0);

      std::swap(new_state, *(*fab_array[n])[lev]);
    }

  t_new[lev] = time;
  t_old[lev] = time - 1.e200;
}

//
// DELETE EVERYTHING
//
// (OVERRIDES PURE VIRTUAL METHOD - DO NOT CHANGE)
//
void
GeneralAMRIntegrator::ClearLevel (int lev)
{
  for (int n = 0; n < number_of_fabs; n++)
    {
      (*fab_array[n])[lev].reset(nullptr);
    }
}

//
//
//

void // CUSTOM METHOD - CHANGEABLE
GeneralAMRIntegrator::RegisterNewFab(amrex::Array<std::unique_ptr<amrex::MultiFab> > &new_fab,
				     GeneralAMRIntegratorPhysBC &new_bc,
				     int ncomp,
				     int nghost,
				     std::string name)
{
  int nlevs_max = maxLevel() + 1;

  new_fab.resize(nlevs_max);
  fab_array.push_back((std::unique_ptr<amrex::Array<std::unique_ptr<amrex::MultiFab> > >)&new_fab);
  physbc_array.push_back((std::unique_ptr<GeneralAMRIntegratorPhysBC>)&new_bc); 
  ncomp_array.push_back(ncomp);
  nghost_array.push_back(nghost);
  if (ncomp > 1) for (int i = 1; i <= ncomp; i++) name_array.push_back(amrex::Concatenate(name, i, 3));
  else name_array.push_back(name);
  number_of_fabs++;
}


long // CUSTOM METHOD - CHANGEABLE
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

void  // CUSTOM METHOD - CHANGEABLE
GeneralAMRIntegrator::FillPatch (int lev, Real time,
				 Array<std::unique_ptr<MultiFab> > &source_mf,
				 MultiFab &destination_mf,
				 GeneralAMRIntegratorPhysBC &physbc, int icomp)
{
  if (lev == 0)
    {
      
      Array<MultiFab*> smf;
      smf.push_back(source_mf[lev].get());
      Array<Real> stime;
      stime.push_back(time);

      // GetData(lev, time, smf, stime);

      physbc.SetLevel(lev);
      amrex::FillPatchSingleLevel(destination_mf,		// Multifab
				  time,                         // time
				  smf,				// Vector<MultiFab*> &smf (CONST)
				  stime,			// Vector<Real> &stime    (CONST)
				  0,				// scomp - Source component 
				  icomp,			// dcomp - Destination component
				  destination_mf.nComp(),	// ncomp - Number of components
				  geom[lev],			// Geometry (CONST)
				  physbc);			// BC
    } 
  else
    {
      Array<MultiFab*> cmf, fmf;
      cmf.push_back(source_mf[lev-1].get());
      fmf.push_back(source_mf[lev].get());
      Array<Real> ctime, ftime;
      ctime.push_back(time);
      ftime.push_back(time);
      // GetData(lev-1, time, cmf, ctime);
      // GetData(lev  , time, fmf, ftime);

      physbc.SetLevel(lev);
      Interpolater* mapper = &cell_cons_interp;

      Array<BCRec> bcs(destination_mf.nComp(), physbc.GetBCRec()); // todo
      amrex::FillPatchTwoLevels(destination_mf, time, cmf, ctime, fmf, ftime,
				0, icomp, destination_mf.nComp(), geom[lev-1], geom[lev],
				physbc, physbc,
				refRatio(lev-1),
				mapper, bcs);
    }
}

void // CUSTOM METHOD - CHANGEABLE
GeneralAMRIntegrator::FillCoarsePatch (int lev, Real time, MultiFab& mf, GeneralAMRIntegratorPhysBC &physbc, int icomp, int ncomp)
{
  BL_ASSERT(lev > 0);

  Array<MultiFab* > cmf(number_of_fabs);
  Array<Real> ctime;
  GetData(lev-1, time, cmf, ctime);
    
  if (cmf.size() != 1) 
    amrex::Abort("FillCoarsePatch: how did this happen?");
  physbc.SetLevel(lev);
  Interpolater* mapper = &cell_cons_interp;
    
  Array<BCRec> bcs(ncomp, physbc.GetBCRec());
  amrex::InterpFromCoarseLevel(mf, time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
			       physbc, physbc,
			       refRatio(lev-1),
			       mapper, bcs);
}
 
void // CUSTOM METHOD - CHANGEABLE
GeneralAMRIntegrator::GetData (const int /*lev*/,
			       const Real /*time*/,
			       Array<MultiFab*>& /*data*/,
			       Array<Real>& /*datatime*/)
{
  // data.clear();
  // datatime.clear();

  // TODO - Fix this
  // const Real teps = (t_new[lev] - t_old[lev]) * 1.e-3;

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
  // {
  //   for (int n = 0; n < number_of_fabs; n++)
  // 	{
  // 	  data.push_back((*fab_array[n])[lev].get());
  // 	}
  //   datatime.push_back(t_old[lev]);
  //   datatime.push_back(t_new[lev]);
  // }
  amrex::Abort("GeneralAMRIntegrator::GetData: Used for time interpolation, which is not currently supported!");
}

void
GeneralAMRIntegrator::ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow)
{
  TagCellsForRefinement(lev,tags,time,ngrow);
}


void
GeneralAMRIntegrator::InitData ()
{
    if (restart_chkfile.empty())
    {
	const Real time = 0.0;
	InitFromScratch(time);
	for (int lev = finest_level-1; lev >= 0; --lev)
	  {
	    for (int n = 0; n < number_of_fabs; n++)
	      amrex::average_down(*(*fab_array[n])[lev+1], *(*fab_array[n])[lev],
				geom[lev+1], geom[lev],
				  0, (*fab_array[n])[lev]->nComp(), refRatio(lev));
	  }

	if (plot_int > 0) {
	    WritePlotFile();
	}
    }
    else
    {
	InitFromCheckpoint();
    }
}

void GeneralAMRIntegrator::MakeNewLevelFromScratch (int lev, Real t, const BoxArray& ba,
						    const DistributionMapping& dm)
{
  for (int n = 0 ; n < number_of_fabs; n++)
    (*fab_array[n])[lev].reset(new MultiFab(ba, dm, ncomp_array[n], nghost_array[n]));
  
  t_new[lev] = t;
  t_old[lev] = t - 1.e200;

  Initialize(lev);
  
  for (int n = 0 ; n < number_of_fabs; n++)
    {
      physbc_array[n]->SetLevel(lev);
      physbc_array[n]->FillBoundary(*(*fab_array[n])[lev],0,0,t);
    }
}

std::vector<std::string>
GeneralAMRIntegrator::PlotFileName (int lev) const
{
  std::vector<std::string> name;
  name.push_back(plot_file+"/");
  name.push_back(amrex::Concatenate("", lev, 5));
  return name;
}

#define STR(a) #a
void
GeneralAMRIntegrator::WriteMetaData() const
{
  std::ofstream metadatafile;
  metadatafile.open(plot_file+"/metadata",std::ios_base::out);
  metadatafile << "COMPILATION DETAILS" << std::endl;
  metadatafile << "===================" << std::endl;
  metadatafile << "Git commit hash:         " << METADATA_GITHASH << std::endl;
  metadatafile << "AMReX version:           " << amrex::Version() << std::endl;
  metadatafile << "Dimension:               " << AMREX_SPACEDIM << std::endl;
  metadatafile << "Compiled by:             " << METADATA_USER  << std::endl;
  metadatafile << "Platform:                " << METADATA_PLATFORM  << std::endl;
  metadatafile << "Compiler:                " << METADATA_COMPILER  << std::endl;
  metadatafile << "Compilation Date:        " << METADATA_DATE  << std::endl;
  metadatafile << "Compilation Time:        " << METADATA_TIME  << std::endl;
  metadatafile << std::endl;

  auto starttime = std::chrono::system_clock::now();
  std::time_t now = std::chrono::system_clock::to_time_t(starttime);
  metadatafile << "RUN DETAILS" << std::endl;
  metadatafile << "===========" << std::endl;
  metadatafile << "Simulation started:      " << std::ctime(&now);
  metadatafile << "Number of processors:    " << amrex::ParallelDescriptor::NProcs() << std::endl;
  metadatafile << std::endl;

  metadatafile << "PARAMETERS" << std::endl;
  metadatafile << "==========" << std::endl;
  amrex::ParmParse pp;
  pp.dumpTable(metadatafile,true);
  
  metadatafile.close();
}


void
GeneralAMRIntegrator::WritePlotFile () const
{
  const int nlevels = finest_level+1;

  int total_components = 0;
  for (int i = 0; i < number_of_fabs; i++)
    total_components += ncomp_array[i];

  amrex::Vector<MultiFab> plotmf(nlevels);
  
  for (int ilev = 0; ilev < nlevels; ++ilev)
    {
      plotmf[ilev].define(grids[ilev], dmap[ilev], total_components, 0);

      int n = 0;
      for (int i = 0; i < number_of_fabs; i++)
	{
	  MultiFab::Copy(plotmf[ilev], *(*fab_array[i])[ilev], 0, n, ncomp_array[i], 0);
	  n += ncomp_array[i];
	}
    }

  const std::vector<std::string>& plotfilename = PlotFileName(istep[0]);
  
  WriteMultiLevelPlotfile(plotfilename[0]+plotfilename[1], nlevels, amrex::GetVecOfConstPtrs(plotmf), name_array,
			  Geom(), t_new[0],istep, refRatio());


  if (ParallelDescriptor::IOProcessor())
    {
      std::ofstream outfile;
      if (istep[0]==0) outfile.open(plot_file+"/output.visit",std::ios_base::out);
      else outfile.open(plot_file+"/output.visit",std::ios_base::app);
      outfile << plotfilename[1] + "/Header" << std::endl;
    }
}

void
GeneralAMRIntegrator::InitFromCheckpoint ()
{
  amrex::Abort("GeneralAMRIntegrator::InitFromCheckpoint: todo");
}

void
GeneralAMRIntegrator::Evolve ()
{
  Real cur_time = t_new[0];
  int last_plot_file_step = 0;

  for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
    {
      if (ParallelDescriptor::IOProcessor()) {
	std::cout << "\nSTEP " << step+1 << " starts ..." << std::endl;
      }
      int lev = 0;
      int iteration = 1;
      TimeStep(lev, cur_time, iteration);
      cur_time += dt[0];

      if (ParallelDescriptor::IOProcessor()) {
	std::cout << "STEP " << step+1 << " ends."
		  << " TIME = " << cur_time << " DT = " << dt[0]
		  << std::endl;
      }

      // sync up time
      for (int lev = 0; lev <= finest_level; ++lev) {
	t_new[lev] = cur_time;
      }

      if (plot_int > 0 && (step+1) % plot_int == 0) {
	last_plot_file_step = step+1;
	WritePlotFile();
      }

      if (cur_time >= stop_time - 1.e-6*dt[0]) break;
    }

  if (plot_int > 0 && istep[0] > last_plot_file_step) {
    WritePlotFile();
  }
}

void
GeneralAMRIntegrator::TimeStep (int lev, Real time, int /*iteration*/)
{

  if (regrid_int > 0)  // We may need to regrid
    {
      static Array<int> last_regrid_step(max_level+1, 0);

      // regrid doesn't change the base level, so we don't regrid on max_level
      if (lev < max_level && istep[lev] > last_regrid_step[lev])
	{
          if (istep[lev] % regrid_int == 0)
	    {
	      int old_finest = finest_level; // regrid changes finest_level
	      regrid(lev, time, false); 
	      for (int k = lev; k <= finest_level; ++k) {
		last_regrid_step[k] = istep[k];
	      }
	      for (int k = old_finest+1; k <= finest_level; ++k) {
		dt[k] = dt[k-1] / MaxRefRatio(k-1);
	      }
  	    }
  	}
    }

  if (Verbose() && ParallelDescriptor::IOProcessor()) {
    std::cout << "[Level " << lev 
	      << " step " << istep[lev]+1 << "] ";
    std::cout << "ADVANCE with dt = "
	      << dt[lev]
	      << std::endl;
  }

  for (int n = 0 ; n < number_of_fabs ; n++)
    FillPatch(lev,t_old[lev],*fab_array[n],*(*fab_array[n])[lev],*physbc_array[n],0);
  Advance(lev, time, dt[lev]);
  ++istep[lev];

  if (Verbose() && ParallelDescriptor::IOProcessor())
    {
      std::cout << "[Level " << lev
		<< " step " << istep[lev] << "] ";
      std::cout << "Advanced "
		<< CountCells(lev)
		<< " cells"
		<< std::endl;
    }

  if (lev < finest_level)
    {
      for (int i = 1; i <= nsubsteps[lev+1]; ++i)
	TimeStep(lev+1, time+(i-1)*dt[lev+1], i);

      for (int n = 0; n < number_of_fabs; n++)
	{
	  amrex::average_down(*(*fab_array[n])[lev+1], *(*fab_array[n])[lev],
			      geom[lev+1], geom[lev],
			      0, (*fab_array[n])[lev]->nComp(), refRatio(lev));
	}
    }
}
