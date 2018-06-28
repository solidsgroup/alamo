///
/// \file Integrator.cpp
/// \brief Compute the volume integral of two multiplied Fourier series
///

#include "Integrator.H"
#include "IO/FileNameParse.H"
#include "Util/Util.H"
#include <numeric>


using namespace amrex;

namespace Integrator
{

Integrator::Integrator ()
{
	{
		ParmParse pp;   // Basic run parameters
		pp.query("max_step", max_step);
		pp.query("stop_time", stop_time);
		pp.query("timestep",timestep);
	}
	{
		ParmParse pp("amr"); // AMR specific parameters
		pp.query("regrid_int", regrid_int);     // ALL processors
		//pp.query("check_file", check_file);   // Not currently used
		//pp.query("check_int", check_int);     // Not currently used
		pp.query("plot_int", plot_int);         // ALL processors
		pp.query("plot_file", plot_file);       // IO Processor only
		//pp.query("restart", restart_chkfile); // Not currently used
		pp.query("plot_file", plot_file);       // IO Processor only

		IO::FileNameParse(plot_file);

		nsubsteps.resize(maxLevel()+1,1);
		int cnt = pp.countval("nsubsteps");
		if (cnt != 0)
			if (cnt == maxLevel()){
				pp.queryarr("nsubsteps",nsubsteps);
				nsubsteps.insert(nsubsteps.begin(),1);
				nsubsteps.pop_back();
			}
			else if (cnt == 1)
			{
				int nsubsteps_all;
				pp.query("nsubsteps",nsubsteps_all);
				for (int lev = 1; lev <= maxLevel(); ++lev) nsubsteps[lev] = nsubsteps_all;
			}
			else
				Util::Abort("number of nsubsteps input must equal either 1 or amr.max_level");
		else
			for (int lev = 1; lev <= maxLevel(); ++lev) 
				nsubsteps[lev] = MaxRefRatio(lev-1);
	}

	int nlevs_max = maxLevel() + 1;

	istep.resize(nlevs_max, 0);

	t_new.resize(nlevs_max, 0.0);
	t_old.resize(nlevs_max, -1.e100);
	dt.resize(nlevs_max, 1.e100);
	dt[0] = timestep;
	for (int i = 1; i < nlevs_max; i++)
		dt[i] = dt[i-1] / (amrex::Real)nsubsteps[i];

	plot_file = Util::GetFileName();
	IO::WriteMetaData(plot_file);
}

///
/// \func  ~Integrator
/// \brief Does nothing -- check here first if there are memory leaks
///
Integrator::~Integrator ()
{
	if (ParallelDescriptor::IOProcessor())
		IO::WriteMetaData(plot_file);
}

/// \fn    Integrator::MakeNewLevelFromCoarse
/// \brief Wrapper to call FillCoarsePatch
/// \note **THIS OVERRIDES A PURE VIRTUAL METHOD - DO NOT CHANGE**
///
void
Integrator::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba, const DistributionMapping& dm)
{
	for (int n = 0; n < number_of_fabs; n++)
	{
		const int ncomp = (*fab_array[n])[lev-1]->nComp();
		const int nghost = (*fab_array[n])[lev-1]->nGrow();

		(*fab_array[n])[lev].reset(new MultiFab(ba, dm, ncomp, nghost));

		t_new[lev] = time;
		t_old[lev] = time - 1.e200;

		if (physbc_array[n]) // only do this if there is a BC object
			FillCoarsePatch(lev, time, *fab_array[n], *physbc_array[n], 0, ncomp);

	}
}


///
/// RESETS ALL MULTIFABS AT A GIVEN LEVEL
///
/// (OVERRIDES PURE VIRTUAL METHOD - DO NOT CHANGE)
///
void
Integrator::RemakeLevel (int lev,       ///<[in] AMR Level
			 Real time,     ///<[in] Simulation time
			 const BoxArray& ba, 
			 const DistributionMapping& dm)
{
	for (int n=0; n < number_of_fabs; n++)
	{

		const int ncomp = (*fab_array[n])[lev]->nComp();
		const int nghost = (*fab_array[n])[lev]->nGrow();

		MultiFab new_state(ba, dm, ncomp, nghost); 

		if (physbc_array[n]) 
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
Integrator::ClearLevel (int lev)
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
Integrator::RegisterNewFab(amrex::Vector<std::unique_ptr<amrex::MultiFab> > &new_fab,
			   BC::BC &new_bc,
			   int ncomp,
			   int nghost,
			   std::string name)
{
	int nlevs_max = maxLevel() + 1;
	new_fab.resize(nlevs_max); 
	fab_array.push_back(&new_fab);
	physbc_array.push_back(&new_bc); 
	ncomp_array.push_back(ncomp);
	nghost_array.push_back(nghost);
	name_array.push_back(name);
	number_of_fabs++;
}

void // CUSTOM METHOD - CHANGEABLE
Integrator::RegisterNewFab(amrex::Vector<std::unique_ptr<amrex::MultiFab> > &new_fab,
			   int ncomp,
			   std::string name)
{
	int nlevs_max = maxLevel() + 1;
	new_fab.resize(nlevs_max); 
	fab_array.push_back(&new_fab);
	physbc_array.push_back(NULL); 
	ncomp_array.push_back(ncomp);
	nghost_array.push_back(0);
	name_array.push_back(name);
	number_of_fabs++;
}


long // CUSTOM METHOD - CHANGEABLE
Integrator::CountCells (int lev)
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
Integrator::FillPatch (int lev, Real time,
		       Vector<std::unique_ptr<MultiFab> > &source_mf,
		       MultiFab &destination_mf,
		       BC::BC &physbc, int icomp)
{
	if (lev == 0)
	{
      
		Vector<MultiFab*> smf;
		smf.push_back(source_mf[lev].get());
		Vector<Real> stime;
		stime.push_back(time);

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
		amrex::Vector<MultiFab*> cmf, fmf;
		cmf.push_back(source_mf[lev-1].get());
		fmf.push_back(source_mf[lev].get());
		amrex::Vector<Real> ctime, ftime;
		ctime.push_back(time);
		ftime.push_back(time);

		physbc.SetLevel(lev);
		Interpolater* mapper = &cell_cons_interp;

		amrex::Vector<BCRec> bcs(destination_mf.nComp(), physbc.GetBCRec()); // todo
		amrex::ParallelDescriptor::Barrier();
      
		amrex::ParallelDescriptor::Barrier();
		amrex::FillPatchTwoLevels(destination_mf, time, cmf, ctime, fmf, ftime,
					  0, icomp, destination_mf.nComp(), geom[lev-1], geom[lev],
					  physbc, physbc,
					  refRatio(lev-1),
					  mapper, bcs);
		amrex::ParallelDescriptor::Barrier();
	}
}

/// \fn    Integrator::FillCoarsePatch
/// \brief Fill a fab at current level with the data from one level up
///
/// \note  This is a custom method and is changeable
void
Integrator::FillCoarsePatch (int lev, ///<[in] AMR level
			     Real time, ///<[in] Simulatinon time
			     amrex::Vector<std::unique_ptr<MultiFab> > &mf, ///<[in] Fab to fill
			     BC::BC &physbc, ///<[in] BC object applying to Fab
			     int icomp, ///<[in] start component
			     int ncomp) ///<[in] end component (i.e. applies to components `icomp`...`ncomp`)
{
	AMREX_ASSERT(lev > 0);

	amrex::Vector<amrex::MultiFab* > cmf;
	cmf.push_back(mf[lev-1].get());
	amrex::Vector<amrex::Real> ctime;
	ctime.push_back(time);
  
	physbc.SetLevel(lev);
	Interpolater* mapper = &cell_cons_interp;
    
	amrex::Vector<BCRec> bcs(ncomp, physbc.GetBCRec());
	amrex::InterpFromCoarseLevel(*mf[lev], time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
				     physbc, physbc,
				     refRatio(lev-1),
				     mapper, bcs);
}
 
void
Integrator::ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow)
{
	TagCellsForRefinement(lev,tags,time,ngrow);
}


void
Integrator::InitData ()
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

void
Integrator::MakeNewLevelFromScratch (int lev, Real t, const BoxArray& ba,
				     const DistributionMapping& dm)
{
	for (int n = 0 ; n < number_of_fabs; n++)
	{
		(*fab_array[n])[lev].reset(new MultiFab(ba, dm, ncomp_array[n], nghost_array[n]));
	}

	t_new[lev] = t;
	t_old[lev] = t - 1.e200;

	Initialize(lev);
  
	for (int n = 0 ; n < number_of_fabs; n++)
	{
		if (physbc_array[n]) // NO PROBLEM
		{
			physbc_array[n]->SetLevel(lev);
			physbc_array[n]->FillBoundary(*(*fab_array[n])[lev],0,0,t);
		}

	}
}

std::vector<std::string>
Integrator::PlotFileName (int lev) const
{
	std::vector<std::string> name;
	name.push_back(plot_file+"/");
	name.push_back(amrex::Concatenate("", lev, 5));
	return name;
}

void
Integrator::WritePlotFile () const
{
	const int nlevels = finest_level+1;

	int total_components = 0;
	amrex::Vector<std::string> complete_name_array;
	for (int i = 0; i < number_of_fabs; i++)
	{
		total_components += ncomp_array[i];
		if (ncomp_array[i] > 1)
			for (int j = 1; j <= ncomp_array[i]; j++)
				complete_name_array.push_back(amrex::Concatenate(name_array[i], j, 3));
		else
			complete_name_array.push_back(name_array[i]);
	}

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
  
	WriteMultiLevelPlotfile(plotfilename[0]+plotfilename[1], nlevels, amrex::GetVecOfConstPtrs(plotmf), complete_name_array,
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
Integrator::InitFromCheckpoint ()
{
	Util::Abort("Integrator::InitFromCheckpoint: todo");
}

void
Integrator::Evolve ()
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
		TimeStepBegin(cur_time,step);
		TimeStep(lev, cur_time, iteration);
		TimeStepComplete(cur_time,step);
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
Integrator::TimeStep (int lev, Real time, int /*iteration*/)
{
	if (regrid_int > 0)  // We may need to regrid
	{
		static Vector<int> last_regrid_step(max_level+1, 0);

		// regrid doesn't change the base level, so we don't regrid on max_level
		if (lev < max_level && istep[lev] > last_regrid_step[lev])
		{
			if (istep[lev] % regrid_int == 0)
			{
				/// \todo Delete this section (except for the regrid call) once it is verified that it is no longer needed
				//int old_finest = finest_level; // regrid changes finest_level
				regrid(lev, time, false); 
				// for (int k = lev; k <= finest_level; ++k) {
				// 	last_regrid_step[k] = istep[k];
				// }
				// for (int k = old_finest+1; k <= finest_level; ++k) {
				// 	dt[k] = dt[k-1] / MaxRefRatio(k-1);
				// }
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
}
