///
/// \file Integrator.cpp
/// \brief Compute the volume integral of two multiplied Fourier series
///

#include "Integrator.H"
#include "IO/FileNameParse.H"
#include "Util/Util.H"
#include "BC/Nothing.H"
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
				Util::Abort(INFO, "number of nsubsteps input must equal either 1 or amr.max_level");
		else
			for (int lev = 1; lev <= maxLevel(); ++lev) 
				nsubsteps[lev] = MaxRefRatio(lev-1);
	}
	{
		ParmParse pp("amr.thermo"); // AMR specific parameters
		pp.query("intvar_int", intvar_int);     // ALL processors
		pp.query("intvar_plot", intvar_plot);         // ALL processors
	}


	int nlevs_max = maxLevel() + 1;

	istep.resize(nlevs_max, 0);

	t_new.resize(nlevs_max, 0.0);
	t_old.resize(nlevs_max, -1.e100);
	SetTimestep(timestep);

	plot_file = Util::GetFileName();
	IO::WriteMetaData(plot_file,IO::Status::Running,0);

}

///
/// \func  ~Integrator
/// \brief Does nothing -- check here first if there are memory leaks
///
Integrator::~Integrator ()
{
	if (ParallelDescriptor::IOProcessor())
	  IO::WriteMetaData(plot_file,IO::Status::Complete);
}

void Integrator::SetTimestep(Set::Scalar _timestep)
{
	int nlevs_max = maxLevel() + 1;
	timestep = _timestep;
	dt.resize(nlevs_max, 1.e100);
	dt[0] = timestep;
	for (int i = 1; i < nlevs_max; i++)
		dt[i] = dt[i-1] / (amrex::Real)nsubsteps[i];
}

/// \fn    Integrator::MakeNewLevelFromCoarse
/// \brief Wrapper to call FillCoarsePatch
/// \note **THIS OVERRIDES A PURE VIRTUAL METHOD - DO NOT CHANGE**
///
void
Integrator::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& cgrids, const DistributionMapping& dm)
{
	t_new[lev] = time;
	t_old[lev] = time - 1.e200;

	for (int n = 0; n < cell.number_of_fabs; n++)
	{
		const int ncomp  = (*cell.fab_array[n])[lev-1]->nComp();
		const int nghost = (*cell.fab_array[n])[lev-1]->nGrow();

		(*cell.fab_array[n])[lev].reset(new MultiFab(cgrids, dm, ncomp, nghost));
		
		(*cell.fab_array[n])[lev]->setVal(0.0);

		FillCoarsePatch(lev, time, *cell.fab_array[n], *cell.physbc_array[n], 0, ncomp);
	}

	amrex::BoxArray ngrids = cgrids;
	ngrids.convert(amrex::IntVect::TheNodeVector());

	for (int n = 0; n < node.number_of_fabs; n++)
	{
		const int ncomp = (*node.fab_array[n])[lev-1]->nComp();
		const int nghost = (*node.fab_array[n])[lev-1]->nGrow();

		(*node.fab_array[n])[lev].reset(new MultiFab(ngrids, dm, ncomp, nghost));
		(*node.fab_array[n])[lev]->setVal(0.0);

		FillCoarsePatch(lev, time, *node.fab_array[n], *node.physbc_array[n], 0, ncomp);
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
			 const BoxArray& cgrids, 
			 const DistributionMapping& dm)
{
	for (int n=0; n < cell.number_of_fabs; n++)
	{
		const int ncomp  = (*cell.fab_array[n])[lev]->nComp();
		const int nghost = (*cell.fab_array[n])[lev]->nGrow();

		MultiFab new_state(cgrids, dm, ncomp, nghost); 

		FillPatch(lev, time, *cell.fab_array[n], new_state, *cell.physbc_array[n], 0);
		std::swap(new_state, *(*cell.fab_array[n])[lev]);
	}

	amrex::BoxArray ngrids = cgrids;
	ngrids.convert(amrex::IntVect::TheNodeVector());

	for (int n=0; n < node.number_of_fabs; n++)
	{
		const int ncomp  = (*node.fab_array[n])[lev]->nComp();
		const int nghost = (*node.fab_array[n])[lev]->nGrow();

		MultiFab new_state(ngrids, dm, ncomp, nghost); 

		FillPatch(lev, time, *node.fab_array[n], new_state, *node.physbc_array[n], 0);
		std::swap(new_state, *(*node.fab_array[n])[lev]);
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
	for (int n = 0; n < cell.number_of_fabs; n++)
	{
		(*cell.fab_array[n])[lev].reset(nullptr);
	}
	for (int n = 0; n < node.number_of_fabs; n++)
	{
		(*node.fab_array[n])[lev].reset(nullptr);
	}
}

//
//
//

void // CUSTOM METHOD - CHANGEABLE
Integrator::RegisterNewFab(amrex::Vector<std::unique_ptr<amrex::MultiFab> > &new_fab,
			   BC::BC *new_bc,
			   int ncomp,
			   int nghost,
			   std::string name)
{
	int nlevs_max = maxLevel() + 1;
	new_fab.resize(nlevs_max); 
	cell.fab_array.push_back(&new_fab);
	cell.physbc_array.push_back(new_bc); 
	cell.ncomp_array.push_back(ncomp);
	cell.nghost_array.push_back(nghost);
	cell.name_array.push_back(name);
	cell.number_of_fabs++;
}

void // CUSTOM METHOD - CHANGEABLE
Integrator::RegisterNewFab(amrex::Vector<std::unique_ptr<amrex::MultiFab> > &new_fab,
			   int ncomp,
			   std::string name)
{
	int nlevs_max = maxLevel() + 1;
	new_fab.resize(nlevs_max); 
	cell.fab_array.push_back(&new_fab);
	cell.physbc_array.push_back(new BC::Nothing); 
	cell.ncomp_array.push_back(ncomp);
	cell.nghost_array.push_back(0);
	cell.name_array.push_back(name);
	cell.number_of_fabs++;
}
void // CUSTOM METHOD - CHANGEABLE
Integrator::RegisterNodalFab(amrex::Vector<std::unique_ptr<amrex::MultiFab> > &new_fab,
			     int ncomp,
			     int nghost,
			     std::string name)
{
	int nlevs_max = maxLevel() + 1;
	new_fab.resize(nlevs_max); 
	node.fab_array.push_back(&new_fab);
	node.physbc_array.push_back(new BC::Nothing); 
	node.ncomp_array.push_back(ncomp);
	node.nghost_array.push_back(nghost);
	node.name_array.push_back(name);
	node.number_of_fabs++;
}


void // CUSTOM METHOD - CHANGEABLE
Integrator::RegisterIntegratedVariable(Set::Scalar *integrated_variable, std::string name)
{
	intvar_array.push_back(integrated_variable);
	intvar_names.push_back(name);
	Util::Message(INFO,intvar_array[0]);
	number_of_intvars++;
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

		physbc.define(geom[lev]);
		amrex::FillPatchSingleLevel(destination_mf,		// Multifab
					    time,                         // time
					    smf,				// Vector<MultiFab*> &smf (CONST)
					    stime,			// Vector<Real> &stime    (CONST)
					    0,				// scomp - Source component 
					    icomp,			// dcomp - Destination component
					    destination_mf.nComp(),	// ncomp - Number of components
					    geom[lev],			// Geometry (CONST)
					    physbc, 0);			// BC
	} 
	else
	{
		amrex::Vector<MultiFab*> cmf, fmf;
		cmf.push_back(source_mf[lev-1].get());
		fmf.push_back(source_mf[lev].get());
		amrex::Vector<Real> ctime, ftime;
		ctime.push_back(time);
		ftime.push_back(time);

		physbc.define(geom[lev]);
		Interpolater* mapper = &cell_cons_interp;

		amrex::Vector<BCRec> bcs(destination_mf.nComp(), physbc.GetBCRec()); // todo
		amrex::ParallelDescriptor::Barrier();
      
		amrex::ParallelDescriptor::Barrier();
		amrex::FillPatchTwoLevels(destination_mf, time, cmf, ctime, fmf, ftime,
					  0, icomp, destination_mf.nComp(), geom[lev-1], geom[lev],
					  physbc, 0,
					  physbc, 0,
					  refRatio(lev-1),
					  mapper, bcs, 0);
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
  
	physbc.define(geom[lev]);
	Interpolater* mapper = &cell_cons_interp;
    
	amrex::Vector<BCRec> bcs(ncomp, physbc.GetBCRec());
	amrex::InterpFromCoarseLevel(*mf[lev], time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
				     physbc, 0,
				     physbc, 0,
				     refRatio(lev-1),
				     mapper, bcs, 0);
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
		for (int n = 0; n < cell.number_of_fabs; n++)
			amrex::average_down(*(*cell.fab_array[n])[lev+1], *(*cell.fab_array[n])[lev],
					    geom[lev+1], geom[lev],
					    0, (*cell.fab_array[n])[lev]->nComp(), refRatio(lev));
		Util::Warning(INFO,"Not averaging down nodal fabs");
		// for (int n = 0; n < node.number_of_fabs; n++)
		// 	amrex::average_down_nodal(*(*node.fab_array[n])[lev+1], *(*node.fab_array[n])[lev], refRatio(lev));
	}
  
	if (plot_int > 0) {
		WritePlotFile();
	}
}

void
Integrator::MakeNewLevelFromScratch (int lev, Real t, const BoxArray& cgrids,
				     const DistributionMapping& dm)
{
	for (int n = 0 ; n < cell.number_of_fabs; n++)
	{
		(*cell.fab_array[n])[lev].reset(new MultiFab(cgrids, dm, cell.ncomp_array[n], cell.nghost_array[n]));
		(*cell.fab_array[n])[lev]->setVal(0.0);
	}

	amrex::BoxArray ngrids = cgrids;
	ngrids.convert(amrex::IntVect::TheNodeVector());
	for (int n = 0 ; n < node.number_of_fabs; n++)
	{
		(*node.fab_array[n])[lev].reset(new MultiFab(ngrids, dm, node.ncomp_array[n], node.nghost_array[n]));
		(*node.fab_array[n])[lev]->setVal(0.0);
	}

	t_new[lev] = t;
	t_old[lev] = t - 1.e200;

	Initialize(lev);
  
	for (int n = 0 ; n < cell.number_of_fabs; n++)
	{
		cell.physbc_array[n]->define(geom[lev]);
		cell.physbc_array[n]->FillBoundary(*(*cell.fab_array[n])[lev],0,0,t,0);
	}

	for (int n = 0 ; n < node.number_of_fabs; n++)
	{
		node.physbc_array[n]->define(geom[lev]);
		node.physbc_array[n]->FillBoundary(*(*node.fab_array[n])[lev],0,0,t,0);
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

	int ccomponents = 0, ncomponents = 0;
	amrex::Vector<std::string> cnames, nnames;
	for (int i = 0; i < cell.number_of_fabs; i++)
	{
		ccomponents += cell.ncomp_array[i];
		if (cell.ncomp_array[i] > 1)
			for (int j = 1; j <= cell.ncomp_array[i]; j++)
				cnames.push_back(amrex::Concatenate(cell.name_array[i], j, 3));
		else
			cnames.push_back(cell.name_array[i]);
	}
	for (int i = 0; i < node.number_of_fabs; i++)
	{
		ncomponents += node.ncomp_array[i];
		if (node.ncomp_array[i] > 1)
			for (int j = 1; j <= node.ncomp_array[i]; j++)
				nnames.push_back(amrex::Concatenate(node.name_array[i], j, 3));
		else
			nnames.push_back(node.name_array[i]);
	}

	amrex::Vector<MultiFab> cplotmf(nlevels), nplotmf(nlevels);
  
	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		cplotmf[ilev].define(grids[ilev], dmap[ilev], ccomponents, 0);
		amrex::BoxArray ngrids = grids[ilev];
		ngrids.convert(amrex::IntVect::TheNodeVector());
		if (ncomponents>0) nplotmf[ilev].define(ngrids, dmap[ilev], ncomponents, 0);

		int n = 0;
		for (int i = 0; i < cell.number_of_fabs; i++)
		{
			MultiFab::Copy(cplotmf[ilev], *(*cell.fab_array[i])[ilev], 0, n, cell.ncomp_array[i], 0);
			n += cell.ncomp_array[i];
		}

		n = 0;
		for (int i = 0; i < node.number_of_fabs; i++)
		{
			MultiFab::Copy(nplotmf[ilev], *(*node.fab_array[i])[ilev], 0, n, node.ncomp_array[i], 0);
			n += node.ncomp_array[i];
		}
	}

	const std::vector<std::string>& plotfilename = PlotFileName(istep[0]);
  
	WriteMultiLevelPlotfile(plotfilename[0]+plotfilename[1]+"cell", nlevels, amrex::GetVecOfConstPtrs(cplotmf), cnames,
				Geom(), t_new[0],istep, refRatio());

	if (ncomponents > 0)
		WriteMultiLevelPlotfile(plotfilename[0]+plotfilename[1]+"node", nlevels, amrex::GetVecOfConstPtrs(nplotmf), nnames,
					Geom(), t_new[0],istep, refRatio());

	if (ParallelDescriptor::IOProcessor())
	{
		std::ofstream coutfile, noutfile;
		if (istep[0]==0)
		{
			coutfile.open(plot_file+"/celloutput.visit",std::ios_base::out);
			if (ncomponents > 0) noutfile.open(plot_file+"/nodeoutput.visit",std::ios_base::out);
		}
		else
		{
			coutfile.open(plot_file+"/celloutput.visit",std::ios_base::app);
			if (ncomponents > 0) noutfile.open(plot_file+"/nodeoutput.visit",std::ios_base::app);
		}
		coutfile << plotfilename[1] + "cell" + "/Header" << std::endl;
		if (ncomponents > 0) noutfile << plotfilename[1] + "node" + "/Header" << std::endl;
	}
}

void
Integrator::InitFromCheckpoint ()
{
	Util::Abort(INFO, "Integrator::InitFromCheckpoint: todo");
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
		IntegrateVariables(cur_time,step);
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
			IO::WriteMetaData(plot_file,IO::Status::Running,(int)(100.0*cur_time/stop_time));
		}

		if (cur_time >= stop_time - 1.e-6*dt[0]) break;
	}
	if (plot_int > 0 && istep[0] > last_plot_file_step) {
		WritePlotFile();
	}
}

void
Integrator::IntegrateVariables (Real time, int step)
{
	if (!number_of_intvars) return;

	if (!(step % intvar_int))
	{
		// Zero out all variables
		for (int i = 0; i < number_of_intvars; i++) *intvar_array[i] = 0; 

		// All levels except the finest
		for (int ilev = 0; ilev < max_level; ilev++)
		{
			const BoxArray& cfba = amrex::coarsen(grids[ilev+1], refRatio(ilev));

			for ( amrex::MFIter mfi(grids[ilev],dmap[ilev],true); mfi.isValid(); ++mfi )
			{
				const amrex::Box& box = mfi.tilebox();
				const::BoxArray & comp = amrex::complementIn(box,cfba);

				for (int i = 0; i < comp.size(); i++)
				{
					Integrate(ilev,time, step,
						  mfi, comp[i]);
				}
			}
		}
		// Now do the finest level
		{
			for ( amrex::MFIter mfi(grids[max_level],dmap[max_level],true); mfi.isValid(); ++mfi )
			{
				const amrex::Box& box = mfi.tilebox();
				Integrate(max_level, time, step, mfi, box);
			}
		}

		// Sum up across all processors
		for (int i = 0; i < number_of_intvars; i++) 
		{
			amrex::ParallelDescriptor::ReduceRealSum(*intvar_array[i]);
		}
	}
	if (!(step % intvar_plot) && 	ParallelDescriptor::IOProcessor())
	{
		std::ofstream outfile;
		if (step==0)
		{
			outfile.open(plot_file+"/thermo.dat",std::ios_base::out);
			outfile << "time";
			for (int i = 0; i < number_of_intvars; i++) 
				outfile << "\t" << intvar_names[i];
			outfile << std::endl;
		}
		else outfile.open(plot_file+"/thermo.dat",std::ios_base::app);
		outfile << time;
		for (int i = 0; i < number_of_intvars; i++)
			outfile << "\t" << *intvar_array[i];
		outfile << std::endl;
		outfile.close();
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

	for (int n = 0 ; n < cell.number_of_fabs ; n++)
		FillPatch(lev,t_old[lev],*cell.fab_array[n],*(*cell.fab_array[n])[lev],*cell.physbc_array[n],0);
	for (int n = 0 ; n < node.number_of_fabs ; n++)
		FillPatch(lev,t_old[lev],*node.fab_array[n],*(*node.fab_array[n])[lev],*node.physbc_array[n],0);

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

		for (int n = 0; n < cell.number_of_fabs; n++)
		{
			amrex::average_down(*(*cell.fab_array[n])[lev+1], *(*cell.fab_array[n])[lev],
					    geom[lev+1], geom[lev],
					    0, (*cell.fab_array[n])[lev]->nComp(), refRatio(lev));
		}
		// for (int n = 0; n < node.number_of_fabs; n++)
		// {
		// 	amrex::average_down(*(*node.fab_array[n])[lev+1], *(*node.fab_array[n])[lev],
		// 			    geom[lev+1], geom[lev],
		// 			    0, (*node.fab_array[n])[lev]->nComp(), refRatio(lev));
		// }
	}
}
}
