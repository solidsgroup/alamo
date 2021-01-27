///
/// \file Integrator.cpp
/// \brief Compute the volume integral of two multiplied Fourier series
///

#include "Integrator.H"
#include "IO/FileNameParse.H"
#include "Util/Util.H"
#include <numeric>



namespace Integrator
{

Integrator::Integrator ()
{
	BL_PROFILE("Integrator::Integrator()");
	{
		amrex::ParmParse pp;   // Basic run parameters
		pp.query("max_step", max_step);
		pp.query("stop_time", stop_time);
		pp.query("timestep",timestep);
		pp.query("restart", restart_file_cell); 
		pp.query("restart_cell", restart_file_cell); 
		pp.query("restart_node", restart_file_node); 
	}
	{
		amrex::ParmParse pp("amr"); // AMR specific parameters
		pp.query("regrid_int", regrid_int);     // ALL processors
		pp.query("plot_int", plot_int);         // ALL processors
		pp.query("plot_dt", plot_dt);         // ALL processors
		pp.query("plot_file", plot_file);       // IO Processor only
		
		pp.query("max_plot_level",max_plot_level);		

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
		amrex::ParmParse pp("amr.thermo"); // AMR specific parameters
		thermo.interval = 1; // Default: integrate every time.
		pp.query("int", thermo.interval);     // ALL processors
		pp.query("plot_int", thermo.plot_int);         // ALL processors
		pp.query("plot_dt", thermo.plot_dt);         // ALL processors
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
	BL_PROFILE("Integrator::~Integrator");
	if (amrex::ParallelDescriptor::IOProcessor())
	  IO::WriteMetaData(plot_file,IO::Status::Complete);
}

void Integrator::SetTimestep(Set::Scalar _timestep)
{
	BL_PROFILE("Integrator::SetTimestep");
	int nlevs_max = maxLevel() + 1;
	timestep = _timestep;
	dt.resize(nlevs_max, 1.e100);
	dt[0] = timestep;
	for (int i = 1; i < nlevs_max; i++)
		dt[i] = dt[i-1] / (amrex::Real)nsubsteps[i];
}
void Integrator::SetPlotInt(int a_plot_int)
{
	BL_PROFILE("Integrator::SetPlotInt");
	plot_int = a_plot_int;
}

/// \fn    Integrator::MakeNewLevelFromCoarse
/// \brief Wrapper to call FillCoarsePatch
/// \note **THIS OVERRIDES A PURE VIRTUAL METHOD - DO NOT CHANGE**
///
void
Integrator::MakeNewLevelFromCoarse (int lev, amrex::Real time, const amrex::BoxArray& cgrids, const amrex::DistributionMapping& dm)
{
	Util::Message(INFO);
	BL_PROFILE("Integrator::MakeNewLevelFromCoarse");

	for (int n = 0; n < cell.number_of_fabs; n++)
	{
		const int ncomp  = (*cell.fab_array[n])[lev-1]->nComp();
		const int nghost = (*cell.fab_array[n])[lev-1]->nGrow();

		(*cell.fab_array[n])[lev].reset(new amrex::MultiFab(cgrids, dm, ncomp, nghost));
		
		(*cell.fab_array[n])[lev]->setVal(0.0);

		FillCoarsePatch(lev, time, *cell.fab_array[n], *cell.physbc_array[n], 0, ncomp);
	}

	amrex::BoxArray ngrids = cgrids;
	ngrids.convert(amrex::IntVect::TheNodeVector());

	for (int n = 0; n < node.number_of_fabs; n++)
	{
		const int ncomp = (*node.fab_array[n])[lev-1]->nComp();
		const int nghost = (*node.fab_array[n])[lev-1]->nGrow();

		(*node.fab_array[n])[lev].reset(new amrex::MultiFab(ngrids, dm, ncomp, nghost));
		(*node.fab_array[n])[lev]->setVal(0.0);

		FillCoarsePatch(lev, time, *node.fab_array[n], *node.physbc_array[n], 0, ncomp);
	}

	for (unsigned int n = 0; n < m_basefields.size(); n++)
	{
		m_basefields[n]->MakeNewLevelFromCoarse(lev,time,cgrids,dm);
	}

}


///
/// RESETS ALL MULTIFABS AT A GIVEN LEVEL
///
/// (OVERRIDES PURE VIRTUAL METHOD - DO NOT CHANGE)
///
void
Integrator::RemakeLevel (int lev,       ///<[in] AMR Level
			 amrex::Real time,     ///<[in] Simulation time
			 const amrex::BoxArray& cgrids, 
			 const amrex::DistributionMapping& dm)
{
	BL_PROFILE("Integrator::RemakeLevel");
	for (int n=0; n < cell.number_of_fabs; n++)
	{
		const int ncomp  = (*cell.fab_array[n])[lev]->nComp();
		const int nghost = (*cell.fab_array[n])[lev]->nGrow();

		amrex::MultiFab new_state(cgrids, dm, ncomp, nghost); 

		new_state.setVal(0.0);
		FillPatch(lev, time, *cell.fab_array[n], new_state, *cell.physbc_array[n], 0);
		std::swap(new_state, *(*cell.fab_array[n])[lev]);
	}

	amrex::BoxArray ngrids = cgrids;
	ngrids.convert(amrex::IntVect::TheNodeVector());

	for (int n=0; n < node.number_of_fabs; n++)
	{
		const int ncomp  = (*node.fab_array[n])[lev]->nComp();
		const int nghost = (*node.fab_array[n])[lev]->nGrow();

		amrex::MultiFab new_state(ngrids, dm, ncomp, nghost); 

		new_state.setVal(0.0);
		FillPatch(lev, time, *node.fab_array[n], new_state, *node.physbc_array[n], 0);
		std::swap(new_state, *(*node.fab_array[n])[lev]);
	}

	for (unsigned int n = 0; n < m_basefields.size(); n++)
	{
		m_basefields[n]->RemakeLevel(lev,time,cgrids,dm);
	}
}

//
// DELETE EVERYTHING
//
// (OVERRIDES PURE VIRTUAL METHOD - DO NOT CHANGE)
//
void
Integrator::ClearLevel (int lev)
{
	BL_PROFILE("Integrator::ClearLevel");
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
Integrator::RegisterNewFab(Set::Field<Set::Scalar> &new_fab,
			   BC::BC<Set::Scalar> *new_bc,
			   int ncomp,
			   int nghost,
			   std::string name,
			   bool writeout)
{
	BL_PROFILE("Integrator::RegisterNewFab_1");
	int nlevs_max = maxLevel() + 1;
	new_fab.resize(nlevs_max); 
	cell.fab_array.push_back(&new_fab);
	cell.physbc_array.push_back(new_bc); 
	cell.ncomp_array.push_back(ncomp);
	cell.nghost_array.push_back(nghost);
	cell.name_array.push_back(name);
	cell.writeout_array.push_back(writeout);
	cell.number_of_fabs++;
}

void // CUSTOM METHOD - CHANGEABLE
Integrator::RegisterNewFab(Set::Field<Set::Scalar> &new_fab,
			   int ncomp,
			   std::string name,
			   bool writeout)
{
	BL_PROFILE("Integrator::RegisterNewFab_2");
	int nlevs_max = maxLevel() + 1;
	new_fab.resize(nlevs_max); 
	cell.fab_array.push_back(&new_fab);
	cell.physbc_array.push_back(&bcnothing); 
	cell.ncomp_array.push_back(ncomp);
	cell.nghost_array.push_back(0);
	cell.name_array.push_back(name);
	cell.writeout_array.push_back(writeout);
	cell.number_of_fabs++;
}
void // CUSTOM METHOD - CHANGEABLE
Integrator::RegisterNodalFab(Set::Field<Set::Scalar> &new_fab,
			   	 BC::BC<Set::Scalar> *new_bc,
			     int ncomp,
			     int nghost,
			     std::string name,
				 bool writeout)
{
	BL_PROFILE("Integrator::RegisterNodalFab");
	int nlevs_max = maxLevel() + 1;
	new_fab.resize(nlevs_max); 
	node.fab_array.push_back(&new_fab);
	node.physbc_array.push_back(new_bc); 
	node.ncomp_array.push_back(ncomp);
	node.nghost_array.push_back(nghost);
	node.name_array.push_back(name);
	node.writeout_array.push_back(writeout);
	node.number_of_fabs++;
}
void // CUSTOM METHOD - CHANGEABLE
Integrator::RegisterNodalFab(Set::Field<Set::Scalar> &new_fab,
			     int ncomp,
			     int nghost,
			     std::string name,
				 bool writeout)
{
	RegisterNodalFab(new_fab,&bcnothing,ncomp,nghost,name,writeout);
}


void // CUSTOM METHOD - CHANGEABLE
Integrator::RegisterIntegratedVariable(Set::Scalar *integrated_variable, std::string name)
{
	BL_PROFILE("Integrator::RegisterIntegratedVariable");
	thermo.vars.push_back(integrated_variable);
	thermo.names.push_back(name);
	thermo.number++;
}

long // CUSTOM METHOD - CHANGEABLE
Integrator::CountCells (int lev)
{
	BL_PROFILE("Integrator::CountCells");
	const int N = grids[lev].size();

	long cnt = 0;

	for (int i = 0; i < N; ++i)
	{
		cnt += grids[lev][i].numPts();
	}

	return cnt;
}

void  // CUSTOM METHOD - CHANGEABLE
Integrator::FillPatch (int lev, amrex::Real time,
		       amrex::Vector<std::unique_ptr<amrex::MultiFab>> &source_mf,
		       amrex::MultiFab &destination_mf,
		       BC::BC<Set::Scalar> &physbc, int icomp)
{
	BL_PROFILE("Integrator::FillPatch");
	if (lev == 0)
	{
      
		amrex::Vector<amrex::MultiFab*> smf;
		smf.push_back(source_mf[lev].get());
		amrex::Vector<amrex::Real> stime;
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
					    physbc,
						0);			// BC
	} 
	else
	{
		amrex::Vector<amrex::MultiFab*> cmf, fmf;
		cmf.push_back(source_mf[lev-1].get());
		fmf.push_back(source_mf[lev].get());
		amrex::Vector<amrex::Real> ctime, ftime;
		ctime.push_back(time);
		ftime.push_back(time);

		physbc.define(geom[lev]);
		
		amrex::Interpolater* mapper;

		if (destination_mf.boxArray().ixType() == amrex::IndexType::TheNodeType())
			mapper = &amrex::node_bilinear_interp;
		else
			mapper = &amrex::cell_cons_interp;

		amrex::Vector<amrex::BCRec> bcs(destination_mf.nComp(), physbc.GetBCRec()); // todo
		amrex::FillPatchTwoLevels(destination_mf, time, cmf, ctime, fmf, ftime,
					  0, icomp, destination_mf.nComp(), geom[lev-1], geom[lev],
					  physbc, 0,
					  physbc, 0,
					  refRatio(lev-1),
					  mapper, bcs, 0);
//		if (destination_mf.contains_nan()) Util::Abort(INFO);
	}
}

/// \fn    Integrator::FillCoarsePatch
/// \brief Fill a fab at current level with the data from one level up
///
/// \note  This is a custom method and is changeable
void
Integrator::FillCoarsePatch (int lev, ///<[in] AMR level
			     amrex::Real time, ///<[in] Simulatinon time
			     Set::Field<Set::Scalar> &mf, ///<[in] Fab to fill
			     BC::BC<Set::Scalar> &physbc, ///<[in] BC object applying to Fab
			     int icomp, ///<[in] start component
			     int ncomp) ///<[in] end component (i.e. applies to components `icomp`...`ncomp`)
{
	BL_PROFILE("Integrator::FillCoarsePatch");
	AMREX_ASSERT(lev > 0);
	amrex::Vector<amrex::MultiFab *> cmf;
	cmf.push_back(mf[lev-1].get());
	amrex::Vector<amrex::Real> ctime;
	ctime.push_back(time);
  
	physbc.define(geom[lev]);

	amrex::Interpolater* mapper;
	if (mf[lev]->boxArray().ixType() == amrex::IndexType::TheNodeType())
		mapper = &amrex::node_bilinear_interp;
	else
		mapper = &amrex::cell_cons_interp;
    
	amrex::Vector<amrex::BCRec> bcs(ncomp, physbc.GetBCRec());
	amrex::InterpFromCoarseLevel(*mf[lev], time, *cmf[0], 0, icomp, ncomp, geom[lev-1], geom[lev],
				     physbc, 0,
				     physbc, 0,
				     refRatio(lev-1),
				     mapper, bcs, 0);
}
 
void
Integrator::ErrorEst (int lev, amrex::TagBoxArray& tags, amrex::Real time, int ngrow)
{
	BL_PROFILE("Integrator::ErrorEst");
	TagCellsForRefinement(lev,tags,time,ngrow);
}


void
Integrator::InitData ()
{
	BL_PROFILE("Integrator::InitData");
	
	if (restart_file_cell == "" && restart_file_node == "")
	{
		const amrex::Real time = 0.0;
		InitFromScratch(time);

		for (int lev = finest_level-1; lev >= 0; --lev)
		{
			if (lev < max_level) regrid(lev,0.0);
			for (int n = 0; n < cell.number_of_fabs; n++)
				amrex::average_down(*(*cell.fab_array[n])[lev+1], *(*cell.fab_array[n])[lev],
						    geom[lev+1], geom[lev],
						    0, (*cell.fab_array[n])[lev]->nComp(), refRatio(lev));
		}
		SetFinestLevel(finest_level);
	}
	if (restart_file_cell != "")
	{
		Restart(restart_file_cell,false);
	}
	if (restart_file_node != "")
	{
		Restart(restart_file_node,true);
	}
	
	if (plot_int > 0 || plot_dt > 0.0) {
		WritePlotFile();
	}
}

void
Integrator::Restart(const std::string dirname, bool a_nodal)
{
	BL_PROFILE("Integrator::Restart");
	
	if ( a_nodal && node.fab_array.size() == 0) 
	{
		Util::Message(INFO,"Nothing here for nodal fabs");
		return;
	}
	if (!a_nodal && cell.fab_array.size() == 0) 
	{
		Util::Message(INFO,"Nothing here for cell-based fabs");
		return;
	}
	
	std::string filename = dirname + "/Header";
	std::string chkptfilename = dirname + "/Checkpoint";
	amrex::VisMF::IO_Buffer io_buffer(amrex::VisMF::GetIOBufferSize());
	amrex::Vector<char> fileCharPtr,chkptfileCharPtr;
	amrex::ParallelDescriptor::ReadAndBcastFile(filename,fileCharPtr);
	amrex::ParallelDescriptor::ReadAndBcastFile(chkptfilename,chkptfileCharPtr);
	std::string fileCharPtrString(fileCharPtr.dataPtr());
	std::string chkptfileCharPtrString(chkptfileCharPtr.dataPtr());
	std::istringstream is(fileCharPtrString,std::istringstream::in);
	std::istringstream chkpt_is(chkptfileCharPtrString,std::istringstream::in);

	std::string line, word;

	// Get version
	std::getline(is,line); 
	Util::Message(INFO,"Version: ", line);

	// Get number of fabs
	int tmp_numfabs;
	std::getline(is,line); tmp_numfabs = std::stoi(line);
	Util::Message(INFO,"number of fabs:",tmp_numfabs);
	std::vector<std::string> tmp_name_array;
	
	for (int i = 0; i < tmp_numfabs; i++)
	{
		std::getline(is,line);
		tmp_name_array.push_back(line);
	}
	
	// Dimension?
	std::getline(is,line); 
	Util::Warning(INFO,"Dimension: " + line);
	
	// Current time
	Set::Scalar tmp_time = 0.0;
	std::getline(is,line); tmp_time = std::stof(line); Util::Message(INFO,"Current time: ", tmp_time);
	for (int i = 0; i < max_level + 1; i++)
	{
		t_new[i] = tmp_time; t_old[i] = tmp_time;
	}
	
	// AMR level
	int tmp_max_level;
	std::getline(is,line); tmp_max_level = std::stoi(line); Util::Message(INFO,"Max AMR level: ", line);
	if (tmp_max_level != max_level)
		Util::Abort(INFO,"The max level specified (",max_level,") does not match the max level in the restart file (",tmp_max_level,")");
	finest_level = tmp_max_level;
	// Geometry ?
	std::getline(is,line); Util::Message(INFO,"Input geometry: ", line);
	std::getline(is,line); Util::Message(INFO,"                ", line);

	// Mesh refinement ratio?
	std::getline(is,line); Util::Message(INFO,"Mesh refinement ratio: ",line);
	
	// Domain
	std::getline(is,line); Util::Warning(INFO,"Domain: ",line);

	// Domain
	std::getline(is,line); 
	std::vector<std::string> tmp_iters = Util::String::Split(line);
	if (tmp_iters.size() != (unsigned int)(max_level+1)) Util::Abort(INFO, "Error reading in interation counts: line = ", line);
	for (int lev = 0; lev <= max_level; lev++) {istep[lev] = std::stoi(tmp_iters[lev]); Util::Message(INFO,"Iter on level " , lev , " = ", istep[lev]);}

	amrex::Vector<amrex::MultiFab> tmpdata(tmp_max_level+1);
	int total_ncomp = 0; 
	
	if (a_nodal) for (unsigned int i = 0; i < node.fab_array.size(); i++) total_ncomp += node.ncomp_array[i];
	else         for (unsigned int i = 0; i < cell.fab_array.size(); i++) total_ncomp += cell.ncomp_array[i];

	int total_nghost = a_nodal ? 0 : cell.nghost_array[0];

	for (int lev = 0; lev <= max_level; lev++)
	{
		amrex::BoxArray tmp_ba;
		tmp_ba.readFrom(chkpt_is);
		SetBoxArray(lev,tmp_ba);
		amrex::DistributionMapping tmp_dm(tmp_ba,amrex::ParallelDescriptor::NProcs());
		SetDistributionMap(lev,tmp_dm);

		if (a_nodal)
		{
			amrex::BoxArray ngrids = grids[lev];
			ngrids.convert(amrex::IntVect::TheNodeVector());
			tmpdata[lev].define(ngrids,dmap[lev],total_ncomp,total_nghost);
		}
		else
		{
			tmpdata[lev].define(grids[lev],dmap[lev],total_ncomp,total_nghost);
		}
		amrex::VisMF::Read( tmpdata[lev],
							amrex::MultiFabFileFullPrefix(lev,dirname,"Level_","Cell"));
							

		if (a_nodal)
			for (int i = 0; i < node.number_of_fabs; i++)
			{
				amrex::BoxArray ngrids = grids[lev];
				ngrids.convert(amrex::IntVect::TheNodeVector());
				(*node.fab_array[i])[lev].reset(new amrex::MultiFab(ngrids,dmap[lev],node.ncomp_array[i],node.nghost_array[i]));
				(*node.fab_array[i])[lev]->setVal(0.);
			}
		else
			for (int i = 0; i < cell.number_of_fabs; i++)
				(*cell.fab_array[i])[lev].reset(new amrex::MultiFab(grids[lev],dmap[lev],cell.ncomp_array[i],cell.nghost_array[i]));
		for (int i = 0; i < tmp_numfabs; i++)
		{
			bool match = false;
			if (a_nodal)
			{
				for (int j = 0; j < node.number_of_fabs; j++)
				{
					if (tmp_name_array[i] == node.name_array[j])
					{
						match = true;
						Util::Message(INFO,"Initializing ", node.name_array[j], "; nghost=",node.nghost_array[j], " with ", tmp_name_array[i] );
						amrex::MultiFab::Copy(*((*node.fab_array[j])[lev]).get(),tmpdata[lev],i,0,1,total_nghost);
					}
					for (int k = 0; k < node.ncomp_array[j]; k++)
					{
						if (tmp_name_array[i] == amrex::Concatenate(node.name_array[j],k+1,3))
						{
							match = true;
							Util::Message(INFO,"Initializing ", node.name_array[j], "[",k,"]; ncomp=", node.ncomp_array[j], "; nghost=",node.nghost_array[j], " with ", tmp_name_array[i] );
							amrex::MultiFab::Copy(*((*node.fab_array[j])[lev]).get(),tmpdata[lev],i,k,1,total_nghost);
						}
					}
					Util::RealFillBoundary(*((*node.fab_array[j])[lev]).get(),geom[lev]);
				}
			}
			else
			{
				for (int j = 0; j < cell.number_of_fabs; j++)
					for (int k = 0; k < cell.ncomp_array[j]; k++)
					{
						if (tmp_name_array[i] == amrex::Concatenate(cell.name_array[j],k+1,3))
						{
							match = true;
							Util::Message(INFO,"Initializing ", cell.name_array[j], "[",k,"]; ncomp=", cell.ncomp_array[j], "; nghost=",cell.nghost_array[j], " with ", tmp_name_array[i] );
							amrex::MultiFab::Copy(*((*cell.fab_array[j])[lev]).get(),tmpdata[lev],i,k,1,cell.nghost_array[j]);
						}
					}
			}
			if (!match) Util::Warning(INFO,"Fab ",tmp_name_array[i]," is in the restart file, but there is no fab with that name here.");
		}		

		for (unsigned int n = 0; n < m_basefields.size(); n++)
		{
			Util::Message(INFO,"n = ", n , " size = ", m_basefields.size());
			m_basefields[n]->MakeNewLevelFromScratch(lev,t_new[lev],grids[lev],dmap[lev]);
		}
	}

	SetFinestLevel(max_level);
}

void
Integrator::MakeNewLevelFromScratch (int lev, amrex::Real t, const amrex::BoxArray& cgrids,
				     const amrex::DistributionMapping& dm)
{
	BL_PROFILE("Integrator::MakeNewLevelFromScratch");
	for (int n = 0 ; n < cell.number_of_fabs; n++)
	{
		(*cell.fab_array[n])[lev].reset(new amrex::MultiFab(cgrids, dm, cell.ncomp_array[n], cell.nghost_array[n]));
		(*cell.fab_array[n])[lev]->setVal(0.0);
	}

	amrex::BoxArray ngrids = cgrids;
	ngrids.convert(amrex::IntVect::TheNodeVector());
	for (int n = 0 ; n < node.number_of_fabs; n++)
	{
		(*node.fab_array[n])[lev].reset(new amrex::MultiFab(ngrids, dm, node.ncomp_array[n], node.nghost_array[n]));
		(*node.fab_array[n])[lev]->setVal(0.0);
	}
	for (unsigned int n = 0; n < m_basefields.size(); n++)
	{
		m_basefields[n]->MakeNewLevelFromScratch(lev,t,cgrids,dm);
	}

	t_new[lev] = t;
	t_old[lev] = t - dt[lev];

	Initialize(lev);
  
	for (int n = 0 ; n < cell.number_of_fabs; n++)
	{
		cell.physbc_array[n]->define(geom[lev]);
		cell.physbc_array[n]->FillBoundary(*(*cell.fab_array[n])[lev],0,0,t,0);
//		Util::RealFillBoundary(*(*cell.fab_array[n])[lev],geom[lev]);
	}

	for (int n = 0 ; n < node.number_of_fabs; n++)
	{
		node.physbc_array[n]->define(geom[lev]);
		node.physbc_array[n]->FillBoundary(*(*node.fab_array[n])[lev],0,0,t,0);
//		Util::RealFillBoundary(*(*node.fab_array[n])[lev],geom[lev]);
	}
}

std::vector<std::string>
Integrator::PlotFileName (int lev,std::string prefix) const
{
	BL_PROFILE("Integrator::PlotFileName");
	std::vector<std::string> name;
	name.push_back(plot_file+"/"+prefix+"/");
	name.push_back(amrex::Concatenate("", lev, 5));
	return name;
}

void
Integrator::WritePlotFile (bool initial) const
{
	WritePlotFile(t_new[0], istep, initial, "");
}
void
Integrator::WritePlotFile (std::string prefix, Set::Scalar time, int step) const
{
	amrex::Vector<int> istep(max_level + 1, step);
	WritePlotFile(time, istep, false, prefix);
}

void
Integrator::WritePlotFile (Set::Scalar time, amrex::Vector<int> iter, bool initial, std::string prefix) const
{
	BL_PROFILE("Integrator::WritePlotFile");
	int nlevels = finest_level+1;
	if (max_plot_level >= 0) nlevels = std::min(nlevels,max_plot_level);

	int ccomponents = 0, ncomponents = 0;
	amrex::Vector<std::string> cnames, nnames;
	for (int i = 0; i < cell.number_of_fabs; i++)
	{
		if (!cell.writeout_array[i]) continue;
		ccomponents += cell.ncomp_array[i];
		if (cell.ncomp_array[i] > 1)
			for (int j = 1; j <= cell.ncomp_array[i]; j++)
				cnames.push_back(amrex::Concatenate(cell.name_array[i], j, 3));
		else
			cnames.push_back(cell.name_array[i]);
	}
	for (int i = 0; i < node.number_of_fabs; i++)
	{
		if (!node.writeout_array[i]) continue;
		ncomponents += node.ncomp_array[i];
		if (node.ncomp_array[i] > 1)
			for (int j = 1; j <= node.ncomp_array[i]; j++)
				nnames.push_back(amrex::Concatenate(node.name_array[i], j, 3));
		else
			nnames.push_back(node.name_array[i]);
	}

	amrex::Vector<amrex::MultiFab> cplotmf(nlevels), nplotmf(nlevels);
  
	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		if (ccomponents > 0)
		{
			cplotmf[ilev].define(grids[ilev], dmap[ilev], ccomponents, 0);

			int n = 0;
			for (int i = 0; i < cell.number_of_fabs; i++)
			{
				if (!cell.writeout_array[i]) continue;
				if ((*cell.fab_array[i])[ilev]->contains_nan()) Util::Abort(INFO,cnames[i]," contains nan (i=",i,")");
				if ((*cell.fab_array[i])[ilev]->contains_inf()) Util::Abort(INFO,cnames[i]," contains inf (i=",i,")");
				amrex::MultiFab::Copy(cplotmf[ilev], *(*cell.fab_array[i])[ilev], 0, n, cell.ncomp_array[i], 0);
				n += cell.ncomp_array[i];
			}
		}

		if (ncomponents > 0)
		{
			amrex::BoxArray ngrids = grids[ilev];
			ngrids.convert(amrex::IntVect::TheNodeVector());
			nplotmf[ilev].define(ngrids, dmap[ilev], ncomponents, 0);
			
			int n = 0;
			for (int i = 0; i < node.number_of_fabs; i++)
			{
				if (!node.writeout_array[i]) continue;
				if ((*node.fab_array[i])[ilev]->contains_nan()) 
				{
					Util::Warning(INFO,nnames[i]," contains nan (i=",i,"). Resetting to zero.");
					(*node.fab_array[i])[ilev]->setVal(0.0);
				}
				if ((*node.fab_array[i])[ilev]->contains_inf()) Util::Abort(INFO,nnames[i]," contains inf (i=",i,")");
				amrex::MultiFab::Copy(nplotmf[ilev], *(*node.fab_array[i])[ilev], 0, n, node.ncomp_array[i], 0);
				n += node.ncomp_array[i];
			}
		}
	}

	std::vector<std::string> plotfilename = PlotFileName(istep[0],prefix);
	if (initial) plotfilename[1] = plotfilename[1] + "init";
  
	if (ccomponents > 0)
	{
		WriteMultiLevelPlotfile(plotfilename[0]+plotfilename[1]+"cell", nlevels, amrex::GetVecOfConstPtrs(cplotmf), cnames,
					Geom(), time, iter, refRatio());
	
		std::ofstream chkptfile;
		chkptfile.open(plotfilename[0]+plotfilename[1]+"cell/Checkpoint");
		for (int i = 0; i <= max_level; i++) boxArray(i).writeOn(chkptfile);
		chkptfile.close();
	}

	if (ncomponents > 0)
	{
		WriteMultiLevelPlotfile(plotfilename[0]+plotfilename[1]+"node", nlevels, amrex::GetVecOfConstPtrs(nplotmf), nnames,
					Geom(), time, iter, refRatio());

		std::ofstream chkptfile;
		chkptfile.open(plotfilename[0]+plotfilename[1]+"node/Checkpoint");
		for (int i = 0; i <= max_level; i++) boxArray(i).writeOn(chkptfile);
		chkptfile.close();
	}

	if (amrex::ParallelDescriptor::IOProcessor())
	{
		std::ofstream coutfile, noutfile;
		if (istep[0]==0)
		{
			if (ccomponents > 0) coutfile.open(plot_file+"/celloutput.visit",std::ios_base::out);
			if (ncomponents > 0) noutfile.open(plot_file+"/nodeoutput.visit",std::ios_base::out);
		}
		else
		{
			if (ccomponents > 0) coutfile.open(plot_file+"/celloutput.visit",std::ios_base::app);
			if (ncomponents > 0) noutfile.open(plot_file+"/nodeoutput.visit",std::ios_base::app);
		}
		coutfile << plotfilename[1] + "cell" + "/Header" << std::endl;
		if (ncomponents > 0) noutfile << plotfilename[1] + "node" + "/Header" << std::endl;
	}
}

void
Integrator::Evolve ()
{
	BL_PROFILE("Integrator::Evolve");
	amrex::Real cur_time = t_new[0];
	int last_plot_file_step = 0;

	for (int step = istep[0]; step < max_step && cur_time < stop_time; ++step)
	{
		if (amrex::ParallelDescriptor::IOProcessor()) {
			std::cout << "\nSTEP " << step+1 << " starts ..." << std::endl;
		}
		int lev = 0;
		int iteration = 1;
		TimeStepBegin(cur_time,step);
		if (integrate_variables_before_advance) IntegrateVariables(cur_time,step);
		TimeStep(lev, cur_time, iteration);
		if (integrate_variables_after_advance) IntegrateVariables(cur_time,step);
		TimeStepComplete(cur_time,step);
		cur_time += dt[0];

		if (amrex::ParallelDescriptor::IOProcessor()) {
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
		else if (std::fabs(std::remainder(cur_time,plot_dt)) < 0.5*dt[0])
		{
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
Integrator::IntegrateVariables (amrex::Real time, int step)
{
	BL_PROFILE("Integrator::IntegrateVariables");
	if (!thermo.number) return;

	if ( (thermo.interval > 0 && (step) % thermo.interval == 0) ||
		 ((thermo.dt > 0.0) && (std::fabs(std::remainder(time,plot_dt)) < 0.5*dt[0])) )
	{
		// Zero out all variables
		for (int i = 0; i < thermo.number; i++) *thermo.vars[i] = 0; 

		// All levels except the finest
		for (int ilev = 0; ilev < max_level; ilev++)
		{
			const amrex::BoxArray& cfba = amrex::coarsen(grids[ilev+1], refRatio(ilev));

#ifdef OMP
			#pragma omp parallel
#endif
			for ( amrex::MFIter mfi(grids[ilev],dmap[ilev],true); mfi.isValid(); ++mfi )
			{
				const amrex::Box& box = mfi.tilebox();
				const amrex::BoxArray & comp = amrex::complementIn(box,cfba);

				for (int i = 0; i < comp.size(); i++)
				{
					Integrate(ilev,time, step,
						  mfi, comp[i]);
				}
			}
		}
		// Now do the finest level
		{
			#ifdef OMP
			#pragma omp parallel
			#endif 
			for ( amrex::MFIter mfi(grids[max_level],dmap[max_level],true); mfi.isValid(); ++mfi )
			{
				const amrex::Box& box = mfi.tilebox();
				Integrate(max_level, time, step, mfi, box);
			}
		}

		// Sum up across all processors
		for (int i = 0; i < thermo.number; i++) 
		{
			amrex::ParallelDescriptor::ReduceRealSum(*thermo.vars[i]);
		}
	}
	if ( amrex::ParallelDescriptor::IOProcessor() &&
		 (
			 (thermo.plot_int > 0 && step % thermo.plot_int == 0) ||
			 (thermo.plot_dt > 0.0 && std::fabs(std::remainder(time,thermo.plot_dt)) < 0.5*dt[0])
		 ))
	{
		std::ofstream outfile;
		if (step==0)
		{
			outfile.open(plot_file+"/thermo.dat",std::ios_base::out);
			outfile << "time";
			for (int i = 0; i < thermo.number; i++) 
				outfile << "\t" << thermo.names[i];
			outfile << std::endl;
		}
		else outfile.open(plot_file+"/thermo.dat",std::ios_base::app);
		outfile << time;
		for (int i = 0; i < thermo.number; i++)
			outfile << "\t" << *thermo.vars[i];
		outfile << std::endl;
		outfile.close();
	}

}


void
Integrator::TimeStep (int lev, amrex::Real time, int /*iteration*/)
{
	BL_PROFILE("Integrator::TimeStep");
	if (regrid_int > 0)  // We may need to regrid
	{
		static amrex::Vector<int> last_regrid_step(max_level+1, 0);

		// regrid doesn't change the base level, so we don't regrid on max_level
		if (lev < max_level && istep[lev] > last_regrid_step[lev])
		{
			if (istep[lev] % regrid_int == 0)
			{
				regrid(lev, time, false); 
			}
		}
	}
	SetFinestLevel(finest_level);

	if (Verbose() && amrex::ParallelDescriptor::IOProcessor()) {
		std::cout << "[Level " << lev 
			  << " step " << istep[lev]+1 << "] ";
		std::cout << "ADVANCE with dt = "
			  << dt[lev]
			  << std::endl;
	}

	for (int n = 0 ; n < cell.number_of_fabs ; n++)
		FillPatch(lev,time,*cell.fab_array[n],*(*cell.fab_array[n])[lev],*cell.physbc_array[n],0);
	for (int n = 0 ; n < node.number_of_fabs ; n++)
		FillPatch(lev,time,*node.fab_array[n],*(*node.fab_array[n])[lev],*node.physbc_array[n],0);

	Advance(lev, time, dt[lev]);
	++istep[lev];

	if (Verbose() && amrex::ParallelDescriptor::IOProcessor())
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
