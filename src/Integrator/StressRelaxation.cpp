#include "StressRelaxation.H"
//#if AMREX_SPACEDIM == 1
namespace Integrator
{
StressRelaxation::StressRelaxation():
	Integrator()
{
	//
	// READ INPUT PARAMETERS
	//

	amrex::Vector<std::string> bc_hi_str({AMREX_D_DECL("EXT_DIR", "EXT_DIR", "EXT_DIR")});
	amrex::Vector<std::string> bc_lo_str({AMREX_D_DECL("EXT_DIR", "EXT_DIR", "EXT_DIR")});
	water_bc = new BC::Constant(bc_hi_str, bc_lo_str
							  ,AMREX_D_DECL({0.}, {0.}, {0.})
							  ,AMREX_D_DECL({0.}, {0.}, {0.})
							  );
	RegisterNewFab(water_conc,water_bc, 1, 1, "Water Concentration",true);
	water_ic = new IC::Constant(geom,{0.0});

	// ---------------------------------------------------------------------
	// --------------------- Material model --------------------------------
	// ---------------------------------------------------------------------
	amrex::ParmParse pp_material("material");
	pp_material.query("model",input_material);
	if(input_material == "isotropic")
	{
		amrex::Vector<Set::Scalar> E_i;
		amrex::Vector<Set::Scalar> tau_i;

		amrex::Real nu = 0.3;
		int prony_terms = 0;
		
		amrex::ParmParse pp_material_isotropic("material.isotropic");
		pp_material_isotropic.query("nu",nu);
		pp_material_isotropic.query("number_of_terms",prony_terms);
		//Sanity check
		if(prony_terms < 0) Util::Abort(INFO,"Number of prony series terms must be non-negative");
		if(prony_terms > 8) Util::Abort(INFO,"Number of prony series terms must not exceed 8");

		pp_material_isotropic.queryarr("E_i", E_i);
		//Sanity check
		if (E_i.size() == 0) Util::Abort(INFO, "E_i must contain at least one value");
		if (E_i.size() != prony_terms+1) Util::Abort(INFO, "E_i must have prony_terms + 1 entries");
		//for (int i = 0; i < E_i.size(); i++)
		//	if(E_i[i] < 0.0)  Util::Abort(INFO, "E_i can not be less than zero");


		pp_material_isotropic.queryarr("tau_i", tau_i);
		//Sanity check
		if(tau_i.size() == 0) Util::Abort(INFO, "tau_i must contain at least one value");
		if(tau_i.size() != prony_terms) Util::Abort(INFO, "tau_i must have prony_terms entries");
		//for (int i = 0; i < tau_i.size(); i++)
		//	if(tau_i[i] < 0.0)	Util::Abort(INFO, "tau_i can not be less than zero");


		std::array<Set::Scalar,9> youngs_modulus;
		std::array<Set::Scalar,8> relaxation_times;

		for(int i = 0; i < 9; i++)
		{
			if(i < prony_terms)
			{
				youngs_modulus[i] = E_i[i];
				relaxation_times[i] = tau_i[i];
			}
			else if(i == prony_terms) youngs_modulus[i] = E_i[i];
			else
			{
				youngs_modulus[i] = 0.0;
				relaxation_times[i] = 0.0;
			}

		}
		
		modeltype = Model::Solid::Viscoelastic::Isotropic(nu, prony_terms, youngs_modulus, relaxation_times);
	}
	else
		Util::Abort(INFO, "Not implemented yet");

	// ---------------------------------------------------------------------
	// --------------------- Elasticity parameters -------------------------
	// ---------------------------------------------------------------------
	amrex::ParmParse pp_elastic("elastic");
	pp_elastic.query("on",elastic_on);
	if(elastic_on)
	{
		pp_elastic.query("int",elastic_int);
		pp_elastic.query("type",elastic_type);
		pp_elastic.query("max_iter",elastic_max_iter);
		pp_elastic.query("max_fmg_iter",elastic_max_fmg_iter);
		pp_elastic.query("verbose",elastic_verbose);
		pp_elastic.query("cgverbose",elastic_cgverbose);
		pp_elastic.query("tol_rel",elastic_tol_rel);
		pp_elastic.query("tol_abs",elastic_tol_abs);
		pp_elastic.query("use_fsmooth",elastic_use_fsmooth);
		pp_elastic.query("agglomeration", agglomeration);
		pp_elastic.query("consolidation", consolidation);

		amrex::ParmParse pp_temp;
		Set::Scalar stop_time, timestep;
		pp_temp.query("timestep",timestep);
		pp_temp.query("stop_time",stop_time);

		elastic_tstart = 0.0;
		elastic_tend = stop_time;

		pp_elastic.query("bottom_solver",bottom_solver);
		pp_elastic.query("linop_maxorder", linop_maxorder);
		pp_elastic.query("max_coarsening_level",max_coarsening_level);
		pp_elastic.query("verbose",elastic_verbose);
		pp_elastic.query("cg_verbose", elastic_cgverbose);
		pp_elastic.query("bottom_max_iter", elastic_bottom_max_iter);
		if (pp_elastic.countval("body_force")) pp_elastic.getarr("body_force",body_force);

		amrex::ParmParse pp_elastic_bc("elastic.bc");
		amrex::Vector<std::string> AMREX_D_DECL(bc_x_lo_str,bc_y_lo_str,bc_z_lo_str);
		amrex::Vector<std::string> AMREX_D_DECL(bc_x_hi_str,bc_y_hi_str,bc_z_hi_str);

		AMREX_D_TERM(	pp_elastic_bc.queryarr("bc_x_lo",bc_x_lo_str);,
						pp_elastic_bc.queryarr("bc_y_lo",bc_y_lo_str);,
						pp_elastic_bc.queryarr("bc_z_lo",bc_z_lo_str););
		AMREX_D_TERM(	pp_elastic_bc.queryarr("bc_x_hi",bc_x_hi_str);,
						pp_elastic_bc.queryarr("bc_y_hi",bc_y_hi_str);,
						pp_elastic_bc.queryarr("bc_z_hi",bc_z_hi_str););

		bc_map["displacement"] = Operator::Elastic<Model::Solid::Viscoelastic::Isotropic>::BC::Displacement;
		bc_map["disp"] = Operator::Elastic<Model::Solid::Viscoelastic::Isotropic>::BC::Displacement;
		bc_map["traction"] = Operator::Elastic<Model::Solid::Viscoelastic::Isotropic>::BC::Traction;
		bc_map["trac"] = Operator::Elastic<Model::Solid::Viscoelastic::Isotropic>::BC::Traction;
		bc_map["periodic"] = Operator::Elastic<Model::Solid::Viscoelastic::Isotropic>::BC::Periodic;

		AMREX_D_TERM(	bc_x_lo = {AMREX_D_DECL(bc_map[bc_x_lo_str[0]], bc_map[bc_x_lo_str[1]], bc_map[bc_x_lo_str[2]])};
						bc_x_hi = {AMREX_D_DECL(bc_map[bc_x_hi_str[0]], bc_map[bc_x_hi_str[1]], bc_map[bc_x_hi_str[2]])};
						,
						bc_y_lo = {AMREX_D_DECL(bc_map[bc_y_lo_str[0]], bc_map[bc_y_lo_str[1]], bc_map[bc_y_lo_str[2]])};
						bc_y_hi = {AMREX_D_DECL(bc_map[bc_y_hi_str[0]], bc_map[bc_y_hi_str[1]], bc_map[bc_y_hi_str[2]])};
						,
						bc_z_lo = {AMREX_D_DECL(bc_map[bc_z_lo_str[0]], bc_map[bc_z_lo_str[1]], bc_map[bc_z_lo_str[2]])};
						bc_z_hi = {AMREX_D_DECL(bc_map[bc_z_hi_str[0]], bc_map[bc_z_hi_str[1]], bc_map[bc_z_hi_str[2]])};);

		amrex::Vector<Set::Scalar> bc_lo_1, bc_hi_1;
		amrex::Vector<Set::Scalar> bc_lo_2, bc_hi_2;
		amrex::Vector<Set::Scalar> bc_lo_3, bc_hi_3;
		amrex::Vector<Set::Scalar> bc_lo_1_t, bc_hi_1_t;
		amrex::Vector<Set::Scalar> bc_lo_2_t, bc_hi_2_t;
		amrex::Vector<Set::Scalar> bc_lo_3_t, bc_hi_3_t;

		if (pp_elastic_bc.countval("left_face")) pp_elastic_bc.getarr("left_face",bc_lo_1);
		if(bc_lo_1.size() % AMREX_SPACEDIM !=0)
			Util::Abort(INFO, "Invalid number of values for left_face displacement");

		if (pp_elastic_bc.countval("right_face")) pp_elastic_bc.getarr("right_face",bc_hi_1);
		if(bc_hi_1.size() % AMREX_SPACEDIM !=0)
			Util::Abort(INFO, "Invalid number of values for right_face displacement");

		if (pp_elastic_bc.countval("bottom_face")) pp_elastic_bc.getarr("bottom_face",bc_lo_2);
		if(bc_lo_2.size() % AMREX_SPACEDIM !=0)
			Util::Abort(INFO, "Invalid number of values for bottom_face displacement");

		if (pp_elastic_bc.countval("top_face")) pp_elastic_bc.getarr("top_face",bc_hi_2);
		if(bc_hi_2.size() % AMREX_SPACEDIM !=0)
			Util::Abort(INFO, "Invalid number of values for top_face displacement");

#if AMREX_SPACEDIM>2
		if (pp_elastic_bc.countval("back_face")) pp_elastic_bc.getarr("back_face",bc_lo_3);
		if(bc_lo_3.size() % AMREX_SPACEDIM !=0)
			Util::Abort(INFO, "Invalid number of values for back_face displacement");

		if (pp_elastic_bc.countval("front_face")) pp_elastic_bc.getarr("front_face",bc_hi_3);
		if(bc_hi_3.size() % AMREX_SPACEDIM !=0)
			Util::Abort(INFO, "Invalid number of values for front_face displacement");
#endif

		int tempSize = bc_lo_1.size()/AMREX_SPACEDIM;
		for (int i = 0; i<tempSize; i++)
			elastic_bc_left.push_back(Set::Vector(AMREX_D_DECL(bc_lo_1[AMREX_SPACEDIM*i],bc_lo_1[AMREX_SPACEDIM*i+1],bc_lo_1[AMREX_SPACEDIM*i+2])));
		if(pp_elastic_bc.countval("left_face_t")) pp_elastic_bc.getarr("left_face_t",bc_lo_1_t);
		if(bc_lo_1_t.size() == tempSize)
			elastic_bc_left_t = bc_lo_1_t;
		else
			for (int j = 0; j < tempSize; j++)
				elastic_bc_left_t.push_back(elastic_tstart + j*(elastic_tend - elastic_tstart)/(tempSize-1.0 != 0.0 ? tempSize-1.0 : 1.0));


		tempSize = bc_hi_1.size()/AMREX_SPACEDIM;
		for (int i = 0; i<tempSize; i++)
			elastic_bc_right.push_back(Set::Vector(AMREX_D_DECL(bc_hi_1[AMREX_SPACEDIM*i],bc_hi_1[AMREX_SPACEDIM*i+1],bc_hi_1[AMREX_SPACEDIM*i+2])));
		if(pp_elastic_bc.countval("right_face_t")) pp_elastic_bc.getarr("right_face_t",bc_hi_1_t);
		if(bc_hi_1_t.size() == tempSize)
			elastic_bc_right_t = bc_hi_1_t;
		else
			for (int j = 0; j < tempSize; j++)
				elastic_bc_right_t.push_back(elastic_tstart + j*(elastic_tend - elastic_tstart)/(tempSize-1.0 != 0.0 ? tempSize-1.0 : 1.0));


		tempSize = bc_lo_2.size()/AMREX_SPACEDIM;
		for (int i = 0; i<tempSize; i++)
			elastic_bc_bottom.push_back(Set::Vector(AMREX_D_DECL(bc_lo_2[AMREX_SPACEDIM*i],bc_lo_2[AMREX_SPACEDIM*i+1],bc_lo_2[AMREX_SPACEDIM*i+2])));
		if(pp_elastic_bc.countval("bottom_face_t")) pp_elastic_bc.getarr("bottom_face_t",bc_lo_2_t);
		if(bc_lo_2_t.size() == tempSize)
			elastic_bc_bottom_t = bc_lo_2_t;
		else
			for (int j = 0; j < tempSize; j++)
				elastic_bc_bottom_t.push_back(elastic_tstart + j*(elastic_tend - elastic_tstart)/(tempSize-1.0 != 0.0 ? tempSize-1.0 : 1.0));

		tempSize = bc_hi_2.size()/AMREX_SPACEDIM;
		for (int i = 0; i<tempSize; i++)
			elastic_bc_top.push_back(Set::Vector(AMREX_D_DECL(bc_hi_2[AMREX_SPACEDIM*i],bc_hi_2[AMREX_SPACEDIM*i+1],bc_hi_2[AMREX_SPACEDIM*i+2])));
		if(pp_elastic_bc.countval("top_face_t")) pp_elastic_bc.getarr("top_face_t",bc_hi_2_t);
		if(bc_hi_2_t.size() == tempSize)
			elastic_bc_top_t = bc_hi_2_t;
		else
			for (int j = 0; j < tempSize; j++)
				elastic_bc_top_t.push_back(elastic_tstart + j*(elastic_tend - elastic_tstart)/(tempSize-1.0 != 0.0 ? tempSize-1.0 : 1.0));

#if AMREX_SPACEDIM > 2
		tempSize = bc_lo_3.size()/AMREX_SPACEDIM;
		for (int i = 0; i<tempSize; i++)
			elastic_bc_back.push_back(Set::Vector(AMREX_D_DECL(bc_lo_3[AMREX_SPACEDIM*i],bc_lo_3[AMREX_SPACEDIM*i+1],bc_lo_3[AMREX_SPACEDIM*i+2])));
		if(pp_elastic_bc.countval("back_face_t")) pp_elastic_bc.getarr("back_face_t",bc_lo_3_t);
		if(bc_lo_3_t.size() == tempSize)
			elastic_bc_back_t = bc_lo_3_t;
		else
			for (int j = 0; j < tempSize; j++)
				elastic_bc_back_t.push_back(elastic_tstart + j*(elastic_tend - elastic_tstart)/(tempSize-1.0 != 0.0 ? tempSize-1.0 : 1.0));

		tempSize = bc_hi_3.size()/AMREX_SPACEDIM;
		for (int i = 0; i<tempSize; i++)
			elastic_bc_front.push_back(Set::Vector(AMREX_D_DECL(bc_hi_3[AMREX_SPACEDIM*i],bc_hi_3[AMREX_SPACEDIM*i+1],bc_hi_3[AMREX_SPACEDIM*i+2])));
		if(pp_elastic_bc.countval("front_face_t")) pp_elastic_bc.getarr("front_face_t",bc_hi_3_t);
		if(bc_hi_3_t.size() == tempSize)
			elastic_bc_front_t = bc_hi_3_t;
		else
			for (int j = 0; j < tempSize; j++)
				elastic_bc_front_t.push_back(elastic_tstart + j*(elastic_tend - elastic_tstart)/(tempSize-1.0 != 0.0 ? tempSize-1.0 : 1.0));
#endif
		
		//----------------------------------------------------------------------
		// The following routine should be replaced by RegisterNewFab in the
		// future. For now, we are manually defining and resizing
		//-----------------------------------------------------------------------
		amrex::ParmParse pp_amr("amr");
		int maxLev = 0, n_cell = 8, max_grid_size = 64;
		int ref_ratio = 2;
		pp_amr.query("max_level",maxLev);
		pp_amr.query("n_cell",n_cell);
		pp_amr.query("max_grid_size",max_grid_size);
		pp_amr.query("ref_ratio",ref_ratio);

		nlevels = maxLev+1;

		ngrids.resize(nlevels);
		
		displacement.resize(nlevels);
		rhs.resize(nlevels);
		strain.resize(nlevels);
		stress.resize(nlevels);
		stress_vm.resize(nlevels);
		energy.resize(nlevels);
		model.resize(nlevels);

		Util::Message(INFO);
	}
}


void
StressRelaxation::Advance (int /*lev*/, amrex::Real /*time*/, amrex::Real /*dt*/)
{
}

void
StressRelaxation::Initialize (int lev)
{
	water_ic->Initialize(lev,water_conc);
}


void
StressRelaxation::TagCellsForRefinement (int /*lev*/, amrex::TagBoxArray& /*tags*/, amrex::Real /*time*/, int /*ngrow*/)
{
}


std::vector<std::string>
StressRelaxation::PlotFileNameNode (std::string plot_file_name, int lev) const
{
	std::vector<std::string> name;
	name.push_back(plot_file_name+"/");
	name.push_back(amrex::Concatenate("", lev, 5));
	return name;
}

void 
StressRelaxation::TimeStepComplete(amrex::Real /*time*/, int /*iter*/)
{
	if (! elastic_on) return;

#if AMREX_SPACEDIM == 1
	const int ncomp = 5;
	Vector<std::string> varname = {"u01", "rhs01", "stress11", "strain11", "energy"};

#elif AMREX_SPACEDIM == 2
	const int ncomp = 11;
	Vector<std::string> varname = {"u01", "u02", "rhs01", "rhs02", "stress11", "stress22", "stress12",
									"strain11", "strain22", "strain12", "energy"};

#elif AMREX_SPACEDIM == 3
	const int ncomp = 19;
	Vector<std::string> varname = {"u01", "u02", "u03", "rhs01", "rhs02", "rhs03","stress11", "stress22",
								 "stress33", "stress23", "stress13", "stress12", 
								 "strain11", "strain22","strain33", "strain23", "strain13", 
								 "strain12", 
								 "energy"};
#endif

	Vector<amrex::MultiFab> plotmf(nlevels);
	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		//plotmf[ilev].define(ngrids[ilev], ndmap[ilev], ncomp, 0);
		plotmf[ilev].define(ngrids[ilev], dmap[ilev], ncomp, 0);
#if AMREX_SPACEDIM == 2
		MultiFab::Copy(plotmf[ilev], displacement      [ilev], 0, 0, 1, 0);
		MultiFab::Copy(plotmf[ilev], displacement      [ilev], 1, 1, 1, 0);
		MultiFab::Copy(plotmf[ilev], rhs    [ilev], 0, 2, 1, 0);
		MultiFab::Copy(plotmf[ilev], rhs    [ilev], 1, 3, 1, 0);
		//MultiFab::Copy(plotmf[ilev], res    [ilev], 0, 4, 1, 0);
		//MultiFab::Copy(plotmf[ilev], res    [ilev], 1, 5, 1, 0);
		MultiFab::Copy(plotmf[ilev], stress [ilev], 0, 4, 1, 0);
		MultiFab::Copy(plotmf[ilev], stress [ilev], 3, 5, 1, 0);
		MultiFab::Copy(plotmf[ilev], stress [ilev], 1, 6, 1, 0);
		MultiFab::Copy(plotmf[ilev], strain [ilev], 0, 7, 1, 0);
		MultiFab::Copy(plotmf[ilev], strain [ilev], 3, 8, 1, 0);
		MultiFab::Copy(plotmf[ilev], strain [ilev], 1, 9, 1, 0);
		MultiFab::Copy(plotmf[ilev], energy [ilev], 0, 10, 1, 0);
#elif AMREX_SPACEDIM == 3
		MultiFab::Copy(plotmf[ilev], displacement      [ilev], 0, 0,  1, 0);
		MultiFab::Copy(plotmf[ilev], displacement      [ilev], 1, 1,  1, 0);
		MultiFab::Copy(plotmf[ilev], displacement      [ilev], 2, 2,  1, 0);
		MultiFab::Copy(plotmf[ilev], rhs    [ilev], 0, 3,  1, 0);
		MultiFab::Copy(plotmf[ilev], rhs    [ilev], 1, 4,  1, 0);
		MultiFab::Copy(plotmf[ilev], rhs    [ilev], 2, 5,  1, 0);
		//MultiFab::Copy(plotmf[ilev], res    [ilev], 0, 6,  1, 0);
		//MultiFab::Copy(plotmf[ilev], res    [ilev], 1, 7,  1, 0);
		//MultiFab::Copy(plotmf[ilev], res    [ilev], 2, 8,  1, 0);
		MultiFab::Copy(plotmf[ilev], stress [ilev], 0, 6,  1, 0);
		MultiFab::Copy(plotmf[ilev], stress [ilev], 4, 7,  1, 0);
		MultiFab::Copy(plotmf[ilev], stress [ilev], 8, 8,  1, 0);
		MultiFab::Copy(plotmf[ilev], stress [ilev], 5, 9,  1, 0);
		MultiFab::Copy(plotmf[ilev], stress [ilev], 2, 10, 1, 0);
		MultiFab::Copy(plotmf[ilev], stress [ilev], 1, 11, 1, 0);
		MultiFab::Copy(plotmf[ilev], strain [ilev], 0, 12,  1, 0);
		MultiFab::Copy(plotmf[ilev], strain [ilev], 4, 13,  1, 0);
		MultiFab::Copy(plotmf[ilev], strain [ilev], 8, 14,  1, 0);
		MultiFab::Copy(plotmf[ilev], strain [ilev], 5, 15,  1, 0);
		MultiFab::Copy(plotmf[ilev], strain [ilev], 2, 16, 1, 0);
		MultiFab::Copy(plotmf[ilev], strain [ilev], 1, 17, 1, 0);
		MultiFab::Copy(plotmf[ilev], energy	[ilev], 0, 18, 1, 0);
#endif 
	}
	//Util::Message(INFO);
	std::string plot_file_node = plot_file+"-node";
	const std::vector<std::string>& plotfilename = PlotFileNameNode(plot_file_node,istep[0]);
	//Util::Message(INFO);
	WriteMultiLevelPlotfile(plotfilename[0]+plotfilename[1], nlevels, amrex::GetVecOfConstPtrs(plotmf), varname,
				Geom(), t_new[0], istep, refRatio());
	//Util::Message(INFO);
	if (ParallelDescriptor::IOProcessor())
	{
		std::ofstream outfile;
		if (istep[0]==0) outfile.open(plot_file_node+"/output.visit",std::ios_base::out);
		else outfile.open(plot_file_node+"/output.visit",std::ios_base::app);
		outfile << plotfilename[1] + "/Header" << std::endl;
	}
}

void 
StressRelaxation::TimeStepBegin(amrex::Real time, int /*iter*/)
{
	if (!elastic_on) return;
	
	ngrids.resize(nlevels);

	displacement.resize(nlevels);
	rhs.resize(nlevels);
	strain.resize(nlevels);
	stress.resize(nlevels);
	stress_vm.resize(nlevels);
	energy.resize(nlevels);
	model.resize(nlevels);
	
	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		ngrids[ilev] = grids[ilev];
		ngrids[ilev].convert(IntVect::TheNodeVector());
		Util::Message(INFO,"ngrid size = ",ngrids[ilev].size(), ". Grids size = ", grids[ilev].size());
	}

	int number_of_stress_components = AMREX_SPACEDIM*AMREX_SPACEDIM;
	int number_of_components = AMREX_SPACEDIM;
	int number_of_ghost_cells = 1;

	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		displacement[ilev].define(ngrids[ilev], dmap[ilev], number_of_components, number_of_ghost_cells);
		rhs     	[ilev].define(ngrids[ilev], dmap[ilev], number_of_components, number_of_ghost_cells);
		stress 		[ilev].define(ngrids[ilev], dmap[ilev], number_of_stress_components, number_of_ghost_cells);
		strain		[ilev].define(ngrids[ilev], dmap[ilev], number_of_stress_components, number_of_ghost_cells);
		energy 		[ilev].define(ngrids[ilev], dmap[ilev], 1, number_of_ghost_cells);
		stress_vm	[ilev].define(ngrids[ilev], dmap[ilev], 1, number_of_ghost_cells);
		model		[ilev].define(ngrids[ilev], dmap[ilev], 1, 1);
	}
	
	info.setAgglomeration(agglomeration);
	info.setConsolidation(consolidation);
	info.setMaxCoarseningLevel(max_coarsening_level);
	
	elastic_operator.define(geom, grids, dmap, info);
	elastic_operator.setMaxOrder(linop_maxorder);
	
	geom[0].isPeriodic(0);
	elastic_operator.SetBC({{AMREX_D_DECL(bc_x_lo,bc_y_lo,bc_z_lo)}},
		     				{{AMREX_D_DECL(bc_x_hi,bc_y_hi,bc_z_hi)}});

	AMREX_D_TERM(	Numeric::Interpolator::Linear<Set::Vector> interpolate_left(elastic_bc_left,elastic_bc_left_t);
					Numeric::Interpolator::Linear<Set::Vector> interpolate_right(elastic_bc_right,elastic_bc_right_t);
					,
					Numeric::Interpolator::Linear<Set::Vector> interpolate_bottom(elastic_bc_bottom,elastic_bc_bottom_t);
					Numeric::Interpolator::Linear<Set::Vector> interpolate_top(elastic_bc_top,elastic_bc_top_t);
					,
					Numeric::Interpolator::Linear<Set::Vector> interpolate_back(elastic_bc_back,elastic_bc_back_t);
					Numeric::Interpolator::Linear<Set::Vector> interpolate_front(elastic_bc_front,elastic_bc_front_t););
	
	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		modeltype.SetTime(time);
		model[ilev].setVal(modeltype);
		//if(iter == 0 || time == elastic_tstart)
		{
			displacement[ilev].setVal(0.0);
			strain[ilev].setVal(0.0);
			stress[ilev].setVal(0.0);
			stress_vm[ilev].setVal(0.0);
			rhs[ilev].setVal(0.0);
			energy[ilev].setVal(0.0);	
		}
	}

	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		const Real* DX = geom[ilev].CellSize();
		Set::Scalar volume = AMREX_D_TERM(DX[0],*DX[1],*DX[2]);

		//Util::Message(INFO);
		AMREX_D_TERM(rhs[ilev].setVal(body_force[0]*volume,0,1);,
					rhs[ilev].setVal(body_force[1]*volume,1,1);,
					rhs[ilev].setVal(body_force[2]*volume,2,1););

		elastic_operator.SetModel(ilev,model[ilev]);
		
		for (amrex::MFIter mfi(rhs[ilev],true); mfi.isValid(); ++mfi)
		{
		 	const amrex::Box& box = mfi.tilebox();
		 	amrex::BaseFab<amrex::Real> &rhsfab = (rhs[ilev])[mfi];

		 	AMREX_D_TERM(for (int i = box.loVect()[0]; i<=box.hiVect()[0]; i++),
		 		     for (int j = box.loVect()[1]; j<=box.hiVect()[1]; j++),
		 		     for (int k = box.loVect()[2]; k<=box.hiVect()[2]; k++))
		 	{
		 		bool AMREX_D_DECL(xmin = false, ymin = false, zmin = false);
				bool AMREX_D_DECL(xmax = false, ymax = false, zmax = false);

		 		AMREX_D_TERM(	xmin = (i == geom[ilev].Domain().loVect()[0]);
		 						xmax = (i == geom[ilev].Domain().hiVect()[0]+1);
		 						,
		 						ymin = (j == geom[ilev].Domain().loVect()[1]);
		 						ymax = (j == geom[ilev].Domain().hiVect()[1]+1);
		 						,
		 						zmin = (k == geom[ilev].Domain().loVect()[2]);
		 						zmax = (k == geom[ilev].Domain().hiVect()[2]+1););

		 		if (false || AMREX_D_TERM(xmin || xmax, || ymin || ymax, || zmin || zmax))
		 		{
		 			AMREX_D_TERM(	rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),0) = 0.0;,
		 							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),1) = 0.0;,
		 							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),2) = 0.0;);
		 		}
				for(int l = 0; l<AMREX_SPACEDIM; l++)
				{
					if(xmin && bc_x_lo[l]==Operator::Elastic<Model::Solid::Viscoelastic::Isotropic>::BC::Displacement)
					{
						if(elastic_bc_left.size() == 1)
							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = elastic_bc_left[0](l);
						else	
							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = interpolate_left(time)[l];
					}
					if(xmax && bc_x_hi[l]==Operator::Elastic<Model::Solid::Viscoelastic::Isotropic>::BC::Displacement)
					{
						if(elastic_bc_right.size() == 1)
							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = elastic_bc_right[0](l);
						else	
							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = interpolate_right(time)[l];
					}
#if AMREX_SPACEDIM > 1
					if(ymin && bc_y_lo[l]==Operator::Elastic<Model::Solid::Viscoelastic::Isotropic>::BC::Displacement)
					{
						if(elastic_bc_bottom.size() == 1)
							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = elastic_bc_bottom[0](l);
						else	
							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = interpolate_bottom(time)[l];
					}
					if(ymax && bc_y_hi[l]==Operator::Elastic<Model::Solid::Viscoelastic::Isotropic>::BC::Displacement)
					{
						if(elastic_bc_top.size() == 1)
							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = elastic_bc_top[0](l);
						else	
							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = interpolate_top(time)[l];
					}
#if AMREX_SPACEDIM > 2
					if(zmin && bc_z_lo[l]==Operator::Elastic<Model::Solid::Viscoelastic::Isotropic>::BC::Displacement)
					{
						if(elastic_bc_back.size() == 1)
							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = elastic_bc_back[0](l);
						else	
							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = interpolate_back(time)[l];
					}
					if(zmax && bc_z_hi[l]==Operator::Elastic<Model::Solid::Viscoelastic::Isotropic>::BC::Displacement)
					{
						if(elastic_bc_front.size() == 1)
							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = elastic_bc_front[0](l);
						else	
							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = interpolate_front(time)[l];
					}
#endif
#endif
				}
		 	}
		}
	}
	
	amrex::MLMG solver(elastic_operator);
	solver.setMaxIter(elastic_max_iter);
	solver.setMaxFmgIter(elastic_max_fmg_iter);
	solver.setVerbose(elastic_verbose);
	solver.setCGVerbose(elastic_cgverbose);
	solver.setBottomMaxIter(elastic_bottom_max_iter);
	if (bottom_solver == "cg")
		solver.setBottomSolver(MLMG::BottomSolver::cg);
	else if (bottom_solver == "bicgstab")
		solver.setBottomSolver(MLMG::BottomSolver::bicgstab);
	if (!elastic_use_fsmooth)// <<< put in to NOT require FSmooth
	{
		solver.setFinalSmooth(0);
		solver.setBottomSmooth(0);
	}
	
	solver.solve(GetVecOfPtrs(displacement),
		GetVecOfConstPtrs(rhs),
		elastic_tol_rel,
		elastic_tol_abs);

	for (int lev = 0; lev < nlevels; lev++)
	{
		elastic_operator.Strain(lev,strain[lev],displacement[lev]);
		elastic_operator.Stress(lev,stress[lev],displacement[lev]);
		elastic_operator.Energy(lev,energy[lev],displacement[lev]);
	}
	/*

#if AMREX_SPACEDIM == 1
	static amrex::IntVect dx(1);
#elif AMREX_SPACEDIM == 2
	static amrex::IntVect dx(1,0), dy(0,1);
#elif AMREX_SPACEDIM == 3
	static amrex::IntVect dx(1,0,0), dy(0,1,0), dz(0,0,1);
#endif
	//std::cout<< __FILE__ << ":" << __LINE__ << ": TimeStepBegin()" << std::endl;

	for (int lev = 0; lev < displacement.size(); lev++)
	{
		const amrex::Real* DX = geom[lev].CellSize();
		for ( amrex::MFIter mfi(displacement[lev],true); mfi.isValid(); ++mfi )
		{
			const Box& bx = mfi.tilebox();

			FArrayBox &ufab  = (displacement[lev])[mfi];
			FArrayBox &epsfab  = (strain[lev])[mfi];
			FArrayBox &sigmafab  = (stress[lev])[mfi];
			FArrayBox &sigmavmfab  = (stress_vm[lev])[mfi];

			elastic_operator.Stress(lev,stress[lev],displacement[lev]);
			elastic_operator.Energy(lev,energy[lev],displacement[lev]);

			for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
				for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
#if AMREX_SPACEDIM > 2
					for (int k = bx.loVect()[2]; k<= bx.hiVect()[2]; k++)
#endif
				{
					amrex::IntVect m(AMREX_D_DECL(i,j,k));

#if AMREX_SPACEDIM==2
					epsfab(m,0) = (ufab(m+dx,0) - ufab(m-dx,0))/(2.0*DX[0]);
					epsfab(m,1) = (ufab(m+dy,0) - ufab(m-dy,0))/(2.0*DX[1]);
					epsfab(m,2) = (ufab(m+dx,1) - ufab(m-dx,1))/(2.0*DX[0]);
					epsfab(m,3) = (ufab(m+dy,1) - ufab(m-dy,1))/(2.0*DX[1]);
					Set::Scalar trace = (sigmafab(m,0) + sigmafab(m,3))/3.0;
					sigmavmfab(m) = sqrt(
								1.5*(
									(sigmafab(m,0)-trace)*(sigmafab(m,0)-trace) +
									(sigmafab(m,1))*(sigmafab(m,1)) +
									(sigmafab(m,2))*(sigmafab(m,2)) +
									((sigmafab(m,3)-trace)*(sigmafab(m,3)-trace))
								)
							);

#elif AMREX_SPACEDIM==3
					epsfab(m,0) = (ufab(m+dx,0) - ufab(m-dx,0))/(2.0*DX[0]);
					epsfab(m,1) = (ufab(m+dy,0) - ufab(m-dy,0))/(2.0*DX[1]);
					epsfab(m,2) = (ufab(m+dz,0) - ufab(m-dz,0))/(2.0*DX[2]);
					epsfab(m,3) = (ufab(m+dx,1) - ufab(m-dx,1))/(2.0*DX[0]);
					epsfab(m,4) = (ufab(m+dy,1) - ufab(m-dy,1))/(2.0*DX[1]);
					epsfab(m,5) = (ufab(m+dz,1) - ufab(m-dz,1))/(2.0*DX[2]);
					epsfab(m,6) = (ufab(m+dx,2) - ufab(m-dx,2))/(2.0*DX[0]);
					epsfab(m,7) = (ufab(m+dy,2) - ufab(m-dy,2))/(2.0*DX[1]);
					epsfab(m,8) = (ufab(m+dz,2) - ufab(m-dz,2))/(2.0*DX[2]);
					Set::Scalar trace = (sigmafab(m,0) + sigmafab(m,4) + sigmafab(m,8))/3.0;
					sigmavmfab(m) = sqrt(
								1.5*(
									(sigmafab(m,0)-trace)*(sigmafab(m,0)-trace) +
									(sigmafab(m,1))*(sigmafab(m,1)) +
									(sigmafab(m,2))*(sigmafab(m,2)) +
									(sigmafab(m,3))*(sigmafab(m,3)) +
									((sigmafab(m,4)-trace)*(sigmafab(m,4)-trace)) +
									(sigmafab(m,5))*(sigmafab(m,5)) +
									(sigmafab(m,6))*(sigmafab(m,6)) +
									(sigmafab(m,7))*(sigmafab(m,7)) +
									((sigmafab(m,8)-trace)*(sigmafab(m,8)-trace))
								)
							);

#endif

				}

			//FArrayBox &energyfab  = (*energy[lev])[mfi];
			//elastic_operator->Energy(energyfab,ufab,lev,mfi);

			//FArrayBox &energiesfab  = (*energies[lev])[mfi];
			//elastic_operator->Energies(energiesfab,ufab,lev,0,mfi);
		}
	}
	*/
}
}
//#endif
