#include "PolymerDegradation.H"
#include "AMReX_MLCGSolver.H"
//#if AMREX_SPACEDIM == 1
namespace Integrator
{
PolymerDegradation::PolymerDegradation():
	Integrator()
{
	//
	// READ INPUT PARAMETERS
	//

	// ---------------------------------------------------------------------
	// --------------------- Water diffusion -------------------------------
	// ---------------------------------------------------------------------
	amrex::ParmParse pp_water("water");
	pp_water.query("on",water_diffusion_on);
	if(water_diffusion_on)
	{
		pp_water.query("diffusivity", water_diffusivity);
		pp_water.query("refinement_threshold", water_refinement_threshold);
		pp_water.query("ic_type", water_ic_type);

		// // Determine initial condition
		if (water_ic_type == "constant")
		{
			amrex::ParmParse pp_water_ic("water.ic");
			amrex::Vector<amrex::Real> value;
			pp_water_ic.queryarr("value",value);
			water_ic = new IC::Constant(geom,value);
		}
		else
			Util::Abort(INFO, "This kind of IC has not been implemented yet");

		amrex::ParmParse pp_water_bc("water.bc");

		amrex::Vector<std::string> bc_hi_str(AMREX_SPACEDIM);
		amrex::Vector<std::string> bc_lo_str(AMREX_SPACEDIM);

		pp_water_bc.queryarr("lo",bc_lo_str,0,AMREX_SPACEDIM);
		pp_water_bc.queryarr("hi",bc_hi_str,0,AMREX_SPACEDIM);
		amrex::Vector<amrex::Real> bc_lo_1, bc_hi_1;
		amrex::Vector<amrex::Real> bc_lo_2, bc_hi_2;
		amrex::Vector<amrex::Real> bc_lo_3, bc_hi_3;

		if (pp_water_bc.countval("lo_1")) pp_water_bc.getarr("lo_1",bc_lo_1);
		if (pp_water_bc.countval("hi_1")) pp_water_bc.getarr("hi_1",bc_hi_1);
#if AMREX_SPACEDIM>1
		if (pp_water_bc.countval("lo_2")) pp_water_bc.getarr("lo_2",bc_lo_2);
		if (pp_water_bc.countval("hi_2")) pp_water_bc.getarr("hi_2",bc_hi_2);

#if AMREX_SPACEDIM>2
		if (pp_water_bc.countval("lo_3")) pp_water_bc.getarr("lo_3",bc_lo_3);
		if (pp_water_bc.countval("hi_3")) pp_water_bc.getarr("hi_3",bc_hi_3);
#endif
#endif
		water_bc = new BC::Constant(bc_hi_str, bc_lo_str
							  ,AMREX_D_DECL(bc_lo_1, bc_lo_2, bc_lo_3)
							  ,AMREX_D_DECL(bc_hi_1, bc_hi_2, bc_hi_3)
							  );

		RegisterNewFab(water_conc,     water_bc, 1, number_of_ghost_cells, "Water Concentration");
		RegisterNewFab(water_conc_old, water_bc, 1, number_of_ghost_cells, "Water Concentration Old");
	}
	// ---------------------------------------------------------------------
	// --------------------- Heat diffusion -------------------------------
	// ---------------------------------------------------------------------
	amrex::ParmParse pp_heat("thermal");
	pp_heat.query("on",heat_diffusion_on);
	if(heat_diffusion_on)
	{
		pp_heat.query("diffusivity", thermal_diffusivity);
		pp_heat.query("refinement_threshold",thermal_refinement_threshold);
		pp_heat.query("ic_type",thermal_ic_type);

		if (thermal_ic_type == "constant")
		{
			amrex::ParmParse pp_heat_ic("thermal.ic");
			amrex::Vector<amrex::Real> T;
			pp_heat_ic.queryarr("value",T);
			thermal_ic = new IC::Constant(geom,T);
		}
		else
			Util::Abort(INFO, "This kind of IC has not been implemented yet");

		amrex::ParmParse pp_heat_bc("thermal.bc");

		amrex::Vector<std::string> bc_hi_str(AMREX_SPACEDIM);
		amrex::Vector<std::string> bc_lo_str(AMREX_SPACEDIM);

		pp_heat_bc.queryarr("lo",bc_lo_str,0,BL_SPACEDIM);
		pp_heat_bc.queryarr("hi",bc_hi_str,0,BL_SPACEDIM);
		amrex::Vector<amrex::Real> bc_lo_1, bc_hi_1;
		amrex::Vector<amrex::Real> bc_lo_2, bc_hi_2;
		amrex::Vector<amrex::Real> bc_lo_3, bc_hi_3;

		if (pp_heat_bc.countval("lo_1")) pp_heat_bc.getarr("lo_1",bc_lo_1);
		if (pp_heat_bc.countval("hi_1")) pp_heat_bc.getarr("hi_1",bc_hi_1);

		if (pp_heat_bc.countval("lo_2")) pp_heat_bc.getarr("lo_2",bc_lo_2);
		if (pp_heat_bc.countval("hi_2")) pp_heat_bc.getarr("hi_2",bc_hi_2);

#if AMREX_SPACEDIM>2
		if (pp_heat_bc.countval("lo_3")) pp_heat_bc.getarr("lo_3",bc_lo_3);
		if (pp_heat_bc.countval("hi_3")) pp_heat_bc.getarr("hi_3",bc_hi_3);
#endif

		thermal_bc = new BC::Constant(bc_hi_str, bc_lo_str
								,AMREX_D_DECL(bc_lo_1, bc_lo_2, bc_lo_3)
							  	,AMREX_D_DECL(bc_hi_1, bc_hi_2, bc_hi_3)
						);
		RegisterNewFab(Temp,     thermal_bc, 1, number_of_ghost_cells, "Temperature");
		RegisterNewFab(Temp_old, thermal_bc, 1, number_of_ghost_cells, "Temperature Old");
	}
	// ---------------------------------------------------------------------
	// --------------------- Material model --------------------------------
	// ---------------------------------------------------------------------
	amrex::ParmParse pp_material("material");
	pp_material.query("model",input_material);
	if(input_material == "isotropic")
	{
		amrex::Real lambda = 1.0;
		amrex::Real mu = 1.0;
		amrex::ParmParse pp_material_isotropic("material.isotropic");
		pp_material_isotropic.query("lambda",lambda);
		pp_material_isotropic.query("mu",mu);
		if(lambda <=0)
		{
			std::cout<< "Warning. Lambda must be positive. Resetting back to default value" <<std::endl;
			lambda = 1.0;
		}
		if(mu <= 0)
		{
			std::cout<< "Warning. Mu must be positive. Resetting back to default value" <<std::endl;
			mu = 1.0;
		}
		modeltype = model_type(lambda,mu);
		//modeltype = model_type(1.0);
		//Model::Solid::Elastic::Degradable::Isotropic modeltype(lambda,mu);
		//modeltype(lambda,mu);
	}
	else if(input_material == "cubic")
	{
		Util::Abort(INFO, "Not implemented yet");
		/*Set::Scalar C11 = 1.0;
		Set::Scalar C12 = 1.0;
		Set::Scalar C44 = 1.0;
		amrex::ParmParse pp_material_isotropic("material.cubic");
		pp_material_isotropic.query("C11",C11);
		pp_material_isotropic.query("C12",C12);
		pp_material_isotropic.query("C44",C44);
		if(C11 <= 0.0)
		{
			std::cout<< "Warning. C11 must be positive. Resetting back to default value" <<std::endl;
			C11 = 1.0;
		}
		if(C12 <= 0.0)
		{
			std::cout<< "Warning. C12 must be positive. Resetting back to default value" <<std::endl;
			C12 = 1.0;
		}
		if(C44 <= 0)
		{
			std::cout<< "Warning. C44 must be positive. Resetting back to default value" <<std::endl;
			C44 = 1.0;
		}
		models.push_back(new Operator::Elastic::Degradation::Cubic(C11,C12,C44,0.0,0.0,0.0));*/
	}
	// ---------------------------------------------------------------------
	// --------------------- Damage model ----------------------------------
	// ---------------------------------------------------------------------
	amrex::ParmParse pp_damage("damage"); // Phase-field model parameters
	pp_damage.query("anisotropy",damage_anisotropy);

	if(damage_anisotropy == 0)
		number_of_eta = 1;
	else
		number_of_eta = AMREX_SPACEDIM;

	pp_damage.query("type",damage_type);
	if(damage_type == "relaxation") //meant for isotropic right now.
	{
		number_of_eta = 1;
		damage_anisotropy = 0;
		pp_damage.query("d_final",d_final);
		pp_damage.query("number_of_terms",number_of_terms);
		pp_damage.queryarr("d_i",d_i);
		pp_damage.queryarr("tau_i",tau_i);
		pp_damage.queryarr("t_start_i",t_start_i);

		if(d_final > 1.0)
		{
			std::cout << "Warning: d_final can not be greater than 1. Resetting it to default" << std::endl;
			d_final = 1.0;
		}

		if(d_i.size() != number_of_terms || tau_i.size() != number_of_terms || t_start_i.size() != number_of_terms)
			Util::Abort(INFO, "missing entries in d_i, tau_i or t_start_i");

		amrex::Real sum = 0;
		for (int temp = 0; temp < d_i.size(); temp++)
		{
			if(d_i[temp] < 0.0 || d_i[temp] > 1.0)
			 	Util::Abort(INFO, "Invalid values for d_i. Must be between 0 and 1");

			sum += d_i[temp];
		}

		if(sum != d_final) //need to replace this in the future
			Util::Abort(INFO, "d_final is not equal to the sum of d_i");
	}
	else
		Util::Abort(INFO, "This kind of damage model has not been implemented yet");

	pp_damage.query("ic_type",eta_ic_type);
	pp_damage.query("refinement_threshold",damage_refinement_threshold);
	if(eta_ic_type == "constant")
	{
		amrex::ParmParse pp_damage_ic("damage.ic");
		amrex::Vector<amrex::Real> eta_init;
		pp_damage_ic.queryarr("value",eta_init);
		eta_ic = new IC::Constant(geom,eta_init);
	}
	else
		Util::Abort(INFO, "This kind of IC has not been implemented yet");

	amrex::ParmParse pp_damage_bc("damage.bc");

	amrex::Vector<std::string> bc_hi_str(AMREX_SPACEDIM);
	amrex::Vector<std::string> bc_lo_str(AMREX_SPACEDIM);

	pp_damage_bc.queryarr("lo",bc_lo_str,0,BL_SPACEDIM);
	pp_damage_bc.queryarr("hi",bc_hi_str,0,BL_SPACEDIM);
	amrex::Vector<amrex::Real> bc_lo_1, bc_hi_1;
	amrex::Vector<amrex::Real> bc_lo_2, bc_hi_2;
	amrex::Vector<amrex::Real> bc_lo_3, bc_hi_3;

	if (pp_damage_bc.countval("lo_1")) pp_damage_bc.getarr("lo_1",bc_lo_1);
	if (pp_damage_bc.countval("hi_1")) pp_damage_bc.getarr("hi_1",bc_hi_1);

	if (pp_damage_bc.countval("lo_2")) pp_damage_bc.getarr("lo_2",bc_lo_2);
	if (pp_damage_bc.countval("hi_2")) pp_damage_bc.getarr("hi_2",bc_hi_2);

	if (pp_damage_bc.countval("lo_3")) pp_damage_bc.getarr("lo_3",bc_lo_3);
	if (pp_damage_bc.countval("hi_3")) pp_damage_bc.getarr("hi_3",bc_hi_3);

	eta_bc = new BC::Constant(bc_hi_str, bc_lo_str
							,AMREX_D_DECL(bc_lo_1, bc_lo_2, bc_lo_3)
							,AMREX_D_DECL(bc_hi_1, bc_hi_2, bc_hi_3)
					);

	RegisterNewFab(eta_new, eta_bc, number_of_eta, number_of_ghost_cells, "Eta");
	RegisterNewFab(eta_old, eta_bc, number_of_eta, number_of_ghost_cells, "Eta old");

	//std::cout << __FILE__ << ": " << __LINE__ << std::endl;
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

		pp_elastic.query("tstart",elastic_tstart);
		if(elastic_tstart < 0.0)
		{
			std::cout << "Warning: Invalid value for elasitc t_start. Setting it to zero" << std::endl;
			elastic_tstart = 0.0;
		}
		else if(elastic_tstart > stop_time)
		{
			std::cout << "Warning: Invalid value for elastic t_start. Setting it to stop_time" <<std::endl;
			elastic_tstart = stop_time;
		}

		pp_elastic.query("tend",elastic_tend);
		if(elastic_tend < elastic_tstart || elastic_tend > stop_time)
		{
			std::cout << "Warning: Invalid value for elastic t_end. Setting it to stop_time" << std::endl;
			elastic_tend = stop_time;
			if(elastic_tstart == stop_time) elastic_tstart = stop_time - timestep;

		}

		pp_elastic.query("bottom_solver",bottom_solver);
		pp_elastic.query("linop_maxorder", linop_maxorder);
		pp_elastic.query("max_coarsening_level",max_coarsening_level);
		pp_elastic.query("verbose",elastic_verbose);
		pp_elastic.query("cg_verbose", elastic_cgverbose);
		pp_elastic.query("bottom_max_iter", elastic_bottom_max_iter);
		pp_elastic.query("max_fixed_iter", elastic_max_fixed_iter);
		pp_elastic.query("bottom_tol", elastic_bottom_tol);

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

		bc_map["displacement"] = Operator::Elastic<model_type>::BC::Displacement;
		bc_map["disp"] = Operator::Elastic<model_type>::BC::Displacement;
		bc_map["traction"] = Operator::Elastic<model_type>::BC::Traction;
		bc_map["trac"] = Operator::Elastic<model_type>::BC::Traction;
		bc_map["periodic"] = Operator::Elastic<model_type>::BC::Periodic;

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

#if AMREX_SPACEDIM > 1
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

#if AMREX_SPACEDIM > 1
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
		ndmap.resize(nlevels);

		const int number_of_stress_components = AMREX_SPACEDIM*AMREX_SPACEDIM;
		//displacement.resize(nlevels);
		RegisterNodalFab (displacement,AMREX_SPACEDIM,2,"displacement");;
		RegisterNodalFab (rhs,AMREX_SPACEDIM,2,"rhs");;
		RegisterNodalFab (strain,number_of_stress_components,2,"strain");;
		RegisterNodalFab (stress,number_of_stress_components,2,"stress");;
		RegisterNodalFab (stress_vm,1,2,"stress_vm");;
		RegisterNodalFab (energy,1,2,"energy");;
		RegisterNodalFab (residual,AMREX_SPACEDIM,2,"residual");;

	}
}


#define ETA(i,j,k,n) eta_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n)

void
PolymerDegradation::Advance (int lev, amrex::Real time, amrex::Real dt)
{
	std::swap(*eta_old[lev], *eta_new[lev]);

	if(water_diffusion_on) std::swap(*water_conc_old[lev],*water_conc[lev]);

	if(heat_diffusion_on) std::swap(*Temp_old[lev], *Temp[lev]);

	static amrex::IntVect AMREX_D_DECL(	dx(AMREX_D_DECL(1,0,0)),
						dy(AMREX_D_DECL(0,1,0)),
						dz(AMREX_D_DECL(0,0,1)));
	const amrex::Real* DX = geom[lev].CellSize();

	if(water_diffusion_on)
	{
		for ( amrex::MFIter mfi(*water_conc[lev],true); mfi.isValid(); ++mfi )
		{
			const amrex::Box& bx = mfi.tilebox();

			amrex::FArrayBox &water_conc_old_box = (*water_conc_old[lev])[mfi];
			amrex::FArrayBox &water_conc_box = (*water_conc[lev])[mfi];

			AMREX_D_TERM(	for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++),
							for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++),
							for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++)
						)
				{
					amrex::IntVect m(AMREX_D_DECL(i,j,k));
					if(std::isnan(water_conc_old_box(m))) Util::Abort(INFO, "Nan found in WATER_OLD(i,j,k)");
					if(std::isinf(water_conc_old_box(m))) Util::Abort(INFO, "Inf found in WATER_OLD(i,j,k)");
					if(water_conc_old_box(m) > 1.0)
					{
						std::cout << "dt = " << dt << ", time = " << time << ", lev = " << lev << std::endl;
						Util::Abort(INFO, "water concentration exceeded 2");
					}
					water_conc_box(m) =
						water_conc_old_box(m)
						+ dt * water_diffusivity * (AMREX_D_TERM((water_conc_old_box(m+dx) + water_conc_old_box(m-dx) - 2.0*water_conc_old_box(m)) / DX[0] / DX[0],
						+ (water_conc_old_box(m+dy) + water_conc_old_box(m-dy) - 2.0*water_conc_old_box(m)) / DX[1] / DX[1],
						+ (water_conc_old_box(m+dz) + water_conc_old_box(m-dz) - 2.0*water_conc_old_box(m)) / DX[2] / DX[2]));
					if(water_conc_box(m) > 1.0)
					{
						std::cout << "Warning: water concentration has exceeded one after computation. Resetting it to one" << std::endl;
						water_conc_box(m) = 1.0;
					}
				}
		}
	}

	if(heat_diffusion_on)
	{
		for ( amrex::MFIter mfi(*Temp[lev],true); mfi.isValid(); ++mfi )
		{
			const amrex::Box& bx = mfi.tilebox();

			amrex::FArrayBox &Temp_old_box = (*Temp_old[lev])[mfi];
			amrex::FArrayBox &Temp_box = (*Temp[lev])[mfi];

			AMREX_D_TERM(	for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++),
							for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++),
							for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++)
						)
				{
					amrex::IntVect m(AMREX_D_DECL(i,j,k));
					Temp_box(m) =
						Temp_old_box(m)
						+ dt * thermal_diffusivity * (AMREX_D_TERM((Temp_old_box(m+dx) + Temp_old_box(m-dx) - 2.0*Temp_old_box(m)) / DX[0] / DX[0],
						+ (Temp_old_box(m+dy) + Temp_old_box(m-dy) - 2.0*Temp_old_box(m)) / DX[1] / DX[1],
						+ (Temp_old_box(m+dz) + Temp_old_box(m-dz) - 2.0*Temp_old_box(m)) / DX[2] / DX[2]));
				}
		}
	}

	for ( amrex::MFIter mfi(*eta_new[lev],true); mfi.isValid(); ++mfi )
	{
		const amrex::Box& bx = mfi.tilebox();
		FArrayBox &eta_new_box     = (*eta_new[lev])[mfi];
		FArrayBox &eta_old_box     = (*eta_old[lev])[mfi];
		amrex::FArrayBox &water_conc_box = (*water_conc[lev])[mfi];

		//FArrayBox &energiesfab = (*energy)

		AMREX_D_TERM(	for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++),
						for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++),
						for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++)
					)
			{
				amrex::IntVect m(AMREX_D_DECL(i,j,k));
				for (int n = 0; n < number_of_eta; n++)
				{
					if(damage_type == "relaxation") // Need to replace this with more sophisticated check
					{
						amrex::Real rhs = 0.0;
						if (water_conc_box(m) >0 && eta_old_box(m,n) < d_final)
						{
							for (int l = 0; l<number_of_terms; l++)
								rhs += d_i[l]*water_conc_box(m)*std::exp(-std::max(0.0,time-t_start_i[l])/tau_i[l])/(tau_i[l]);
						}
						eta_new_box(m,n) = eta_old_box(m,n) + rhs*dt;
						if(eta_new_box(m,n) > d_final)
						{
							Util::Abort(INFO, "eta exceeded ",d_final, ". Rhs = ", rhs, ", Water = ", water_conc_box(m,n));
						}
					}

					//
					// ELASTIC DRIVING FORCE
					//
					//if (elastic_on)
					//{
					//	FArrayBox &energiesfab     = (*energies[lev])[mfi];

					//	eta_new_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),m)
					//		+= M*dt*( energiesfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),m));
					//}

				}
			}
	}
}

void
PolymerDegradation::Initialize (int lev)
{
	if(water_diffusion_on)
	{
		water_ic->Initialize(lev,water_conc);
		water_ic->Initialize(lev,water_conc_old);
	}

	if(heat_diffusion_on)
	{
		thermal_ic->Initialize(lev,Temp);
		thermal_ic->Initialize(lev,Temp_old);
	}


	eta_ic->Initialize(lev,eta_new);
	eta_ic->Initialize(lev,eta_old);
	//Util::Message(INFO,"Grids size = ",grids[0].size());
	//Util::Abort(INFO);
}


void
PolymerDegradation::TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/)
{
	const amrex::Real* dx      = geom[lev].CellSize();

	amrex::Vector<int>  itags;
	if(water_diffusion_on)
	{
		for (amrex::MFIter mfi(*water_conc[lev],true); mfi.isValid(); ++mfi)
		{
			const amrex::Box&  bx  = mfi.tilebox();
			amrex::TagBox&     tag  = tags[mfi];

			amrex::FArrayBox &water_conc_box = (*water_conc[lev])[mfi];

			AMREX_D_TERM(	for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++),
							for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++),
							for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++)
						)
				{
					AMREX_D_TERM(	amrex::Real grad1 = (WATER(i+1,j,k) - WATER(i-1,j,k))/(2.*dx[0]);,
									amrex::Real grad2 = (WATER(i,j+1,k) - WATER(i,j-1,k))/(2.*dx[1]);,
									amrex::Real grad3 = (WATER(i,j,k+1) - WATER(i,j,k-1))/(2.*dx[2]););
					amrex::Real grad = sqrt(AMREX_D_TERM(grad1*grad1,
						+ grad2*grad2,
						+ grad3*grad3));

					amrex::Real dr = sqrt(AMREX_D_TERM(dx[0]*dx[0],
						+ dx[1]*dx[1],
						+ dx[2]*dx[2]));

					if (grad*dr > water_refinement_threshold)
						tag(amrex::IntVect(AMREX_D_DECL(i,j,k))) = amrex::TagBox::SET;
				}
		}
	}

	if(heat_diffusion_on)
	{
		for (amrex::MFIter mfi(*Temp[lev],true); mfi.isValid(); ++mfi)
		{
			const amrex::Box&  bx  = mfi.tilebox();
			amrex::TagBox&     tag  = tags[mfi];

			amrex::FArrayBox &Temp_box = (*Temp[lev])[mfi];

			AMREX_D_TERM(	for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++),
							for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++),
							for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++)
						)
				{
					AMREX_D_TERM(	amrex::Real grad1 = (TEMP(i+1,j,k) - TEMP(i-1,j,k))/(2.*dx[0]);,
									amrex::Real grad2 = (TEMP(i,j+1,k) - TEMP(i,j-1,k))/(2.*dx[1]);,
									amrex::Real grad3 = (TEMP(i,j,k+1) - TEMP(i,j,k-1))/(2.*dx[2]);
								)
					amrex::Real grad = sqrt(AMREX_D_TERM(grad1*grad1,
						+ grad2*grad2,
						+ grad3*grad3));

					amrex::Real dr = sqrt(AMREX_D_TERM(dx[0]*dx[0],
						+ dx[1]*dx[1],
						+ dx[2]*dx[2]));

					if (grad*dr > thermal_refinement_threshold)
						tag(amrex::IntVect(AMREX_D_DECL(i,j,k))) = amrex::TagBox::SET;
				}
		}
	}


	for (amrex::MFIter mfi(*eta_new[lev],true); mfi.isValid(); ++mfi)
	{
		const amrex::Box&  bx  = mfi.tilebox();
		amrex::TagBox&     tag  = tags[mfi];
		amrex::BaseFab<amrex::Real> &eta_new_box = (*eta_new[lev])[mfi];

		AMREX_D_TERM(		for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++),
							for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++),
							for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++)
					)
			{
				for (int m = 0; m < number_of_eta; m++)
				{
					AMREX_D_TERM(	Set::Scalar gradx = (ETA_NEW(i+1,j,k,m) - ETA_NEW(i-1,j,k,m))/(2.*dx[0]);,
									Set::Scalar grady = (ETA_NEW(i,j+1,k,m) - ETA_NEW(i,j-1,k,m))/(2.*dx[1]);,
									Set::Scalar gradz = (ETA_NEW(i,j,k+1,m) - ETA_NEW(i,j,k-1,m))/(2.*dx[2]);
					)
					Set::Scalar grad = sqrt(AMREX_D_TERM(gradx*gradx, + grady*grady, + gradz*gradz));
					Set::Scalar dr = sqrt(AMREX_D_TERM(dx[0]*dx[0], + dx[1]*dx[1], + dx[2]*dx[2]));

					if(grad*dr > damage_refinement_threshold)
						tag(amrex::IntVect(AMREX_D_DECL(i,j,k))) = amrex::TagBox::SET;
				}
			}

	}
}

void
PolymerDegradation::DegradeMaterial(int lev)
{
	/*
	This function is supposed to degrade material parameters based on certain
	damage model.
	For now we are just using isotropic degradation.
	*/
	if(damage_anisotropy)
		Util::Abort(__FILE__,"DegradeModulus",__LINE__,"Not implemented yet");

	static std::array<amrex::IntVect,AMREX_SPACEDIM> dx = {AMREX_D_DECL(amrex::IntVect(AMREX_D_DECL(1,0,0)),
																		amrex::IntVect(AMREX_D_DECL(0,1,0)),
																		amrex::IntVect(AMREX_D_DECL(0,0,1)))};

	//bool isMFIterSafe  = (model[lev].DistributionMap() == (*eta_new[lev]).DistributionMap());
	//Util::Message(INFO, "isMFIterSafe = ",isMFIterSafe);

	for (amrex::MFIter mfi(model[lev],true); mfi.isValid(); ++mfi)
	{
		const amrex::Box& box = mfi.tilebox();
	 	amrex::BaseFab<model_type> &modelfab = (model[lev])[mfi];
		amrex::BaseFab<amrex::Real> &etafab = (*eta_new[lev])[mfi];

	 	AMREX_D_TERM(for (int i = box.loVect()[0]; i<=box.hiVect()[0]; i++),
		 		     for (int j = box.loVect()[1]; j<=box.hiVect()[1]; j++),
		 		     for (int k = box.loVect()[2]; k<=box.hiVect()[2]; k++))
	 	{
			amrex::IntVect m(AMREX_D_DECL(i,j,k));
			Set::Scalar mul = 1.0/(AMREX_D_TERM(2.0,+2.0,+4.0));
			Set::Scalar temp = mul*(AMREX_D_TERM(	etafab(m) 	+ etafab(m-dx[0])
													,
													+ etafab(m-dx[1]) + etafab(m-dx[0]-dx[1])
													,
													+ etafab(m-dx[2])	+ etafab(m-dx[0]-dx[2])
													+ etafab(m-dx[1]-dx[2]) + etafab(m-dx[0]-dx[1]-dx[2])
												));
			//if (false || AMREX_D_TERM( i == geom[lev].Domain().loVect()[0] || i == geom[lev].Domain().hiVect()[0]+1,
			//						|| j == geom[lev].Domain().loVect()[1] || j == geom[lev].Domain().hiVect()[1]+1,
			//						|| k == geom[lev].Domain().loVect()[2] || k == geom[lev].Domain().hiVect()[2]+1))
			//	continue;
			if(temp > d_final || std::isnan(temp) || std::isinf(temp))
			{
				Util::Message(INFO,"Invalid value of temp = ", temp);
				Util::Message(INFO," mul = ", mul);
				Util::Message(INFO,"etafab(m) = ", etafab(m));
				Util::Message(INFO,"etafab(m-dx) = ", etafab(m-dx[0]));
				Util::Message(INFO,"etafab(m-dy) = ", etafab(m-dx[1]));
				Util::Message(INFO,"etafab(m-dz) = ", etafab(m-dx[2]));
				Util::Message(INFO,"etafab(m-dx-dy) = ", etafab(m-dx[0]-dx[1]));
				Util::Message(INFO,"etafab(m-dx-dz) = ", etafab(m-dx[0]-dx[2]));
				Util::Message(INFO,"etafab(m-dy-dz) = ", etafab(m-dx[2]-dx[1]));
				Util::Message(INFO,"etafab(m-dx-dy-dz) = ", etafab(m-dx[0]-dx[1]-dx[2]));
			}
			modelfab(m).DegradeModulus(temp);
		}
		/*AMREX_D_TERM(for (int i = box.loVect()[0]-1; i<=box.hiVect()[0]+1; i++),
		 		     for (int j = box.loVect()[1]-1; j<=box.hiVect()[1]+1; j++),
		 		     for (int k = box.loVect()[2]-1; k<=box.hiVect()[2]+1; k++))
	 	{
	 		amrex::IntVect m(AMREX_D_DECL(i,j,k));
	 		if(i == box.loVect()[0]-1)
	 			modelfab(m) = modelfab(m+dx[0]);
	 		if(i == box.hiVect()[0]+1)
	 			modelfab(m) = modelfab(m-dx[0]);
	 		#if AMREX_SPACEDIM > 1
	 		if(j == box.loVect()[1]-1)
	 			modelfab(m) = modelfab(m+dx[1]);
	 		if(j == box.hiVect()[1]+1)
	 			modelfab(m) = modelfab(m-dx[1]);
	 		#if AMREX_SPACEDIM > 2
	 		if(k == box.loVect()[2]-1)
	 			modelfab(m) = modelfab(m+dx[2]);
	 		if(k == box.hiVect()[2]+1)
	 			modelfab(m) = modelfab(m-dx[2]);
	 		#endif
	 		#endif
	 	}*/
	}
}

std::vector<std::string>
PolymerDegradation::PlotFileNameNode (std::string plot_file_name, int lev) const
{
	std::vector<std::string> name;
	name.push_back(plot_file_name+"/");
	name.push_back(amrex::Concatenate("", lev, 5));
	return name;
}

void 
PolymerDegradation::TimeStepComplete(amrex::Real time, int iter)
{
	if (! elastic_on) return;
	if (iter % elastic_int) return;
	if (time < elastic_tstart) return;
	if (time > elastic_tend) return;

}

void 
PolymerDegradation::TimeStepBegin(amrex::Real time, int iter)
{
	if (!elastic_on) return;
	if (iter%elastic_int) return;
	if (time < elastic_tstart) return;
	if (time > elastic_tend) return;

	int number_of_stress_components = AMREX_SPACEDIM*AMREX_SPACEDIM;
	int number_of_components = AMREX_SPACEDIM;
	int number_of_ghost_cells = 2;

	if (iter==0 || time == elastic_tstart)
	{
		model.resize(nlevels);
		for (int ilev = 0; ilev < nlevels; ++ilev)
		{
			model[ilev].define(displacement[ilev]->boxArray(), displacement[ilev]->DistributionMap(), 1, number_of_ghost_cells);
		}
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
		const Real* DX = geom[ilev].CellSize();
		Set::Scalar volume = AMREX_D_TERM(DX[0],*DX[1],*DX[2]);

		Util::Message(INFO);
		if (iter==0 || time == elastic_tstart)
		{
			displacement[ilev]->setVal(0.0);
			strain[ilev]->setVal(0.0);
			stress[ilev]->setVal(0.0);
			stress_vm[ilev]->setVal(0.0);
			rhs[ilev]->setVal(0.0);
			energy[ilev]->setVal(0.0);
			residual[ilev]->setVal(0.0);
		}
		model[ilev].setVal(modeltype);

		AMREX_D_TERM(rhs[ilev]->setVal(body_force[0]*volume,0,1);,
					rhs[ilev]->setVal(body_force[1]*volume,1,1);,
					rhs[ilev]->setVal(body_force[2]*volume,2,1););

		//Util::Message(INFO);
		DegradeMaterial(ilev);
		//Util::Message(INFO);
		elastic_operator.SetModel(ilev,model[ilev]);
		//Util::Message(INFO);

		for (amrex::MFIter mfi(*rhs[ilev],true); mfi.isValid(); ++mfi)
		{
		 	const amrex::Box& box = mfi.tilebox();

		 	amrex::BaseFab<amrex::Real> &rhsfab = (*rhs[ilev])[mfi];

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

		 		if (	false
						|| xmin || xmax
#if AMREX_SPACEDIM > 1
						|| ymin || ymax
#if AMREX_SPACEDIM > 2
						|| zmin || zmax
#endif
#endif
						)
		 		{
		 			AMREX_D_TERM(	rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),0) = 0.0;,
		 							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),1) = 0.0;,
		 							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),2) = 0.0;);
		 		}
				for(int l = 0; l<AMREX_SPACEDIM; l++)
				{
					if(xmin && bc_x_lo[l]==Operator::Elastic<model_type>::BC::Displacement)
					{
						if(elastic_bc_left.size() == 1)
							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = elastic_bc_left[0](l);
						else	
							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = interpolate_left(time)[l];
					}
					if(xmax && bc_x_hi[l]==Operator::Elastic<model_type>::BC::Displacement)
					{
						if(elastic_bc_right.size() == 1)
							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = elastic_bc_right[0](l);
						else	
							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = interpolate_right(time)[l];
					}
#if AMREX_SPACEDIM > 1
					if(ymin && bc_y_lo[l]==Operator::Elastic<model_type>::BC::Displacement)
					{
						if(elastic_bc_bottom.size() == 1)
							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = elastic_bc_bottom[0](l);
						else	
							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = interpolate_bottom(time)[l];
					}
					if(ymax && bc_y_hi[l]==Operator::Elastic<model_type>::BC::Displacement)
					{
						if(elastic_bc_top.size() == 1)
							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = elastic_bc_top[0](l);
						else	
							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = interpolate_top(time)[l];
					}
#if AMREX_SPACEDIM > 2
					if(zmin && bc_z_lo[l]==Operator::Elastic<model_type>::BC::Displacement)
					{
						if(elastic_bc_back.size() == 1)
							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = elastic_bc_back[0](l);
						else	
							rhsfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),l) = interpolate_back(time)[l];
					}
					if(zmax && bc_z_hi[l]==Operator::Elastic<model_type>::BC::Displacement)
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
	Util::Message(INFO);
	amrex::MLMG solver(elastic_operator);
	solver.setMaxIter(elastic_max_iter);
	solver.setMaxFmgIter(elastic_max_fmg_iter);
	solver.setFixedIter(elastic_max_fixed_iter);
	solver.setVerbose(elastic_verbose);
	solver.setCGVerbose(elastic_cgverbose);
	solver.setBottomMaxIter(elastic_bottom_max_iter);
	solver.setBottomTolerance(elastic_bottom_tol) ;
	if (bottom_solver == "cg")
		solver.setBottomSolver(MLMG::BottomSolver::cg);
	else if (bottom_solver == "bicgstab")
		solver.setBottomSolver(MLMG::BottomSolver::bicgstab);
	if (!elastic_use_fsmooth)// <<< put in to NOT require FSmooth
	{
		solver.setFinalSmooth(0);
		solver.setBottomSmooth(0);
	}


	amrex::MLCGSolver mlcg(&solver,elastic_operator);
	elastic_operator.PrepareForSolve();
	mlcg.setVerbose(4);
	mlcg.setMaxIter(elastic_bottom_max_iter);
	int ret = mlcg.solve(*displacement[0],*rhs[0],elastic_tol_rel,elastic_tol_abs);
	if (ParallelDescriptor::IOProcessor() &&  ret) Util::Abort(INFO,"Solver did not converge");

	// solver.solve(GetVecOfPtrs(displacement),
	// 	GetVecOfConstPtrs(rhs),
	// 	elastic_tol_rel,
	// 	elastic_tol_abs);
	// solver.compResidual(GetVecOfPtrs(residual),GetVecOfPtrs(displacement),GetVecOfConstPtrs(rhs));

	//solver.compResidual(GetVecOfPtrs(residual),GetVecOfPtrs(displacement),GetVecOfConstPtrs(rhs));

	for (int lev = 0; lev < nlevels; lev++)
	{
		elastic_operator.Strain(lev,*strain[lev],*displacement[lev]);
		elastic_operator.Stress(lev,*stress[lev],*displacement[lev]);
		elastic_operator.Energy(lev,*energy[lev],*displacement[lev]);
	}
}
}
//#endif
