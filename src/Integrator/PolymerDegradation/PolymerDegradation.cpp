#include "PolymerDegradation.H"
#if AMREX_SPACEDIM == 1
namespace Integrator
{
PolymerDegradation::PolymerDegradation():
	Integrator()
{

	//
	// READ INPUT PARAMETERS
	//
	amrex::ParmParse pp_water("water"); // Water diffusion parameters
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
		
		if (pp_water_bc.countval("lo_2")) pp_water_bc.getarr("lo_2",bc_lo_2);
		if (pp_water_bc.countval("hi_2")) pp_water_bc.getarr("hi_2",bc_hi_2);
		
#if AMREX_SPACEDIM>2
		if (pp_water_bc.countval("lo_3")) pp_water_bc.getarr("lo_3",bc_lo_3);
		if (pp_water_bc.countval("hi_3")) pp_water_bc.getarr("hi_3",bc_hi_3);
#endif

		water_bc = new BC::Constant(geom, bc_hi_str, bc_lo_str
							  ,bc_lo_1, bc_hi_1
							  ,bc_lo_2, bc_hi_2
#if AMREX_SPACEDIM > 2
							  ,bc_lo_3, bc_hi_3
#endif
							  );

		RegisterNewFab(water_conc,     water_bc, 1, number_of_ghost_cells, "Water Concentration");
		RegisterNewFab(water_conc_old, water_bc, 1, number_of_ghost_cells, "Water Concentration Old");
	}

	amrex::ParmParse pp_heat("thermal"); //Heat diffusion parameters
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

		thermal_bc = new BC::Constant(geom, bc_hi_str, bc_lo_str
								,bc_lo_1,bc_hi_1
								,bc_lo_2,bc_hi_2
#if AMREX_SPACEDIM > 2
								,bc_lo_3,bc_hi_3
#endif
					);
		RegisterNewFab(Temp,     thermal_bc, 1, number_of_ghost_cells, "Temperature");
		RegisterNewFab(Temp_old, thermal_bc, 1, number_of_ghost_cells, "Temperature Old");
	}
	//
	// Material model
	//
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
		models.push_back(new Model::Solid::Elastic::Degradable::Isotropic(lambda,mu));
	}
	else if(input_material == "cubic")
	{
		Util::Abort("Not implemented yet");
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

	//
	// Damage model
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
			amrex::Abort("missing entries in d_i, tau_i or t_start_i");

		amrex::Real sum = 0;		
		for (int temp = 0; temp < d_i.size(); temp++)
		{
			if(d_i[temp] < 0.0 || d_i[temp] > 1.0)
			 	amrex::Abort("Invalid values for d_i. Must be between 0 and 1");

			sum += d_i[temp];
		}

		if(sum != d_final) //need to replace this in the future
			amrex::Abort("d_final is not equal to the sum of d_i");
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

	eta_bc = new BC::Constant(geom, bc_hi_str, bc_lo_str
						,bc_lo_1, bc_hi_1
						,bc_lo_2, bc_hi_2
#if AMREX_SPACEDIM>2
						,bc_lo_3, bc_hi_3
#endif
						);

	RegisterNewFab(eta_new, eta_bc, number_of_eta, number_of_ghost_cells, "Eta");
	RegisterNewFab(eta_old, eta_bc, number_of_eta, number_of_ghost_cells, "Eta old");

  
	// Elasticity
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
		if (pp_elastic.countval("body_force")) pp_elastic.getarr("body_force",body_force);

		amrex::ParmParse pp_elastic_bc("elastic.bc");
		AMREX_D_TERM(	amrex::Vector<std::string> bc_x_lo_str(AMREX_SPACEDIM);
						amrex::Vector<std::string> bx_x_hi_str(AMREX_SPACEDIM);
						,
						amrex::Vector<std::string> bc_y_lo_str(AMREX_SPACEDIM);
						amrex::Vector<std::string> bx_y_hi_str(AMREX_SPACEDIM);
						,
						amrex::Vector<std::string> bc_z_lo_str(AMREX_SPACEDIM);
						amrex::Vector<std::string> bx_z_hi_str(AMREX_SPACEDIM););

		//amrex::Vector<std::string> bc_hi_str(AMREX_SPACEDIM);
		//amrex::Vector<std::string> bc_lo_str(AMREX_SPACEDIM);
		AMREX_D_TERM(	pp_elastic_bc.queryarr("bc_x_lo",bc_x_lo_str);,
						pp_elastic_bc.queryarr("bc_y_lo",bc_y_lo_str);,
						pp_elastic_bc.queryarr("bc_z_lo",bc_z_lo_str););
		AMREX_D_TERM(	pp_elastic_bc.queryarr("bc_x_hi",bc_x_hi_str);,
						pp_elastic_bc.queryarr("bc_y_hi",bc_y_hi_str);,
						pp_elastic_bc.queryarr("bc_z_hi",bc_z_hi_str););
		bc_map["displacement"] = Operator::Elastic::BC::Displacement;
		bc_map["disp"] = Operator::Elastic::BC::Displacement;
		bc_map["traction"] = Operator::Elastic::BC::Traction;
		bc_map["trac"] = Operator::Elastic::BC::Traction;
		bc_map["periodic"] = Operator::Elastic::BC::Periodic;
		//pp_elastic_bc.queryarr("type_lo",bc_lo_str,0,AMREX_SPACEDIM);
		//pp_elastic_bc.queryarr("type_hi",bc_hi_str,0,AMREX_SPACEDIM);

		AMREX_D_TERM(	bc_x_lo = {AMREX_D_DECL(bc[bc_x_lo_str[0]], bc[bc_x_lo_str[1]], bc[bc_x_lo_str[2]])};
						bc_x_hi = {AMREX_D_DECL(bc[bc_x_hi_str[0]], bc[bc_x_hi_str[1]], bc[bc_x_hi_str[2]])};
						,
						bc_y_lo = {AMREX_D_DECL(bc[bc_y_lo_str[0]], bc[bc_y_lo_str[1]], bc[bc_y_lo_str[2]])};
						bc_y_hi = {AMREX_D_DECL(bc[bc_y_hi_str[0]], bc[bc_y_hi_str[1]], bc[bc_y_hi_str[2]])};
						,
						bc_z_lo = {AMREX_D_DECL(bc[bc_z_lo_str[0]], bc[bc_z_lo_str[1]], bc[bc_z_lo_str[2]])};
						bc_z_hi = {AMREX_D_DECL(bc[bc_z_hi_str[0]], bc[bc_z_hi_str[1]], bc[bc_z_hi_str[2]])};);
	
		amrex::Vector<Set::Scalar> bc_lo_1, bc_hi_1;
		amrex::Vector<Set::Scalar> bc_lo_2, bc_hi_2;
		amrex::Vector<Set::Scalar> bc_lo_3, bc_hi_3;
		amrex::Vector<Set::Scalar> bc_lo_1_t, bc_hi_1_t;
		amrex::Vector<Set::Scalar> bc_lo_2_t, bc_hi_2_t;
		amrex::Vector<Set::Scalar> bc_lo_3_t, bc_hi_3_t;

		if (pp_elastic_bc.countval("left_face")) pp_elastic_bc.getarr("left_face",bc_lo_1);
		if(bc_lo_1.size() % AMREX_SPACEDIM !=0)
			amrex::Abort("Invalid number of values for left_face displacement");

		if (pp_elastic_bc.countval("right_face")) pp_elastic_bc.getarr("right_face",bc_hi_1);
		if(bc_hi_1.size() % AMREX_SPACEDIM !=0)
			amrex::Abort("Invalid number of values for right_face displacement");
	
		if (pp_elastic_bc.countval("bottom_face")) pp_elastic_bc.getarr("bottom_face",bc_lo_2);
		if(bc_lo_2.size() % AMREX_SPACEDIM !=0)
			amrex::Abort("Invalid number of values for bottom_face displacement");

		if (pp_elastic_bc.countval("top_face")) pp_elastic_bc.getarr("top_face",bc_hi_2);
		if(bc_hi_2.size() % AMREX_SPACEDIM !=0)
			amrex::Abort("Invalid number of values for top_face displacement");
		
#if AMREX_SPACEDIM>2
		if (pp_elastic_bc.countval("back_face")) pp_elastic_bc.getarr("back_face",bc_lo_3);
		if(bc_lo_3.size() % AMREX_SPACEDIM !=0)
			amrex::Abort("Invalid number of values for back_face displacement");

		if (pp_elastic_bc.countval("front_face")) pp_elastic_bc.getarr("front_face",bc_hi_3);
		if(bc_hi_3.size() % AMREX_SPACEDIM !=0)
			amrex::Abort("Invalid number of values for front_face displacement");
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
		/*elastic_bc = new BC::Changeable(geom, bc_hi_str, bc_lo_str
						,elastic_bc_left[0], elastic_bc_right[0]
						,elastic_bc_bottom[0], elastic_bc_top[0]
#if AMREX_SPACEDIM>2
						,elastic_bc_back[0], elastic_bc_front[0]
#endif
						);*/
		/*BC::BC *temp_elastic_bc;
		amrex::Vector<std::string> temp_str = {AMREX_D_DECL("dirichlet","dirichlet","dirichlet")
							,AMREX_D_DECL("dirichlet","dirichlet","dirichlet")
#if AMREX_SPACEDIM>2
							,AMREX_D_DECL("dirichlet","dirichlet","dirichlet")
#endif
							};
		amrex::Vector<Set::Scalar> temp_bc_val = {AMREX_D_DECL(0.0,0.0,0.0)
							,AMREX_D_DECL(0.0,0.0,0.0)
#if AMREX_SPACEDIM>2
							,AMREX_D_DECL(0.0,0.0,0.0)
#endif
							};
		temp_elastic_bc = new BC::Constant(geom,temp_str,temp_str
						,temp_bc_val,temp_bc_val
						,temp_bc_val,temp_bc_val
#if AMREX_SPACEDIM>2
						,temp_bc_val,temp_bc_val
#endif
						);*/

		int number_of_stress_components = AMREX_SPACEDIM*AMREX_SPACEDIM;
		RegisterNewFab(displacement, AMREX_SPACEDIM, "disp");
		RegisterNewFab(rhs, AMREX_SPACEDIM, "rhs");
		RegisterNewFab(strain, number_of_stress_components, "eps");
		RegisterNewFab(stress, number_of_stress_components, "sig");
		RegisterNewFab(stress_vm, 1, "sig_VM");
		RegisterNewFab(energy, 1, "W");
		RegisterNewFab(energies, number_of_eta, "W_1");
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

			for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
				for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
#if AMREX_SPACEDIM>2
				for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++)
#endif
				{
					amrex::IntVect m(AMREX_D_DECL(i,j,k));
					if(std::isnan(water_conc_old_box(m))) amrex::Abort("Nan found in WATER_OLD(i,j,k)");
					if(std::isinf(water_conc_old_box(m))) amrex::Abort("Inf found in WATER_OLD(i,j,k)");
					if(water_conc_old_box(m) > 2.0) 
					{
						std::cout << "dt = " << dt << ", time = " << time << ", lev = " << lev << std::endl;
						amrex::Abort("water concentration exceeded 2");
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

			for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
				for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
#if AMREX_SPACEDIM>2
				for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++)
#endif
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

		for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
			for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
#if AMREX_SPACEDIM > 2
				for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++)
#endif
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
	
	if(elastic_on)
	{  
		displacement[lev].get()->setVal(0.0);
		strain[lev].get()->setVal(0.0); 
		stress[lev].get()->setVal(0.0); 
		stress_vm[lev].get()->setVal(0.0);
		rhs[lev].get()->setVal(0.0);
		energy[lev].get()->setVal(0.0); 
		energies[lev].get()->setVal(0.0); 
	}
	eta_ic->Initialize(lev,eta_new);
	eta_ic->Initialize(lev,eta_old);

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

			for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
				for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
#if AMREX_SPACEDIM>2
					for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++)
#endif
				{
					amrex::Real grad1 = (WATER(i+1,j,k) - WATER(i-1,j,k))/(2.*dx[0]);
					amrex::Real grad2 = (WATER(i,j+1,k) - WATER(i,j-1,k))/(2.*dx[1]);
#if AMREX_SPACEDIM>2
					amrex::Real grad3 = (WATER(i,j,k+1) - WATER(i,j,k-1))/(2.*dx[2]);
#endif
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

			for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
				for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
#if AMREX_SPACEDIM>2
					for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++)
#endif
				{
					amrex::Real grad1 = (TEMP(i+1,j,k) - TEMP(i-1,j,k))/(2.*dx[0]);
					amrex::Real grad2 = (TEMP(i,j+1,k) - TEMP(i,j-1,k))/(2.*dx[1]);
#if AMREX_SPACEDIM>2
					amrex::Real grad3 = (TEMP(i,j,k+1) - TEMP(i,j,k-1))/(2.*dx[2]);
#endif
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

		for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
			for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
#if AMREX_SPACEDIM > 2
				for (int k = bx.loVect()[2]; k <= bx.hiVect()[2]; k++)
#endif
			{
				for (int m = 0; m < number_of_eta; m++)
				{
					Set::Scalar gradx = (ETA_NEW(i+1,j,k,m) - ETA_NEW(i-1,j,k,m))/(2.*dx[0]);
					Set::Scalar grady = (ETA_NEW(i,j+1,k,m) - ETA_NEW(i,j-1,k,m))/(2.*dx[1]);
#if AMREX_SPACEDIM >2
					Set::Scalar gradz = (ETA_NEW(i,j,k+1,m) - ETA_NEW(i,j,k-1,m))/(2.*dx[2]);
#endif
					Set::Scalar grad = sqrt(AMREX_D_TERM(gradx*gradx, + grady*grady, + gradz*gradz));
					Set::Scalar dr = sqrt(AMREX_D_TERM(dx[0]*dx[0], + dx[1]*dx[1], + dx[2]*dx[2]));

					if(grad*dr > damage_refinement_threshold) 
						tag(amrex::IntVect(AMREX_D_DECL(i,j,k))) = amrex::TagBox::SET;
				}
			}

	}
}

void PolymerDegradation::TimeStepBegin(amrex::Real time, int iter)
{
	if (!elastic_on) return;
	if (iter%elastic_int) return;
	if (time < elastic_tstart) return;
	if (time > elastic_tend) return;
	LPInfo info;
	info.setAgglomeration(true);
	info.setConsolidation(true);

	elastic_operator = new Operator::Elastic::Degradation::Degradation(0.0,damage_anisotropy,damage_type);
  
	geom[0].isPeriodic(0);
	elastic_operator->define(geom,grids,dmap,info);
	elastic_operator->setMaxOrder(2);
	elastic_operator->SetBC({{AMREX_D_DECL(bc_x_lo,bc_y_lo,bc_z_lo)}},
		     				{{AMREX_D_DECL(bc_x_hi,bc_y_hi,bc_z_hi)}});
	//elastic_operator->setDomainBC(((BC::Changeable *)elastic_bc)->GetBCTypes<amrex::LinOpBCType>()[0],
	//				((BC::Changeable *)elastic_bc)->GetBCTypes<amrex::LinOpBCType>()[1]);

	elastic_operator->SetEta(eta_new,*eta_bc,models);
	Numeric::Interpolator::Linear::Linear<Set::Vector> interpolate_left(elastic_bc_left,elastic_bc_left_t);
	//((BC::Changeable *)elastic_bc)->SetBC(interpolate_left(time),1);

	Numeric::Interpolator::Linear::Linear<Set::Vector> interpolate_right(elastic_bc_right,elastic_bc_right_t);
	//((BC::Changeable *)elastic_bc)->SetBC(interpolate_right(time),2);

	Numeric::Interpolator::Linear::Linear<Set::Vector> interpolate_bottom(elastic_bc_bottom,elastic_bc_bottom_t);
	//((BC::Changeable *)elastic_bc)->SetBC(interpolate_bottom(time),3);

	Numeric::Interpolator::Linear::Linear<Set::Vector> interpolate_top(elastic_bc_top,elastic_bc_top_t);
	//((BC::Changeable *)elastic_bc)->SetBC(interpolate_top(time),4);

#if AMREX_SPACEDIM > 2
	Numeric::Interpolator::Linear::Linear<Set::Vector> interpolate_back(elastic_bc_back,elastic_bc_back_t);
	//((BC::Changeable *)elastic_bc)->SetBC(interpolate_back(time),5);

	Numeric::Interpolator::Linear::Linear<Set::Vector> interpolate_front(elastic_bc_front,elastic_bc_front_t);
	//((BC::Changeable *)elastic_bc)->SetBC(interpolate_front(time),6);
#endif

	for (int ilev = 0; ilev < displacement.size(); ++ilev)
	{
		//elastic_bc->FillBoundary(*displacement[ilev],0,0,0.0);
		//elastic_operator->setLevelBC(ilev,displacement[ilev].get());

		const Real* DX = geom[ilev].CellSize();
		Set::Scalar volume = AMREX_D_TERM(DX[0],*DX[1],*DX[2]);

		AMREX_D_TERM(rhs[ilev]->setVal(body_force[0]*volume,0,1);,
					rhs[ilev]->setVal(body_force[1]*volume,1,1);,
					rhs[ilev]->setVal(body_force[2]*volume,2,1););

		if (iter==0)
		{
			displacement[ilev]->setVal(0.0);
		}
	}

	amrex::MLMG solver(*elastic_operator);
	solver.setMaxIter(elastic_max_iter);
	solver.setMaxFmgIter(elastic_max_fmg_iter);
	solver.setVerbose(elastic_verbose);
	solver.setCGVerbose(elastic_cgverbose);
	solver.setFinalFillBC(true);
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
		for ( amrex::MFIter mfi(*displacement[lev],true); mfi.isValid(); ++mfi )
		{
			const Box& bx = mfi.tilebox();

			FArrayBox &ufab  = (*displacement[lev])[mfi];
			FArrayBox &epsfab  = (*strain[lev])[mfi];
			FArrayBox &sigmafab  = (*stress[lev])[mfi];
			FArrayBox &sigmavmfab  = (*stress_vm[lev])[mfi];

			elastic_operator->Stress(sigmafab,ufab,lev,mfi);

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

			FArrayBox &energyfab  = (*energy[lev])[mfi];
			elastic_operator->Energy(energyfab,ufab,lev,mfi);

			FArrayBox &energiesfab  = (*energies[lev])[mfi];
			elastic_operator->Energies(energiesfab,ufab,lev,0,mfi);
		}
	}
}
}
#endif
