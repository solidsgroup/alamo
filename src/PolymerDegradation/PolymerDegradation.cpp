#include "PolymerDegradation/PolymerDegradation.H"

#if BL_SPACEDIM == 2

PolymerDegradation::PolymerDegradation() :
	Integrator::Integrator(),
	mybc(geom)
{

	//
	// READ INPUT PARAMETERS
	//

	amrex::ParmParse pp("water"); // Water diffusion parameters
	pp.query("on",water_diffusion_on);
	if(water_diffusion_on)
	{
		pp.query("diffusivity", water_diffusivity);
		pp.query("refinement_threshold", water_refinement_threshold);
		pp.query("ic_type", water_ic_type);

		// // Determine initial condition
		if (water_ic_type == "constant")
		{
			amrex::ParmParse pp("water.ic");
			amrex::Vector<amrex::Real> value;
			pp.query("value",value);
			water_ic = new IC::Constant(geom,value);
		}
		else
			amrex::Abort("This kind of IC has not been implemented yet");
		

		amrex::ParmParse pp("water.bc");

		amrex::Vector<std::string> bc_hi_str(AMREX_SPACEDIM);
		amrex::Vector<std::string> bc_lo_str(AMREX_SPACEDIM);

		pp.queryarr("lo",bc_lo_str,0,BL_SPACEDIM);
		pp.queryarr("hi",bc_hi_str,0,BL_SPACEDIM);
		amrex::Vector<amrex::Real> bc_lo_1, bc_hi_1;
		amrex::Vector<amrex::Real> bc_lo_2, bc_hi_2;
		amrex::Vector<amrex::Real> bc_lo_3, bc_hi_3;

		if (pp.countval("lo_1")) pp.getarr("lo_1",bc_lo_1);
		if (pp.countval("hi_1")) pp.getarr("hi_1",bc_hi_1);
		
		if (pp.countval("lo_2")) pp.getarr("lo_2",bc_lo_2);
		if (pp.countval("hi_2")) pp.getarr("hi_2",bc_hi_2);
		
		if (pp.countval("lo_3")) pp.getarr("lo_3",bc_lo_3);
		if (pp.countval("hi_3")) pp.getarr("hi_3",bc_hi_3);

		water_bc = new BC::BC(geom, bc_hi_str, bc_lo_str,
		      bc_lo_1, bc_hi_1, bc_lo_2, bc_hi_2,
		      bc_lo_3, bc_hi_3);

		RegisterNewFab(water_conc,     water_bc, 1, number_of_ghost_cells, "Water Concentration");
		RegisterNewFab(water_conc_old, water_bc, 1, number_of_ghost_cells, "Water Concentration Old");
	}

	amrex::ParmParse pp("thermal"); //Heat diffusion parameters
	pp.query("on",heat_diffusion_on);
	if(heat_diffusion_on)
	{
		pp.query("diffusivity", thermal_diffusivity);
		pp.query("refinement_threshold",thermal_refinement_threshold);
		pp.query("ic_type",thermal_ic_type);
	
		if (thermal_ic_type == "constant")
		{
			amrex::ParmParse pp("thermal.ic");
			amrex::Vector<amrex::Real> T;
			pp.query("value",T);
			thermal_ic = new IC::Constant(geom,T);
		}
		else
			amrex::Abort("This kind of IC has not been implemented yet");

		amrex::ParmParse pp("thermal.bc");

		amrex::Vector<std::string> bc_hi_str(AMREX_SPACEDIM);
		amrex::Vector<std::string> bc_lo_str(AMREX_SPACEDIM);

		pp.queryarr("lo",bc_lo_str,0,BL_SPACEDIM);
		pp.queryarr("hi",bc_hi_str,0,BL_SPACEDIM);
		amrex::Vector<amrex::Real> bc_lo_1, bc_hi_1;
		amrex::Vector<amrex::Real> bc_lo_2, bc_hi_2;
		amrex::Vector<amrex::Real> bc_lo_3, bc_hi_3;

		if (pp.countval("lo_1")) pp.getarr("lo_1",bc_lo_1);
		if (pp.countval("hi_1")) pp.getarr("hi_1",bc_hi_1);
		
		if (pp.countval("lo_2")) pp.getarr("lo_2",bc_lo_2);
		if (pp.countval("hi_2")) pp.getarr("hi_2",bc_hi_2);
		
		if (pp.countval("lo_3")) pp.getarr("lo_3",bc_lo_3);
		if (pp.countval("hi_3")) pp.getarr("hi_3",bc_hi_3);

		thermal_bc = new BC::BC(geom, bc_hi_str, bc_lo_str,
		      bc_lo_1, bc_hi_1, bc_lo_2, bc_hi_2,
		      bc_lo_3, bc_hi_3);

		RegisterNewFab(Temp,     thermal_bc, 1, number_of_ghost_cells, "Temperature");
		RegisterNewFab(Temp_old, thermal_bc, 1, number_of_ghost_cells, "Temperature Old");
	}

	amrex::ParmParse pp("damage"); // Phase-field model parameters
	pp.query("anisotropy",damage_anisotropy);

	if(damage_anisotropy == 0)
		number_of_eta = 1;
	else
		number_of_eta = BL_SPACEDIM;

	pp.query("type",damage_type);
	if(damage_type == "relaxation") //meant for isotropic right now.
	{
		pp.qeury("E0",E0); // to be replaced to get the value from elastic operator
		pp.query("number_of_terms",number_of_terms);
		pp.queryarr("E_i",E_i);
		pp.queryarr("tau_i",tau_i);
		pp.queryarr("t_start_i",t_start_i);

		if(E_i.size() != number_of_terms || tau_i.size() != number_of_terms || t_start_i.size() != number_of_terms)
			amrex::Abort("missing entries in E_i, tau_i or t_start_i");

		amrex::Real sum = 0;		
		for (int temp = 0; temp < E_i.size(); temp++)
			sum += E_i[temp];

		if(sum != E0) //need to replace this in the future
			amrex::Abort("E0 is not equal to the sum of E_i");
	}
	else
		amrex::Abort("This kind of damage model has not been implemented yet");

	pp.query("ic_type",eta_ic_type)
	if(eta_ic_type == "constant")
	{
		amrex::ParmParse pp("damage.ic");
		amrex::Vector<amrex::Real> eta_init;
		pp.qeury("value",eta_init);
		eta_ic = new IC::Constant(geom,eta_init);
	}
	else
		amrex::Abort("This kind of IC has not been implemented yet");
	
	amrex::ParmParse pp("damage.bc");

	amrex::Vector<std::string> bc_hi_str(AMREX_SPACEDIM);
	amrex::Vector<std::string> bc_lo_str(AMREX_SPACEDIM);

	pp.queryarr("lo",bc_lo_str,0,BL_SPACEDIM);
	pp.queryarr("hi",bc_hi_str,0,BL_SPACEDIM);
	amrex::Vector<amrex::Real> bc_lo_1, bc_hi_1;
	amrex::Vector<amrex::Real> bc_lo_2, bc_hi_2;
	amrex::Vector<amrex::Real> bc_lo_3, bc_hi_3;

	if (pp.countval("lo_1")) pp.getarr("lo_1",bc_lo_1);
	if (pp.countval("hi_1")) pp.getarr("hi_1",bc_hi_1);
	
	if (pp.countval("lo_2")) pp.getarr("lo_2",bc_lo_2);
	if (pp.countval("hi_2")) pp.getarr("hi_2",bc_hi_2);
		
	if (pp.countval("lo_3")) pp.getarr("lo_3",bc_lo_3);
	if (pp.countval("hi_3")) pp.getarr("hi_3",bc_hi_3);

	damage_bc = new BC::BC(geom, bc_hi_str, bc_lo_str,
	      bc_lo_1, bc_hi_1, bc_lo_2, bc_hi_2,
	      bc_lo_3, bc_hi_3);

	RegisterNewFab(eta_new, damage_bc, number_of_eta, number_of_ghost_cells, "Eta");
	RegisterNewFab(eta_old, damage_bc, number_of_eta, number_of_ghost_cells, "Eta old");

  
	// Elasticity
	amrex::ParmParse pp("elastic");
	pp.query("on",elastic_on);
	pp.query("int",elastic_int);
	pp.query("type",elastic_type);
	pp.query("max_iter",elastic_max_iter);
	pp.query("max_fmg_iter",elastic_max_fmg_iter);
	pp.query("verbose",elastic_verbose);
	pp.query("cgverbose",elastic_cgverbose);
	pp.query("tol_rel",elastic_tol_rel);
	pp.query("tol_abs",elastic_tol_abs);
	pp.query("tstart",elastic_tstart);

	pp.queryarr("load_t",elastic_load_t);
	pp.queryarr("load_disp",elastic_load_disp);
	if (elastic_load_t.size() != elastic_load_disp.size())
		amrex::Abort("load_t and load_disp must have the same number of entries");

	if (elastic_on)
	{
		RegisterNewFab(displacement, mybc, AMREX_SPACEDIM, 1, "u");
		RegisterNewFab(body_force, mybc, AMREX_SPACEDIM, 1, "b");
		RegisterNewFab(strain, mybc, 3, 1, "eps");
		RegisterNewFab(stress, mybc, 3, 1, "sig");
		RegisterNewFab(stress_vm, mybc, 1, 1, "sig_VM");
		RegisterNewFab(energy, mybc, 1, 1, "W");
		RegisterNewFab(energies, mybc, number_of_grains, 1, "W");
	}
}


#define ETA(i,j,k,n) eta_old_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),n)

void
PolymerDegradation::Advance (int lev, amrex::Real time, amrex::Real dt)
{
	std::swap(eta_old[lev], eta_new[lev]);

	if(water_diffusion_on) std::swap(water_conc_old[lev],water_conc[lev]);
	if(heat_diffusion_on) std::swap(Temp_old[lev], Temp[lev]);

	const amrex::Real* dx = geom[lev].CellSize();

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
					WATER(i,j,k) =
						WATER_OLD(i,j,k)
						+ dt * water_diffusivity * (AMREX_D_TERM((WATER_OLD(i+1,j,k) + WATER_OLD(i-1,j,k) - 2.0*WATER_OLD(i,j,k)) / dx[0] / dx[0],
						+ (WATER_OLD(i,j+1,k) + WATER_OLD(i,j-1,k) - 2.0*WATER_OLD(i,j,k)) / dx[1] / dx[1],
						+ (WATER_OLD(i,j,k+1) + WATER_OLD(i,j,k-1) - 2.0*WATER_OLD(i,j,k)) / dx[2] / dx[2]));
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
					TEMP(i,j,k) =
						TEMP_OLD(i,j,k)
						+ dt * thermal_diffusivity * (AMREX_D_TERM((TEMP_OLD(i+1,j,k) + TEMP_OLD(i-1,j,k) - 2.0*TEMP_OLD(i,j,k)) / dx[0] / dx[0],
						+ (TEMP_OLD(i,j+1,k) + TEMP_OLD(i,j-1,k) - 2.0*TEMP_OLD(i,j,k)) / dx[1] / dx[1],
						+ (TEMP_OLD(i,j,k+1) + TEMP_OLD(i,j,k-1) - 2.0*TEMP_OLD(i,j,k)) / dx[2] / dx[2]));
				}
		}
	}

	for ( amrex::MFIter mfi(*eta_new[lev],true); mfi.isValid(); ++mfi )
	{
		const amrex::Box& bx = mfi.tilebox();
		FArrayBox &eta_new_box     = (*eta_new[lev])[mfi];
		FArrayBox &eta_old_box     = (*eta_old[lev])[mfi];

		//FArrayBox &energiesfab = (*energy)

		for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
			for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
#if AMREX_SPACEDIM > 2
				for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++)
#endif
			{
				for (int m = 0; m < number_of_eta; m++)
				{
					if(damage_type == "relaxation") // Need to replace this with more sophisticated check
					{
						amrex::Real rhs = 0.0;
						for (int l = 0; l<number_of_terms; l++)
							rhs += E_i[l]*std::exp(-std::max(0.0,time-t_start_i[l])/tau_i[l])/(E0*tau_i[l]);

						eta_new_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),m) =
							ETA(i,j,k,m) +
							dt*rhs;
					}

					//
					// ELASTIC DRIVING FORCE
					//
					if (elastic_on)
					{
						FArrayBox &energiesfab     = (*energies[lev])[mfi];

						eta_new_box(amrex::IntVect(AMREX_D_DECL(i,j,k)),m)
							+= M*dt*( energiesfab(amrex::IntVect(AMREX_D_DECL(i,j,k)),m));
					}

				}
			}
	}
}

void
PolymerDegradation::Initialize (int lev)
{
	ic->Initialize(lev,eta_new);
	ic->Initialize(lev,eta_old);
  
	displacement[lev].get()->setVal(0.0);
	strain[lev].get()->setVal(0.0); 
	stress[lev].get()->setVal(0.0); 
	stress_vm[lev].get()->setVal(0.0);
	body_force[lev].get()->setVal(0.0);
	energy[lev].get()->setVal(0.0); 
	energies[lev].get()->setVal(0.0); 
}


void
PolymerDegradation::TagCellsForRefinement (int lev, amrex::TagBoxArray& tags, amrex::Real /*time*/, int /*ngrow*/)
{
	const amrex::Real* dx      = geom[lev].CellSize();

	amrex::Vector<int>  itags;

	for (amrex::MFIter mfi(*eta_new[lev],true); mfi.isValid(); ++mfi)
	{
		const amrex::Box&  bx  = mfi.tilebox();
		amrex::TagBox&     tag  = tags[mfi];
		amrex::BaseFab<amrex::Real> &eta_new_box = (*eta_new[lev])[mfi];

#if BL_SPACEDIM==2
		for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
			for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
			{
				for (int n = 0; n < number_of_grains; n++)
				{
					amrex::Real gradx = (eta_new_box(amrex::IntVect(i+1,j),n) - eta_new_box(amrex::IntVect(i-1,j),n))/(2.*dx[0]);
					amrex::Real grady = (eta_new_box(amrex::IntVect(i,j+1),n) - eta_new_box(amrex::IntVect(i,j-1),n))/(2.*dx[1]);
					if (dx[0]*sqrt(gradx*gradx + grady*grady)>0.1) tag(amrex::IntVect(i,j)) = amrex::TagBox::SET;
				}
			}

#elif BL_SPACEDIM==3

		for (int i = bx.loVect()[0]; i<=bx.hiVect()[0]; i++)
			for (int j = bx.loVect()[1]; j<=bx.hiVect()[1]; j++)
				for (int k = bx.loVect()[2]; k<=bx.hiVect()[2]; k++)
				{
					for (int n = 0; n < number_of_grains; n++)
					{
						amrex::Real gradx = (eta_new_box(amrex::IntVect(i+1,j,k),n) 
							- eta_new_box(amrex::IntVect(i-1,j,k),n))/(2.*dx[0]);
						amrex::Real grady = (eta_new_box(amrex::IntVect(i,j+1,k),n) 
							- eta_new_box(amrex::IntVect(i,j-1,k),n))/(2.*dx[1]);
						amrex::Real gradz = (eta_new_box(amrex::IntVect(i,j,k+1),n) 
							- eta_new_box(amrex::IntVect(i,j,k-1),n))/(2.*dx[2]);
						if (dx[0]*sqrt(gradx*gradx + grady*grady + gradz*gradz)>0.1) tag(amrex::IntVect(i,j,k)) = amrex::TagBox::SET;
					}
				}
#endif

	}
}

void PolymerDegradation::TimeStepBegin(amrex::Real time, int iter)
{
	if (!elastic_on) return;
	if (iter%elastic_int) return;
	if (time < elastic_tstart) return;

	LPInfo info;
	info.setAgglomeration(true);
	info.setConsolidation(true);

	elastic_operator = new Operator::Elastic::PolyCrystal::PolyCrystal();
  
	geom[0].isPeriodic(0);
	elastic_operator->define(geom,grids,dmap,info);
	elastic_operator->setMaxOrder(2);
	elastic_operator->setDomainBC(
		{AMREX_D_DECL(geom[0].isPeriodic(0) ? LinOpBCType::Periodic : LinOpBCType::Dirichlet ,
			geom[0].isPeriodic(1) ? LinOpBCType::Periodic : LinOpBCType::Dirichlet,
			geom[0].isPeriodic(2) ? LinOpBCType::Periodic : LinOpBCType::Dirichlet)},
		{AMREX_D_DECL(geom[0].isPeriodic(0) ? LinOpBCType::Periodic : LinOpBCType::Dirichlet,
			geom[0].isPeriodic(1) ? LinOpBCType::Periodic : LinOpBCType::Dirichlet,
			geom[0].isPeriodic(2) ? LinOpBCType::Periodic : LinOpBCType::Dirichlet)});

	std::vector<Operator::Elastic::PolyCrystal::PolyCrystalModel *> models;
	
	// for (int n = 0; n <  number_of_grains; n++) {} <<<<<< need to replace with this!
	models.push_back(new Operator::Elastic::PolyCrystal::Cubic(107.3, 60.9, 28.30,
							     2.49, 2.49, 4.328));
	models.push_back(new Operator::Elastic::PolyCrystal::Cubic(107.3, 60.9, 28.30,
							     0.99, 0.511, 1.39));
	//models.push_back(&g2);

	elastic_operator->SetEta(eta_new,mybc,models);

	amrex::Real ushear = 0.0;

	for (int i = 0; i<elastic_load_t.size(); i++)
	{
		if (time < elastic_load_t[0])
		{
			ushear = elastic_load_disp[0];
			break;
		}
		if (i < elastic_load_t.size()-1 &&  elastic_load_t[i] < time && time < elastic_load_t[i+1])
		{
			ushear = elastic_load_disp[i] +
				(time - elastic_load_t[i]) * (elastic_load_disp[i+1] - elastic_load_disp[i]) /
				(elastic_load_t[i+1] - elastic_load_t[i]);
			break;
		}
		else
		{
			ushear = elastic_load_disp[elastic_load_t.size()-1];
			break;
		}
	}
	std::cout << "ushear = " << ushear << std::endl;


	for (int ilev = 0; ilev < displacement.size(); ++ilev)
	{
		amrex::Box domain(geom[ilev].Domain());
		for (amrex::MFIter mfi(*displacement[ilev], true); mfi.isValid(); ++mfi)
		{
			const amrex::Box& box = mfi.tilebox();
			amrex::BaseFab<amrex::Real> &disp_box = (*displacement[ilev])[mfi];

			for (int i = box.loVect()[0] - displacement[ilev]->nGrow(); i<=box.hiVect()[0] + displacement[ilev]->nGrow(); i++)
				for (int j = box.loVect()[1] - displacement[ilev]->nGrow(); j<=box.hiVect()[1] + displacement[ilev]->nGrow(); j++)
				{
					if (j > domain.hiVect()[1]) // Top boundary
					{
						disp_box(amrex::IntVect(i,j),0) = ushear;
						disp_box(amrex::IntVect(i,j),1) = 0.0;
					}
					else if (i > domain.hiVect()[0]) // Right boundary
					{
						disp_box(amrex::IntVect(i,j),0) = 0.0;
						disp_box(amrex::IntVect(i,j),1) = 0.0;
					}
					else if (j < domain.loVect()[1]) // Bottom
					{
						disp_box(amrex::IntVect(i,j),0) = 0.0;
						disp_box(amrex::IntVect(i,j),1) = 0.0;
					}
					else if (i < domain.loVect()[0]) // Left boundary
					{
						disp_box(amrex::IntVect(i,j),0) = 0.0;
						disp_box(amrex::IntVect(i,j),1) = 0.0;
					}
				}
		}
		elastic_operator->setLevelBC(ilev,displacement[ilev].get());

		/// \todo Replace with proper driving force initialization
		body_force[ilev]->setVal(0.0,0);
		body_force[ilev]->setVal(0.0,1);

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

	solver.solve(GetVecOfPtrs(displacement),
		GetVecOfConstPtrs(body_force),
		elastic_tol_rel,
		elastic_tol_abs);

	for (int lev = 0; lev < displacement.size(); lev++)
	{
		const amrex::Real* dx = geom[lev].CellSize();
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
				{
					epsfab(amrex::IntVect(i,j),0) = (ufab(amrex::IntVect(i+1,j),0) - ufab(amrex::IntVect(i-1,j),0))/(2.0*dx[0]);
					epsfab(amrex::IntVect(i,j),1) = (ufab(amrex::IntVect(i,j+1),1) - ufab(amrex::IntVect(i,j-1),1))/(2.0*dx[1]);
					epsfab(amrex::IntVect(i,j),2) = 0.5*(ufab(amrex::IntVect(i+1,j),1) - ufab(amrex::IntVect(i-1,j),1))/(2.0*dx[0]) +
						0.5*(ufab(amrex::IntVect(i,j+1),0) - ufab(amrex::IntVect(i,j-1),0))/(2.0*dx[1]);

					sigmavmfab(amrex::IntVect(i,j)) =
						sqrt(0.5*((sigmafab(amrex::IntVect(i,j),0) - sigmafab(amrex::IntVect(i,j),1)
						*(sigmafab(amrex::IntVect(i,j),0) - sigmafab(amrex::IntVect(i,j),1))
						+ sigmafab(amrex::IntVect(i,j),0)*sigmafab(amrex::IntVect(i,j),0)
						+ sigmafab(amrex::IntVect(i,j),1)*sigmafab(amrex::IntVect(i,j),1)
						+ 6.0*sigmafab(amrex::IntVect(i,j),2)*sigmafab(amrex::IntVect(i,j),2))));
				}

			FArrayBox &energyfab  = (*energy[lev])[mfi];
			elastic_operator->Energy(energyfab,ufab,lev,mfi);

			FArrayBox &energiesfab  = (*energies[lev])[mfi];
			elastic_operator->Energies(energiesfab,ufab,lev,mfi);
		}
	}
}


#endif
