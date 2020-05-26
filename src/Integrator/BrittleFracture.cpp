#include "BrittleFracture.H"

namespace Integrator
{
BrittleFracture::BrittleFracture() :
	Integrator()
{
	// Crack model
	IO::ParmParse pp_crack("crack");
	std::string crack_type;
	pp_crack.query("type",crack_type);
	pp_crack.query("modulus_scaling_max",scaleModulusMax);
	pp_crack.query("refinement_threshold",refinement_threshold);

	if(crack_type=="constant")
	{
		Model::Interface::Crack::Constant *tmpbdy = new Model::Interface::Crack::Constant();
		pp_crack.queryclass("constant",*tmpbdy);
		crack.boundary = tmpbdy;
	}
	else if(crack_type == "sin")
	{
		Model::Interface::Crack::Sin *tmpbdy = new Model::Interface::Crack::Sin();
		pp_crack.queryclass("sin", *tmpbdy);
		crack.boundary = tmpbdy;
	}
	else
		Util::Abort(INFO,"This crack model hasn't been implemented yet");

	// Anistropy
	IO::ParmParse pp_anisotropy("crack.anisotropy");
	pp_anisotropy.query("on", anisotropy.on);
	pp_anisotropy.query("tstart",anisotropy.tstart);
	anisotropy.timestep = timestep;
	pp_anisotropy.query("timestep", anisotropy.timestep);
	anisotropy.plot_int = plot_int;
	pp_anisotropy.query("plot_int", anisotropy.plot_int);
	anisotropy.plot_dt = plot_dt;
	pp_anisotropy.query("plot_dt", anisotropy.plot_dt);
	pp_anisotropy.query("beta", anisotropy.beta);

	// ICs
	IO::ParmParse pp("ic"); // Phase-field model parameters
	pp.query("type", crack.ic_type);

	if(crack.ic_type == "ellipsoid")
	{
		IC::Ellipsoid *tmpic = new IC::Ellipsoid(geom);
		pp.queryclass("ellipsoid", *tmpic);
		crack.ic = tmpic;
	}
	else if(crack.ic_type == "notch")
	{
		IC::Notch *tmpic = new IC::Notch(geom);
		pp.queryclass("notch",*tmpic);
		crack.ic = tmpic;
	}
		
	else
		Util::Abort(INFO,"This type of IC hasn't been implemented yet");
	
	pp_crack.query("tol_crack",tol_crack);
	pp_crack.query("tol_step",tol_step);

	// BCs
	// Crack field should have a zero neumann BC. So we just code it up here.
	// In case this needs to change, we can add options to read it from input.
	// Below are conditions for full simulation.  If this doesn't work, we can
	// try symmetric simulation.
	crack.mybc = new BC::Constant(1);
	crack.mybcdf = new BC::Constant(4);
	pp_crack.queryclass("bc",*static_cast<BC::Constant *>(crack.mybc));
	pp_crack.queryclass("bc_df",*static_cast<BC::Constant *>(crack.mybcdf));

	RegisterNewFab(m_c,     crack.mybc, 1, number_of_ghost_cells+1, "c",		true);
	RegisterNewFab(m_c_old, crack.mybc, 1, number_of_ghost_cells+1, "c_old",	true);
	RegisterNewFab(m_driving_force, crack.mybcdf, 4, number_of_ghost_cells+1, "driving_force",true);
	
	crack_err_norm = 0.; crack_err_temp_norm = 0.;
	crack_err_norm_init = 1.e4; crack_err_temp_norm_init = 1.e4;

	disp_err_norm = 0.; disp_err_norm_init = 1.e4;

	RegisterIntegratedVariable(&crack_err_norm, "crack_err_norm");
	RegisterIntegratedVariable(&c_new_norm,"c_new_norm");
	
	// Material input
	IO::ParmParse pp_material("material");
	pp_material.query("model",material.input_material);
	if(material.input_material == "isotropic")
	{
		pp_material.queryclass("isotropic",material.modeltype);
	}
	else
		Util::Abort(INFO,"This model has not been implemented yet.");

	// Elasticity properties
	IO::ParmParse pp_elastic("elastic");
	pp_elastic.query("int",				elastic.interval);
	pp_elastic.query("type",			elastic.type);
	pp_elastic.query("max_iter",		elastic.max_iter);
	pp_elastic.query("max_fmg_iter",	elastic.max_fmg_iter);
	pp_elastic.query("verbose",			elastic.verbose);
	pp_elastic.query("cgverbose",		elastic.cgverbose);
	pp_elastic.query("tol_rel",			elastic.tol_rel);
	pp_elastic.query("tol_abs",			elastic.tol_abs);
	pp_elastic.query("cg_tol_rel",		elastic.cg_tol_rel);
	pp_elastic.query("cg_tol_abs",		elastic.cg_tol_abs);
	pp_elastic.query("use_fsmooth",		elastic.use_fsmooth);
	pp_elastic.query("agglomeration", 	elastic.agglomeration);
	pp_elastic.query("consolidation", 	elastic.consolidation);

	pp_elastic.query("bottom_solver",elastic.bottom_solver);
	pp_elastic.query("linop_maxorder", elastic.linop_maxorder);
	pp_elastic.query("max_coarsening_level",elastic.max_coarsening_level);
	pp_elastic.query("verbose",elastic.verbose);
	pp_elastic.query("cg_verbose", elastic.cgverbose);
	pp_elastic.query("bottom_max_iter", elastic.bottom_max_iter);
	pp_elastic.query("max_fixed_iter", elastic.max_fixed_iter);
	pp_elastic.query("bottom_tol", elastic.bottom_tol);
	pp_elastic.query("df_mult", elastic.df_mult);

	if (pp_elastic.countval("body_force")) pp_elastic.queryarr("body_force",elastic.body_force);

	pp_elastic.queryclass("bc",elastic.bc);
	IO::ParmParse pp_elastic_bc("elastic.bc");
	pp_elastic_bc.query("disp_step",elastic.test_rate);
	pp_elastic_bc.query("disp_init",elastic.test_init);
	pp_elastic_bc.query("max_disp",elastic.test_max);
	pp_elastic_bc.query("crackStressTest",crack.crackStressTest);

	elastic.bc_top = elastic.test_init;

	if(elastic.test_rate < 0.) { Util::Warning(INFO,"Rate can't be less than zero. Resetting to 1.0"); elastic.test_rate = 0.1; }
	if(elastic.test_max < 0. ||  elastic.test_max < elastic.test_rate) {Util::Warning(INFO,"Max can't be less than load step. Resetting to load step"); elastic.test_max = elastic.test_rate;}

	const int number_of_stress_components = AMREX_SPACEDIM*AMREX_SPACEDIM;
	
	RegisterNodalFab (m_disp, 		AMREX_SPACEDIM, 				number_of_ghost_nodes, "Disp",true);
	RegisterNodalFab (m_rhs,  		AMREX_SPACEDIM, 				number_of_ghost_nodes, "RHS",true);
	RegisterNodalFab (m_strain,		number_of_stress_components,	number_of_ghost_nodes,	"strain",true);
	RegisterNodalFab (m_stress,		number_of_stress_components,	number_of_ghost_nodes,	"stress",true);
	RegisterNodalFab (m_stressvm,	1,								number_of_ghost_nodes,	"stress_vm",true);
	RegisterNodalFab (m_energy,		1,								number_of_ghost_nodes,	"energy",true);
	RegisterNodalFab (m_energy_pristine,		1,					number_of_ghost_nodes,	"energyP",true);
	RegisterNodalFab (m_energy_pristine_old,	1,					number_of_ghost_nodes,	"energyPOld",true);
	RegisterNodalFab (m_residual,	AMREX_SPACEDIM,					number_of_ghost_nodes,	"residual",true);

	RegisterGeneralFab(material.model, 1, 2);

	nlevels = maxLevel() + 1;
}

BrittleFracture::~BrittleFracture()
{
}

void
BrittleFracture::Initialize (int lev)
{
	crack.ic->Initialize(lev,m_c);
	crack.ic->Initialize(lev,m_c_old);
	
	m_driving_force[lev]->setVal(0.0);
	
	m_disp[lev]->setVal(0.0);
	m_strain[lev]->setVal(0.0);
	m_stress[lev]->setVal(0.0);
	m_stressvm[lev]->setVal(0.0);
	m_rhs[lev]->setVal(0.0);
	m_energy[lev]->setVal(0.0);
	m_residual[lev]->setVal(0.0);
	m_energy_pristine[lev] -> setVal(0.);
	m_energy_pristine_old[lev] -> setVal(0.);
	
	material.model[lev]->setVal(material.modeltype);
}

void 
BrittleFracture::TimeStepBegin(amrex::Real time, int iter)
{
	if (anisotropy.on && time >= anisotropy.tstart)
	{
		SetTimestep(anisotropy.timestep);
		if (anisotropy.elastic_int > 0) 
			if (iter % anisotropy.elastic_int) return;
	}
	if(iter%elastic.interval) return;

	elastic.bc_top = elastic.test_init + ((double)elastic.test_step)*elastic.test_rate;
	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		std::swap(*m_energy_pristine_old[ilev], *m_energy_pristine[ilev]);
		m_c[ilev]->FillBoundary();

		static amrex::IntVect AMREX_D_DECL(dx(AMREX_D_DECL(1,0,0)),
					   dy(AMREX_D_DECL(0,1,0)),
					   dz(AMREX_D_DECL(0,0,1)));

		for (amrex::MFIter mfi(*(material.model)[ilev],true); mfi.isValid(); ++mfi)
		{
			amrex::Box box = mfi.growntilebox(1);
			amrex::Array4<const amrex::Real> const& c_new = (*m_c[ilev]).array(mfi);
			amrex::Array4<fracture_model_type> const& modelfab = (material.model)[ilev]->array(mfi);

			amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k){
				Set::Scalar _temp = 0;
				Set::Scalar mul = AMREX_D_PICK(0.5,0.25,0.125);
				_temp =  mul*(AMREX_D_TERM(	
									crack.boundary->g_phi(c_new(i,j,k,0),0.) + crack.boundary->g_phi(c_new(i-1,j,k,0),0.)
									, 
									+ crack.boundary->g_phi(c_new(i,j-1,k,0),0.) + crack.boundary->g_phi(c_new(i-1,j-1,k,0),0.)
									, 
									+ crack.boundary->g_phi(c_new(i,j,k-1,0),0.) + crack.boundary->g_phi(c_new(i-1,j,k-1,0),0.)
									+ crack.boundary->g_phi(c_new(i,j-1,k-1,0),0.) + crack.boundary->g_phi(c_new(i-1,j-1,k-1,0),0.))
									);
				if (std::isnan(_temp)) Util::Abort(INFO);
				if(_temp < 0.0) _temp = 0.;
				if(_temp > 1.0) _temp = 1.0;
				modelfab(i,j,k,0).DegradeModulus(std::min(1.-_temp,1.-scaleModulusMax));
			});
		}
	}

	for (int ilev=0; ilev < nlevels; ++ilev) Util::RealFillBoundary(*material.model[ilev],geom[ilev]);

	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		m_disp[ilev]->setVal(0.0);
		const Real* DX = geom[ilev].CellSize();
		Set::Scalar volume = AMREX_D_TERM(DX[0],*DX[1],*DX[2]);

		AMREX_D_TERM(m_rhs[ilev]->setVal(elastic.body_force(0)*volume,0,1);,
			     m_rhs[ilev]->setVal(elastic.body_force(1)*volume,1,1);,
			     m_rhs[ilev]->setVal(elastic.body_force(2)*volume,2,1););
	}

	// Mode I - displacement
	elastic.bc.Set(elastic.bc.Face::YHI, elastic.bc.Direction::Y, BC::Operator::Elastic<fracture_model_type>::Type::Displacement, elastic.bc_top);
	elastic.bc.Set(elastic.bc.Face::XLO_YHI, elastic.bc.Direction::Y, BC::Operator::Elastic<fracture_model_type>::Type::Displacement, elastic.bc_top);
	elastic.bc.Set(elastic.bc.Face::XHI_YHI, elastic.bc.Direction::Y, BC::Operator::Elastic<fracture_model_type>::Type::Displacement, elastic.bc_top);

	// Mode I - force
	//elastic.bc.Set(elastic.bc.Face::YHI, elastic.bc.Direction::Y, BC::Operator::Elastic<fracture_model_type>::Type::Traction, elastic.bc_top);
	//elastic.bc.Set(elastic.bc.Face::XLO_YHI, elastic.bc.Direction::Y, BC::Operator::Elastic<fracture_model_type>::Type::Traction, elastic.bc_top);
	//elastic.bc.Set(elastic.bc.Face::XHI_YHI, elastic.bc.Direction::Y, BC::Operator::Elastic<fracture_model_type>::Type::Traction, elastic.bc_top);

	//Mode II
	//elastic.bc.Set(elastic.bc.Face::YHI, elastic.bc.Direction::X, BC::Operator::Elastic<fracture_model_type>::Type::Displacement, elastic.bc_top);
	//elastic.bc.Set(elastic.bc.Face::XLO_YHI, elastic.bc.Direction::X, BC::Operator::Elastic<fracture_model_type>::Type::Displacement, elastic.bc_top);
	//elastic.bc.Set(elastic.bc.Face::XHI_YHI, elastic.bc.Direction::X, BC::Operator::Elastic<fracture_model_type>::Type::Displacement, elastic.bc_top);
	
	elastic.bc.Init(m_rhs,geom);
	
	LPInfo info;
	info.setAgglomeration(elastic.agglomeration);
	info.setConsolidation(elastic.consolidation);
	info.setMaxCoarseningLevel(elastic.max_coarsening_level);

	Operator::Elastic<fracture_model_type> elastic_op;
	elastic_op.define(geom, grids, dmap, info);
	elastic_op.setMaxOrder(elastic.linop_maxorder);

	elastic_op.SetBC(&(elastic.bc));

	Solver::Nonlocal::Newton<fracture_model_type>  solver(elastic_op);
	solver.setMaxIter(elastic.max_iter);
	solver.setMaxFmgIter(elastic.max_fmg_iter);
	solver.setFixedIter(elastic.max_fixed_iter);
	solver.setVerbose(elastic.verbose);
	solver.setCGVerbose(elastic.cgverbose);
	solver.setBottomMaxIter(elastic.bottom_max_iter);
	solver.setBottomTolerance(elastic.cg_tol_rel) ;
	solver.setBottomToleranceAbs(elastic.cg_tol_abs) ;

	for (int ilev = 0; ilev < nlevels; ilev++) if (m_disp[ilev]->contains_nan()) Util::Warning(INFO);

	if (elastic.bottom_solver == "cg") solver.setBottomSolver(MLMG::BottomSolver::cg);
	else if (elastic.bottom_solver == "bicgstab") solver.setBottomSolver(MLMG::BottomSolver::bicgstab);

	solver.solve(m_disp, m_rhs, material.model, elastic.tol_rel, elastic.tol_abs);
	solver.compResidual(m_residual,m_disp,m_rhs,material.model);
	
	for (int lev = 0; lev < nlevels; lev++)
	{
		elastic_op.Strain(lev,*m_strain[lev],*m_disp[lev]);
		elastic_op.Stress(lev,*m_stress[lev],*m_disp[lev]);
		elastic_op.Energy(lev,*m_energy[lev],*m_disp[lev]);
	}
	for (int lev = 0; lev < nlevels; lev++)
	{
		m_strain[lev]->FillBoundary();
		m_energy_pristine_old[lev]->FillBoundary();
		
		for (amrex::MFIter mfi(*m_strain[lev],true); mfi.isValid(); ++mfi)
		{
			const amrex::Box& box = mfi.growntilebox(number_of_ghost_nodes);//validbox();
			amrex::Array4<const Set::Scalar> const& strain_box 	= (*m_strain[lev]).array(mfi);
			amrex::Array4<Set::Scalar> const& energy_box 		= (*m_energy_pristine[lev]).array(mfi);
			amrex::Array4<Set::Scalar> const& energy_box_old 	= (*m_energy_pristine_old[lev]).array(mfi);
			amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k){
				Set::Matrix eps = Numeric::FieldToMatrix(strain_box,i,j,k);
				energy_box(i,j,k,0) = material.modeltype.W(eps);
				energy_box(i,j,k,0) = energy_box(i,j,k,0) > energy_box_old(i,j,k,0) ? energy_box(i,j,k,0) : energy_box_old(i,j,k,0);
			});
		}
		m_energy_pristine[lev]->FillBoundary();
	}
}


void
BrittleFracture::Advance (int lev, amrex::Real time, amrex::Real dt)
{
	if(crack.crackStressTest) return;
	std::swap(*m_c_old[lev], *m_c[lev]);

	static amrex::IntVect AMREX_D_DECL(	dx(AMREX_D_DECL(1,0,0)), dy(AMREX_D_DECL(0,1,0)), dz(AMREX_D_DECL(0,0,1)));
	const amrex::Real* DX = geom[lev].CellSize();
	
	for ( amrex::MFIter mfi(*m_c[lev],true); mfi.isValid(); ++mfi )
	{
		const amrex::Box& bx = mfi.validbox();
		amrex::Array4<amrex::Real> const& c_old = (*m_c_old[lev]).array(mfi);
		amrex::Array4<amrex::Real> const& df = (*m_driving_force[lev]).array(mfi);
		amrex::Array4<amrex::Real> const& c_new = (*m_c[lev]).array(mfi);
		amrex::Array4<const Set::Scalar> const& energy_box = (*m_energy_pristine[lev]).array(mfi);

		amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k){

			Set::Scalar rhs = 0.0;

			// Elastic component of the driving force
			Set::Scalar en_cell = Numeric::Interpolate::NodeToCellAverage(energy_box,i,j,k,0);
			df(i,j,k,0) = crack.boundary->Dg_phi(c_old(i,j,k,0),0.)*en_cell*elastic.df_mult;
			rhs += crack.boundary->Dg_phi(c_old(i,j,k,0),0.)*en_cell*elastic.df_mult;

			// Boundary terms
			Set::Vector Dc = Numeric::Gradient(c_old, i, j, k, 0, DX);
			Set::Scalar Theta = atan2(Dc(1),Dc(0));

			Set::Scalar normgrad = Dc.lpNorm<2>();
			if (normgrad < 1E-5) // testing to speed things up a bit.
			{
				df(i,j,k,1) = 0.0;
				df(i,j,k,2) = 0.0;
				rhs = 0.0;
			}
			else
			{
			Set::Matrix DDc = Numeric::Hessian(c_old, i, j, k, 0, DX);
			Set::Scalar laplacian = DDc.trace();

			if (!anisotropy.on || time < anisotropy.tstart)
			{
				df(i,j,k,1) = crack.boundary->Gc(Theta)*crack.boundary->Dw_phi(c_old(i,j,k,0),0.)/(4.0*crack.boundary->Zeta(Theta));
				df(i,j,k,2) = 2.0*crack.boundary->Zeta(Theta)*laplacian;

				rhs += crack.boundary->Gc(Theta)*crack.boundary->Dw_phi(c_old(i,j,k,0),0.)/(4.0*crack.boundary->Zeta(Theta));
				rhs -= 2.0*crack.boundary->Zeta(Theta)*laplacian;
			}
			else
			{
#if AMREX_SPACEDIM == 1
				Util::Abort(INFO, "Anisotropy is not enabled in 1D.");
#elif AMREX_SPACEDIM == 2
				Set::Vector normal = Dc / normgrad;
				Set::Vector tangent(normal[1],-normal[0]);

				Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Full> DDDDc = Numeric::DoubleHessian<AMREX_SPACEDIM>(c_old, i, j, k, 0, DX);
				Set::Scalar sinTheta = sin(Theta);
				Set::Scalar cosTheta = cos(Theta);
				Set::Scalar sin2Theta = sin(2.0*Theta);
				Set::Scalar cos2Theta = cos(2.0*Theta);

				Set::Scalar zeta = crack.boundary->Zeta(Theta);
				Set::Scalar ws = crack.boundary->w_phi(c_old(i,j,k,0),0.0)/(4.0*zeta*normgrad*normgrad);
				
				if( std::isnan(ws)) Util::Abort(INFO, "nan at m=",i,",",j,",",k);
				if( std::isinf(ws)) ws = 1.0E6;

				Set::Scalar Boundary_term = 0.;
				Boundary_term += crack.boundary->Gc(Theta)*crack.boundary->Dw_phi(c_old(i,j,k,0),0.)/(4.0*crack.boundary->Zeta(Theta));
				Boundary_term -= 2.0*crack.boundary->Gc(Theta)*crack.boundary->Zeta(Theta)*laplacian;

				Boundary_term += crack.boundary->DGc(Theta)
								* (zeta - ws)
								* (sin2Theta*(DDc(0,0)-DDc(1,1)) - 2.0*cos2Theta*DDc(0,1));
				Boundary_term += crack.boundary->DDGc(Theta)
								* (-zeta - ws)
								* (sinTheta*sinTheta*DDc(0,0) + cosTheta*cosTheta*DDc(1,1) - sin2Theta*DDc(0,1));

				if(std::isnan(Boundary_term)) Util::Abort(INFO,"nan at m=",i,",",j,",",k);
				df(i,j,k,1) = Boundary_term;
				rhs += Boundary_term;


				Set::Scalar Curvature_term =
						DDDDc(0,0,0,0)*(    sinTheta*sinTheta*sinTheta*sinTheta) +
						DDDDc(0,0,0,1)*(4.0*sinTheta*sinTheta*sinTheta*cosTheta) +
						DDDDc(0,0,1,1)*(6.0*sinTheta*sinTheta*cosTheta*cosTheta) +
						DDDDc(0,1,1,1)*(4.0*sinTheta*cosTheta*cosTheta*cosTheta) +
						DDDDc(1,1,1,1)*(    cosTheta*cosTheta*cosTheta*cosTheta);

				if(std::isnan(Curvature_term)) Util::Abort(INFO,"nan at m=",i,",",j,",",k);
				df(i,j,k,2) = anisotropy.beta*Curvature_term;
				rhs += anisotropy.beta*Curvature_term;
				
				
#elif AMREX_SPACEDIM == 3
				Util::Abort(INFO, "3D model hasn't been implemented yet");
#endif
			}
			} // disable this line to remove normgrad 1E-4 check
			df(i,j,k,3) = std::max(0.,rhs - crack.boundary->DrivingForceThreshold(c_old(i,j,k,0)));
			
			if(std::isnan(rhs)) Util::Abort(INFO, "Dwphi = ", crack.boundary->Dw_phi(c_old(i,j,k,0),0.),". c_old(i,j,k,0) = ",c_old(i,j,k,0));
			c_new(i,j,k,0) = c_old(i,j,k,0) - dt*std::max(0., rhs - crack.boundary->DrivingForceThreshold(c_old(i,j,k,0)))*crack.boundary->Mobility(c_old(i,j,k,0));

			if(c_new(i,j,k,0) > 1.0) {Util::Message(INFO, "cnew = ", c_new(i,j,k,0) ,", resetting to 1.0"); c_new(i,j,k,0) = 1.;}
			if(c_new(i,j,k,0) < 0.0) {Util::Message(INFO, "cnew = ", c_new(i,j,k,0) ,", resetting to 0.0"); c_new(i,j,k,0) = 0.;}
		});
	}
}

void
BrittleFracture::TimeStepComplete(amrex::Real time,int iter)
{
	if(crack.crackStressTest) 
	{
		amrex::Vector<Set::Scalar> plottime;
		amrex::Vector<int> plotstep;
		std::string plotfolder = "crack";

		plottime.resize(nlevels);
		plotstep.resize(nlevels);
		for (int lev = 0; lev < nlevels; lev++) {plottime[lev] = (double)elastic.test_step; plotstep[lev]=elastic.test_step;}
		WritePlotFile(plotfolder,plottime,plotstep);

		SetStopTime(time-0.01);return;
	}
	IntegrateVariables(time,iter);
	
	Util::Message(INFO, "crack_err_norm = ", crack_err_norm);
	Util::Message(INFO, "c_new_norm = ", c_new_norm);
	Util::Message(INFO, "relative error = ", crack_err_norm/c_new_norm);
	
	if(crack_err_norm/c_new_norm > tol_crack) return;
	
	amrex::Vector<Set::Scalar> plottime;
	amrex::Vector<int> plotstep;
	std::string plotfolder = "crack";

	plottime.resize(nlevels);
	plotstep.resize(nlevels);
	for (int lev = 0; lev < nlevels; lev++) {plottime[lev] = (double)elastic.test_step; plotstep[lev]=elastic.test_step;}
	WritePlotFile(plotfolder,plottime,plotstep);

	crack.newCrackProblem = true;
	elastic.test_step++;
	if(elastic.bc_top >= elastic.test_max) SetStopTime(time-0.01);
	
	crack_err_norm = 0.; c_new_norm = 0.;
}

void
BrittleFracture::TagCellsForRefinement (int lev, amrex::TagBoxArray &a_tags, amrex::Real /*time*/, int /*ngrow*/)
{
	const amrex::Real *DX = geom[lev].CellSize();
	const Set::Vector dx(DX);
	const Set::Scalar dxnorm = dx.lpNorm<2>();

	for (amrex::MFIter mfi(*m_c[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		const amrex::Box 							&bx 	= mfi.tilebox();
		amrex::Array4<char> const 					&tags 	= a_tags.array(mfi);
		amrex::Array4<const Set::Scalar> const 		&c_new 	= (*m_c[lev]).array(mfi);
		
		amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
			Set::Vector grad = Numeric::Gradient(c_new, i, j, k, 0, DX);
			if (dxnorm * grad.lpNorm<2>() > refinement_threshold)
				tags(i, j, k) = amrex::TagBox::SET;
		});
	}
}

void 
BrittleFracture::Integrate(int amrlev, Set::Scalar /*time*/, int /*step*/,const amrex::MFIter &mfi, const amrex::Box &box)
{
	const amrex::Real* DX = geom[amrlev].CellSize();

	amrex::Array4<const Set::Scalar> const &c_new = (*m_c[amrlev]).array(mfi);
	amrex::Array4<const Set::Scalar> const &c_old = (*m_c_old[amrlev]).array(mfi);

	amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
		crack_err_norm += ((c_new(i,j,k,0)-c_old(i,j,k,0))*(c_new(i,j,k,0)-c_old(i,j,k,0)))*(AMREX_D_TERM(DX[0],*DX[1],*DX[2]));
		c_new_norm += c_new(i,j,k,0)*c_new(i,j,k,0)*(AMREX_D_TERM(DX[0],*DX[1],*DX[2]));
	});
}
}
