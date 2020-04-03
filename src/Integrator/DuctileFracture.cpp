#include "DuctileFracture.H"

namespace Integrator
{
DuctileFracture::DuctileFracture() :
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
	
	pp_crack.query("tol_crack",tol_crack);
	pp_crack.query("tol_ep",tol_ep);

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


	// BCs
	// Crack field should have a zero neumann BC. So we just code it up here.
	// In case this needs to change, we can add options to read it from input.
	crack.mybc = new BC::Constant(1);
	crack.mybcdf = new BC::Constant(6);
	pp_crack.queryclass("bc",*static_cast<BC::Constant *>(crack.mybc));
	pp_crack.queryclass("bc_df",*static_cast<BC::Constant *>(crack.mybcdf));

	RegisterNewFab(m_c,     		crack.mybc, 	1, number_of_ghost_cells+1, "c",		true);
	RegisterNewFab(m_c_old, 		crack.mybc, 	1, number_of_ghost_cells+1, "c_old",	true);
	RegisterNewFab(m_driving_force, crack.mybcdf, 	6, number_of_ghost_cells+1, "driving_force",true);

	RegisterIntegratedVariable(&crack_err_norm, "crack_err_norm");
	RegisterIntegratedVariable(&c_new_norm,"c_new_norm");
	RegisterIntegratedVariable(&ep_norm, "Ep_norm");
	RegisterIntegratedVariable(&ep_err_norm, "Ep_err_norm");
	
	// Material input
	IO::ParmParse pp_material("material");
	pp_material.query("model",material.input_material);
	if(material.input_material == "isotropicj2plastic")
	{
		pp_material.queryclass("isotropicj2plastic",material.modeltype);
	}
	else if(material.input_material == "cubiccrystalplastic")
	{
		pp_material.queryclass("cubiccrystalplastic",material.modeltype);
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

	if (pp_elastic.countval("body_force")) pp_elastic.queryarr("body_force",elastic.body_force);

	/* Need to replace this later with a proper specification of boundary.
	   Right now we are hard-coding the tensile test and just requesting
	   rate of pulling. */
	pp_elastic.queryclass("bc",elastic.bc);
	
	pp_elastic.query("disp_step",elastic.test_rate);
	pp_elastic.query("disp_init",elastic.test_init);
	pp_elastic.query("max_disp",elastic.test_max);
	pp_elastic.query("crackStressTest",crack.crackStressTest);

	elastic.bc_top = elastic.test_init;

	if(elastic.test_rate < 0.) { Util::Warning(INFO,"Rate can't be less than zero. Resetting to 1.0"); elastic.test_rate = 0.1; }
	if(elastic.test_max < 0. ||  elastic.test_max < elastic.test_rate) {Util::Warning(INFO,"Max can't be less than load step. Resetting to load step"); elastic.test_max = elastic.test_rate;}

	const int number_of_stress_components = AMREX_SPACEDIM*AMREX_SPACEDIM;
	
	RegisterNodalFab (m_disp, 		AMREX_SPACEDIM, 				number_of_ghost_nodes, 	"Disp",			true);
	RegisterNodalFab (m_rhs,  		AMREX_SPACEDIM, 				number_of_ghost_nodes, 	"RHS",			true);
	RegisterNodalFab (m_strain,		number_of_stress_components,	number_of_ghost_nodes,	"strain",		true);
	RegisterNodalFab (m_stress,		number_of_stress_components,	number_of_ghost_nodes,	"stress",		true);
	RegisterNodalFab (m_stressdev,	number_of_stress_components,	number_of_ghost_nodes,	"stress_dev",	true);
	RegisterNodalFab (m_energy,		1,								number_of_ghost_nodes,	"energy",		true);
	RegisterNodalFab (m_energy_pristine,		1,					number_of_ghost_nodes,	"energyP",		true);
	RegisterNodalFab (m_energy_pristine_old,		1,				number_of_ghost_nodes,	"energyPOld",		true);
	RegisterNodalFab (m_residual,	AMREX_SPACEDIM,					number_of_ghost_nodes,	"residual",		true);

	RegisterNodalFab (m_strain_p,	number_of_stress_components,	number_of_ghost_nodes,	"strainp",		true);
	RegisterNodalFab (m_strain_p_old,	number_of_stress_components,	number_of_ghost_nodes,	"strainp_old",	true);
	RegisterNodalFab (m_beta,		number_of_stress_components,	number_of_ghost_nodes,	"beta",			true);
	RegisterNodalFab (m_alpha,		1,								number_of_ghost_nodes,	"strain",		true);

	nlevels = maxLevel() + 1;

	material.model.resize(nlevels);
}

DuctileFracture::~DuctileFracture(){}

void
DuctileFracture::Initialize (int lev)
{
	// Initialize crack fiel
	crack.ic->Initialize(lev,m_c);
	crack.ic->Initialize(lev,m_c_old);
	m_driving_force[lev]->setVal(0.0);

	// Initialize elastic fields
	m_disp[lev]->setVal(0.0);
	m_strain[lev]->setVal(0.0);
	m_energy[lev]->setVal(0.0);
	m_energy_pristine[lev] -> setVal(0.);
	m_energy_pristine_old[lev] -> setVal(0.);

	// Initialize stress fields
	m_stress[lev]->setVal(0.0);
	m_stressdev[lev]->setVal(0.0);
	m_rhs[lev]->setVal(0.0);
	m_residual[lev]->setVal(0.0);

	// Initialize plastic fields
	m_strain_p[lev] -> setVal(0.0);
	m_strain_p_old[lev] -> setVal(0.0);
	m_beta[lev] -> setVal(0.0);
	m_alpha[lev] -> setVal(0.0);
	
	//material.model[lev].reset(new amrex::FabArray<amrex::BaseFab<ductile_fracture_model_type>>(m_disp[lev]->boxArray(), m_disp[lev]->DistributionMap(), 1, number_of_ghost_cells));
	//material.model[lev]->setVal(*(material.modeltype));
}

void
DuctileFracture::ScaledModulus(int lev, amrex::FabArray<amrex::BaseFab<ductile_fracture_model_type> > &model)
{
	/*
	  This function is supposed to degrade material parameters based on certain
	  fracture model.
	  For now we are just using isotropic degradation.
	*/
	static amrex::IntVect AMREX_D_DECL(	dx(AMREX_D_DECL(1,0,0)), dy(AMREX_D_DECL(0,1,0)), dz(AMREX_D_DECL(0,0,1)));
	//Util::Message(INFO);

	m_c[lev]->FillBoundary();
	//Util::Message(INFO);

	for (amrex::MFIter mfi(model,amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		//amrex::Box box = mfi.growntilebox(2);
		amrex::Box box = mfi.tilebox();
		box.grow(1);
		amrex::Array4<const amrex::Real> 			const& c_new	= (*m_c[lev]).array(mfi);
		amrex::Array4<ductile_fracture_model_type> 	const& modelfab	= model.array(mfi);
		amrex::Array4<const Set::Scalar> 			const& p_box	= (*m_alpha[lev]).array(mfi);

		amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k){
			Set::Scalar mul = AMREX_D_PICK(0.5, 0.25, 0.125);
			amrex::Vector<Set::Scalar> _temp;
			_temp.push_back( mul*(AMREX_D_TERM(	
								crack.boundary->g_phi(c_new(i,j,k,0),Numeric::Interpolate::NodeToCellAverage(p_box,i,j,k,0)) 
								+ crack.boundary->g_phi(c_new(i-1,j,k,0),Numeric::Interpolate::NodeToCellAverage(p_box,i-1,j,k,0))
								, 
								+ crack.boundary->g_phi(c_new(i,j-1,k,0),Numeric::Interpolate::NodeToCellAverage(p_box,i,j-1,k,0)) 
								+ crack.boundary->g_phi(c_new(i-1,j-1,k,0),Numeric::Interpolate::NodeToCellAverage(p_box,i-1,j-1,k,0))
								, 
								+ crack.boundary->g_phi(c_new(i,j,k-1,0),Numeric::Interpolate::NodeToCellAverage(p_box,i,j,k-1,0)) 
								+ crack.boundary->g_phi(c_new(i-1,j,k-1,0),Numeric::Interpolate::NodeToCellAverage(p_box,i-1,j,k-1,0))
								+ crack.boundary->g_phi(c_new(i,j-1,k-1,0),Numeric::Interpolate::NodeToCellAverage(p_box,i,j-1,k-1,0)) 
								+ crack.boundary->g_phi(c_new(i-1,j-1,k-1,0),Numeric::Interpolate::NodeToCellAverage(p_box,i-1,j-1,k-1,0))
								)));
			if(_temp[0] < 0.0) _temp[0] = 0.;
			if(_temp[0] > 1.0) _temp[0] = 1.0;
			modelfab(i,j,k,0).DegradeModulus(std::min(1.-_temp[0],1.-scaleModulusMax));
			modelfab(i,j,k,0).DegradeYieldSurface(std::min(1.-_temp[0],1.-scaleModulusMax));
		});
	//Util::Message(INFO);

	}
	amrex::Geometry tmp_geom = geom[lev];
    for (int i = 0; i < 2; i++)
    {
		amrex::FabArray<amrex::BaseFab<ductile_fracture_model_type>> &mf = model;
        mf.FillBoundary(tmp_geom.periodicity());
        mf.FillBoundary();
        const int ncomp = mf.nComp();
        const int ng1 = 1;
        const int ng2 = 2;
        amrex::FabArray<amrex::BaseFab<ductile_fracture_model_type>> tmpmf(mf.boxArray(), mf.DistributionMap(), ncomp, ng1);
        amrex::Copy(tmpmf, mf, 0, 0, ncomp, ng1);
        mf.ParallelCopy(tmpmf, 0, 0, ncomp, ng1, ng2, tmp_geom.periodicity());
    }
	//Util::Message(INFO);
}

void 
DuctileFracture::TimeStepBegin(amrex::Real /*time*/, int /*iter*/)
{
	Util::Message(INFO);
	if(crack.newCrackProblem)
	{
		Util::Message(INFO);
		elastic.bc_top = elastic.test_init + ((double)elastic.test_step)*elastic.test_rate;
		crack.newCrackProblem = false;
	}
	//elastic.bc_top = elastic.test_init + ((double)elastic.test_step)*elastic.test_rate;

	material.model.resize(nlevels);
	material.modeltype.UpdateF0(Set::Matrix::Zero());

	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		material.model[ilev].reset(new amrex::FabArray<amrex::BaseFab<ductile_fracture_model_type>>(m_disp[ilev]->boxArray(), m_disp[ilev]->DistributionMap(), 1, 2));
		material.model[ilev]->setVal((material.modeltype));

		std::swap(*m_strain_p_old[ilev], *m_strain_p[ilev]);
		std::swap(*m_energy_pristine_old[ilev], *m_energy_pristine[ilev]);
		
		m_energy_pristine[ilev]->setVal(0.);
		m_strain_p[ilev]->setVal(0.);
		ScaledModulus(ilev,*(material.model)[ilev]);
	}
	
	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		m_strain_p[ilev]->FillBoundary();
		m_alpha[ilev]->FillBoundary();
		m_beta[ilev]->FillBoundary();
		Set::Vector DX(geom[ilev].CellSize());

		for (MFIter mfi(*(material.model)[ilev], false); mfi.isValid(); ++mfi)
		{
			amrex::Box bx = mfi.growntilebox(2);
			amrex::Array4<ductile_fracture_model_type>	const &model_box	= material.model[ilev]->array(mfi);
			amrex::Array4<const Set::Scalar>			const &strain_p_box	= (*m_strain_p[ilev]).array(mfi);
			amrex::Array4<const Set::Scalar>			const &beta_box		= (*m_beta[ilev]).array(mfi);
			amrex::Array4<const Set::Scalar>			const &alpha_box	= (*m_alpha[ilev]).array(mfi);

			amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
			{
				Set::Matrix eps, beta;
				AMREX_D_PICK(
					eps(0,0) = strain_p_box(i,j,k,0);
					beta(0,0) = beta_box(i,j,k,0);
					,
					eps(0,0) = strain_p_box(i,j,k,0); eps(0,1) = strain_p_box(i,j,k,1);
					eps(1,0) = strain_p_box(i,j,k,2); eps(1,1) = strain_p_box(i,j,k,3);
					beta(0,0) = beta_box(i,j,k,0); beta(0,1) = beta_box(i,j,k,1);
					beta(1,0) = beta_box(i,j,k,2); beta(1,1) = beta_box(i,j,k,3);
					,
					eps(0,0) = strain_p_box(i,j,k,0); eps(0,1) = strain_p_box(i,j,k,1); eps(0,2) = strain_p_box(i,j,k,2);
					eps(1,0) = strain_p_box(i,j,k,3); eps(1,1) = strain_p_box(i,j,k,4); eps(1,2) = strain_p_box(i,j,k,5);
					eps(2,0) = strain_p_box(i,j,k,6); eps(2,1) = strain_p_box(i,j,k,7); eps(2,2) = strain_p_box(i,j,k,8);
					beta(0,0) = beta_box(i,j,k,0); beta(0,1) = beta_box(i,j,k,1); beta(0,2) = beta_box(i,j,k,2);
					beta(1,0) = beta_box(i,j,k,3); beta(1,1) = beta_box(i,j,k,4); beta(1,2) = beta_box(i,j,k,5);
					beta(2,0) = beta_box(i,j,k,6); beta(2,1) = beta_box(i,j,k,7); beta(2,2) = beta_box(i,j,k,8);
				);
				Model::Solid::Affine::PlasticState state;
				state.epsp = eps; state.beta = beta; state.alpha = alpha_box(i,j,k,0);
				model_box(i, j, k, 0).SetPlasticStrains(state);
            });
		}
	}
	//Util::Message(INFO);

	for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		const Real* DX = geom[ilev].CellSize();
		Set::Scalar volume = AMREX_D_TERM(DX[0],*DX[1],*DX[2]);

		AMREX_D_TERM(m_rhs[ilev]->setVal(elastic.body_force[0]*volume,0,1);,
			     m_rhs[ilev]->setVal(elastic.body_force[1]*volume,1,1);,
			     m_rhs[ilev]->setVal(elastic.body_force[2]*volume,2,1););
	}
	
	elastic.bc.Set(elastic.bc.Face::YHI, elastic.bc.Direction::Y, BC::Operator::Elastic<ductile_fracture_model_type>::Type::Displacement, elastic.bc_top);
	elastic.bc.Set(elastic.bc.Face::XLO_YHI, elastic.bc.Direction::Y, BC::Operator::Elastic<ductile_fracture_model_type>::Type::Displacement, elastic.bc_top);
	elastic.bc.Set(elastic.bc.Face::XHI_YHI, elastic.bc.Direction::Y, BC::Operator::Elastic<ductile_fracture_model_type>::Type::Displacement, elastic.bc_top);
	elastic.bc.Init(m_rhs,geom);
	
	LPInfo info;
	info.setAgglomeration(elastic.agglomeration);
	info.setConsolidation(elastic.consolidation);
	info.setMaxCoarseningLevel(elastic.max_coarsening_level);

	Operator::Elastic<ductile_fracture_model_type> elastic_op;
	elastic_op.define(geom, grids, dmap, info);
	elastic_op.setMaxOrder(elastic.linop_maxorder);
	elastic_op.SetBC(&(elastic.bc));
	
	Solver::Nonlocal::Newton<ductile_fracture_model_type> solver(elastic_op);

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

	solver.solve(m_disp,m_rhs,material.model,elastic.tol_rel,elastic.tol_abs);

	for (int lev = 0; lev < nlevels; lev++)
	{
		elastic_op.Strain(lev,*m_strain[lev],*m_disp[lev]);
		elastic_op.Stress(lev,*m_stress[lev],*m_disp[lev]);
		elastic_op.Energy(lev,*m_energy[lev],*m_disp[lev]);

		for (amrex::MFIter mfi(*m_strain[lev],true); mfi.isValid(); ++mfi)
		{
			const amrex::Box& box = mfi.validbox();
			amrex::Array4<Set::Scalar>					const& sig_box 		= (*m_stress[lev]).array(mfi);
			amrex::Array4<const Set::Scalar>			const& eps_box 		= (*m_strain[lev]).array(mfi);
			amrex::Array4<Set::Scalar>					const& strainp_box 	= (*m_strain_p[lev]).array(mfi);
			amrex::Array4<Set::Scalar>					const& beta_box 	= (*m_beta[lev]).array(mfi);
			amrex::Array4<Set::Scalar>					const& alpha_box 	= (*m_alpha[lev]).array(mfi);
			amrex::Array4<Set::Scalar>					const& sigdev_box 	= (*m_stressdev[lev]).array(mfi);
			amrex::Array4<Set::Scalar>					const& energy_box 	= (*m_energy_pristine[lev]).array(mfi);
			//amrex::Array4<Set::Scalar>					const& energy_box_old 	= (*m_energy_pristine_old[lev]).array(mfi);
			amrex::Array4<ductile_fracture_model_type>	const& model_box 	= material.model[lev]->array(mfi);

			amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k){
				Set::Matrix eps, sig, sigdev;
				AMREX_D_PICK(
					sig(0,0) = sig_box(i,j,k,0);
					eps(0,0) = eps_box(i,j,k,0);
					,
					sig(0,0) = sig_box(i,j,k,0); sig(0,1) = sig_box(i,j,k,1);
					sig(1,0) = sig_box(i,j,k,2); sig(1,1) = sig_box(i,j,k,3);
					eps(0,0) = eps_box(i,j,k,0); eps(0,1) = eps_box(i,j,k,1);
					eps(1,0) = eps_box(i,j,k,2); eps(1,1) = eps_box(i,j,k,3);
					,
					sig(0,0) = sig_box(i,j,k,0); sig(0,1) = sig_box(i,j,k,1); sig(0,2) = sig_box(i,j,k,2);
					sig(1,0) = sig_box(i,j,k,3); sig(1,1) = sig_box(i,j,k,4); sig(1,2) = sig_box(i,j,k,5);
					sig(2,0) = sig_box(i,j,k,6); sig(2,1) = sig_box(i,j,k,7); sig(2,2) = sig_box(i,j,k,8);
					eps(0,0) = eps_box(i,j,k,0); eps(0,1) = eps_box(i,j,k,1); eps(0,2) = eps_box(i,j,k,2);
					eps(1,0) = eps_box(i,j,k,3); eps(1,1) = eps_box(i,j,k,4); eps(1,2) = eps_box(i,j,k,5);
					eps(2,0) = eps_box(i,j,k,6); eps(2,1) = eps_box(i,j,k,7); eps(2,2) = eps_box(i,j,k,8);
				);
				
				sigdev = sig - (1.0/((double)AMREX_SPACEDIM))*sig.trace()*Set::Matrix::Identity();
				
				AMREX_D_PICK(
					sigdev_box(i,j,k,0) = sigdev(0,0);
					,
					sigdev_box(i,j,k,0) = sigdev(0,0); sigdev_box(i,j,k,1) = sigdev(0,1);
					sigdev_box(i,j,k,2) = sigdev(1,0); sigdev_box(i,j,k,3) = sigdev(1,1);
					,
					sigdev_box(i,j,k,0) = sigdev(0,0); sigdev_box(i,j,k,1) = sigdev(0,1); sigdev_box(i,j,k,2) = sigdev(0,2);
					sigdev_box(i,j,k,3) = sigdev(1,0); sigdev_box(i,j,k,4) = sigdev(1,1); sigdev_box(i,j,k,5) = sigdev(1,2);
					sigdev_box(i,j,k,6) = sigdev(2,0); sigdev_box(i,j,k,7) = sigdev(2,1); sigdev_box(i,j,k,8) = sigdev(2,2);
				);
				//energy_box(i,j,k,0) > energy_box_old(i,j,k,0) ? energy_box(i,j,k,0) : energy_box_old(i,j,k,0);

				model_box(i,j,k,0).EvolvePlasticStrain(sig,eps,0);
				sig = model_box(i,j,k,0).DW(eps);
				
				AMREX_D_PICK(
					strainp_box(i,j,k,0) = model_box(i,j,k,0).curr.epsp(0,0);
					beta_box(i,j,k,0) = model_box(i,j,k,0).curr.beta(0,0);
					sig_box(i,j,k,0) = sig(0,0);
					,
					strainp_box(i,j,k,0) = model_box(i,j,k,0).curr.epsp(0,0); strainp_box(i,j,k,1) = model_box(i,j,k,0).curr.epsp(0,1);
					strainp_box(i,j,k,2) = model_box(i,j,k,0).curr.epsp(1,0); strainp_box(i,j,k,3) = model_box(i,j,k,0).curr.epsp(1,1);
					beta_box(i,j,k,0) = model_box(i,j,k,0).curr.beta(0,0); beta_box(i,j,k,1) = model_box(i,j,k,0).curr.beta(0,1);
					beta_box(i,j,k,2) = model_box(i,j,k,0).curr.beta(1,0); beta_box(i,j,k,3) = model_box(i,j,k,0).curr.beta(1,1);
					sig_box(i,j,k,0) = sig(0,0); sig_box(i,j,k,1) = sig(0,1);
					sig_box(i,j,k,2) = sig(1,0); sig_box(i,j,k,3) = sig(1,1);
					,
					strainp_box(i,j,k,0) = model_box(i,j,k,0).curr.epsp(0,0); strainp_box(i,j,k,1) = model_box(i,j,k,0).curr.epsp(0,1); strainp_box(i,j,k,2) = model_box(i,j,k,0).curr.epsp(0,2);
					strainp_box(i,j,k,3) = model_box(i,j,k,0).curr.epsp(1,0); strainp_box(i,j,k,4) = model_box(i,j,k,0).curr.epsp(1,1); strainp_box(i,j,k,5) = model_box(i,j,k,0).curr.epsp(1,2);
					strainp_box(i,j,k,6) = model_box(i,j,k,0).curr.epsp(2,0); strainp_box(i,j,k,7) = model_box(i,j,k,0).curr.epsp(2,1); strainp_box(i,j,k,8) = model_box(i,j,k,0).curr.epsp(2,2);
					beta_box(i,j,k,0) = model_box(i,j,k,0).curr.beta(0,0); beta_box(i,j,k,1) = model_box(i,j,k,0).curr.beta(0,1); beta_box(i,j,k,2) = model_box(i,j,k,0).curr.beta(0,2);
					beta_box(i,j,k,3) = model_box(i,j,k,0).curr.beta(1,0); beta_box(i,j,k,4) = model_box(i,j,k,0).curr.beta(1,1); beta_box(i,j,k,5) = model_box(i,j,k,0).curr.beta(1,2);
					beta_box(i,j,k,6) = model_box(i,j,k,0).curr.beta(2,0); beta_box(i,j,k,7) = model_box(i,j,k,0).curr.beta(2,1); beta_box(i,j,k,8) = model_box(i,j,k,0).curr.beta(2,2);
					sig_box(i,j,k,0) = sig(0,0); sig_box(i,j,k,1) = sig(0,1); sig_box(i,j,k,2) = sig(0,2);
					sig_box(i,j,k,3) = sig(1,0); sig_box(i,j,k,4) = sig(1,1); sig_box(i,j,k,5) = sig(1,2);
					sig_box(i,j,k,6) = sig(2,0); sig_box(i,j,k,7) = sig(2,1); sig_box(i,j,k,8) = sig(2,2);
				);
				alpha_box(i,j,k,0) = model_box(i,j,k,0).curr.alpha;

				material.modeltype.UpdateF0(model_box(i,j,k).curr.epsp);
				energy_box(i,j,k,0) = material.modeltype.W(eps);
				//energy_box(i,j,k,0) = energy_box(i,j,k,0) > energy_box_old(i,j,k,0) ? energy_box(i,j,k,0) : energy_box_old(i,j,k,0);
			});
		}

		elastic_op.Strain(lev,*m_strain[lev],*m_disp[lev]);
		//elastic_op.Stress(lev,*m_stress[lev],*m_disp[lev]);
		elastic_op.Energy(lev,*m_energy[lev],*m_disp[lev]);
	}
	Util::Message(INFO);
}

void
DuctileFracture::Advance (int lev, amrex::Real /*time*/, amrex::Real dt)
{
	static amrex::IntVect AMREX_D_DECL(	dx(AMREX_D_DECL(1,0,0)), dy(AMREX_D_DECL(0,1,0)), dz(AMREX_D_DECL(0,0,1)));
	const amrex::Real* DX = geom[lev].CellSize();

	std::swap(*m_c[lev],*m_c_old[lev]);

	for ( amrex::MFIter mfi(*m_c[lev],true); mfi.isValid(); ++mfi )
	{
		const amrex::Box& box = mfi.validbox();
		amrex::Array4<Set::Scalar>			const& c_new		= (*m_c[lev]).array(mfi);
		amrex::Array4<const Set::Scalar>	const& c_old		= (*m_c_old[lev]).array(mfi);
		amrex::Array4<const Set::Scalar> 	const& p_box		= (*m_alpha[lev]).array(mfi);
		amrex::Array4<Set::Scalar>			const& df			= (*m_driving_force[lev]).array(mfi);
		amrex::Array4<const Set::Scalar>	const& energy_box	= (*m_energy_pristine[lev]).array(mfi);
		
		amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k){
			
			Set::Scalar laplacian = Numeric::Laplacian(c_old,i,j,k,0,DX);
			Set::Scalar rhs = 0.;	
			Set::Scalar en_cell = Numeric::Interpolate::NodeToCellAverage(energy_box,i,j,k,0);

			df(i,j,k,0) = crack.boundary->Dg_phi(c_old(i,j,k,0),Numeric::Interpolate::NodeToCellAverage(p_box,i,j,k,0))*en_cell;
			df(i,j,k,1) = crack.boundary->Epc(c_old(i,j,k,0))*crack.boundary->Dw_phi(c_old(i,j,k,0),Numeric::Interpolate::NodeToCellAverage(p_box,i,j,k,0));
			df(i,j,k,2) = crack.boundary->kappa(c_old(i,j,k,0))*laplacian;
			df(i,j,k,3) = crack.boundary->Dg_phi(c_old(i,j,k,0),Numeric::Interpolate::NodeToCellAverage(p_box,i,j,k,0))*(material.modeltype.OriginalPlasticEnergy(Numeric::Interpolate::NodeToCellAverage(p_box,i,j,k,0)));

			rhs += crack.boundary->Dg_phi(c_old(i,j,k,0),Numeric::Interpolate::NodeToCellAverage(p_box,i,j,k,0))*en_cell;
			rhs += crack.boundary->Epc(c_old(i,j,k,0))*crack.boundary->Dw_phi(c_old(i,j,k,0),Numeric::Interpolate::NodeToCellAverage(p_box,i,j,k,0));
			rhs -= crack.boundary->kappa(c_old(i,j,k,0))*laplacian;
			rhs += crack.boundary->Dg_phi(c_old(i,j,k,0),Numeric::Interpolate::NodeToCellAverage(p_box,i,j,k,0))*(material.modeltype.OriginalPlasticEnergy(Numeric::Interpolate::NodeToCellAverage(p_box,i,j,k,0)));
			
			df(i,j,k,4) = rhs;
			df(i,j,k,5) = max(0.,rhs-crack.boundary->DrivingForceThreshold(c_old(i,j,k,0)));
			
			c_new(i,j,k,0) = c_old(i,j,k,0) - dt*std::max(0.,(rhs-crack.boundary->DrivingForceThreshold(c_old(i,j,k,0))))*crack.boundary->GetMobility(c_old(i,j,k,0));

			if(c_new(i,j,k,0) > 1.0) {Util::Warning(INFO, "cnew exceeded 1.0, resetting to 1.0"); c_new(i,j,k,0) = 1.;}
			if(c_new(i,j,k,0) < 0.0) {Util::Warning(INFO, "cnew is below 0.0, resetting to 0.0"); c_new(i,j,k,0) = 0.;}
		});
	}
}

void
DuctileFracture::TimeStepComplete(amrex::Real time,int iter)
{
	//Util::Message(INFO);
	IntegrateVariables(time,iter);
	
	Set::Scalar rel_c_err = crack_err_norm/c_new_norm;
	Set::Scalar rel_ep_err = ep_err_norm/ep_norm;
	Util::Message(INFO, "crack_err_norm = ", crack_err_norm);
	Util::Message(INFO, "c_new_norm = ", c_new_norm);
	Util::Message(INFO, "plastic strain norm = ", ep_norm);
	Util::Message(INFO, "plastic strain error norm = ", ep_err_norm);
	Util::Message(INFO, "crack relative error= ", rel_c_err);
	Util::Message(INFO, "plastic strain relative error= ", rel_ep_err);
	
	if(rel_c_err > tol_crack) return;
	if(rel_ep_err > tol_ep) return;

	crack_err_norm = 0.; c_new_norm = 0.;
	ep_err_norm = 0.; ep_norm = 0.;
	
	/*amrex::Vector<Set::Scalar> plottime;
	amrex::Vector<int> plotstep;
	std::string plotfolder = "crack";

	plottime.resize(nlevels);
	plotstep.resize(nlevels);
	for (int lev = 0; lev < nlevels; lev++) {plottime[lev] = (double)elastic.test_step; plotstep[lev]=elastic.test_step;}
	WritePlotFile(plotfolder,plottime,plotstep);*/

	crack.newCrackProblem = true;
	elastic.test_step++;
	if(elastic.bc_top >= elastic.test_max) SetStopTime(time-0.01);
	
	//crack_err_norm = 0.; c_new_norm = 0.;
}

void
DuctileFracture::TagCellsForRefinement (int lev, amrex::TagBoxArray& a_tags, amrex::Real /*time*/, int /*ngrow*/)
{
	Set::Vector DX(geom[lev].CellSize());
    Set::Scalar DXnorm = DX.lpNorm<2>();
	a_tags.setVal(amrex::TagBox::CLEAR);

	for (amrex::MFIter mfi(*m_c[lev],TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		const amrex::Box&  bx  = mfi.tilebox();
		amrex::Array4<char> const &tags = a_tags.array(mfi);
		amrex::Array4<Set::Scalar> const &c_new = (*m_c[lev]).array(mfi);;

		amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
			Set::Vector grad = Numeric::Gradient(c_new, i, j, k, 0, DX.data());
            if (grad.lpNorm<2>() * DXnorm > 0.01)
				tags(i, j, k) = amrex::TagBox::SET;
		});
	}
}

void 
DuctileFracture::Integrate(int amrlev, Set::Scalar /*time*/, int /*step*/,const amrex::MFIter &mfi, const amrex::Box &box)
{
	const amrex::Real* DX = geom[amrlev].CellSize();

	amrex::Array4<const Set::Scalar> const &c_new = (*m_c[amrlev]).array(mfi);
	amrex::Array4<const Set::Scalar> const &c_old = (*m_c_old[amrlev]).array(mfi);
	amrex::Array4<const Set::Scalar> const &ep_old = (*m_strain_p[amrlev]).array(mfi);
	amrex::Array4<const Set::Scalar> const &ep_new = (*m_strain_p_old[amrlev]).array(mfi);

	amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
		crack_err_norm += ((c_new(i,j,k,0)-c_old(i,j,k,0))*(c_new(i,j,k,0)-c_old(i,j,k,0)))*(AMREX_D_TERM(DX[0],*DX[1],*DX[2]));
		c_new_norm += c_new(i,j,k,0)*c_new(i,j,k,0)*(AMREX_D_TERM(DX[0],*DX[1],*DX[2]));

		for (int n=0; n<AMREX_SPACEDIM*AMREX_SPACEDIM; n++)
		{
			ep_norm += ep_new(i,j,k,n)*ep_new(i,j,k,n)*(AMREX_D_TERM(DX[0],*DX[1],*DX[2]));
			ep_err_norm += (ep_new(i,j,k,n)-ep_old(i,j,k,n))*(ep_new(i,j,k,n)-ep_old(i,j,k,n))*(AMREX_D_TERM(DX[0],*DX[1],*DX[2]));
			//ep_norm += (Numeric::Interpolate::NodeToCellAverage(ep_new,i,j,k,n)*Numeric::Interpolate::NodeToCellAverage(ep_new,i,j,k,n))*(AMREX_D_TERM(DX[0],*DX[1],*DX[2]));
			//ep_err_norm += 	(Numeric::Interpolate::NodeToCellAverage(ep_new,i,j,k,n)
			//				- Numeric::Interpolate::NodeToCellAverage(ep_old,i,j,k,n))
			//				*(Numeric::Interpolate::NodeToCellAverage(ep_new,i,j,k,n)
			//				- Numeric::Interpolate::NodeToCellAverage(ep_old,i,j,k,n))
			//				*(AMREX_D_TERM(DX[0],*DX[1],*DX[2]));
		}
		
	});
}

}