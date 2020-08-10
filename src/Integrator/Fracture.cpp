#include "Fracture.H"

namespace Integrator
{
Fracture::Fracture() :
	Integrator()
{
    //==================================================
    // Problem type
    {
        std::string fracture_problem;
        IO::ParmParse pp("fracture");
        pp.query("problem_type", fracture_problem);

        if (fracture_problem == "brittle") fracture_type = FractureType::Brittle;
        if (fracture_problem == "ductile") fracture_type = FractureType::Ductile;
        
    }
    //==================================================

    //==================================================
    // Crack model
    {
        IO::ParmParse pp_crack("crack");
        std::string crack_type;
        pp_crack.query("type",crack_type);
        pp_crack.query("modulus_scaling_max",crack.scaleModulusMax);
        pp_crack.query("refinement_threshold",crack.refinement_threshold);
        pp_crack.query("tol_rel", crack.tol_rel);
        pp_crack.query("tol_abs", crack.tol_abs);

        if(crack_type=="constant")
        {
            Model::Interface::Crack::Constant *tmpbdy = new Model::Interface::Crack::Constant();
            pp_crack.queryclass("constant",*tmpbdy);
            crack.cracktype = tmpbdy;
        }
        else if(crack_type == "sin")
        {
            Model::Interface::Crack::Sin *tmpbdy = new Model::Interface::Crack::Sin();
            pp_crack.queryclass("sin", *tmpbdy);
            crack.cracktype = tmpbdy;
        }
        else
            Util::Abort(INFO,"This crack model hasn't been implemented yet");

        IO::ParmParse pp_crack_df("crack.df");
        pp_crack_df.query("mult_Gc", crack.mult_df_Gc);
        pp_crack_df.query("mult_Lap", crack.mult_df_lap);
    }
    //==================================================

    //==================================================
    // Anistropy
    {
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
    }
    //==================================================

    //==================================================
    // ICs
    {
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
    }
    //==================================================

    //==================================================
    // Crack BCs
    {
        // Crack field should have a zero neumann BC. So we just code it up here.
        // In case this needs to change, we can add options to read it from input.
        // Below are conditions for full simulation.  If this doesn't work, we can
        // try symmetric simulation.
        crack.bc = new BC::Constant(1);
        crack.bcdf = new BC::Constant(4);
        IO::ParmParse pp_crack("crack");
        pp_crack.queryclass("bc",*static_cast<BC::Constant *>(crack.bc));
        pp_crack.queryclass("bc_df",*static_cast<BC::Constant *>(crack.bcdf));

        RegisterNewFab(crack.field,     crack.bc, 1, number_of_ghost_cells, "c",		true);
        RegisterNewFab(crack.field_old, crack.bc, 1, number_of_ghost_cells, "c_old",	true);
        switch (fracture_type)
        {
            case FractureType::Brittle: 
                RegisterNewFab(crack.driving_force, crack.bcdf, 4, number_of_ghost_cells, "driving_force",true);
                break;
            case FractureType::Ductile:
                RegisterNewFab(crack.driving_force, crack.bcdf, 6, number_of_ghost_cells, "driving_force",true);
                break;
            default:
                break;
        }
        
        RegisterIntegratedVariable(&crack.error_norm, "crack_err_norm");
	    RegisterIntegratedVariable(&crack.norm,"c_new_norm");
    }
    //==================================================

    //==================================================
    // Material model
    {
        IO::ParmParse pp_material("material");
        pp_material.query("model",material.input_material);

        switch (fracture_type)
        {
            case FractureType::Brittle:
                if(material.input_material == "isotropic") pp_material.queryclass("isotropic",material.brittlemodeltype);
                else Util::Abort(INFO,"This model has not been implemented yet.");
                break;
            case FractureType::Ductile:
                if(material.input_material == "isotropicj2plastic") pp_material.queryclass("isotropicj2plastic",material.ductilemodeltype);
                else if(material.input_material == "cubiccrystalplastic")pp_material.queryclass("cubiccrystalplastic",material.ductilemodeltype);
                else Util::Abort(INFO,"This model has not been implemented yet.");
            default:
                break;
        }       
    }
    //==================================================

    //==================================================
    // Solver properties
    {
        IO::ParmParse pp_elastic("solver");
        pp_elastic.query("int",				sol.interval);
        pp_elastic.query("type",			sol.type);
        pp_elastic.query("max_iter",		sol.max_iter);
        pp_elastic.query("max_fmg_iter",	sol.max_fmg_iter);
        pp_elastic.query("verbose",			sol.verbose);
        pp_elastic.query("cgverbose",		sol.cgverbose);
        pp_elastic.query("tol_rel",			sol.tol_rel);
        pp_elastic.query("tol_abs",			sol.tol_abs);
        pp_elastic.query("cg_tol_rel",		sol.cg_tol_rel);
        pp_elastic.query("cg_tol_abs",		sol.cg_tol_abs);
        pp_elastic.query("use_fsmooth",		sol.use_fsmooth);
        pp_elastic.query("agglomeration", 	sol.agglomeration);
        pp_elastic.query("consolidation", 	sol.consolidation);

        pp_elastic.query("bottom_solver",       sol.bottom_solver);
        pp_elastic.query("linop_maxorder",      sol.linop_maxorder);
        pp_elastic.query("max_coarsening_level",sol.max_coarsening_level);
        pp_elastic.query("verbose",             sol.verbose);
        pp_elastic.query("cg_verbose",          sol.cgverbose);
        pp_elastic.query("bottom_max_iter",     sol.bottom_max_iter);
        pp_elastic.query("max_fixed_iter",      sol.max_fixed_iter);
        pp_elastic.query("bottom_tol",          sol.bottom_tol);
    }
    //==================================================

    //==================================================
    // Loading conditions and multipliers
    {
        IO::ParmParse pp_elastic("elastic");
        pp_elastic.query("df_mult", elastic.df_mult);

        std::string loading_mode, loading_type;
        IO::ParmParse pp_load("loading");
        if (pp_load.countval("body_force")) pp_load.queryarr("body_force",loading.body_force);
        pp_load.query("mode",loading_mode);
        pp_load.query("type",loading_type);
        pp_load.query("init", loading.init);
        pp_load.query("max", loading.max);
        pp_load.query("rate", loading.rate);
        loading.val = loading.init;

        if(loading_mode == "modeI" || loading_mode == "mode1") loading.mode = ModeType::ModeI;
        if(loading_mode == "modeII" || loading_mode == "mode2") loading.mode = ModeType::ModeII;
        if(loading_mode == "modeIII" || loading_mode == "mode3") loading.mode = ModeType::ModeIII;

        if(loading_type == "force") loading.load = LoadType::Force;
        if(loading_type == "disp" || loading_type == "displacement") loading.load = LoadType::Displacement;

        if( fracture_type == FractureType::Brittle ) pp_elastic.queryclass("bc",elastic.brittlebc);
        else if (fracture_type == FractureType::Ductile) pp_elastic.queryclass("bc",elastic.ductilebc);
    }
    //==================================================

    //==================================================
    // Plastic variables if needed
    {
        IO::ParmParse pp_plastic("plastic");
        pp_plastic.query("tol_rel", plastic.tol_rel);
        pp_plastic.query("tol_abs", plastic.tol_abs);
    }

    //==================================================
    // Registering fabs now
    {
        nlevels = maxLevel() + 1;
        RegisterNodalFab(elastic.disp,  AMREX_SPACEDIM, number_of_ghost_nodes, "disp", true);
        RegisterNodalFab(elastic.rhs,  AMREX_SPACEDIM, number_of_ghost_nodes, "rhs", true);
        RegisterNodalFab(elastic.residual,  AMREX_SPACEDIM, number_of_ghost_nodes, "res", true);
        RegisterNodalFab(elastic.strain,  AMREX_SPACEDIM*AMREX_SPACEDIM, number_of_ghost_nodes, "strain", true);
        RegisterNodalFab(elastic.stress,  AMREX_SPACEDIM*AMREX_SPACEDIM, number_of_ghost_nodes, "stress", true);
        RegisterNodalFab(elastic.energy, 1, number_of_ghost_nodes, "energy", true);
        RegisterNodalFab(elastic.energy_pristine, 1, number_of_ghost_nodes, "energy_pristine", true);
        RegisterNodalFab(elastic.energy_pristine_old, 1, number_of_ghost_nodes, "energy_pristine_old", true);

        if (fracture_type == FractureType::Brittle)
        {
            RegisterGeneralFab(material.brittlemodel, 1, number_of_ghost_nodes);
            material.brittlemodel.resize(nlevels);
        }
        else if (fracture_type == FractureType::Ductile)
        {
            RegisterNodalFab(plastic.strain, AMREX_SPACEDIM*AMREX_SPACEDIM, number_of_ghost_nodes, "strainp", true);
            RegisterNodalFab(plastic.strain_old, AMREX_SPACEDIM*AMREX_SPACEDIM, number_of_ghost_nodes, "strainp_old", true);
            RegisterGeneralFab(material.ductilemodel, 1, number_of_ghost_nodes);
            RegisterIntegratedVariable(&plastic.error_norm, "plastic_err_norm");
            RegisterIntegratedVariable(&plastic.norm, "plastic_norm");
            material.ductilemodel.resize(nlevels);
        }
    }
    //==================================================
}

Fracture::~Fracture()
{
}

void
Fracture::Initialize (int ilev)
{
    //==================================================
    // Initialization of crack fields
    {
        crack.ic->Initialize(ilev,crack.field);
        crack.ic->Initialize(ilev,crack.field_old);
        
        crack.driving_force[ilev]->setVal(0.0);
    }
    //==================================================
	
    //==================================================
    // Initialization of elastic fields
    {
        elastic.disp[ilev]->setVal(0.0);
        elastic.strain[ilev]->setVal(0.0);
        elastic.stress[ilev]->setVal(0.0);
        elastic.rhs[ilev]->setVal(0.0);
        elastic.energy[ilev]->setVal(0.0);
        elastic.residual[ilev]->setVal(0.0);
        elastic.energy_pristine[ilev] -> setVal(0.);
        elastic.energy_pristine_old[ilev] -> setVal(0.);
    }
	//==================================================

    //==================================================
    // Initialization of brittle and ductile specific fields
    {
        if (fracture_type == FractureType::Brittle)
            material.brittlemodel[ilev]->setVal(material.brittlemodeltype);
        
        else if (fracture_type == FractureType::Ductile)
        {
            plastic.strain[ilev] -> setVal(0.);
            plastic.strain_old[ilev] -> setVal(0.);
            material.ductilemodel[ilev] -> setVal(material.ductilemodeltype);
        }
    }
    //==================================================
}

void
Fracture::TimeStepBegin(amrex::Real time, int iter)
{
    if (anisotropy.on && time >= anisotropy.tstart)
	{
		SetTimestep(anisotropy.timestep);
		if (anisotropy.elastic_int > 0) 
			if (iter % anisotropy.elastic_int) return;
	}
    if(iter%sol.interval) return;
    loading.val = loading.init + ((double)loading.step)*loading.rate;

    for (int ilev = 0; ilev < nlevels; ++ilev)
	{
		std::swap(*elastic.energy_pristine_old[ilev], *elastic.energy_pristine[ilev]);
        if (fracture_type == FractureType::Ductile) std::swap(*plastic.strain[ilev], *plastic.strain_old[ilev]);
		crack.field[ilev]->FillBoundary();
        if (fracture_type == FractureType::Brittle) material.brittlemodel[ilev]->setVal(material.brittlemodeltype);
        //else material.ductilemodel[ilev]->setVal(material.ductilemodeltype);
    }

    //==================================================
    // Scaling the modulus
    {
        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            for (amrex::MFIter mfi(*elastic.disp[ilev],true); mfi.isValid(); ++mfi)
            {
                amrex::Box box = mfi.growntilebox(2);
                amrex::Array4<const amrex::Real> const& c_new = (*crack.field[ilev]).array(mfi);

                amrex::Array4<brittle_fracture_model_type> modelfab_b;
                amrex::Array4<ductile_fracture_model_type> modelfab_d;

                if (fracture_type == FractureType::Brittle)
                    modelfab_b = (material.brittlemodel)[ilev]->array(mfi);
                else if (fracture_type == FractureType::Ductile)
                    modelfab_d = (material.ductilemodel)[ilev]->array(mfi);
                
                amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k){
                    Set::Scalar _temp = 0;
                    Set::Scalar mul = AMREX_D_PICK(0.5,0.25,0.125);
                    if (fracture_type == FractureType::Brittle)
                    {
                        _temp =  mul*(AMREX_D_TERM(	
                                            crack.cracktype->g_phi(c_new(i,j,k,0),0.) + crack.cracktype->g_phi(c_new(i-1,j,k,0),0.)
                                            , 
                                            + crack.cracktype->g_phi(c_new(i,j-1,k,0),0.) + crack.cracktype->g_phi(c_new(i-1,j-1,k,0),0.)
                                            , 
                                            + crack.cracktype->g_phi(c_new(i,j,k-1,0),0.) + crack.cracktype->g_phi(c_new(i-1,j,k-1,0),0.)
                                            + crack.cracktype->g_phi(c_new(i,j-1,k-1,0),0.) + crack.cracktype->g_phi(c_new(i-1,j-1,k-1,0),0.))
                                            );
                    }
                    else if (fracture_type == FractureType::Ductile)
                    {
                        Set::Scalar p = modelfab_d(i,j,k).curr.alpha;
                        _temp =  mul*(AMREX_D_TERM(	
                                            crack.cracktype->g_phi(c_new(i,j,k,0),p) + crack.cracktype->g_phi(c_new(i-1,j,k,0),p)
                                            , 
                                            + crack.cracktype->g_phi(c_new(i,j-1,k,0),p) + crack.cracktype->g_phi(c_new(i-1,j-1,k,0),p)
                                            , 
                                            + crack.cracktype->g_phi(c_new(i,j,k-1,0),p) + crack.cracktype->g_phi(c_new(i-1,j,k-1,0),p)
                                            + crack.cracktype->g_phi(c_new(i,j-1,k-1,0),p) + crack.cracktype->g_phi(c_new(i-1,j-1,k-1,0),p))
                                            );
                    }
                    
                    if (std::isnan(_temp)) Util::Abort(INFO);
                    if(_temp < 0.0) _temp = 0.;
                    if(_temp > 1.0) _temp = 1.0;

                    if (fracture_type == FractureType::Brittle)
                        modelfab_b(i,j,k,0).DegradeModulus(std::min(1.-_temp,1.-crack.scaleModulusMax));
                    else if (fracture_type == FractureType::Ductile)
                    {
                        modelfab_d(i,j,k,0).DegradeModulus(std::min(1.-_temp,1.-crack.scaleModulusMax));
                        modelfab_d(i,j,k,0).DegradeYieldSurface(std::min(1.-_temp,1.-crack.scaleModulusMax));
                    }
                });
            }
            if (fracture_type == FractureType::Brittle) Util::RealFillBoundary(*material.brittlemodel[ilev],geom[ilev]);
            else Util::RealFillBoundary(*material.ductilemodel[ilev],geom[ilev]);
        }
    }
    //==================================================

    //==================================================
    // Setting the body forces
    {
        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            const Real* DX = geom[ilev].CellSize();
            Set::Scalar volume = AMREX_D_TERM(DX[0],*DX[1],*DX[2]);

            AMREX_D_TERM(elastic.rhs[ilev]->setVal(loading.body_force(0)*volume,0,1);,
                    elastic.rhs[ilev]->setVal(loading.body_force(1)*volume,1,1);,
                    elastic.rhs[ilev]->setVal(loading.body_force(2)*volume,2,1););
        }
    }
    //==================================================

    //==================================================
    // Setting the elastic boundary conditions
    {
        BC::Operator::Elastic::Constant::Type bctype_d = BC::Operator::Elastic::Constant::Type::Displacement;
        BC::Operator::Elastic::Constant::Type bctype_b = BC::Operator::Elastic::Constant::Type::Displacement;
        if (loading.load == LoadType::Force)
        {
            bctype_d = BC::Operator::Elastic::Constant::Type::Traction;
            bctype_b = BC::Operator::Elastic::Constant::Type::Traction;
        }
        
        if (loading.mode == ModeType::ModeI)
        {
            if (fracture_type == FractureType::Brittle)
            {
                elastic.brittlebc.Set(elastic.brittlebc.Face::YHI, elastic.brittlebc.Direction::Y, bctype_b, loading.val);
                elastic.brittlebc.Set(elastic.brittlebc.Face::XLO_YHI, elastic.brittlebc.Direction::Y, bctype_b, loading.val);
                elastic.brittlebc.Set(elastic.brittlebc.Face::XHI_YHI, elastic.brittlebc.Direction::Y, bctype_b, loading.val);
                elastic.brittlebc.Init(elastic.rhs,geom);
            }
            else
            {
                elastic.ductilebc.Set(elastic.ductilebc.Face::YHI, elastic.ductilebc.Direction::Y, bctype_d, loading.val);
                elastic.ductilebc.Set(elastic.ductilebc.Face::XLO_YHI, elastic.ductilebc.Direction::Y, bctype_d, loading.val);
                elastic.ductilebc.Set(elastic.ductilebc.Face::XHI_YHI, elastic.ductilebc.Direction::Y, bctype_d, loading.val);
                elastic.ductilebc.Init(elastic.rhs,geom);
            }
            
        }
        else if (loading.mode == ModeType::ModeII)
        {
            if (fracture_type == FractureType::Brittle)
            {
                elastic.brittlebc.Set(elastic.brittlebc.Face::YHI, elastic.brittlebc.Direction::X, bctype_b, loading.val);
                elastic.brittlebc.Set(elastic.brittlebc.Face::XLO_YHI, elastic.brittlebc.Direction::X, bctype_b, loading.val);
                elastic.brittlebc.Set(elastic.brittlebc.Face::XHI_YHI, elastic.brittlebc.Direction::X, bctype_b, loading.val);
                elastic.brittlebc.Init(elastic.rhs,geom);
            }
            else
            {
                elastic.ductilebc.Set(elastic.ductilebc.Face::YHI, elastic.ductilebc.Direction::X, bctype_d, loading.val);
                elastic.ductilebc.Set(elastic.ductilebc.Face::XLO_YHI, elastic.ductilebc.Direction::X, bctype_d, loading.val);
                elastic.ductilebc.Set(elastic.ductilebc.Face::XHI_YHI, elastic.ductilebc.Direction::X, bctype_d, loading.val);
                elastic.ductilebc.Init(elastic.rhs,geom);
            }
        }
        else
            Util::Abort(INFO, "Loading mode not implemented yet");
    }
    //==================================================

    //==================================================
    // Setting up the solver parameters
    Operator::Elastic<brittle_fracture_model_type::sym> op_b;
    Operator::Elastic<ductile_fracture_model_type::sym> op_d;
    {
        LPInfo info;
        info.setAgglomeration(sol.agglomeration);
        info.setConsolidation(sol.consolidation);
        info.setMaxCoarseningLevel(sol.max_coarsening_level);

        for (int ilev = 0; ilev < nlevels; ilev++) if (elastic.disp[ilev]->contains_nan()) Util::Warning(INFO);

        if (fracture_type == FractureType::Brittle)
        {
            op_b.define(geom, grids, dmap, info);
            op_b.setMaxOrder(sol.linop_maxorder);
            op_b.SetBC(&elastic.brittlebc);
            Solver::Nonlocal::Newton<brittle_fracture_model_type>  solver(op_b);
            solver.setMaxIter(sol.max_iter);
            solver.setMaxFmgIter(sol.max_fmg_iter);
            solver.setFixedIter(sol.max_fixed_iter);
            solver.setVerbose(sol.verbose);
            solver.setCGVerbose(sol.cgverbose);
            solver.setBottomMaxIter(sol.bottom_max_iter);
            solver.setBottomTolerance(sol.cg_tol_rel) ;
            solver.setBottomToleranceAbs(sol.cg_tol_abs) ;
            if (sol.bottom_solver == "cg") solver.setBottomSolver(MLMG::BottomSolver::cg);
            else if (sol.bottom_solver == "bicgstab") solver.setBottomSolver(MLMG::BottomSolver::bicgstab);
            solver.solve(elastic.disp, elastic.rhs, material.brittlemodel, sol.tol_rel, sol.tol_abs);
            solver.compResidual(elastic.residual,elastic.disp,elastic.rhs,material.brittlemodel);
        }
        else
        {
            op_d.define(geom, grids, dmap, info);
            op_d.setMaxOrder(sol.linop_maxorder);
            op_d.SetBC(&elastic.ductilebc);
            Solver::Nonlocal::Newton<ductile_fracture_model_type>  solver(op_d);
            solver.setMaxIter(sol.max_iter);
            solver.setMaxFmgIter(sol.max_fmg_iter);
            solver.setFixedIter(sol.max_fixed_iter);
            solver.setVerbose(sol.verbose);
            solver.setCGVerbose(sol.cgverbose);
            solver.setBottomMaxIter(sol.bottom_max_iter);
            solver.setBottomTolerance(sol.cg_tol_rel) ;
            solver.setBottomToleranceAbs(sol.cg_tol_abs) ;
            if (sol.bottom_solver == "cg") solver.setBottomSolver(MLMG::BottomSolver::cg);
            else if (sol.bottom_solver == "bicgstab") solver.setBottomSolver(MLMG::BottomSolver::bicgstab);
            solver.solve(elastic.disp, elastic.rhs, material.ductilemodel, sol.tol_rel, sol.tol_abs);
            solver.compResidual(elastic.residual,elastic.disp,elastic.rhs,material.ductilemodel);
        }
    }
    //==================================================

    //==================================================
    // Computing new stresses, strains and energies
    {
        for (int ilev = 0; ilev < nlevels; ilev++)
        {
            if (fracture_type == FractureType::Brittle)
            {
                op_b.Strain(ilev,*elastic.strain[ilev],*elastic.disp[ilev]);
                op_b.Stress(ilev,*elastic.stress[ilev],*elastic.disp[ilev]);
                op_b.Energy(ilev,*elastic.energy[ilev],*elastic.disp[ilev]);
            }
            else
            {
                op_d.Strain(ilev,*elastic.strain[ilev],*elastic.disp[ilev]);
                op_d.Stress(ilev,*elastic.stress[ilev],*elastic.disp[ilev]);
                op_d.Energy(ilev,*elastic.energy[ilev],*elastic.disp[ilev]);
            }
        }
        for (int ilev = 0; ilev < nlevels; ilev++)
        {
            elastic.strain[ilev]->FillBoundary();
            elastic.energy_pristine[ilev]->setVal(0.0);
            elastic.energy_pristine_old[ilev]->FillBoundary();
            
            for (amrex::MFIter mfi(*elastic.strain[ilev],true); mfi.isValid(); ++mfi)
            {
                const amrex::Box& box = mfi.validbox();
                amrex::Array4<Set::Scalar>	const& sig_box 		    = (*elastic.stress[ilev]).array(mfi);
                amrex::Array4<const Set::Scalar> const& strain_box 	= (*elastic.strain[ilev]).array(mfi);
                amrex::Array4<Set::Scalar> const& energy_box 		= (*elastic.energy_pristine[ilev]).array(mfi);
                //amrex::Array4<Set::Scalar> const& energy_box_old 	= (*elastic.energy_pristine_old[ilev]).array(mfi);

                amrex::Array4<ductile_fracture_model_type>	model_box;
                amrex::Array4<Set::Scalar>	strainp_box;

                if (fracture_type == FractureType::Ductile)
                {
                    model_box 	= material.ductilemodel[ilev]->array(mfi);
                    strainp_box 	= (*plastic.strain[ilev]).array(mfi);
                }

                amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k){
                    Set::Matrix eps = Numeric::FieldToMatrix(strain_box,i,j,k);

                    if( eps != eps ) Util::Abort(INFO,"Nans in eps, eps = ", eps);

                    if(fracture_type == FractureType::Ductile)
                    {
                        Set::Matrix sig = Numeric::FieldToMatrix(sig_box,i,j,k);
                        model_box(i,j,k,0).EvolvePlasticStrain(sig,eps,0);
                        sig = model_box(i,j,k,0).DW(eps);
                        Numeric::MatrixToField(strainp_box,i,j,k,model_box(i,j,k,0).curr.epsp);
                        material.ductilemodeltype.SetF0(model_box(i,j,k).curr.epsp);
                    }

                    energy_box(i,j,k,0) = (fracture_type == FractureType::Brittle) ? material.brittlemodeltype.W(eps) : material.ductilemodeltype.W(eps);
                    if (std::isnan(energy_box(i,j,k,0))) Util::Abort(INFO, "Nans detected in energy_box. material.brittlemodeltype.W(eps) = ", material.brittlemodeltype.W(eps));
                    //energy_box(i,j,k,0) = energy_box(i,j,k,0) > energy_box_old(i,j,k,0) ? energy_box(i,j,k,0) : energy_box_old(i,j,k,0);
                });
            }
            elastic.energy_pristine[ilev]->FillBoundary();
            if (fracture_type == FractureType::Ductile)
            {
                plastic.strain[ilev]->FillBoundary();
                Util::RealFillBoundary(*material.ductilemodel[ilev],geom[ilev]);
            }
        }
    }
    //==================================================
}

void 
Fracture::Advance (int lev, Set::Scalar time, Set::Scalar dt)
{
    std::swap(*crack.field_old[lev], *crack.field[lev]);

	static amrex::IntVect AMREX_D_DECL(	dx(AMREX_D_DECL(1,0,0)), dy(AMREX_D_DECL(0,1,0)), dz(AMREX_D_DECL(0,0,1)));
	const Set::Scalar* DX = geom[lev].CellSize();

    for ( amrex::MFIter mfi(*crack.field[lev],true); mfi.isValid(); ++mfi )
	{
		const amrex::Box& bx = mfi.validbox();
		amrex::Array4<const Set::Scalar> const& c_old = (*crack.field_old[lev]).array(mfi);
		amrex::Array4<Set::Scalar> const& df = (*crack.driving_force[lev]).array(mfi);
		amrex::Array4<Set::Scalar> const& c_new = (*crack.field[lev]).array(mfi);
		amrex::Array4<const Set::Scalar> const& energy_box = (*elastic.energy_pristine[lev]).array(mfi);
        amrex::Array4<ductile_fracture_model_type>	model_box;

        if (fracture_type == FractureType::Ductile)
            model_box 	= material.ductilemodel[lev]->array(mfi);

        amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k){

			Set::Scalar rhs = 0.0;
            Set::Scalar p = (fracture_type == FractureType::Brittle) ? 0.0 : model_box(i,j,k).curr.alpha;

			// Elastic component of the driving force
			Set::Scalar en_cell = Numeric::Interpolate::NodeToCellAverage(energy_box,i,j,k,0);
            if (std::isnan(en_cell)) Util::Abort(INFO, "Nans detected in en_cell. energy_box(i,j,k,0) = ", energy_box(i,j,k,0));
			df(i,j,k,0) = crack.cracktype->Dg_phi(c_old(i,j,k,0),p)*en_cell*elastic.df_mult;
			rhs += crack.cracktype->Dg_phi(c_old(i,j,k,0),p)*en_cell*elastic.df_mult;

			// Boundary terms
			Set::Vector Dc = Numeric::Gradient(c_old, i, j, k, 0, DX);
			Set::Scalar Theta = atan2(Dc(1),Dc(0));

			Set::Matrix DDc = Numeric::Hessian(c_old, i, j, k, 0, DX);
			Set::Scalar laplacian = DDc.trace();

			if (!anisotropy.on || time < anisotropy.tstart)
			{
				df(i,j,k,1) = crack.cracktype->Gc(Theta)*crack.cracktype->Dw_phi(c_old(i,j,k,0),p)/(4.0*crack.cracktype->Zeta(Theta))*crack.mult_df_Gc;
				df(i,j,k,2) = 2.0*crack.cracktype->Zeta(Theta)*crack.cracktype->Gc(Theta)*laplacian*crack.mult_df_lap;

				rhs += crack.cracktype->Gc(Theta)*crack.cracktype->Dw_phi(c_old(i,j,k,0),p)/(4.0*crack.cracktype->Zeta(Theta))*crack.mult_df_Gc;
				rhs -= 2.0*crack.cracktype->Zeta(Theta)*crack.cracktype->Gc(Theta)*laplacian*crack.mult_df_lap;
			}
			else
			{
#if AMREX_SPACEDIM == 1
				Util::Abort(INFO, "Anisotropy is not enabled in 1D.");
#elif AMREX_SPACEDIM == 2
			    Set::Scalar normgrad = Dc.lpNorm<2>();
				Set::Vector normal = Dc / normgrad;
				Set::Vector tangent(normal[1],-normal[0]);

				Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Full> DDDDc = Numeric::DoubleHessian<AMREX_SPACEDIM>(c_old, i, j, k, 0, DX);
				Set::Scalar sinTheta = sin(Theta);
				Set::Scalar cosTheta = cos(Theta);
				Set::Scalar sin2Theta = sin(2.0*Theta);
				Set::Scalar cos2Theta = cos(2.0*Theta);

				Set::Scalar zeta = crack.cracktype->Zeta(Theta);
				Set::Scalar ws = crack.cracktype->w_phi(c_old(i,j,k,0),p)/(4.0*zeta*normgrad*normgrad);
				
				if( std::isnan(ws)) Util::Abort(INFO, "nan at m=",i,",",j,",",k);
				if( std::isinf(ws)) ws = 1.0E6;

				Set::Scalar Boundary_term = 0.;
				Boundary_term += crack.cracktype->Gc(Theta)*crack.cracktype->Dw_phi(c_old(i,j,k,0),p)/(4.0*crack.cracktype->Zeta(Theta));
				Boundary_term -= 2.0*crack.cracktype->Gc(Theta)*crack.cracktype->Zeta(Theta)*laplacian;

				Boundary_term += crack.cracktype->DGc(Theta)
								* (zeta - ws)
								* (sin2Theta*(DDc(0,0)-DDc(1,1)) - 2.0*cos2Theta*DDc(0,1));
				Boundary_term += crack.cracktype->DDGc(Theta)
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
			//} // disable this line to remove normgrad 1E-4 check

            if (fracture_type == FractureType::Brittle)
			    df(i,j,k,3) = std::max(0.,rhs - crack.cracktype->DrivingForceThreshold(c_old(i,j,k,0)));

            else
            {
                df(i,j,k,3) = crack.cracktype->Dg_phi(c_old(i,j,k,0),p)*(material.ductilemodeltype.OriginalPlasticEnergy());
                rhs += crack.cracktype->Dg_phi(c_old(i,j,k,0),p)*(material.ductilemodeltype.OriginalPlasticEnergy());

                df(i,j,k,4) = rhs;
			    df(i,j,k,5) = max(0.,rhs-crack.cracktype->DrivingForceThreshold(c_old(i,j,k,0)));
            }
            
			if(std::isnan(rhs)) Util::Abort(INFO, "Dwphi = ", crack.cracktype->Dw_phi(c_old(i,j,k,0),p),". c_old(i,j,k,0) = ",c_old(i,j,k,0));
			c_new(i,j,k,0) = c_old(i,j,k,0) - dt*std::max(0., rhs - crack.cracktype->DrivingForceThreshold(c_old(i,j,k,0)))*crack.cracktype->Mobility(c_old(i,j,k,0));

			if(c_new(i,j,k,0) > 1.0) {Util::Message(INFO, "cnew = ", c_new(i,j,k,0) ,", resetting to 1.0"); c_new(i,j,k,0) = 1.;}
			if(c_new(i,j,k,0) < 0.0) {Util::Message(INFO, "cnew = ", c_new(i,j,k,0) ,", resetting to 0.0"); c_new(i,j,k,0) = 0.;}
		});

    }
}

void 
Fracture::Integrate(int amrlev, Set::Scalar /*time*/, int /*step*/,const amrex::MFIter &mfi, const amrex::Box &box)
{
	const amrex::Real* DX = geom[amrlev].CellSize();

	amrex::Array4<const Set::Scalar> const &c_new = (*crack.field[amrlev]).array(mfi);
	amrex::Array4<const Set::Scalar> const &c_old = (*crack.field_old[amrlev]).array(mfi);
    amrex::Array4<const Set::Scalar> ep_new;
    amrex::Array4<const Set::Scalar> ep_old;

    if(fracture_type == FractureType::Ductile)
    {
        ep_new = (*plastic.strain[amrlev]).array(mfi);
        ep_old = (*plastic.strain_old[amrlev]).array(mfi);
    }

    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
    {
		crack.error_norm += ((c_new(i,j,k,0)-c_old(i,j,k,0))*(c_new(i,j,k,0)-c_old(i,j,k,0)))*(AMREX_D_TERM(DX[0],*DX[1],*DX[2]));
		crack.norm += c_new(i,j,k,0)*c_new(i,j,k,0)*(AMREX_D_TERM(DX[0],*DX[1],*DX[2]));
        
        if(fracture_type == FractureType::Ductile)
        {
            for (int n = 0; n < AMREX_SPACEDIM*AMREX_SPACEDIM; n++)
            {
                plastic.norm += ep_new(i,j,k,n)*ep_new(i,j,k,n)*(AMREX_D_TERM(DX[0],*DX[1],*DX[2]));
                plastic.error_norm += (ep_new(i,j,k,n)-ep_old(i,j,k,n))*(ep_new(i,j,k,n)-ep_old(i,j,k,n))*(AMREX_D_TERM(DX[0],*DX[1],*DX[2]));
            }
            plastic.norm *= crack.cracktype->g_phi(c_new(i,j,k,0),0.0);
		    plastic.error_norm *= crack.cracktype->g_phi(c_new(i,j,k,0),0.);
        }
	});
}

void
Fracture::TagCellsForRefinement (int lev, amrex::TagBoxArray &a_tags, amrex::Real /*time*/, int /*ngrow*/)
{
	const amrex::Real *DX = geom[lev].CellSize();
	const Set::Vector dx(DX);
	const Set::Scalar dxnorm = dx.lpNorm<2>();

	for (amrex::MFIter mfi(*crack.field[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		const amrex::Box 							&bx 	= mfi.tilebox();
		amrex::Array4<char> const 					&tags 	= a_tags.array(mfi);
		amrex::Array4<const Set::Scalar> const 		&c_new 	= (*crack.field[lev]).array(mfi);
		
		amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
			Set::Vector grad = Numeric::Gradient(c_new, i, j, k, 0, DX);
			if (dxnorm * grad.lpNorm<2>() > crack.refinement_threshold)
				tags(i, j, k) = amrex::TagBox::SET;
		});
	}
}

void
Fracture::TimeStepComplete(amrex::Real time,int iter)
{
	IntegrateVariables(time,iter);

    Set::Scalar rel_c_err = crack.error_norm/crack.norm;
    Util::Message(INFO, "crack_err_norm = ", crack.error_norm);
	Util::Message(INFO, "c_new_norm = ", crack.norm);
    Util::Message(INFO, "crack relative error= ", rel_c_err);

    if (fracture_type == FractureType::Ductile)
    {
        Set::Scalar rel_ep_err = plastic.error_norm/plastic.norm;
        Util::Message(INFO, "plastic strain norm = ", plastic.norm);
        Util::Message(INFO, "plastic strain error norm = ", plastic.error_norm);
        Util::Message(INFO, "plastic strain relative error= ", rel_ep_err);
    }

    if (crack.error_norm > crack.tol_abs || rel_c_err > crack.tol_rel) return;
    if (fracture_type == FractureType::Ductile && (plastic.norm > plastic.tol_abs || plastic.error_norm > plastic.tol_rel)) return;

	WritePlotFile("crack",(double)loading.step,loading.step);

	loading.step++;
	if(loading.val >= loading.max) SetStopTime(time-0.01);
}

}
