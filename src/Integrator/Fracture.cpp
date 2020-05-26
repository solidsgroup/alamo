#include "Fracture.H"

namespace Integrator
{
Fracture::Fracture() :
	Integrator()
{
    //==================================================
    // Problem type
    {
        IO::ParmParse pp("fracture");
        pp.query("problem_type", fracture_problem);
        if (!(fracture_problem == "brittle" || fracture_problem == "ductile"))
            Util::Abort(INFO, "Incorrect fracture problem");
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
        pp_crack.queryclass("bc",*static_cast<BC::Constant *>(crack.mybc));
        pp_crack.queryclass("bc_df",*static_cast<BC::Constant *>(crack.mybcdf));

        RegisterNewFab(crack.field,     crack.bc, 1, number_of_ghost_cells, "c",		true);
        RegisterNewFab(crack.field_old, crack.bc, 1, number_of_ghost_cells, "c_old",	true);
        if (fracture_problem == "brittle")
            RegisterNewFab(crack.driving_force, crack.bcdf, 4, number_of_ghost_cells, "driving_force",true);
        else
            RegisterNewFab(crack.driving_force, crack.bcdf, 6, number_of_ghost_cells, "driving_force",true);
    }
    //==================================================

    //==================================================
    // Material model
    {
        IO::ParmParse pp_material("material");
        pp_material.query("model",material.input_material);

        if (fracture_problem == "brittle")
        {
            if(material.input_material == "isotropic") pp_material.queryclass("isotropic",material.brittlemodeltype);
            else Util::Abort(INFO,"This model has not been implemented yet.");
        }
        else
        {
            if(material.input_material == "isotropicj2plastic") pp_material.queryclass("isotropicj2plastic",material.ductilemodeltype);
            else if(material.input_material == "cubiccrystalplastic")pp_material.queryclass("cubiccrystalplastic",material.ductilemodeltype);
            else Util::Abort(INFO,"This model has not been implemented yet.");
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

        IO::ParmParse pp_load("loading");
        if (pp_load.countval("body_force")) pp_load.queryarr("body_force",loading.body_force);
        pp_load.query("mode",loading.mode);
        pp_load.query("type",loading.type);
        pp_load.query("init", loading.init);
        pp_load.query("max", loading.max);
        pp_load.query("rate", loading.rate);
        loading.val = loading.init;

        if(fracture_problem == "brittle") pp_elastic.queryclass("bc",elastic.brittlebc);
        else pp_elastic.queryclass("bc",elastic.ductilebc);
    }
    //==================================================

    //==================================================
    // Registering fabs now
    {
        nlevels = maxLevel() + 1;
        RegisterNodalFab(elastic.disp,  AMREX_SPACEDIM, number_of_ghost_nodes, "disp", True);
        RegisterNodalFab(elastic.rhs,  AMREX_SPACEDIM, number_of_ghost_nodes, "rhs", True);
        RegisterNodalFab(elastic.residual,  AMREX_SPACEDIM, number_of_ghost_nodes, "res", True);
        RegisterNodalFab(elastic.strain,  AMREX_SPACEDIM*AMREX_SPACEDIM, number_of_ghost_nodes, "strain", True);
        RegisterNodalFab(elastic.stress,  AMREX_SPACEDIM*AMREX_SPACEDIM, number_of_ghost_nodes, "stress", True);
        RegisterNodalFab(elastic.energy, 1, number_of_ghost_nodes, "energy", True);
        RegisterNodalFab(elastic.energy_pristine, 1, number_of_ghost_nodes, "energy_pristine", True);
        RegisterNodalFab(elastic.energy_pristine_old, 1, number_of_ghost_nodes, "energy_pristine_old", True);
        if (fracture_problem == "brittle") 
        {
            RegisterGeneralFab(material.brittlemodel, 1, number_of_ghost_nodes);
            material.brittlemodel.resize(nlevels);
        }
        else
        {
            RegisterNodalFab(plastic.strain, AMREX_SPACEDIM*AMREX_SPACEDIM, number_of_ghost_nodes, "strainp", True);
            RegisterGeneralFab(material.ductilemodel, 1, number_of_ghost_nodes);
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
        if(fracture_problem == "brittle")
            material.brittlemodel[ilev]->setVal(material.brittlemodeltype);
        else
        {
            plastic.strain[ilev] -> setVal(0.);
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
		crack.field[ilev]->FillBoundary();
    }

    //==================================================
    // Scaling the modulus
    {
        for (int ilev = 0; ilev < nlevels; ++ilev)
        {
            for (amrex::MFIter mfi(*(fracture_problem == "brittle" ? material.brittlemodel : material.ductilemodel)[ilev],true); mfi.isValid(); ++mfi)
            {
                amrex::Box box = mfi.growntilebox(1);
                amrex::Array4<const amrex::Real> const& c_new = (*crack.field[ilev]).array(mfi);

                if (fracture_problem == "brittle")
                    amrex::Array4<brittle_fracture_model_type> const& modelfab = (material.brittlemodel)[ilev]->array(mfi);
                else
                    amrex::Array4<ductile_fracture_model_type> const& modelfab = (material.ductilemodel)[ilev]->array(mfi);

                amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k){
                    Set::Scalar _temp = 0;
                    Set::Scalar mul = AMREX_D_PICK(0.5,0.25,0.125);
                    if (fracture_problem == "brittle")
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
                    else
                    {
                        Set::Scalar p = modelfab(i,j,k).curr.alpha;
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

                    modelfab(i,j,k,0).DegradeModulus(std::min(1.-_temp,1.-crack.scaleModulusMax));
                    if (fracture_problem == "ductile") modelfab(i,j,k,0).DegradeYieldSurface(std::min(1.-_temp,1.-crack.scaleModulusMax));
                });
            }
            if (fracture_problem == "brittle") Util::RealFillBoundary(*material.brittlemodel[ilev],geom[ilev]);
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
        BC::Operator::Elastic<brittle_fracture_model_type>::Type brittletype = BC::Operator::Elastic<brittle_fracture_model_type>::Type::Displacement;
        BC::Operator::Elastic<ductile_fracture_model_type>::Type ductiletype = BC::Operator::Elastic<ductile_fracture_model_type>::Type::Displacement;
        if (loading.type == "force")
        {
            brittletype = BC::Operator::Elastic<brittle_fracture_model_type>::Type::Traction;
            ductiletype = BC::Operator::Elastic<ductile_fracture_model_type>::Type::Traction;
        }
        if (loading.mode == "mode1")
        {
            if (fracture_problem == "brittle")
            {
                elastic.brittlebc.Set(elastic.brittlebc.Face::YHI, elastic.brittlebc.Direction::Y, brittletype, loading.val);
                elastic.brittlebc.Set(elastic.brittlebc.Face::XLO_YHI, elastic.brittlebc.Direction::Y, brittletype, loading.val);
                elastic.brittlebc.Set(elastic.brittlebc.Face::XHI_YHI, elastic.brittlebc.Direction::Y, brittletype, loading.val);
                elastic.brittlebc.Init(elastic.rhs,geom);
            }
            else
            {
                elastic.ductilebc.Set(elastic.ductilebc.Face::YHI, elastic.ductilebc.Direction::Y, ductiletype, loading.val);
                elastic.ductilebc.Set(elastic.ductilebc.Face::XLO_YHI, elastic.ductilebc.Direction::Y, ductiletype, loading.val);
                elastic.ductilebc.Set(elastic.ductilebc.Face::XHI_YHI, elastic.ductilebc.Direction::Y, ductiletype, loading.val);
                elastic.ductilebc.Init(elastic.rhs,geom);
            }
            
        }
        else if (loading.mode == "mode2")
        {
            if (fracture_problem == "brittle")
            {
                elastic.brittlebc.Set(elastic.brittlebc.Face::YHI, elastic.brittlebc.Direction::X, brittletype, loading.val);
                elastic.brittlebc.Set(elastic.brittlebc.Face::XLO_YHI, elastic.brittlebc.Direction::X, brittletype, loading.val);
                elastic.brittlebc.Set(elastic.brittlebc.Face::XHI_YHI, elastic.brittlebc.Direction::X, brittletype, loading.val);
                elastic.brittlebc.Init(elastic.rhs,geom);
            }
            else
            {
                elastic.ductilebc.Set(elastic.ductilebc.Face::YHI, elastic.ductilebc.Direction::X, ductiletype, loading.val);
                elastic.ductilebc.Set(elastic.ductilebc.Face::XLO_YHI, elastic.ductilebc.Direction::X, ductiletype, loading.val);
                elastic.ductilebc.Set(elastic.ductilebc.Face::XHI_YHI, elastic.ductilebc.Direction::X, ductiletype, loading.val);
                elastic.ductilebc.Init(elastic.rhs,geom);
            }
        }
        else
            Util::Abort(INFO, "Loading mode not implemented yet");
    }
    //==================================================

    //==================================================
    // Setting up the solver parameters
    Operator::Elastic<(fracture_model == "brittle" ? brittle_fracture_model_type ? ductile_fracture_model_type)> elastic_op;
    {
        LPInfo info;
        info.setAgglomeration(sol.agglomeration);
        info.setConsolidation(sol.consolidation);
        info.setMaxCoarseningLevel(sol.max_coarsening_level);

        elastic_op.define(geom, grids, dmap, info);
        elastic_op.setMaxOrder(sol.linop_maxorder);

        elastic_op.SetBC(&(fracture_model == "brittle" ? elastic.brittlebc : elastic.ductilebc));

        Solver::Nonlocal::Newton<(fracture_model == "brittle" ? brittle_fracture_model_type ? ductile_fracture_model_type)>  solver(elastic_op);
        solver.setMaxIter(sol.max_iter);
        solver.setMaxFmgIter(sol.max_fmg_iter);
        solver.setFixedIter(sol.max_fixed_iter);
        solver.setVerbose(sol.verbose);
        solver.setCGVerbose(sol.cgverbose);
        solver.setBottomMaxIter(sol.bottom_max_iter);
        solver.setBottomTolerance(sol.cg_tol_rel) ;
        solver.setBottomToleranceAbs(sol.cg_tol_abs) ;

        for (int ilev = 0; ilev < nlevels; ilev++) if (elastic.disp[ilev]->contains_nan()) Util::Warning(INFO);

        if (sol.bottom_solver == "cg") solver.setBottomSolver(MLMG::BottomSolver::cg);
        else if (sol.bottom_solver == "bicgstab") solver.setBottomSolver(MLMG::BottomSolver::bicgstab);

        solver.solve(elastic.disp, elastic.rhs, fracture_problem == "brittle" ? material.brittlemodel : material.ductilemodel, sol.tol_rel, sol.tol_abs);
        solver.compResidual(elastic.residual,elastic.disp,elastic.rhs,fracture_problem == "brittle" ? material.brittlemodel : material.ductilemodel);
    }
    //==================================================

    //==================================================
    // Computing new stresses, strains and energies
    {
        for (int ilev = 0; ilev < nlevels; ilev++)
        {
            elastic_op.Strain(ilev,*elastic.strain[ilev],*elastic.disp[ilev]);
            elastic_op.Stress(ilev,*elastic.stress[ilev],*elastic.disp[ilev]);
            elastic_op.Energy(ilev,*elastic.energy[ilev],*elastic.disp[ilev]);
        }
        for (int ilev = 0; ilev < nlevels; ilev++)
        {
            elastic.strain[ilev]->FillBoundary();
            elastic.energy_pristine_old[ilev]->FillBoundary();
            
            for (amrex::MFIter mfi(*elastic.strain[ilev],true); mfi.isValid(); ++mfi)
            {
                const amrex::Box& box = mfi.growntilebox(number_of_ghost_nodes);//validbox();
                amrex::Array4<Set::Scalar>	const& sig_box 		    = (*elastic.stress[ilev]).array(mfi);
                amrex::Array4<const Set::Scalar> const& strain_box 	= (*elastic.strain[ilev]).array(mfi);
                amrex::Array4<Set::Scalar> const& energy_box 		= (*elastic.energy_pristine[ilev]).array(mfi);
                amrex::Array4<Set::Scalar> const& energy_box_old 	= (*elastic.energy_pristine_old[ilev]).array(mfi);

                if (fracture_problem == "ductile")
                {
                    amrex::Array4<ductile_fracture_model_type>	const& model_box 	= material.ductilemodel[ilev]->array(mfi);
                    amrex::Array4<Set::Scalar>	const& strainp_box 	= (*plastic.strain[ilev]).array(mfi);
                }

                amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k){
                    Set::Matrix eps = Numeric::FieldToMatrix(strain_box,i,j,k);

                    if(fracture_problem == "ductile")
                    {
                        Set::Matrix sig = Numeric::FieldToMatrix(sig_box,i,j,k);
                        model_box(i,j,k,0).EvolvePlasticStrain(sig,eps,0);
                        sig = model_box(i,j,k,0).DW(eps);
                    }

                    energy_box(i,j,k,0) = (fracture_problem == "brittle") ? material.brittlemodeltype.W(eps) : material.ductilemodeltype.W(eps);
                    energy_box(i,j,k,0) = energy_box(i,j,k,0) > energy_box_old(i,j,k,0) ? energy_box(i,j,k,0) : energy_box_old(i,j,k,0);
                });
            }
            m_energy_pristine[lev]->FillBoundary();
        }
    }
    //==================================================
}

}