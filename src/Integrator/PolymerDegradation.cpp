#include "PolymerDegradation.H"
#include "Solver/Nonlocal/Newton.H"

//#if AMREX_SPACEDIM == 1
namespace Integrator
{
PolymerDegradation::PolymerDegradation():
    Integrator()
{
    Util::Message(INFO);
    //
    // READ INPUT PARAMETERS
    //

    // ---------------------------------------------------------------------
    // --------------------- Water diffusion -------------------------------
    // ---------------------------------------------------------------------
    IO::ParmParse pp_water("water");
    pp_water.query("on",water.on);
    if(water.on)
    {
        pp_water.query("diffusivity", water.diffusivity);
        pp_water.query("refinement_threshold", water.refinement_threshold);
        pp_water.query("ic_type", water.ic_type);

        // // Determine initial condition
        if (water.ic_type == "constant")
        {
            IO::ParmParse pp_water_ic("water.ic");
            std::vector<amrex::Real> value;
            pp_water_ic.queryarr("value",value);
            water.ic = new IC::Constant(geom,value);
        }
        else
            Util::Abort(INFO, "This kind of IC has not been implemented yet");

        water.bc = new BC::Constant(1);
        pp_water.queryclass("bc",*static_cast<BC::Constant *>(water.bc));

        RegisterNewFab(water_conc,     water.bc, 1, number_of_ghost_cells, "Water Concentration",true);
        RegisterNewFab(water_conc_old, water.bc, 1, number_of_ghost_cells, "Water Concentration Old",false);
    }

    Util::Message(INFO);
    // ---------------------------------------------------------------------
    // --------------------- Heat diffusion -------------------------------
    // ---------------------------------------------------------------------
    IO::ParmParse pp_heat("thermal");
    pp_heat.query("on",thermal.on);
    if(thermal.on)
    {
        pp_heat.query("diffusivity", thermal.diffusivity);
        pp_heat.query("refinement_threshold",thermal.refinement_threshold);
        pp_heat.query("ic_type",thermal.ic_type);

        if (thermal.ic_type == "constant")
        {
            amrex::ParmParse pp_heat_ic("thermal.ic");
            amrex::Vector<amrex::Real> T;
            pp_heat_ic.queryarr("value",T);
            thermal.ic = new IC::Constant(geom,T);
        }
        else
            Util::Abort(INFO, "This kind of IC has not been implemented yet");

        thermal.bc = new BC::Constant(1);
        pp_heat.queryclass("bc",*static_cast<BC::Constant *>(thermal.bc));
    }
    else
        thermal.bc = new BC::Nothing();

    RegisterNewFab(Temp,     thermal.bc, 1, number_of_ghost_cells, "Temperature",true);
    RegisterNewFab(Temp_old, thermal.bc, 1, number_of_ghost_cells, "Temperature Old",false);

    // ---------------------------------------------------------------------
    // --------------------- Material model --------------------------------
    // ---------------------------------------------------------------------
    Util::Message(INFO);
    std::string input_material;
    IO::ParmParse pp_material("material");
    pp_material.query("model",input_material);
    Util::Message(INFO);
    if(input_material == "isotropic2")
    {
        pp_material.queryclass("isotropic2",material.modeltype);
        pp_material.query("yield_strength",material.yieldstrength);
    }
    
    else if(input_material == "isotropic")
    { Util::Message(INFO);
        pp_material.queryclass("isotropic",material.modeltype);
        pp_material.query("yield_strength",material.yieldstrength);
    }
    
    else
        Util::Abort(INFO, "Not implemented yet");

    // ---------------------------------------------------------------------
    // --------------------- Damage model ----------------------------------
    // ---------------------------------------------------------------------
    Util::Message(INFO);
    IO::ParmParse pp_damage("damage"); // Phase-field model parameters
    pp_damage.query("anisotropy",damage.anisotropy);

    if(damage.anisotropy == 0)
        damage.number_of_eta = 1;
    else
        damage.number_of_eta = AMREX_SPACEDIM;

    pp_damage.query("type",damage.type);

    if(damage.type == "water" || damage.type == "water2") 
    {
        if(damage.type == "water") damage.number_of_eta = 2;
        else damage.number_of_eta = 4;

        damage.anisotropy = 0;

        damage.d_final.resize(damage.number_of_eta);
        damage.d_i.resize(damage.number_of_eta);
        damage.tau_i.resize(damage.number_of_eta);
        damage.t_start_i.resize(damage.number_of_eta);
        damage.number_of_terms.resize(damage.number_of_eta);

        amrex::Vector<Set::Scalar> dfinal;
        pp_damage.queryarr("d_final",dfinal);

        if(dfinal.size()!=damage.number_of_eta) Util::Abort(INFO, "Incorrect size of d_final");
        for (int i = 0; i < damage.number_of_eta; i++)
        {
            if(dfinal[i] <0 || dfinal[i] > 1.0) Util::Abort(INFO,"Incorrect value of d_final");
            damage.d_final[i] = dfinal[i];
        }

        int totalnumberofterms = 0, tempIndex = 0, tempIndex2 = 0;
        amrex::Vector<int> numberofterms;
        pp_damage.queryarr("number_of_terms",numberofterms);
        if(numberofterms.size()!=damage.number_of_eta) Util::Abort(INFO, "Incorrect size of number_of_terms");
        for (int i = 0; i<damage.number_of_eta; i++)
        {
            if(numberofterms[i] < 1) {Util::Abort(INFO, "number_of_terms can not be less than 1. Resetting"); numberofterms[i]=1;}
            damage.number_of_terms[i] = numberofterms[i];
            totalnumberofterms += numberofterms[i];
        }

        amrex::Vector<Set::Scalar> di, taui, tstarti;
        Set::Scalar sum = 0.;
        pp_damage.queryarr("d_i",di);
        pp_damage.queryarr("tau_i",taui);
        pp_damage.queryarr("t_start_i",tstarti);
        if(di.size()!=totalnumberofterms || taui.size()!=totalnumberofterms || tstarti.size()!=totalnumberofterms) Util::Abort(INFO,"Incorrect number of terms in di, taui or tstarti");
        for (int i=0; i<totalnumberofterms; i++)
        {
            if(tempIndex2 == damage.number_of_terms[tempIndex]) 
            {
                if(sum != 1.0) Util::Abort(INFO, "d_i's don't add to 1");
                tempIndex2 = 0; tempIndex+=1; sum = 0.;
            }
            
            if(di[i] < 0. || di[i] > 1.) Util::Abort(INFO, "Invalid values of d_i. Must be between 0 and 1");
            if(taui[i] < 0.) Util::Abort(INFO,"Invalid values of tau");

            sum += di[i];
            damage.d_i[tempIndex].push_back(di[i]);
            damage.tau_i[tempIndex].push_back(taui[i]);
            damage.t_start_i[tempIndex].push_back(tstarti[i]);

            tempIndex2 += 1;
        }
    }
    else
        Util::Abort(INFO, "This kind of damage model has not been implemented yet");

    pp_damage.query("ic_type",damage.ic_type);
    pp_damage.query("refinement_threshold",damage.refinement_threshold);
    if(damage.ic_type == "constant")
    {
        amrex::ParmParse pp_damage_ic("damage.ic");
        amrex::Vector<amrex::Real> eta_init;
        pp_damage_ic.queryarr("value",eta_init);
        damage.ic = new IC::Constant(geom,eta_init);
    }
    else
        Util::Abort(INFO, "This kind of IC has not been implemented yet");

    damage.bc = new BC::Constant(damage.number_of_eta);
    pp_damage.queryclass("bc",*static_cast<BC::Constant *>(damage.bc));

    damage.bc_time = new BC::Constant(1);
    pp_damage.queryclass("bc_time",*static_cast<BC::Constant *>(damage.bc_time));

    RegisterNewFab(eta_new, damage.bc, damage.number_of_eta, number_of_ghost_cells, "Eta",true);
    RegisterNewFab(eta_old, damage.bc, damage.number_of_eta, number_of_ghost_cells, "Eta old",true);
    RegisterNewFab(damage_start_time,damage.bc_time,1,number_of_ghost_cells,"Start time",true);

    Util::Message(INFO);
    // ---------------------------------------------------------------------
    // --------------------- Elasticity parameters -------------------------
    // ---------------------------------------------------------------------
    IO::ParmParse pp_elastic("elastic");
    pp_elastic.query("on",elastic.on);
    if(elastic.on)
    {
        pp_elastic.query("int",                elastic.interval);
        pp_elastic.query("type",            elastic.type);
        pp_elastic.query("max_iter",        elastic.max_iter);
        pp_elastic.query("max_fmg_iter",    elastic.max_fmg_iter);
        pp_elastic.query("verbose",            elastic.verbose);
        pp_elastic.query("cgverbose",        elastic.cgverbose);
        pp_elastic.query("tol_rel",            elastic.tol_rel);
        pp_elastic.query("tol_abs",            elastic.tol_abs);
        pp_elastic.query("cg_tol_rel",        elastic.cg_tol_rel);
        pp_elastic.query("cg_tol_abs",        elastic.cg_tol_abs);
        pp_elastic.query("use_fsmooth",        elastic.use_fsmooth);
        pp_elastic.query("agglomeration",     elastic.agglomeration);
        pp_elastic.query("consolidation",     elastic.consolidation);

        pp_elastic.query("bottom_solver",elastic.bottom_solver);
        pp_elastic.query("linop_maxorder", elastic.linop_maxorder);
        pp_elastic.query("max_coarsening_level",elastic.max_coarsening_level);
        pp_elastic.query("verbose",elastic.verbose);
        pp_elastic.query("cg_verbose", elastic.cgverbose);
        pp_elastic.query("bottom_max_iter", elastic.bottom_max_iter);
        pp_elastic.query("max_fixed_iter", elastic.max_fixed_iter);
        pp_elastic.query("bottom_tol", elastic.bottom_tol);

        if (pp_elastic.countval("body_force")) pp_elastic.getarr("body_force",elastic.body_force);

        IO::ParmParse pp_elastic_bc("elastic.bc");
        amrex::Vector<std::string> AMREX_D_DECL(bc_x_lo_str,bc_y_lo_str,bc_z_lo_str);
        amrex::Vector<std::string> AMREX_D_DECL(bc_x_hi_str,bc_y_hi_str,bc_z_hi_str);

        
        elastic.bc_map["displacement"]     = BC::Operator::Elastic::Constant::Type::Displacement;
        elastic.bc_map["disp"]             = BC::Operator::Elastic::Constant::Type::Displacement;
        elastic.bc_map["traction"]         = BC::Operator::Elastic::Constant::Type::Traction;
        elastic.bc_map["trac"]             = BC::Operator::Elastic::Constant::Type::Traction;
        elastic.bc_map["neumann"]         = BC::Operator::Elastic::Constant::Type::Neumann;
        elastic.bc_map["periodic"]         = BC::Operator::Elastic::Constant::Type::Periodic;

        
        amrex::ParmParse pp_temp;
        Set::Scalar stop_time, timestep;
        pp_temp.query("timestep",timestep);
        pp_temp.query("stop_time",stop_time);

        if(elastic.type == "tensile_test" || elastic.type == "tensile")
        {
            pp_elastic.queryclass("bc",elastic.bc);
            pp_elastic.query("rate",elastic.test_rate);
            if(elastic.test_rate < 0.) { Util::Warning(INFO,"Rate can't be less than zero. Resetting to 1.0"); elastic.test_rate = 1.0; }
        
            pp_elastic.queryarr("test_time",elastic.test_time);
            if(elastic.test_time.size() < 1){Util::Warning(INFO,"No test time provided. Providing default value"); elastic.test_time={0.,stop_time};}
            std::sort(elastic.test_time.begin(),elastic.test_time.end());
            pp_elastic.query("test_duration",elastic.test_duration);
            if(elastic.test_duration < 0.) { Util::Warning(INFO,"Test duration must be positive. Resetting it to 2.0"); elastic.test_duration = 1.0; }
            pp_elastic.query("test_dt",elastic.test_dt);
            if(elastic.test_dt < 0.) { Util::Warning(INFO,"Test dt must be positive. Resetting it to 0.01"); elastic.test_duration = 0.01; }
        }

        else if(elastic.type == "tensile_single" || elastic.type == "single")
        {
            pp_elastic.queryclass("bc",elastic.bc);
            pp_elastic.query("tstart",elastic.tstart);
            if(elastic.tstart < 0.0)
            {
                Util::Warning(INFO,"Invalid value for elasitc t_start (",elastic.tstart,"). Setting it to zero");
                elastic.tstart = 0.0;
            }
            else if(elastic.tstart > stop_time)
            {
                Util::Warning(INFO,"Invalid value for elastic t_start (",elastic.tstart,"). Setting it to stop_time");
                elastic.tstart = stop_time;
            }

            pp_elastic.query("tend",elastic.tend);
            if(elastic.tend < elastic.tstart || elastic.tend > stop_time)
            {            
                Util::Warning(INFO,"Invalid value for elastic t_end (",elastic.tend,"). Setting it to stop_time");
                elastic.tend = stop_time;
                if(elastic.tstart == stop_time) elastic.tstart = stop_time - timestep;
            }
        }
        else
            Util::Abort(INFO, "Not implemented this type of test yet");
    
        //----------------------------------------------------------------------
        // The following routine should be replaced by RegisterNewFab in the
        // future. For now, we are manually defining and resizing
        //-----------------------------------------------------------------------

        const int number_of_stress_components = AMREX_SPACEDIM*AMREX_SPACEDIM;
        RegisterNodalFab (displacement,    AMREX_SPACEDIM,                    2,    "displacement",true);;
        RegisterNodalFab (rhs,            AMREX_SPACEDIM,                    2,    "rhs",true);;
        RegisterNodalFab (strain,        number_of_stress_components,    2,    "strain",true);;
        RegisterNodalFab (stress,        number_of_stress_components,    2,    "stress",true);;
        RegisterNodalFab (stress_vm,    1,                                2,    "stress_vm",true);;
        RegisterNodalFab (energy,        1,                                2,    "energy",true);;
        RegisterNodalFab (residual,        AMREX_SPACEDIM,                    2,    "residual",true);;

    }
    RegisterGeneralFab(material.model, 1, 2);
    nlevels = maxLevel() + 1;
}


void
PolymerDegradation::Advance (int lev, amrex::Real time, amrex::Real dt)
{
    Util::Message(INFO, "Enter");
    std::swap(*eta_old[lev],     *eta_new[lev]);
    //    if (elastic.on) if (rhs[lev]->contains_nan()) Util::Abort(INFO);

    if(water.on) std::swap(*water_conc_old[lev],*water_conc[lev]);
    if(thermal.on) std::swap(*Temp_old[lev], *Temp[lev]);

    static amrex::IntVect AMREX_D_DECL(    dx(AMREX_D_DECL(1,0,0)),
                        dy(AMREX_D_DECL(0,1,0)),
                        dz(AMREX_D_DECL(0,0,1)));
    const amrex::Real* DX = geom[lev].CellSize();

    if(water.on)
    {
        Util::Message(INFO);
        for ( amrex::MFIter mfi(*water_conc[lev],true); mfi.isValid(); ++mfi )
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<amrex::Real> const& water_old_box = (*water_conc_old[lev]).array(mfi);
            amrex::Array4<amrex::Real> const& water_box = (*water_conc[lev]).array(mfi);
            amrex::Array4<amrex::Real> const& time_box = (*damage_start_time[lev]).array(mfi);

            amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k){
                if(std::isnan(water_old_box(i,j,k,0))) Util::Abort(INFO, "Nan found in WATER_OLD(i,j,k)");
                if(std::isinf(water_old_box(i,j,k,0))) Util::Abort(INFO, "Nan found in WATER_OLD(i,j,k)");
                
                if(water_old_box(i,j,k,0) > 1.0)
                {
                    Util::Warning(INFO,"Water concentration exceeded 1 at (", i, ",", j, ",", "k) and lev = ", lev, " Resetting");
                    water_old_box(i,j,k,0) = 1.0;
                }
                
                water_box(i,j,k,0) = water_old_box(i,j,k,0) + dt * water.diffusivity * Numeric::Hessian(water_old_box,i,j,k,0,DX).trace();
                
                if(water_box(i,j,k,0) > 1.0)
                {
                    Util::Warning(INFO, "Water concentration has exceeded one after computation. Resetting it to one");
                    water_box(i,j,k,0) = 1.0;
                }
                if(water_old_box(i,j,k,0) < 1.E-2 && water_box(i,j,k,0) > 1.E-2)
                    time_box(i,j,k,0) = time;
            });
        }
        Util::Message(INFO);
        //water_conc[lev]->FillBoundary();
    }

    if(thermal.on)
    {
        for ( amrex::MFIter mfi(*Temp[lev],true); mfi.isValid(); ++mfi )
        {
            const amrex::Box& bx = mfi.validbox();
            amrex::Array4<const amrex::Real> const& Temp_old_box = (*Temp_old[lev]).array(mfi);
            amrex::Array4<amrex::Real> const& Temp_box = (*Temp[lev]).array(mfi);
            
            amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k){
                Temp_box(i,j,k,0) = Temp_old_box(i,j,k,0) + dt * thermal.diffusivity * Numeric::Hessian(Temp_old_box,i,j,k,0,DX).trace();
            });
        }
    }
    Util::Message(INFO);
    for ( amrex::MFIter mfi(*eta_new[lev],true); mfi.isValid(); ++mfi )
    {
        const amrex::Box& bx = mfi.growntilebox(1);
        amrex::Array4<amrex::Real> const& eta_new_box             = (*eta_new[lev]).array(mfi);
        amrex::Array4<const amrex::Real> const& eta_old_box     = (*eta_old[lev]).array(mfi);
        amrex::Array4<const amrex::Real> const& water_box         = (*water_conc[lev]).array(mfi);
        amrex::Array4<const amrex::Real> const& time_box         = (*damage_start_time[lev]).array(mfi);
        
        amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k){
            for (int n = 0; n < damage.number_of_eta; n++)
            {
                if(damage.type == "water" || damage.type == "water2")
                {
                    //Util::Message(INFO);
                    Set::Scalar temp1 = 0.0;
                    if(water_box(i,j,k,0) > 0.0 && eta_old_box(i,j,k,n) < damage.d_final[n])
                    {
                        for (int l = 0; l < damage.number_of_terms[n]; l++)
                                temp1 += damage.d_final[n]*damage.d_i[n][l]*water_box(i,j,k,0)*std::exp(-std::max(0.0,time-time_box(i,j,k,0)-damage.t_start_i[n][l])/damage.tau_i[n][l])/(damage.tau_i[n][l]);
                    }
                    eta_new_box(i,j,k,n) = eta_old_box(i,j,k,n) + temp1*dt;
                    if(eta_new_box(i,j,k,n) > damage.d_final[n])
                    {
                        Util::Warning(INFO, "eta exceeded ",damage.d_final[n], ". Rhs = ", temp1, ", Water = ", water_box(i,j,k,n));
                        eta_new_box(i,j,k,n) = damage.d_final[n];
                    }
                }
                else
                    Util::Abort(INFO, "Damage model not implemented yet");
            }
        });
    }
    //if(elastic.on)    if (rhs[lev]->contains_nan()) Util::Abort(INFO);
    Util::Message(INFO,"Exit");
}

void
PolymerDegradation::Initialize (int lev)
{
    Util::Message(INFO);
    if(water.on)
    {
        water.ic->Initialize(lev,water_conc);
        water.ic->Initialize(lev,water_conc_old);
        IC::IC *time_ic;
        time_ic = new IC::Constant(geom,{0.0});
        time_ic->Initialize(lev,damage_start_time);
    }

    if(thermal.on)
    {
        thermal.ic->Initialize(lev,Temp);
        thermal.ic->Initialize(lev,Temp_old);
    }

    Util::Message(INFO);
    damage.ic->Initialize(lev,eta_new);
    damage.ic->Initialize(lev,eta_old);
    
    if (elastic.on)
    {
        displacement[lev]->setVal(0.0);
        strain[lev]->setVal(0.0);
        stress[lev]->setVal(0.0);
        stress_vm[lev]->setVal(0.0);
        rhs[lev]->setVal(0.0);
        energy[lev]->setVal(0.0);
        residual[lev]->setVal(0.0);
    }
    material.model[lev]->setVal(material.modeltype);
    Util::Message(INFO);
}

PolymerDegradation::~PolymerDegradation()
{
    delete water.ic;
    delete damage.ic;
    delete water.bc;
    delete damage.bc;
}

void
PolymerDegradation::TagCellsForRefinement (int lev, amrex::TagBoxArray& a_tags, amrex::Real /*time*/, int /*ngrow*/)
{
    const amrex::Real *DX = geom[lev].CellSize();
    const Set::Vector dx(DX);
    const Set::Scalar dxnorm = dx.lpNorm<2>();

    if(water.on)
    {
        for (amrex::MFIter mfi(*water_conc[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box                             &bx         = mfi.tilebox();
            amrex::Array4<char> const                     &tags         = a_tags.array(mfi);
            amrex::Array4<const Set::Scalar> const         &water_box     = (*water_conc[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Vector grad = Numeric::Gradient(water_box, i, j, k, 0, DX);
                if (dxnorm * grad.lpNorm<2>() > water.refinement_threshold)
                    tags(i, j, k) = amrex::TagBox::SET;
            });
        }
    }

    if(thermal.on)
    {
        for (amrex::MFIter mfi(*Temp[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box                             &bx         = mfi.tilebox();
            amrex::Array4<char> const                     &tags         = a_tags.array(mfi);
            amrex::Array4<const Set::Scalar> const         &Temp_box     = (*Temp[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                Set::Vector grad = Numeric::Gradient(Temp_box, i, j, k, 0, DX);
                if (dxnorm * grad.lpNorm<2>() > thermal.refinement_threshold)
                    tags(i, j, k) = amrex::TagBox::SET;
            });
        }
    }

    for (amrex::MFIter mfi(*eta_new[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box                             &bx         = mfi.tilebox();
        amrex::Array4<char> const                     &tags         = a_tags.array(mfi);
        amrex::Array4<const Set::Scalar> const         &eta_box     = (*eta_new[lev]).array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Set::Vector grad = Numeric::Gradient(eta_box, i, j, k, 0, DX);
            if (dxnorm * grad.lpNorm<2>() > damage.refinement_threshold)
                tags(i, j, k) = amrex::TagBox::SET;
        });
    }
    Util::Message(INFO);
}

void
PolymerDegradation::DegradeMaterial(int lev, amrex::FabArray<amrex::BaseFab<pd_model_type> > &model)
{
    Util::Message(INFO);
    /*
    This function is supposed to degrade material parameters based on certain
    damage model.
    For now we are just using isotropic degradation.
    */
    if(damage.anisotropy) Util::Abort(__FILE__,"DegradeModulus",__LINE__,"Not implemented yet");

    static amrex::IntVect AMREX_D_DECL(dx(AMREX_D_DECL(1,0,0)),
                        dy(AMREX_D_DECL(0,1,0)),
                        dz(AMREX_D_DECL(0,0,1)));
    eta_new[lev]->FillBoundary();

    for (amrex::MFIter mfi(model,true); mfi.isValid(); ++mfi)
    {
        amrex::Box box = mfi.growntilebox(2);
        //box.grow(1);
        amrex::Array4<const amrex::Real> const& eta_box = (*eta_new[lev]).array(mfi);
        amrex::Array4<pd_model_type> const& modelfab = model.array(mfi);

        amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k){
            amrex::Vector<Set::Scalar> _temp;
            for(int n=0; n<damage.number_of_eta; n++)
            {
                Set::Scalar _temp2 = Numeric::Interpolate::CellToNodeAverage(eta_box,i,j,k,n);
                _temp.push_back(std::min(damage.d_final[n],std::max(0.,_temp2)));
            }
            if(damage.type == "water") modelfab(i,j,k,0).DegradeModulus(_temp[0]);
            else if (damage.type == "water2") modelfab(i,j,k,0).DegradeModulus(_temp);
            else Util::Abort(INFO, "Damage model not implemented yet");
        });
    }

    amrex::Geometry tmp_geom = geom[lev];
    for (int i = 0; i < 2; i++)
    {
        Util::Message(INFO);
        amrex::FabArray<amrex::BaseFab<pd_model_type>> &mf = model;
        mf.FillBoundary(tmp_geom.periodicity());
        mf.FillBoundary();
        const int ncomp = mf.nComp();
        const int ng1 = 1;
        const int ng2 = 2;
        amrex::FabArray<amrex::BaseFab<pd_model_type>> tmpmf(mf.boxArray(), mf.DistributionMap(), ncomp, ng1);
        amrex::Copy(tmpmf, mf, 0, 0, ncomp, ng1);
        mf.ParallelCopy(tmpmf, 0, 0, ncomp, ng1, ng2, tmp_geom.periodicity());
    }
    //Util::Message(INFO, "Exit");
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
    if (! elastic.on) return;
    if (iter % elastic.interval) return;
    if (time < elastic.tstart) return;
    if (time > elastic.tend) return;

}

void 
PolymerDegradation::TimeStepBegin(amrex::Real time, int iter)
{
    if (!elastic.on) return;
    //if (time < elastic.tstart) return;
    //if (time > elastic.tend) return;

    if ((elastic.type == "tensile.single" || elastic.type == "single") && iter%elastic.interval) return;
    if ((elastic.type == "tensile.single" || elastic.type == "single") && time < elastic.tstart) return;
    if ((elastic.type == "tensile.single" || elastic.type == "single") && time > elastic.tend) return;
    if ((elastic.type == "tensile_test" || elastic.type == "tensile") && std::abs(time-elastic.test_time[elastic.current_test]) > 1.e-4) return;

    Util::Message(INFO);

    LPInfo info;
    info.setAgglomeration(elastic.agglomeration);
    info.setConsolidation(elastic.consolidation);
    info.setMaxCoarseningLevel(elastic.max_coarsening_level);

    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        
        eta_new[ilev]->FillBoundary();

        for (amrex::MFIter mfi(*material.model[ilev],true); mfi.isValid(); ++mfi)
        {
            amrex::Box box = mfi.growntilebox(2);
            amrex::Array4<const amrex::Real> const& eta_box = (*eta_new[ilev]).array(mfi);
            amrex::Array4<pd_model_type> const& modelfab = material.model[ilev]->array(mfi);

            amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k){
                amrex::Vector<Set::Scalar> _temp;
                for(int n=0; n<damage.number_of_eta; n++)
                {
                    Set::Scalar _temp2 = Numeric::Interpolate::CellToNodeAverage(eta_box,i,j,k,n);
                    _temp.push_back(std::min(damage.d_final[n],std::max(0.,_temp2)));
                }
                if(damage.type == "water") modelfab(i,j,k,0).DegradeModulus(_temp[0]);
                else if (damage.type == "water2") modelfab(i,j,k,0).DegradeModulus(_temp);
                else Util::Abort(INFO, "Damage model not implemented yet");
            });
        }
        Util::RealFillBoundary(*material.model[ilev],geom[ilev]);

        rhs[ilev]->setVal(0.0);
    }

    Operator::Elastic<pd_model_type::sym> elastic_op;
    elastic_op.define(geom, grids, dmap, info);

    elastic_op.setMaxOrder(elastic.linop_maxorder);
    
    for (int ilev = 0; ilev < nlevels; ++ilev)
    {
        //displacement[ilev]->setVal(0.0);
        const Real* DX = geom[ilev].CellSize();
        Set::Scalar volume = AMREX_D_TERM(DX[0],*DX[1],*DX[2]);

        AMREX_D_TERM(    rhs[ilev]->setVal(elastic.body_force[0]*volume,0,1);,
                        rhs[ilev]->setVal(elastic.body_force[1]*volume,1,1);,
                        rhs[ilev]->setVal(elastic.body_force[2]*volume,2,1););
    }


    if (elastic.type == "single" || elastic.type == "tensile_single")
    {
        //elastic.bc.SetTime(0.0);
        elastic.bc.Init(rhs,geom);

        elastic_op.SetBC(&(elastic.bc));

        //Util::Message(INFO);
        Solver::Nonlocal::Newton<pd_model_type> solver(elastic_op);
        solver.setMaxIter(elastic.max_iter);
        solver.setMaxFmgIter(elastic.max_fmg_iter);
        solver.setFixedIter(elastic.max_fixed_iter);
        solver.setVerbose(elastic.verbose);
        solver.setBottomVerbose(elastic.cgverbose);
        solver.setBottomMaxIter(elastic.bottom_max_iter);
        solver.setBottomTolerance(elastic.cg_tol_rel) ;
        solver.setBottomToleranceAbs(elastic.cg_tol_abs) ;
        
        for (int ilev = 0; ilev < nlevels; ilev++) if (displacement[ilev]->contains_nan()) Util::Warning(INFO);

        if (elastic.bottom_solver == "cg") solver.setBottomSolver(MLMG::BottomSolver::cg);
        else if (elastic.bottom_solver == "bicgstab") solver.setBottomSolver(MLMG::BottomSolver::bicgstab);
        solver.solve(displacement,rhs,material.model,elastic.tol_rel,elastic.tol_abs);
        solver.compResidual(residual,displacement,rhs,material.model);
        
        for (int lev = 0; lev < nlevels; lev++)
        {
            elastic_op.Strain(lev,*strain[lev],*displacement[lev]);
            elastic_op.Stress(lev,*stress[lev],*displacement[lev]);
            elastic_op.Energy(lev,*energy[lev],*displacement[lev]);
        }
        for (int lev = 0; lev < nlevels; lev++)
        {
            for (amrex::MFIter mfi(*stress[lev],true); mfi.isValid(); ++mfi)
            {
                const amrex::Box& box = mfi.validbox();
                amrex::Array4<const Set::Scalar> const& stress_box = (*stress[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const& stress_vm_box = (*stress_vm[lev]).array(mfi);
                amrex::Array4<const Set::Scalar> const& eta_box = (*eta_new[lev]).array(mfi);
                amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k){
                    Set::Matrix sigma = Numeric::FieldToMatrix(stress_box,i,j,k);
                    Set::Scalar temp = Numeric::Interpolate::CellToNodeAverage(eta_box,i,j,k,damage.number_of_eta-1);
                    Set::Matrix sigmadev = sigma - sigma.trace()/((double) AMREX_SPACEDIM)*Set::Matrix::Identity();
                    Set::Scalar temp2 = std::sqrt(1.5*sigmadev.squaredNorm());
                    stress_vm_box(i,j,k,0) =  temp2 < (1.-temp)*material.yieldstrength ? temp2 : (1.-temp)*material.yieldstrength;
                });
            }
        }
    }
    else if (elastic.type == "tensile" || elastic.type == "tensile_test")
    {
        Set::Scalar test_t = 0.;
        int countstep = 0;
        Util::Message(INFO, "Performing tensile test at t = ",elastic.test_time[elastic.current_test]);
        
        std::string plotfolder = "elastic_"+ std::to_string(elastic.current_test);

        while (test_t < elastic.test_duration)
        {
            Util::Message(INFO, "test time = ", test_t);
            test_t += elastic.test_dt; countstep++;

            for (int lev = 0; lev < nlevels; lev++) 
            {
                //displacement[lev]->setVal(0.0);
                const Real* DX = geom[lev].CellSize();
                Set::Scalar volume = AMREX_D_TERM(DX[0],*DX[1],*DX[2]);

                AMREX_D_TERM(    rhs[lev]->setVal(elastic.body_force[0]*volume,0,1);,
                                rhs[lev]->setVal(elastic.body_force[1]*volume,1,1);,
                                rhs[lev]->setVal(elastic.body_force[2]*volume,2,1););
            }
            
            elastic.bc.SetTime(test_t);
            elastic.bc.Init(rhs,geom);
            elastic_op.SetBC(&(elastic.bc));

            //Util::Message(INFO);
            Solver::Nonlocal::Newton<pd_model_type> solver(elastic_op);
            solver.setMaxIter(elastic.max_iter);
            solver.setMaxFmgIter(elastic.max_fmg_iter);
            solver.setFixedIter(elastic.max_fixed_iter);
            solver.setVerbose(elastic.verbose);
            solver.setBottomVerbose(elastic.cgverbose);
            solver.setBottomMaxIter(elastic.bottom_max_iter);
            solver.setBottomTolerance(elastic.cg_tol_rel) ;
            solver.setBottomToleranceAbs(elastic.cg_tol_abs) ;
            for (int ilev = 0; ilev < nlevels; ilev++) if (displacement[ilev]->contains_nan()) Util::Warning(INFO);

            if (elastic.bottom_solver == "cg") solver.setBottomSolver(MLMG::BottomSolver::cg);
            else if (elastic.bottom_solver == "bicgstab") solver.setBottomSolver(MLMG::BottomSolver::bicgstab);
            solver.solve(displacement, rhs, material.model, elastic.tol_rel, elastic.tol_abs);
            //solver.solve(GetVecOfPtrs(displacement), GetVecOfConstPtrs(rhs), elastic.tol_rel, elastic.tol_abs);
            //solver.compResidual(GetVecOfPtrs(residual),GetVecOfPtrs(displacement),GetVecOfConstPtrs(rhs));
            for (int lev = 0; lev < nlevels; lev++)
            {
                elastic_op.Strain(lev,*strain[lev],*displacement[lev]);
                elastic_op.Stress(lev,*stress[lev],*displacement[lev]);
                elastic_op.Energy(lev,*energy[lev],*displacement[lev]);
            }
            for (int lev = 0; lev < nlevels; lev++)
            {
                for (amrex::MFIter mfi(*stress[lev],true); mfi.isValid(); ++mfi)
                {
                    const amrex::Box& box = mfi.validbox();
                    amrex::Array4<const Set::Scalar> const& stress_box = (*stress[lev]).array(mfi);
                    amrex::Array4<Set::Scalar> const& stress_vm_box = (*stress_vm[lev]).array(mfi);
                    amrex::Array4<const Set::Scalar> const& eta_box = (*eta_new[lev]).array(mfi);
                    amrex::ParallelFor (box,[=] AMREX_GPU_DEVICE(int i, int j, int k){
                        Set::Matrix sigma = Numeric::FieldToMatrix(stress_box,i,j,k);
                        Set::Scalar temp = Numeric::Interpolate::CellToNodeAverage(eta_box,i,j,k,damage.number_of_eta-1);
                        Set::Matrix sigmadev = sigma - sigma.trace()/((double) AMREX_SPACEDIM)*Set::Matrix::Identity();
                        Set::Scalar temp2 = std::sqrt(1.5*sigmadev.squaredNorm());
                        stress_vm_box(i,j,k,0) =  temp2 < (1.-temp)*material.yieldstrength ? temp2 : (1.-temp)*material.yieldstrength;
                    });
                }
            }
            WritePlotFile(plotfolder,test_t,countstep-1);
        }
        elastic.current_test++;
    }
    else
        return;
}

}
//#endif
