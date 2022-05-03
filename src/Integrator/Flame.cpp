#include "Flame.H"
#include "IO/ParmParse.H"
#include "BC/Constant.H"
#include "Numeric/Stencil.H"
#include "IC/Laminate.H"
#include "IC/Constant.H"
#include "IC/PSRead.H"

namespace Integrator
{

    Flame::Flame() : Integrator()
    {
        BL_PROFILE("Integrator::Flame::Flame()");
        {
            // These are the phase field method parameters
            // that you use to inform the phase field method.
            IO::ParmParse pp("pf"); 
           // pp.query("M", pf.M); // Mobility parameter
            pp.query("eps", pf.eps); // Burn width thickness
            pp.query("kappa", pf.kappa); // Interface energy param
            pp.query("gamma",pf.gamma); // Scaling factor for mobility
            pp.query("lambda",pf.lambda); // Chemical potential multiplier
            pp.query("w1", pf.w1); // Unburned rest energy
            pp.query("w12", pf.w12);  // Barrier energy
            pp.query("w0", pf.w0);    // Burned rest energy
            pp.query("P", pf.P);             // Pressure [UNITS?]
            pp.query("r_ap", pf.r_ap);       // AP Power law multiplier
            pp.query("n_ap", pf.n_ap);       // AP Power law exponent
            pp.query("r_htpb", pf.r_htpb);   // HTPB Power law multiplier
            pp.query("n_htpb", pf.n_htpb);   // HTPB Power law exponent
            pp.query("r_comb", pf.r_comb);   // Combination power law multiplier
            pp.query("n_comb", pf.n_comb);   // Combination power law exponent

            EtaBC = new BC::Constant(1);
            pp.queryclass("eta.bc", *static_cast<BC::Constant *>(EtaBC)); // See :ref:`BC::Constant`

            RegisterNewFab(Eta_mf, EtaBC, 1, 1, "Eta", true);
            RegisterNewFab(Eta_old_mf, EtaBC, 1, 1, "Eta_old", false);

            RegisterNewFab(mdot_mf, EtaBC, 1, 1, "mdot", true);

            RegisterNewFab(heat_mf, EtaBC, 1, 1, "heat", true);
        }
        

        {
            IO::ParmParse pp("thermal");
            pp.query("on",thermal.on); // Whether to use the Thermal Transport Model
            //Util::Message(INFO,"thermal.on =" , thermal.on);
            pp.query("rho_ap",thermal.rho_ap); // AP Density
            pp.query("rho_htpb", thermal.rho_htpb); // HTPB Density
            pp.query("k_ap", thermal.k_ap); // AP Thermal Conductivity
            pp.query("k_htpb",thermal.k_htpb); // HTPB Thermal Conductivity
            pp.query("k0", thermal.k0); // Thermal conductivity 
            pp.query("k_comb", thermal.k_comb); // Combined Thermal Conductivity
            pp.query("cp_ap", thermal.cp_ap); // AP Specific Heat
            pp.query("cp_htpb", thermal.cp_htpb); //HTPB Specific Heat
            pp.query("q_ap", thermal.q_ap); // AP  Thermal Flux
            pp.query("q_htpb", thermal.q_htpb); // HTPB Thermal Flux
            pp.query("q_comb" , thermal.q_comb); // Interface heat flux
            pp.query("ae_ap", thermal.ae_ap); // AP Activation Energy
            pp.query("ae_htpb", thermal.ae_htpb); // HTPB Activation Energy
            pp.query("ae_comb", thermal.ae_comb);
            pp.query("bound", thermal.bound);
            pp.query("addtemp", thermal.addtemp);
            pp.query("A_ap", thermal.A_ap);
            pp.query("A_htpb", thermal.A_htpb);
            pp.query("A_comb", thermal.A_comb);
            pp.query("beta1", thermal.beta1);
            pp.query("beta2", thermal.beta2);
            pp.query("MinTemp", thermal.MinTemp );

            pp.query("temperature_delay", thermal.temperature_delay); 

            if (thermal.on)
            {
            TempBC = new BC::Constant(1);
            pp.queryclass("temp.bc", *static_cast<BC::Constant *>(TempBC));

            RegisterNewFab(Temp_mf, TempBC, 1, 1, "Temp", true);
            RegisterNewFab(Temp_old_mf, TempBC, 1, 1, "Temp_old", false);
                        
            RegisterNewFab(Mob_mf, TempBC, 1, 1, "Mob", true);
            RegisterNewFab(Mob_old_mf, TempBC, 1, 1, "Mob_old", false);            

            }
        }



        {
            // Additional AMR parameters
            IO::ParmParse pp ("amr");
            pp.query("refinement_criterion", m_refinement_criterion); // Refinement criterion for eta field
            pp.query("refinement_criterion_temp", t_refinement_criterion); // Refinement criterion for temperature field
        }

        std::vector<Set::Scalar> fs(fs_number, 1);
        
        {
            // The material field is referred to as :math:`\phi(\mathbf{x})` and is 
            // specified using these parameters. 
            IO::ParmParse pp("phi.ic");
            std::string type = "packedspheres";
            pp.query("type", type); // IC type - [packedspheres,laminate] - see classes for more information
            if (type == "psread")
            {
                                
                PhiIC = new IC::PSRead(geom);
                pp.queryclass("psread", *static_cast<IC::PSRead *>(PhiIC)); // See :ref:`IC::PackedSpheres`
            }
            else if (type == "laminate")
            {
                PhiIC = new IC::Laminate(geom);
                pp.queryclass("laminate", *static_cast<IC::Laminate *>(PhiIC)); // See :ref:`IC::Laminate`
            }
            else if (type == "constant")
            {
                PhiIC = new IC::Constant(geom);
                pp.queryclass("constant",*static_cast<IC::Constant *>(PhiIC)); // See :ref:`IC::Constant`
            }
            
            RegisterNewFab(phi_mf, EtaBC, 1, 1, "phi", true);
        }

        //
        // This code is mostly temporary.
        //
        Util::Message(INFO,"Checking for unused inputs...");
        int unused_inputs = IO::ParmParse::AllUnusedInputs();
        cout << IO::ParmParse::AllUnusedInputs();
        if (unused_inputs > 0)
        {
            Util::Warning(INFO,"There are a lot of input parameters in the ");
            Util::Warning(INFO,"input file that are no longer valid. Many of them have been");
            Util::Warning(INFO,"renamed, removed, or replaced; the input file must be changed");
            Util::Warning(INFO,"to match.");
            Util::Warning(INFO,"");
            Util::Warning(INFO,"Note: you can run");
            Util::Warning(INFO,"   >   make docs");
            Util::Warning(INFO,"   >   xdg-open docs/build/html/Inputs.html");
            Util::Warning(INFO,"to look at the parameter list.\n");
            Util::Abort(INFO, "Terminating for now until the input decks are brought up to speed.");
        }
    }

    void Flame::Initialize(int lev)
    {
        BL_PROFILE("Integrator::Flame::Initialize");

        if (thermal.on) Temp_mf[lev]->setVal(thermal.bound);
        if (thermal.on) Temp_old_mf[lev]->setVal(thermal.bound);

        Eta_mf[lev]->setVal(1.0);
        Eta_old_mf[lev]->setVal(1.0);
        
        mdot_mf[lev] -> setVal(0.0);

        PhiIC->Initialize(lev, phi_mf);
	  
        heat_mf[lev] -> setVal(0.0);

        if (thermal.on) Mob_mf[lev]->setVal(0.0);
        if (thermal.on) Mob_old_mf[lev]->setVal(0.0);
    }

    
    void Flame::TimeStepBegin(Set::Scalar /*a_time*/, int a_iter)
    {
        BL_PROFILE("Integrator::Flame::TimeStepBegin");
        if (!elastic.interval)
            return;
        if (a_iter % elastic.interval)
            return;

        for (int lev = 0; lev <= finest_level; ++lev)
        {
            Util::RealFillBoundary(*Eta_mf[lev], geom[lev]);
            Util::RealFillBoundary(*phi_mf[lev], geom[lev]);
            Util::RealFillBoundary(*Temp_mf[lev], geom[lev]);

            elastic.rhs_mf[lev]->setVal(0.0);
            elastic.disp_mf[lev]->setVal(0.0);
            elastic.model_mf[lev]->setVal(elastic.model_ap);

            Set::Vector DX(geom[lev].CellSize());

            for (MFIter mfi(*elastic.model_mf[lev], true); mfi.isValid(); ++mfi)
            {

                amrex::Box bx = mfi.nodaltilebox();
                bx.grow(1);
                amrex::Array4<model_type> const &model = elastic.model_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const &eta = Eta_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const &phi = phi_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const &temp = Temp_mf[lev]->array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                                   {
                Set::Scalar phi_avg = Numeric::Interpolate::CellToNodeAverage(phi,i,j,k,0);
                Set::Scalar eta_avg = Numeric::Interpolate::CellToNodeAverage(eta,i,j,k,0);
                Set::Scalar temp_avg = Numeric::Interpolate::CellToNodeAverage(temp,i,j,k,0);

                model_type model_ap = elastic.model_ap;
                model_ap.F0 *= temp_avg;
                model_type model_htpb = elastic.model_htpb;
                model_htpb.F0 *= temp_avg;

                model_type solid = model_ap*phi_avg + model_htpb*(1.-phi_avg);

                model(i,j,k) = solid*eta_avg + elastic.model_void*(1.0-eta_avg); });
            }

            Util::RealFillBoundary(*elastic.model_mf[lev], geom[lev]);
        }
        elastic.bc.Init(elastic.rhs_mf, geom);

        amrex::LPInfo info;
        Operator::Elastic<Model::Solid::Affine::Isotropic::sym> elastic_op(Geom(0, finest_level), grids, DistributionMap(0, finest_level), info);
        elastic_op.SetUniform(false);
        elastic_op.SetBC(&elastic.bc);
        elastic_op.SetAverageDownCoeffs(true);

        Set::Scalar tol_rel = 1E-8, tol_abs = 1E-8;
        // Parameters for the elastic solver (when used with elasticity)
        IO::ParmParse pp("elastic");
        elastic.solver = new Solver::Nonlocal::Newton<Model::Solid::Affine::Isotropic>(elastic_op);
        pp.queryclass("solver", *elastic.solver); // See :ref:`Solver::Nonlocal::Newton`

        elastic.solver->solve(elastic.disp_mf, elastic.rhs_mf, elastic.model_mf, tol_rel, tol_abs);
        elastic.solver->compResidual(elastic.res_mf, elastic.disp_mf, elastic.rhs_mf, elastic.model_mf);

        for (int lev = 0; lev <= elastic.disp_mf.finest_level; lev++)
        {
            const amrex::Real *DX = geom[lev].CellSize();
            for (MFIter mfi(*elastic.disp_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.nodaltilebox();
                amrex::Array4<model_type> const &model = elastic.model_mf[lev]->array(mfi);
                amrex::Array4<Set::Scalar> const &stress = elastic.stress_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const &disp = elastic.disp_mf[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                                   {
                                    std::array<Numeric::StencilType, AMREX_SPACEDIM>
                                        sten = Numeric::GetStencil(i, j, k, bx);
                                    if (model(i, j, k).kinvar == Model::Solid::KinematicVariable::F)
                                    {
                                        Set::Matrix F = Set::Matrix::Identity() + Numeric::Gradient(disp, i, j, k, DX, sten);
                                        Set::Matrix P = model(i, j, k).DW(F);
                                        Numeric::MatrixToField(stress, i, j, k, P);
                                    }
                                    else
                                    {
                                        Set::Matrix gradu = Numeric::Gradient(disp, i, j, k, DX, sten);
                                        Set::Matrix sigma = model(i, j, k).DW(gradu);
                                        Numeric::MatrixToField(stress, i, j, k, sigma);
                                    } });
            }
        }
    }


    void Flame::Advance(int lev, amrex::Real time, amrex::Real dt)
    {
        BL_PROFILE("Integrador::Flame::Advance");
        const amrex::Real *DX = geom[lev].CellSize();

        if (lev == finest_level)
        {
            std::swap(Eta_old_mf[lev], Eta_mf[lev]);
            std::swap(Temp_old_mf[lev], Temp_mf[lev]);
            std::swap(Mob_old_mf[lev], Mob_mf[lev]);

            Set::Scalar a0 = pf.w0, 
                        a1 = 0.0,  
                        a2 = -5.0 * pf.w1 + 16.0 * pf.w12 - 11.0 * a0, 
                        a3 = 14.0 * pf.w1 - 32.0 * pf.w12 + 18.0 * a0, 
                        a4 = -8.0 * pf.w1 + 16.0 * pf.w12 -  8.0 * a0;
                        
            for (amrex::MFIter mfi(*Eta_mf[lev], true); mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.tilebox();
                amrex::Array4<Set::Scalar> const &Eta = (*Eta_mf[lev]).array(mfi);
                amrex::Array4<const Set::Scalar> const &Eta_old = (*Eta_old_mf[lev]).array(mfi);
                amrex::Array4<const Set::Scalar> const &phi = (*phi_mf[lev]).array(mfi);

                amrex::Array4<Set::Scalar> const  &Mob = (*Mob_mf[lev]).array(mfi);
                amrex::Array4<const Set::Scalar> const &Mob_old = (*Mob_old_mf[lev]).array(mfi);
                
                amrex::Array4<Set::Scalar> const &Temp = (*Temp_mf[lev]).array(mfi);
                amrex::Array4<const Set::Scalar> const &Temp_old = (*Temp_old_mf[lev]).array(mfi);

                amrex::Array4<Set::Scalar> const &mdot = (*mdot_mf[lev]).array(mfi);

                amrex::Array4<Set::Scalar> const &heat = (*heat_mf[lev]).array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar eta_lap = Numeric::Laplacian(Eta_old, i, j, k, 0, DX);
                    amrex::Real rho = (thermal.rho_ap - thermal.rho_htpb) * Eta_old(i,j,k) + thermal.rho_htpb;

                    // --------------------
                    // Eta computation
                    // --------------------

		    if (Eta_old(i,j,k) < 0.01) 
                    {
                        Eta(i,j,k) = 0.0;
                    }
		    else
                    {   
                        Eta(i, j, k) = 
                            Eta_old(i, j, k) 
                            - Mob_old(i,j,k) * dt * (
                            (pf.lambda/pf.eps)*(a1 + 2.0 * a2 * Eta_old(i, j, k) + 3.0 * a3 * Eta_old(i, j, k) * Eta_old(i, j, k) + 4 * a4 * Eta_old(i, j, k) * Eta_old(i, j, k) * Eta_old(i, j, k))
                            - pf.eps * pf.kappa * eta_lap);
		            }

                    Set::Vector eta_grad = Numeric::Gradient(Eta_old, i, j, k, 0, DX);

                    // --------------------
                    // Mass flux
                    // --------------------

                    {
			            mdot(i,j,k) = ((thermal.rho_ap - thermal.rho_htpb) * Eta_old(i,j,k) + thermal.rho_htpb) * Numeric::Laplacian(Eta_old, i, j, k, 0, DX);
                    }

                    
                    Set::Vector temp_grad = Numeric::Gradient(Temp_old, i, j, k, 0, DX);
                    Set::Scalar temp_lap = Numeric::Laplacian(Temp_old, i, j, k, 0, DX);
                    Set::Scalar eta_grad_mag = eta_grad.lpNorm<2>();
                    Set::Vector normvec = eta_grad/Eta_old(i,j,k);
                    
                    amrex::Real K_ap = (thermal.k_ap - thermal.k0) * Eta_old(i,j,k) + thermal.k0;
                    amrex::Real K_htpb = (thermal.k_htpb - thermal.k0) * Eta_old(i,j,k) + thermal.k0;
                    amrex::Real K_comb = (thermal.k_comb - thermal.k0) * Eta_old(i,j,k) + thermal.k0;

                    amrex::Real K = K_ap * phi(i,j,k) + K_htpb * (1 - phi(i,j,k)) + 4.0 * phi(i,j,k) * (1 - phi(i,j,k)) * K_comb;  

                    amrex::Real cp = (thermal.cp_ap - thermal.cp_htpb) * Eta_old(i,j,k) + thermal.cp_htpb;

                    Set::Scalar test = normvec.dot(temp_grad);

                    Set::Scalar Bn = thermal.q_ap * phi(i,j,k) + thermal.q_htpb * (1 - phi(i,j,k)) + 4.0 * phi(i,j,k) * (1 - phi(i,j,k)) * thermal.q_comb;


                    // --------------------
                    // Temperature computation
                    // --------------------
                    if (Eta_old(i,j,k) > 0.001)
			        { 
                        Temp(i,j,k) = Temp_old(i,j,k) + dt*(K/cp/rho) * (temp_lap + test + eta_grad_mag/Eta_old(i,j,k) * Bn / K );
			        }
                    else
                    {
                        Temp(i,j,k) = 0.0; 
                    }

	
		    if(Temp(i,60,k) >= 850.0){
		      Util::Message(INFO, "Temp = ", Temp(i,60,k));
		      Util::Message(INFO, "Eta = ", Eta(i,60,k));
		      Util::Message(INFO, "Mob = ", Mob(i,60,k));
		    }

                    // --------------------
                    // Mobility computation
                    // --------------------  
		    Set::Scalar R = 8.314;
		    Set::Scalar T0 = 273.15;

		    Set::Scalar e1 = exp(-thermal.ae_ap   / pf.P / R / (Temp_old(i,j,k) - T0) );
		    Set::Scalar e2 = exp(-thermal.ae_htpb / pf.P / R / (Temp_old(i,j,k) - T0) );
		    Set::Scalar e3 = exp(-thermal.ae_comb / pf.P / R / (Temp_old(i,j,k) - T0) );

		    Mob(i,j,k) = ( 
				  thermal.A_ap   * e1 * phi(i,j,k) + 
				  thermal.A_htpb * e2 * (1.0 - phi(i,j,k)) + 
				  thermal.A_comb * e3 * 4.0 * phi(i,j,k) * (1.0 - phi(i,j,k))
				 ) / pf.gamma / (pf.w1 - pf.w0);
            
                    // -------------------
                    // Fluid Heat
                    // -------------------

		    
                    
                });
                
            }
        }
    }


    void Flame::TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, amrex::Real /*time*/, int /*ngrow*/)
    {
        BL_PROFILE("Integrator::Flame::TagCellsForRefinement");
        
        const amrex::Real *DX = geom[lev].CellSize();
        Set::Scalar dr = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));

        // Eta criterion for refinement
        for (amrex::MFIter mfi(*Eta_mf[lev], true); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.tilebox();
            amrex::Array4<char> const &tags = a_tags.array(mfi);
            amrex::Array4<const Set::Scalar> const &Eta = (*Eta_mf[lev]).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                               {
                Set::Vector gradeta = Numeric::Gradient(Eta, i, j, k, 0, DX);
                if (gradeta.lpNorm<2>() * dr * 2 > m_refinement_criterion)
                    tags(i, j, k) = amrex::TagBox::SET; });
        }

        // Thermal criterion for refinement
        
        if (thermal.on)
        {
            for (amrex::MFIter mfi(*Temp_mf[lev], true); mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.tilebox();
                amrex::Array4<char> const &tags = a_tags.array(mfi);
                amrex::Array4<const Set::Scalar> const &Temp = (*Temp_mf[lev]).array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                                   {
                    Set::Vector tempgrad = Numeric::Gradient(Temp, i, j, k, 0, DX);
                    if (tempgrad.lpNorm<2>() * dr  > t_refinement_criterion){
                        tags(i, j, k) = amrex::TagBox::SET;
                        //Util::Message(INFO, "Refinament = ", tempgrad.lpNorm<2>() * dr)
                        ;}
                        });
            }
        } 
         
    }
    void Flame::Regrid(int lev, Set::Scalar /* time */)
    {
        BL_PROFILE("Integrator::Flame::Regrid");
        if (lev < finest_level) return;
        phi_mf[lev]->setVal(0.0);
        PhiIC->Initialize(lev, phi_mf);
        Util::Message(INFO, "Regridding on level ", lev);
    } 
} // namespace Integrator


