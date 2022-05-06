#include "Flame.H"
#include "IO/ParmParse.H"
#include "BC/Constant.H"
#include "Numeric/Stencil.H"
#include "IC/Laminate.H"
#include "IC/Constant.H"
#include "IC/PSRead.H"
#include "Numeric/Function.H"

namespace Integrator
{

    Flame::Flame() : MechanicsBase<Model::Solid::Affine::Isotropic>() {}

    Flame::Flame(IO::ParmParse &pp) : MechanicsBase<Model::Solid::Affine::Isotropic>() 
    {pp.queryclass(*this);}

    void 
    Flame::Parse(Flame &value, IO::ParmParse &pp)
    {
        BL_PROFILE("Integrator::Flame::Flame()");
        {
            // These are the phase field method parameters
            // that you use to inform the phase field method.
            pp.query("pf.eps", value.pf.eps); // Burn width thickness
            pp.query("pf.kappa", value.pf.kappa); // Interface energy param
            pp.query("pf.gamma",value.pf.gamma); // Scaling factor for mobility
            pp.query("pf.lambda",value.pf.lambda); // Chemical potential multiplier
            pp.query("pf.w1", value.pf.w1); // Unburned rest energy
            pp.query("pf.w12", value.pf.w12);  // Barrier energy
            pp.query("pf.w0", value.pf.w0);    // Burned rest energy
            pp.query("pf.P", value.pf.P);             // Pressure [UNITS?]
            pp.query("pf.r_ap", value.pf.r_ap);       // AP Power law multiplier
            pp.query("pf.n_ap", value.pf.n_ap);       // AP Power law exponent
            pp.query("pf.r_htpb", value.pf.r_htpb);   // HTPB Power law multiplier
            pp.query("pf.n_htpb", value.pf.n_htpb);   // HTPB Power law exponent
            pp.query("pf.r_comb", value.pf.r_comb);   // Combination power law multiplier
            pp.query("pf.n_comb", value.pf.n_comb);   // Combination power law exponent
            pp.query("pf.m0",value.pf.m0);
            pp.query("pf.T0",value.pf.T0);
            pp.query("pf.Ea",value.pf.Ea);

            value.bc_eta = new BC::Constant(1);
            pp.queryclass("pf.eta.bc", *static_cast<BC::Constant *>(value.bc_eta)); // See :ref:`BC::Constant`

            value.RegisterNewFab(value.eta_mf,     value.bc_eta, 1, 1, "eta", true);
            value.RegisterNewFab(value.eta_old_mf, value.bc_eta, 1, 1, "eta_old", false);
            value.RegisterNewFab(value.mdot_mf,    value.bc_eta, 1, 1, "mdot", true);
        }
        

        {
            //IO::ParmParse pp("thermal");
            pp.query("thermal.on",value.thermal.on); // Whether to use the Thermal Transport Model
            //Util::Message(INFO,"thermal.on =" , thermal.on);
            pp.query("thermal.rho_ap",value.thermal.rho_ap); // AP Density
            pp.query("thermal.rho_htpb", value.thermal.rho_htpb); // HTPB Density
            pp.query("thermal.k_ap", value.thermal.k_ap); // AP Thermal Conductivity
            pp.query("thermal.k_htpb",value.thermal.k_htpb); // HTPB Thermal Conductivity
            pp.query("thermal.k0", value.thermal.k0); // Thermal conductivity 
            pp.query("thermal.k_comb", value.thermal.k_comb); // Combined Thermal Conductivity
            pp.query("thermal.cp_ap", value.thermal.cp_ap); // AP Specific Heat
            pp.query("thermal.cp_htpb", value.thermal.cp_htpb); //HTPB Specific Heat
            pp.query("thermal.q_ap", value.thermal.q_ap); // AP  Thermal Flux
            pp.query("thermal.q_htpb", value.thermal.q_htpb); // HTPB Thermal Flux
            pp.query("thermal.q_comb" , value.thermal.q_comb); // Interface heat flux
            pp.query("thermal.q0",value.thermal.q0); // Baseline heat flux
            pp.query("thermal.ae_ap", value.thermal.ae_ap); // AP Activation Energy
            pp.query("thermal.ae_htpb", value.thermal.ae_htpb); // HTPB Activation Energy
            pp.query("thermal.ae_comb", value.thermal.ae_comb);
            pp.query("thermal.bound", value.thermal.bound);
            pp.query("thermal.addtemp", value.thermal.addtemp);
            pp.query("thermal.A_ap", value.thermal.A_ap);
            pp.query("thermal.A_htpb", value.thermal.A_htpb);
            pp.query("thermal.A_comb", value.thermal.A_comb);
            pp.query("thermal.beta1", value.thermal.beta1);
            pp.query("thermal.beta2", value.thermal.beta2);
            pp.query("thermal.MinTemp", value.thermal.MinTemp );

            pp.query("thermal.temperature_delay", value.thermal.temperature_delay); 

            if (value.thermal.on)
            {
                value.bc_temp = new BC::Constant(1);
                pp.queryclass("thermal.temp.bc", *static_cast<BC::Constant *>(value.bc_temp));
                value.RegisterNewFab(value.temp_mf, value.bc_temp, 1, 1, "temp", true);
                value.RegisterNewFab(value.temp_old_mf, value.bc_temp, 1, 1, "temp_old", false);
                value.RegisterNewFab(value.Mob_mf, value.bc_temp, 1, 1, "Mob", true);
                value.RegisterNewFab(value.alpha_mf,value.bc_temp,1,1,"alpha",true);
            }
        }


        // Refinement criterion for eta field
        pp.query("amr.refinement_criterion", value.m_refinement_criterion); 
        // Refinement criterion for temperature field
        pp.query("amr.refinement_criterion_temp", value.t_refinement_criterion); 


        {
            // The material field is referred to as :math:`\phi(\mathbf{x})` and is 
            // specified using these parameters. 
            //IO::ParmParse pp("phi.ic");
            std::string type = "packedspheres";
            pp.query("phi.ic.type", type); // IC type (psread, laminate, constant)
            if      (type == "psread")   value.ic_phi = new IC::PSRead(value.geom, pp, "phi.is.psread");
            else if (type == "laminate") value.ic_phi = new IC::Laminate(value.geom,pp,"phi.ic.laminate");
            else if (type == "constant") value.ic_phi = new IC::Constant(value.geom,pp,"phi.ic.constant");
            else Util::Abort(INFO,"Invalid IC type ",type);
            
            value.RegisterNewFab(value.phi_mf, value.bc_eta, 1, 1, "phi", true);
        }

        pp.queryclass<MechanicsBase<Model::Solid::Affine::Isotropic>>("elastic",value);
        if (value.m_type  != Type::Disable)
        {
            pp.queryclass("model_ap",value.elastic.model_ap);
            Util::Message(INFO,value.elastic.model_ap);
            pp.queryclass("model_htpb",value.elastic.model_htpb);
            Util::Message(INFO,value.elastic.model_htpb);
            pp.queryclass("model_void",value.elastic.model_void);
            Util::Message(INFO,value.elastic.model_void);
        }
    }

    void Flame::Initialize(int lev)
    {
        BL_PROFILE("Integrator::Flame::Initialize");
        Util::Message(INFO,m_type);
        MechanicsBase<Model::Solid::Affine::Isotropic>::Initialize(lev);

        if (thermal.on)
        {
            temp_mf[lev]->setVal(thermal.bound);
            temp_old_mf[lev]->setVal(thermal.bound);
            alpha_mf[lev]->setVal(0.0);
            //heat_mf[lev] -> setVal(0.0);
            Mob_mf[lev]->setVal(0.0);
        } 
        
        eta_mf[lev]->setVal(1.0);
        eta_old_mf[lev]->setVal(1.0);
        
        mdot_mf[lev] -> setVal(0.0);

        ic_phi->Initialize(lev, phi_mf);
    }

    
    void Flame::UpdateModel(int a_step)
    {
        if (a_step % m_interval) return;

        for (int lev = 0; lev <= finest_level; ++lev)
        {
            phi_mf[lev]->FillBoundary();
            eta_mf[lev]->FillBoundary();
            temp_mf[lev]->FillBoundary();
            //Util::RealFillBoundary(*phi_mf[lev],geom[lev]);
            //Util::RealFillBoundary(*eta_mf[lev],geom[lev]);   
            //Util::RealFillBoundary(*temp_mf[lev],geom[lev]);   
            
            for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.nodaltilebox();
                //bx.grow(1);
                amrex::Array4<model_type>        const &model = model_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const &eta = eta_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const &phi = phi_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const &temp = temp_mf[lev]->array(mfi);

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
                    model(i,j,k) = solid*eta_avg + elastic.model_void*(1.0-eta_avg); 
                });
            }

            Util::RealFillBoundary(*model_mf[lev],geom[lev]);
        }
    }
    
    void Flame::TimeStepBegin(Set::Scalar a_time, int a_iter)
    {
        BL_PROFILE("Integrator::Flame::TimeStepBegin");
        MechanicsBase<Model::Solid::Affine::Isotropic>::TimeStepBegin(a_time,a_iter);
    }


    void Flame::Advance(int lev, Set::Scalar time, Set::Scalar dt)
    {
        BL_PROFILE("Integrador::Flame::Advance");
        MechanicsBase<Model::Solid::Affine::Isotropic>::Advance(lev,time,dt);

        const Set::Scalar *DX = geom[lev].CellSize();
        const Set::Scalar small = 1E-8;

        if (true) // (lev == finest_level)
        {
            std::swap(eta_old_mf[lev], eta_mf[lev]);
            std::swap(temp_old_mf[lev], temp_mf[lev]);

            Numeric::Function::Polynomial<4> w(pf.w0, 
                                            0.0, 
                                            -5.0 * pf.w1 + 16.0 * pf.w12 - 11.0 * pf.w0,
                                            14.0 * pf.w1 - 32.0 * pf.w12 + 18.0 * pf.w0,
                                            -8.0 * pf.w1 + 16.0 * pf.w12 -  8.0 * pf.w0 );
            Numeric::Function::Polynomial<3> dw = w.D();

            for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.tilebox();
                // Phase fields
                amrex::Array4<Set::Scalar> const &etanew = (*eta_mf[lev]).array(mfi);
                amrex::Array4<const Set::Scalar> const &eta = (*eta_old_mf[lev]).array(mfi);
                amrex::Array4<const Set::Scalar> const &phi = (*phi_mf[lev]).array(mfi);

                // Heat transfer fields
                amrex::Array4<Set::Scalar>       const &tempnew = (*temp_mf[lev]).array(mfi);
                amrex::Array4<const Set::Scalar> const &temp    = (*temp_old_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar>       const &alpha   = (*alpha_mf[lev]).array(mfi);

                // Diagnostic fields
                amrex::Array4<Set::Scalar> const  &Mob = (*Mob_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const &mdot = (*mdot_mf[lev]).array(mfi);
                //amrex::Array4<Set::Scalar> const &heat = (*heat_mf[lev]).array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar eta_lap = Numeric::Laplacian(eta, i, j, k, 0, DX);

                    // =============== TODO ==================
                    // This part is probably ok for now. Howver I want to look into
                    // splitting the mobility into two sections so that curvature
                    // is not temperature-dependent. 
                    etanew(i, j, k) = 
                            eta(i, j, k) 
                            - Mob(i,j,k) * dt * ( 
                            (pf.lambda/pf.eps)*dw(eta(i,j,k)) - pf.eps * pf.kappa * eta_lap);
                    
                    // Calculate effective thermal conductivity
                    // No special interface mixure rule is needed here.
                    Set::Scalar K   = thermal.k_ap   * phi(i,j,k) + thermal.k_htpb   * (1.0 - phi(i,j,k));
                    Set::Scalar rho = thermal.rho_ap * phi(i,j,k) + thermal.rho_htpb * (1.0 - phi(i,j,k));
                    Set::Scalar cp  = thermal.cp_ap  * phi(i,j,k) + thermal.cp_htpb  * (1.0 - phi(i,j,k));
                    //K = (K - thermal.k0) * eta(i,j,k) + thermal.k0;

                    // Calculate thermal diffusivity and store in field
                    alpha(i,j,k) = eta(i,j,k) * K / rho / cp;

                    // Calculate mass flux
                    mdot(i,j,k) = - rho * (etanew(i,j,k) - eta(i,j,k)) / dt;
                });

                    
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Vector grad_eta = Numeric::Gradient(eta, i, j, k, 0, DX);
                    Set::Vector grad_temp = Numeric::Gradient(temp, i, j, k, 0, DX);
                    Set::Scalar lap_temp = Numeric::Laplacian(temp, i, j, k, 0, DX);
                    Set::Scalar grad_eta_mag = grad_eta.lpNorm<2>();
                    Set::Vector grad_alpha = Numeric::Gradient(alpha,i,j,k,0,DX);

                    // =============== TODO ==================
                    // We need to get heat flux from mass flux HERE
                    // This is a primitive preliminary implementation.
                    // Note: "thermal.q0" is an initiation heat flux - think of it
                    // like a laser that is heating up the interface. 
                    Set::Scalar qdot = 
                        thermal.q_ap * phi(i,j,k) + thermal.q_htpb * (1-phi(i,j,k)) + 4.0 * thermal.q_comb * phi(i,j,k) * (1.-phi(i,j,k))
                        +
                        thermal.q0;
		    //Set::Scalar qdot_ap = (pf.P * 0.11 + 0.34) * phi(i,j,k) / thermal.k_ap;
		    //Set::Scalar qdot_htpb = (pf.P * 0.05 + 0.09) * (1.0 - phi(i,j,k)) / thermal.k_htpb;
		    //Set::Scalar qdot_comb = (pf.P * 0.17719 + 2.0) * 4.0 * phi(i,j,k) * (1.0 - phi(i,j,k)) / thermal.k_comb;

		    //qdot = qdot_ap + qdot_htpb + qdot_comb + thermal.q0;

		    // Note: This qdot equations work for Pressure in atm units and return qdot in kW/cmˆ2 units. I am probably going to change this equation.


                    //
                    // Evolve temperature with the qdot flux term in place
                    //
                    // Calculate modified spatial derivative of temperature
                    Set::Scalar dTdt = 0.0;
                    dTdt += grad_eta.dot(alpha(i,j,k) * grad_temp) / (eta(i,j,k) + small);
                    dTdt += grad_alpha.dot(grad_temp);
                    dTdt += alpha(i,j,k) * lap_temp;        
                    // Calculate the source term
                    dTdt += grad_eta_mag * alpha(i,j,k) * qdot / (eta(i,j,k) + small);
                    // Now, evolve temperature with explicit forward Euler
                    tempnew(i,j,k) = temp(i,j,k) + dt * dTdt;
                    
                    
                    // =============== TODO ==================
                    // We need to more accurately calculate our effective mobility here.
                    // Right now it uses a simple three-parameter exponential model, however,
                    // it does not take material heterogeneity into account. This model
                    // should account for the differing mobilities for AP and HTPB
                    // (but not necessarily for the combination region).
                    Mob(i,j,k) = pf.m0 * exp(- pf.Ea / (pf.T0 + tempnew(i,j,k)));//grad_eta_mag;




                    //// --------------------
                    //// Temperature computation
                    //// --------------------
                    //if (eta(i,j,k) > 0.001)
                    //{ 
                    //    tempnew(i,j,k) = temp(i,j,k) + dt*(K/cp/rho) * (temp_lap + test + eta_grad_mag/eta(i,j,k) * Bn / K );
                    //}
                    //else
                    //{
                    //    tempnew(i,j,k) = 0.0; 
                    //}

                    /*
                    // --------------------
                    // Mobility computation
                    // --------------------

		    if (Eta_old(i,j,k) < 0.01) 
                    {
                        Mob(i,j,k) = 0.0;
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
                        Mob(i,j,k) = 0.01;
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
                        Set::Scalar R = 8.314;
                        Set::Scalar T0 = 399.0;
                        
                        Set::Scalar e1 = exp(-thermal.ae_ap   / pf.P / R / (temp(i,j,k) - T0) );
                        Set::Scalar e2 = exp(-thermal.ae_htpb / pf.P / R / (temp(i,j,k) - T0) );
                        Set::Scalar e3 = exp(-thermal.ae_comb / pf.P / R / (temp(i,j,k) - T0) );

	
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

                    if (eta(i,j,k) > 0.001)
                    {
                        heat(i,j,k) = mdot(i,j,k) * thermal.cp_ap * temp(i,j,k);
                    }
                    else
                    {
                        heat(i,j,k) = mdot(i,j,k) * thermal.cp_ap * temp(i,j,k);
                    }
                    */

                });
                
            }
        }
    }


    void Flame::TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, Set::Scalar time, int ngrow)
    {
        BL_PROFILE("Integrator::Flame::TagCellsForRefinement");
        MechanicsBase<Model::Solid::Affine::Isotropic>::TagCellsForRefinement(lev,a_tags,time,ngrow);
        
        const Set::Scalar *DX = geom[lev].CellSize();
        Set::Scalar dr = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));

        // Eta criterion for refinement
        for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.tilebox();
            amrex::Array4<char> const &tags = a_tags.array(mfi);
            amrex::Array4<const Set::Scalar> const &Eta = (*eta_mf[lev]).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                               {
                Set::Vector gradeta = Numeric::Gradient(Eta, i, j, k, 0, DX);
                if (gradeta.lpNorm<2>() * dr * 2 > m_refinement_criterion)
                    tags(i, j, k) = amrex::TagBox::SET; });
        }

        // Thermal criterion for refinement
        
        if (thermal.on)
        {
            for (amrex::MFIter mfi(*temp_mf[lev], true); mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.tilebox();
                amrex::Array4<char> const &tags = a_tags.array(mfi);
                amrex::Array4<const Set::Scalar> const &temp = (*temp_mf[lev]).array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Vector tempgrad = Numeric::Gradient(temp, i, j, k, 0, DX);
                    if (tempgrad.lpNorm<2>() * dr  > t_refinement_criterion)
                        tags(i, j, k) = amrex::TagBox::SET;
                });
            }
        } 
         
    }
    void Flame::Regrid(int lev, Set::Scalar /* time */)
    {
        BL_PROFILE("Integrator::Flame::Regrid");
        //if (lev < finest_level) return;
        phi_mf[lev]->setVal(0.0);
        ic_phi->Initialize(lev, phi_mf);
        Util::Message(INFO, "Regridding on level ", lev);
    } 
} // namespace Integrator


