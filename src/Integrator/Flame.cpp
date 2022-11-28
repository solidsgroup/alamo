#include "Flame.H"
#include "IO/ParmParse.H"
#include "BC/Constant.H"
#include "Numeric/Stencil.H"
#include "IC/Laminate.H"
#include "IC/Constant.H"
#include "IC/PSRead.H"
#include "Numeric/Function.H"

#include <cmath>

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
            pp.query("timestep", value.base_time);
            // These are the phase field method parameters
            // that you use to inform the phase field method.
            pp.query("pf.eps", value.pf.eps); // Burn width thickness
            pp.query("pf.kappa", value.pf.kappa); // Interface energy param
            pp.query("pf.gamma",value.pf.gamma); // Scaling factor for mobility
            pp.query("pf.lambda",value.pf.lambda); // Chemical potential multiplier
            pp.query("pf.w1", value.pf.w1); // Unburned rest energy
            pp.query("pf.w12", value.pf.w12);  // Barrier energy
            pp.query("pf.w0", value.pf.w0);    // Burned rest energy
            pp.query("pf.min_eta", value.pf.min_eta);
                
            pp.query("pf.time_control", value.pf.time_control);
            
            value.bc_eta = new BC::Constant(1);
            pp.queryclass("pf.eta.bc", *static_cast<BC::Constant *>(value.bc_eta)); // See :ref:`BC::Constant`
            value.RegisterNewFab(value.eta_mf,     value.bc_eta, 1, 1, "eta", true);
            value.RegisterNewFab(value.eta_old_mf, value.bc_eta, 1, 1, "eta_old", false);
            value.RegisterNewFab(value.mdot_mf, 1, "mdot", true);

            std::string eta_bc_str = "constant";
            pp.query("pf.eta.ic.type",eta_bc_str);
            if (eta_bc_str == "constant") value.ic_eta = new IC::Constant(value.geom,pp,"pf.eta.ic.constant");
            else if (eta_bc_str == "expression") value.ic_eta = new IC::Expression(value.geom,pp,"pf.eta.ic.expression");

            std::string eta_ic_type = "constant";
            pp.query("eta.ic.type", eta_ic_type);
            if (eta_ic_type == "laminate") value.ic_eta = new IC::Laminate(value.geom, pp, "eta.ic.laminate");
            else if (eta_ic_type == "constant") value.ic_eta = new IC::Constant(value.geom, pp, "eta.ic.consntant");
            else Util::Abort(INFO, "Invalid eta IC type", eta_ic_type);

        }

        {
            pp.query("pressure.P", value.pressure.P);
            pp.query("pressure.a1", value.pressure.a1);
            pp.query("pressure.a2", value.pressure.a2);
            pp.query("pressure.a3", value.pressure.a3);
            pp.query("pressure.b1", value.pressure.b1);
            pp.query("pressure.b2", value.pressure.b2);
            pp.query("pressure.b3", value.pressure.b3);
            pp.query("pressure.c1", value.pressure.c1);
            pp.query("pressure.E1", value.pressure.E1);
            pp.query("pressure.E2", value.pressure.E2);
        
        }

        {
            pp.query("mass.on", value.mass.on);          
            pp.query("mass.ref_htpb", value.mass.ref_htpb);
            pp.query("mass.a_ap", value.mass.a_ap);
            pp.query("mass.b_ap", value.mass.b_ap);

        }

        {// IO::ParmParse pp(conditional);
          pp.query("conditional.evolve", value.conditional.evolve);
        }

        
        {
            //IO::ParmParse pp("thermal");
            pp.query("thermal.on",value.thermal.on); // Whether to use the Thermal Transport Model

            pp.query("thermal.rho_ap",value.thermal.rho_ap); // AP Density
            pp.query("thermal.rho_htpb", value.thermal.rho_htpb); // HTPB Density
            pp.query("thermal.k_ap", value.thermal.k_ap); // AP Thermal Conductivity
            pp.query("thermal.k_htpb",value.thermal.k_htpb); // HTPB Thermal Conductivity
	    pp.query("thermal.k0", value.thermal.k0); //Numerical Correction variable
            pp.query("thermal.cp_ap", value.thermal.cp_ap); // AP Specific Heat
            pp.query("thermal.cp_htpb", value.thermal.cp_htpb); //HTPB Specific Heat

            pp.query("thermal.q0",value.thermal.q0); // Baseline heat flux
	    pp.query("thermal.q1", value.thermal.q1);
            pp.query("thermal.q_htpb", value.thermal.q_htpb);
            pp.query("thermal.q_ap", value.thermal.q_ap);     

            pp.query("thermal.bound", value.thermal.bound); // System Initial Temperature

            pp.query("thermal.m_ap", value.thermal.m_ap); // AP Pre-exponential factor for Arrhenius Law
            pp.query("thermal.m_htpb", value.thermal.m_htpb); // HTPB Pre-exponential factor for Arrhenius Law
            pp.query("thermal.m_comb", value.thermal.m_comb); // AP/HTPB Pre-exponential factor for Arrhenius Law
            pp.query("thermal.E_ap", value.thermal.E_ap); // AP Activation Energy for Arrhenius Law
            pp.query("thermal.E_htpb", value.thermal.E_htpb); // HTPB Activation Energy for Arrhenius Law
            pp.query("thermal.E_comb", value.thermal.E_comb); // AP/HTPB Activation Energy for Arrhenius Law

            pp.query("thermal.bd", value.thermal.bd);

            pp.query("thermal.r_ap",value.thermal.r_ap);
            pp.query("thermal.r_htpb", value.thermal.r_htpb);
            pp.query("thermal.r_comb", value.thermal.r_comb);
            pp.query("thermal.n_ap", value.thermal.n_ap);

            pp.query("thermal.cut_off", value.thermal.cut_off);
            pp.query("thermal.hc", value.thermal.hc);
	    pp.query("thermal.w1", value.thermal.w1);
	    pp.query("thermal.w2", value.thermal.w2);
            
            pp.query("thermal.qlimit", value.thermal.qlimit);
	    pp.query("thermal.mlocal_ap", value.thermal.mlocal_ap);
	    pp.query("thermal.mlocal_htpb", value.thermal.mlocal_htpb);

	    pp.query("thermal.T_fluid", value.thermal.T_fluid);
	    
	    pp.query("thermal.Wd", value.Wd);
	    pp.query("thermal.mobcap", value.thermal.mobcap);
            value.bc_temp = new BC::Constant(1);
            pp.queryclass("thermal.temp.bc", *static_cast<BC::Constant *>(value.bc_temp));
            value.RegisterNewFab(value.temp_mf,     value.bc_temp, 1, 2, "temp", true);
            value.RegisterNewFab(value.temp_old_mf, value.bc_temp, 1, 2, "temp_old", false);
            value.RegisterNewFab(value.temps_mf, value.bc_temp, 1, 2, "temps", true);
	    value.RegisterNewFab(value.temps_old_mf, value.bc_temp, 1, 2, "temps_old", false);
	    
            value.RegisterNewFab(value.mob_mf, 1, "mob", true);
            value.RegisterNewFab(value.alpha_mf, 1,"alpha",true);  
            value.RegisterNewFab(value.heatflux_mf, 1, "heatflux", true);
        }


        // Refinement criterion for eta field
        pp.query("amr.refinement_criterion", value.m_refinement_criterion); 
        // Refinement criterion for temperature field
        pp.query("amr.refinement_criterion_temp", value.t_refinement_criterion); 
        // Eta value to restrict the refinament for the temperature field
        pp.query("amr.refinament_restriction", value.t_refinement_restriction);
        // small value
        pp.query("small", value.small);

        
        
        {
            // The material field is referred to as :math:`\phi(\mathbf{x})` and is 
            // specified using these parameters. 
            //IO::ParmParse pp("phi.ic");
            std::string type = "packedspheres";
            pp.query("phi.ic.type", type); // IC type (psread, laminate, constant)
            if (type == "psread"){
                value.ic_phi = new IC::PSRead(value.geom, pp, "phi.ic.psread");
                pp.query("phi.ic.psread.eps", value.zeta);
                pp.query("phi.zeta_0", value.zeta_0);
                }
            else if (type == "laminate"){
                value.ic_phi = new IC::Laminate(value.geom,pp,"phi.ic.laminate");
                pp.query("phi.ic.laminate.eps", value.zeta);
                pp.query("phi.zeta_0", value.zeta_0);
                }
            else if (type == "constant")  value.ic_phi = new IC::Constant(value.geom,pp,"phi.ic.constant");
            else Util::Abort(INFO,"Invalid IC type ",type);
            
            value.RegisterNewFab(value.phi_mf, value.bc_eta, 1, 1, "phi", true);
        }

        pp.queryclass<MechanicsBase<Model::Solid::Affine::Isotropic>>("elastic",value);
        if (value.m_type  != Type::Disable)
        {
            pp.queryclass("model_ap",value.elastic.model_ap);
            pp.queryclass("model_htpb",value.elastic.model_htpb);
            pp.queryclass("model_void",value.elastic.model_void);
        }
    }

    void Flame::Initialize(int lev)
    {
        BL_PROFILE("Integrator::Flame::Initialize");
        Util::Message(INFO,m_type);
        MechanicsBase<Model::Solid::Affine::Isotropic>::Initialize(lev);

        temp_mf[lev]->setVal(thermal.bound);
        temp_old_mf[lev]->setVal(thermal.bound);
	temps_mf[lev] ->setVal(thermal.bound);
	temps_old_mf[lev] -> setVal(thermal.bound);

	
        alpha_mf[lev]->setVal(0.0);
        mob_mf[lev] -> setVal(0.0);

        //eta_mf[lev]->setVal(1.0);
        //eta_old_mf[lev]->setVal(1.0);
        ic_eta->Initialize(lev,eta_mf);
        ic_eta->Initialize(lev,eta_old_mf);

        mdot_mf[lev]->setVal(0.0);

        heatflux_mf[lev] -> setVal(0.0);

        ic_phi->Initialize(lev, phi_mf);
	thermal.T_fluid = thermal.bound;

	thermal.w1 = 0.2 * pressure.P + 0.9;

	//for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi){
	//const amrex::Box &bx = mfi.tilebox();
	//amrex::Array4<Set::Scalar> const &temp = (*temp_mf[lev]).array(mfi);
	//amrex::Array4<Set::Scalar> const &phi = (*eta_mf[lev]).array(mfi);
	//amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
	//{
	//  temp(i,j,k) = thermal.bound * (phi(i,j,k));    
	//});}
	
    }

    
    void Flame::UpdateModel(int a_step)
    {
        if (a_step % m_interval) return;

        for (int lev = 0; lev <= finest_level; ++lev)
        {
            phi_mf[lev]->FillBoundary();
            eta_mf[lev]->FillBoundary();
            temp_mf[lev]->FillBoundary();
            
            for (MFIter mfi(*model_mf[lev], true); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.nodaltilebox();
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
        //const Set::Scalar small = 1E-8;

        if (lev == finest_level) //(true)
        {
            std::swap(eta_old_mf[lev], eta_mf[lev]);
            std::swap(temp_old_mf[lev], temp_mf[lev]);
	    std::swap(temps_old_mf[lev], temps_mf[lev]);

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
                //amrex::Array4<Set::Scalar>       const &temph   = (*temph_mf[lev]).array(mfi);
		amrex::Array4<Set::Scalar>       const &tempsnew= (*temps_mf[lev]).array(mfi);
		amrex::Array4<Set::Scalar>       const &temps   = (*temps_old_mf[lev]).array(mfi);

                // Diagnostic fields
                amrex::Array4<Set::Scalar> const  &mob = (*mob_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const &mdot = (*mdot_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const &heatflux = (*heatflux_mf[lev]).array(mfi);


		// Constants
		Set::Scalar k1 = pressure.a1 * pressure.P + pressure.b1 - zeta_0 / zeta;
		Set::Scalar k2 = pressure.a2 * pressure.P + pressure.b2 - zeta_0 / zeta;
		Set::Scalar k3 = 4.0 * log((pressure.c1 * pressure.P * pressure.P + pressure.a3 * pressure.P + pressure.b3) - k1 / 2.0 - k2 / 2.0);
                Wn = 1.0 - Wd;
		
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar eta_lap = Numeric::Laplacian(eta, i, j, k, 0, DX);  
                    Set::Scalar K   = thermal.k0 * ( thermal.k_ap   * phi(i,j,k) + thermal.k_htpb   * (1.0 - phi(i,j,k)) ); // Calculate effective thermal conductivity
                    Set::Scalar rho = thermal.rho_ap * phi(i,j,k) + thermal.rho_htpb * (1.0 - phi(i,j,k)); // No special interface mixure rule is needed here.
                    Set::Scalar cp  = thermal.cp_ap  * phi(i,j,k) + thermal.cp_htpb  * (1.0 - phi(i,j,k));
            
                    etanew(i, j, k) = eta(i, j, k) - mob(i,j,k) * dt * ( (pf.lambda/pf.eps) * dw( eta(i,j,k) ) - pf.eps * pf.kappa * eta_lap );
		    if (eta(i,j,k) < 0.1) etanew(i,j,k) = 0.0;
                   
                    alpha(i,j,k) = thermal.w1 *  K / rho / cp; // Calculate thermal diffusivity and store in fiel
		});
		
		amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
	       {
		 Set::Scalar density = thermal.rho_ap * phi(i,j,k) + thermal.rho_htpb * (1.0 - phi(i,j,k));
		 mdot(i,j,k) = - density * (etanew(i,j,k) - eta(i,j,k)) / dt ;
		});

		amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
		{
		    Set::Scalar K = thermal.k0 * (thermal.k_ap * phi(i,j,k) + thermal.k_htpb * (1.0 - phi(i,j,k)));
		    Set::Scalar qflux = k1 * phi(i,j,k) +
		                        k2 * (1.0 - phi(i,j,k)) +
		                        (zeta_0 / zeta) * exp(k3 * phi(i,j,k) * (1.0 - phi(i,j,k)));
		    
		    Set::Scalar mlocal = thermal.mlocal_ap * phi(i,j,k) + thermal.mlocal_htpb * (1.0 - phi(i,j,k));
		    Set::Scalar mbase;

		    if (mdot(i,j,k) > mlocal) mbase = 1.0;
		    else mbase = mdot(i,j,k) / mlocal; 
		     
		    heatflux(i,j,k) = (thermal.hc * mbase * qflux + thermal.q0 ) / K;
		});
                     
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
		    auto sten = Numeric::GetStencil(i,j,k,bx);
                    Set::Vector grad_eta = Numeric::Gradient(eta, i, j, k, 0, DX);
                    Set::Vector grad_temp = Numeric::Gradient(temp, i, j, k, 0, DX);
                    Set::Scalar lap_temp = Numeric::Laplacian(temp, i, j, k, 0, DX);
                    Set::Scalar grad_eta_mag = grad_eta.lpNorm<2>();
                    Set::Vector grad_alpha = Numeric::Gradient(alpha,i,j,k,0,DX, sten);
                   		    
                    Set::Scalar dTdt = 0.0;

		    dTdt += grad_eta.dot(grad_temp * alpha(i,j,k));
		    dTdt += grad_alpha.dot(eta(i,j,k) * grad_temp);
		    dTdt += eta(i,j,k) * alpha(i,j,k) * lap_temp;

		    dTdt += alpha(i,j,k) * heatflux(i,j,k) * grad_eta_mag;

		    Set::Scalar Tsolid;

		    Tsolid = dTdt + temps(i,j,k) * (etanew(i,j,k) - eta(i,j,k)) / dt; 

		    tempsnew(i,j,k) = temps(i,j,k) + dt * Tsolid;

		    tempnew(i,j,k) = etanew(i,j,k) * tempsnew(i,j,k) + (1.0 - etanew(i,j,k)) * thermal.T_fluid;

                });


                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k){
                    mob(i,j,k) = thermal.m_ap * pressure.P * exp(-thermal.E_ap / tempnew(i,j,k) ) * phi(i,j,k);
		    mob(i,j,k) += thermal.m_htpb * exp(-thermal.E_htpb / tempnew(i,j,k) ) * (1.0 - phi(i,j,k));
		    mob(i,j,k) *= conditional.evolve;

		    if (mob(i,j,k) > thermal.mobcap) mob(i,j,k) = thermal.mobcap;
                });
            } // MFi For loop
        }// For loop
    } //Function


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
            amrex::Array4<const Set::Scalar> const &eta = (*eta_mf[lev]).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                               {
                Set::Vector gradeta = Numeric::Gradient(eta, i, j, k, 0, DX);
                if (gradeta.lpNorm<2>() * dr * 2 > m_refinement_criterion)
                    tags(i, j, k) = amrex::TagBox::SET; });
        }

        // Thermal criterion for refinement
               
            
         for (amrex::MFIter mfi(*temp_mf[lev], true); mfi.isValid(); ++mfi)
         {
             const amrex::Box &bx = mfi.tilebox();
             amrex::Array4<char> const &tags = a_tags.array(mfi);
             amrex::Array4<const Set::Scalar> const &temp = (*temp_mf[lev]).array(mfi);
             amrex::Array4<const Set::Scalar> const &eta  = (*eta_mf[lev]).array(mfi);
             amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
             {
                Set::Vector tempgrad = Numeric::Gradient(temp, i, j, k, 0, DX);
                        if (tempgrad.lpNorm<2>() * dr  > t_refinement_criterion && eta(i,j,k) >= 0.1)
                                tags(i, j, k) = amrex::TagBox::SET;
                });
            }
        
         
    }
    void Flame::Regrid(int lev, Set::Scalar /* time */)
    {
        BL_PROFILE("Integrator::Flame::Regrid");
        if (lev < finest_level) return;
        phi_mf[lev]->setVal(0.0);
        ic_phi->Initialize(lev, phi_mf);
        Util::Message(INFO, "Regridding on level ", lev);
    } 
} // namespace Integrator


