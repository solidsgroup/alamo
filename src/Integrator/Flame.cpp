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
            // These are the phase field method parameters
            // that you use to inform the phase field method.
            pp.query("pf.eps", value.pf.eps); // Burn width thickness
            pp.query("pf.kappa", value.pf.kappa); // Interface energy param
            pp.query("pf.gamma",value.pf.gamma); // Scaling factor for mobility
            pp.query("pf.lambda",value.pf.lambda); // Chemical potential multiplier
            pp.query("pf.w1", value.pf.w1); // Unburned rest energy
            pp.query("pf.w12", value.pf.w12);  // Barrier energy
            pp.query("pf.w0", value.pf.w0);    // Burned rest energy

            value.bc_eta = new BC::Constant(1);
            pp.queryclass("pf.eta.bc", *static_cast<BC::Constant *>(value.bc_eta)); // See :ref:`BC::Constant`

            value.RegisterNewFab(value.eta_mf,     value.bc_eta, 1, 1, "eta", true);
            value.RegisterNewFab(value.eta_old_mf, value.bc_eta, 1, 1, "eta_old", false);
            value.RegisterNewFab(value.mdot_mf,    value.bc_eta, 1, 1, "mdot", true);

            pp.query("pf.i_check", value.pf.i_check);
            pp.query("pf.j_check", value.pf.j_check);
            pp.query("pf.k_check", value.pf.k_check);
          
        }

	    {
            pp.query("pressure.P", value.pressure.P);
	        pp.query("pressure.a_ap", value.pressure.a_ap);
	        pp.query("pressure.a_htpb", value.pressure.a_htpb);
	        pp.query("pressure.a_comb", value.pressure.a_comb);
	        pp.query("pressure.b_ap", value.pressure.b_ap);
	        pp.query("pressure.b_htpb", value.pressure.b_htpb);
	        pp.query("pressure.b_comb", value.pressure.b_comb);
	    }

        {
            //IO::ParmParse pp("thermal");
            pp.query("thermal.on",value.thermal.on); // Whether to use the Thermal Transport Model
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
            
            pp.query("thermal.bound", value.thermal.bound);
            pp.query("thermal.m_ap", value.thermal.m_ap);
            pp.query("thermal.m_htpb", value.thermal.m_htpb);
            pp.query("thermal.m_comb", value.thermal.m_comb);
            pp.query("thermal.E_ap", value.thermal.E_ap);
            pp.query("thermal.E_htpb", value.thermal.E_htpb);
            pp.query("thermal.E_comb", value.thermal.E_comb);
            pp.query("thermal.correction_factor", value.thermal.correction_factor);
            pp.query("thermal.temperature_delay", value.thermal.temperature_delay); 

            if (value.thermal.on)
            {
                value.bc_temp = new BC::Constant(1);
                pp.queryclass("thermal.temp.bc", *static_cast<BC::Constant *>(value.bc_temp));
                value.RegisterNewFab(value.temp_mf, value.bc_temp, 1, 1, "temp", true);
                value.RegisterNewFab(value.temp_old_mf, value.bc_temp, 1, 1, "temp_old", false);
                value.RegisterNewFab(value.mob_mf, value.bc_temp, 1, 1, "mob", true);
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
            pp.queryclass("model_htpb",value.elastic.model_htpb);
            pp.queryclass("model_void",value.elastic.model_void);
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
            mob_mf[lev]->setVal(0.0);
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
                amrex::Array4<Set::Scalar> const  &mob = (*mob_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const &mdot = (*mdot_mf[lev]).array(mfi);


                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Scalar eta_lap = Numeric::Laplacian(eta, i, j, k, 0, DX);

                    // =============== TODO ==================
                    // This part is probably ok for now. Howver I want to look into
                    // splitting the mobility into two sections so that curvature
                    // is not temperature-dependent. 

                    etanew(i, j, k) = eta(i, j, k) - mob(i,j,k) * dt * ( (pf.lambda/pf.eps) * dw( eta(i,j,k) ) - pf.eps * pf.kappa * eta_lap );
                    //if (etanew(i,j,k) < 0.001) etanew(i,j,k) = 0.0;
                  
                    if (etanew(i,j,k) != etanew(i,j,k)){
                    Util::ParallelMessage(INFO, "eta: ", eta(i,j,k));
                    Util::ParallelMessage(INFO, "mob: ", mob(i,j,k));
                    Util::ParallelMessage(INFO, "alpha: ", alpha(i,j,k));
                    //Util::ParallelMessage(INFO, "dw: ", dw);
                    Util::ParallelMessage(INFO, "lamda/eps: ", (pf.lambda/pf.eps));
                    Util::ParallelAbort(INFO, "eta", etanew(i,j,k) == etanew(i,j,k) );

                    }
                    // Calculate effective thermal conductivity
                    // No special interface mixure rule is needed here.
                    Set::Scalar K   = thermal.k_ap   * phi(i,j,k) + thermal.k_htpb   * (1.0 - phi(i,j,k));
                    Set::Scalar rho = thermal.rho_ap * phi(i,j,k) + thermal.rho_htpb * (1.0 - phi(i,j,k));
                    Set::Scalar cp  = thermal.cp_ap  * phi(i,j,k) + thermal.cp_htpb  * (1.0 - phi(i,j,k));
                    //K = (K - thermal.k0) * eta(i,j,k) + thermal.k0;

                    // Calculate thermal diffusivity and store in field
                    alpha(i,j,k) = eta(i,j,k) * K / rho / cp;
                    thermal.alp = K / rho / cp; 

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

                    // =============== TODO ==================  DONE
                    // We need to get heat flux from mass flux HERE
                    // This is a primitive preliminary implementation.
                    // Note: "thermal.q0" is an initiation heat flux - think of it
                    // like a laser that is heating up the interface. 

                    Set::Scalar qdot = 0.0; // Set to work with SI Units. Pressure should be in MPa. qdot is in units of W/m^2 
                    qdot += (pressure.P * pressure.a_ap   + pressure.b_ap)   * phi(i,j,k) / thermal.k_ap; // AP Portion 
                    qdot += (pressure.P * pressure.a_htpb + pressure.b_htpb) * (1.0 - phi(i,j,k)) / thermal.k_htpb; // HTPB Portion 
                    qdot += (pressure.P * pressure.a_comb + pressure.b_comb) * 4.0 * phi(i,j,k) * (1.0 - phi(i,j,k)) / thermal.k_comb; // AP/HTPB Portion
                    qdot += thermal.q0; // initiation heat flux - think of it like a laser that is heating up the interface.

                    //
                    // Evolve temperature with the qdot flux term in place
                    //
                    // Calculate modified spatial derivative of temperature
                    Set::Scalar dTdt = 0.0;
                    
                    dTdt += grad_eta.dot(alpha(i,j,k) * grad_temp) / (eta(i,j,k) + small);                    
                    dTdt += grad_alpha.dot(grad_temp);
                    dTdt += alpha(i,j,k) * lap_temp;                            
                    // Calculate the source term
                    dTdt += grad_eta_mag * alpha(i,j,k) * qdot / thermal.correction_factor / (eta(i,j,k) + small);
                    // Now, evolve temperature with explicit forward Euler
                    //dTdt = thermal.alp * (lap_temp + grad_eta.dot(grad_temp) + grad_eta_mag * qdot / (eta(i,j,k)+small) / thermal.correction_factor);
                    tempnew(i,j,k) = temp(i,j,k) + dt * dTdt;
  
                    // =============== TODO ================== DONE 
                    // We need to more accurately calculate our effective mobility here.
                    // Right now it uses a simple three-parameter exponential model, however,
                    // it does not take material heterogeneity into account. This model
                    // should account for the differing mobilities for AP and HTPB
                    // (but not necessarily for the combination region).
                    //if (tempnew(i,j,k) < 1800.0){
                    thermal.exp_val = -1.0 * thermal.E_ap / tempnew(i,j,k);
                    //}
                    //else{
                    //thermal.exp_val = -1.0 * thermal.E_ap / 1800.0;
                    //}
		            mob(i,j,k)  = (small + thermal.m_ap * exp(thermal.exp_val) ) * phi(i,j,k);

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
        if (lev < finest_level) return;
        phi_mf[lev]->setVal(0.0);
        ic_phi->Initialize(lev, phi_mf);
        Util::Message(INFO, "Regridding on level ", lev);
    } 
} // namespace Integrator


