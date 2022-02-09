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
        }
        
        /*{
            // done TODO: rename variables and assign density/conductivity/spec.heat/thermal flux/ for AP, and HTPB, and comb.
            //       - stick to same naming convention as above
            //       - these params are defined in Flame.H - can change the variable names as well
            //       - you'll need to work on the Advance method as well, which is where these are used.

            // TODO: find real-world actual values for these properties and put them in an input file.

            // TODO: we are going to have to guess what the heat flux values are. That will be done by
            //       running sims and comparing results.

            // These parameters are for the **Thermal transport model**
            IO::ParmParse pp("thermal");
            pp.query("on",thermal.on);       // Whether to use the thermal model
            pp.query("rho1", thermal.rho1); // Density (before)
            pp.query("rho0", thermal.rho0); // Density (after)
            pp.query("ka", thermal.ka); // Thermal conductivity (before and after)
            pp.query("kh", thermal.kh); // Thermal conductivity (before and after)
            pp.query("k0", thermal.k0); // Thermal conductivity (before and after)
            pp.query("cp1", thermal.cp1); // Specific heat (before and after)
            pp.query("cp0", thermal.cp0); // Specific heat (before and after)
            pp.query("delA", thermal.delA); // Thermal flux of each material
            pp.query("delH", thermal.delH); // Thermal flux of each material
            
            
            if (thermal.on)
            {
                TempBC = new BC::Constant(1);
                pp.queryclass("bc", *static_cast<BC::Constant *>(TempBC)); // See :ref:`BC::Constant`
                RegisterNewFab(Temp_mf, TempBC, 1, 1, "Temp", true);
                RegisterNewFab(Temp_old_mf, TempBC, 1, 1, "Temp_old", false);
            }
        }*/

        {
            IO::ParmParse pp("thermal");
            pp.query("on",thermal.on); // Whether to use the Thermal Transport Model
            Util::Message(INFO,"thermal.on =" , thermal.on);
            pp.query("rho_ap",thermal.rho_ap); // AP Density
            pp.query("rho_htpb", thermal.rho_htpb); // HTPB Density
            pp.query("k_ap", thermal.k_ap); // AP Thermal Conductivity
            pp.query("k_htpb",thermal.k_htpb); // HTPB Thermal Conductivity
            pp.query("k0", thermal.k0); // Thermal conductivity 
            pp.query("cp_ap", thermal.cp_ap); // AP Specific Heat
            pp.query("cp_htpb", thermal.cp_htpb); //HTPB Specific Heat
            pp.query("q_ap", thermal.q_ap); // AP  Thermal Flux
            pp.query("q_htpb", thermal.q_htpb); // HTPB Thermal Flux
            pp.query("thermal.q_interface" , thermal.q_interface); // Interface heat flux

            if (thermal.on)
            {
            TempBC = new BC::Constant(1);
            pp.queryclass("bc", *static_cast<BC::Constant *>(TempBC));
            RegisterNewFab(Temp_mf, TempBC, 1, 1, "Temp", true);
            RegisterNewFab(Temp_old_mf, TempBC, 1, 1, "Temp_old", false);            
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
            if (type == "packedspheres")
            {
                PhiIC = new IC::PackedSpheres(geom);
                pp.queryclass("packedspheres", *static_cast<IC::PackedSpheres *>(PhiIC)); // See :ref:`IC::PackedSpheres`
            }
            else if (type == "psread")
            {
                // TODO: run the code in full 3D. Right now everything is in 2D, let's run in 3D.
                //       This method reads in sphere locations and radii from an input file.
                //       We need to update so that we are working with 3D packings.
                //       In the paper repository, look at the packings called 3d.xyzr.
                //       (Note: in the packedspheres input file you can find which packing file is used by
                //       looking in the input file. Swap out the 2D for 3D case.)
                //       Note: you might want to do this before updating the heat flux model.
                
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
        if (thermal.on) Temp_mf[lev]->setVal(0.0);
        if (thermal.on) Temp_old_mf[lev]->setVal(0.0);

        Eta_mf[lev]->setVal(1.0);
        Eta_old_mf[lev]->setVal(1.0);

        PhiIC->Initialize(lev, phi_mf);
    }

    
    void Flame::TimeStepBegin(Set::Scalar /*a_time*/, int a_iter)
    {
        if (!elastic.interval)
            return;
        if (a_iter % elastic.interval)
            return;

        /*{
          IO::ParmParse pp("elastic");
          pp.query("elastic.tol_rel", elastic.tol_rel);
          pp.query("elastic.tol_abs", elastic.tol_abs);
          int nmodels = 1;
          pp.query("nmodels", nmodels);
          for (int i = 0; i<nmodels; i++)
          {
            std::string name = "model" + std::to_string(i+1);
            MODEL tmp_model;
            pp.queryclass(name.data(),tmp_model);
            elastic.models.push_back(tmp_model);
          }
          //Util::Assert(INFO,TEST(elastic.models.size() > 0));
          
    
        }*/
        // TODO: This is where the elastic solve happens. Right now we are just using hard coded values (bad!!)
        //       We want to use actual values that are read in from an input file. Actual values for elastic moduli,
        //       etc. So here, update this code so that the elastic models are read in appropriately.
        //
        //       For reference, consult TensionTest.H. This integrator does full elastic solves with multiple
        //       materials. You can duplicate the process here. Look in there to see how model params are read in
        //       then do the same here.
        //

        // TODO: we need to use actual parameters for AP and HTPB (look them up and put them in the input file)

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
        const amrex::Real *DX = geom[lev].CellSize();


        //
        // Phase field evolution
        //

        if (lev == finest_level)
        {
            std::swap(Eta_old_mf[lev], Eta_mf[lev]);
            
            Set::Scalar 
                a0 = pf.w0, 
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

                // This is the pressure power law stuff.
                Set::Scalar fmod_ap   = pf.r_ap * pow(pf.P, pf.n_ap);
                Set::Scalar fmod_htpb = pf.r_htpb * pow(pf.P, pf.n_htpb);
                Set::Scalar fmod_comb = pf.r_comb * pow(pf.P, pf.n_comb);
                   
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {

                    Set:: Scalar fs_actual;

                    // TODO: note: "fs_acual" and "fmod" variables are actual flame speeds that get converted to
                    //       mobility L below. The relationship between flame speed and mobility is derived in the
                    //       Kanagarajan manuscript. 
                    fs_actual = 
                            fmod_ap * phi(i, j, k) 
                            + fmod_htpb * (1.0 - phi(i, j, k))
                            + 4.0*fmod_comb*phi(i,j,k)*(1.0-phi(i,j,k));
                            
                    // TODO: This parameter, fs_act, is the effective mobility. Notice that it depends on fs_actual (above)
                    //       which is constructed based on "fmod" variables that are defined above, and are where the
                    //       pressure power law is defined
                    //       We want to REPLACE this so that instead of depending on pressure, L depends on the temperature
                    //       and the material properties.
                    //
                    //       BTW: check the references in the Kanagarajan manuscript on level set modeling (Stewart etc)
                    //       They can point you to existing work that similiarly regresses the interface based on temp.
                    Set::Scalar L = fs_actual / pf.gamma / (pf.w1 - pf.w0);

                    Set::Scalar eta_lap = Numeric::Laplacian(Eta_old, i, j, k, 0, DX);

                    Eta(i, j, k) = 
                        Eta_old(i, j, k) 
                        - L * dt * (
                            (pf.lambda/pf.eps)*(a1 + 2.0 * a2 * Eta_old(i, j, k) + 3.0 * a3 * Eta_old(i, j, k) * Eta_old(i, j, k) + 4 * a4 * Eta_old(i, j, k) * Eta_old(i, j, k) * Eta_old(i, j, k))
                            - pf.eps * pf.kappa * eta_lap);
                });
            }
        }

        //
        // Temperature evolution
        //

        //done TODO: replace all the variables here that you redefined in the .H file and the constructor.
        //       This is where they all get used. 

        Set::Scalar temperature_delay = 0.01; // hard coded for now, need to make input
        if (thermal.on && time >= temperature_delay)
        {
            std::swap(Temp_old_mf[lev], Temp_mf[lev]);
            for (amrex::MFIter mfi(*Temp_mf[lev], true); mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.tilebox();

                amrex::Array4<const Set::Scalar> const &Eta_old = (*Eta_old_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const &Temp = (*Temp_mf[lev]).array(mfi);
                amrex::Array4<const Set::Scalar> const &Temp_old = (*Temp_old_mf[lev]).array(mfi);
                amrex::Array4<const Set::Scalar> const &phi = (*phi_mf[lev]).array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Vector eta_grad = Numeric::Gradient(Eta_old, i, j, k, 0, DX);
                    Set::Vector temp_grad = Numeric::Gradient(Temp_old, i, j, k, 0, DX);
                    Set::Scalar temp_lap = Numeric::Laplacian(Temp_old, i, j, k, 0, DX);
                    Set::Scalar eta_grad_mag = eta_grad.lpNorm<2>();
                    Set::Vector normvec = eta_grad/Eta_old(i,j,k);

                    amrex::Real rho = (thermal.rho_ap - thermal.rho_htpb) * Eta_old(i,j,k) + thermal.rho_htpb;
                    amrex::Real K_ap = (thermal.k_ap - thermal.k0) * Eta_old(i,j,k) + thermal.k0;
                    amrex::Real K_htpb = (thermal.k_htpb -thermal.k0) * Eta_old(i,j,k) + thermal.k0;
                    Set:: Scalar K = K_ap*phi(i,j,k) + K_htpb*(1-phi(i,j,k));

                    amrex::Real cp = (thermal.cp_ap - thermal.cp_htpb) * Eta_old(i,j,k) + thermal.cp_htpb;

                    Set::Scalar test = normvec.dot(temp_grad);

                    // TODO: This is where the heat flux gets calculated. But we want to upgrade this
                    //       to reflect heat fluxes at AP, HTPB, AND IN THE INTERFACE. (Right now it is AP or HTPB only.)
                    //       Update this so that we include the interface term as well as the individual species. 
                    Set:: Scalar neumbound = thermal.q_ap*phi(i,j,k) + thermal.q_htpb*(1-phi(i,j,k));
                    
                    // Set::Scalar neumbound = thermal.q_ap * phi(i,j,k) + thermal.q_htpb * (1 - phi(i,j,k)) + 4.0 * phi(i,j,k) * (1 - phi(i,j,k)) * thermal.q_interface;


                    if (Eta_old(i,j,k) > 0.001 && Eta_old(i,j,k)<1)
                    { 
                        Temp(i,j,k) = Temp_old(i,j,k) + dt*(K/cp/rho) * (test + temp_lap + eta_grad_mag/Eta_old(i,j,k) * neumbound);
                    }
                    else
                    {
                        if(Eta_old(i,j,k)<=0.001)
                        {
                            Temp(i,j,k)= 0;
                        }
                    else
                    {
                        Temp(i,j,k) = Temp_old(i,j,k)+ dt*(K/cp/rho) * temp_lap;
                    }
                } });
            }
        }
    }

    void Flame::TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, amrex::Real /*time*/, int /*ngrow*/)
    {
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
                    if (tempgrad.lpNorm<2>() * dr  > t_refinement_criterion)
                        tags(i, j, k) = amrex::TagBox::SET; });
            }
        }
    }
    void Flame::Regrid(int lev, Set::Scalar /* time */)
    {
        if (lev < finest_level) return;
        phi_mf[lev]->setVal(0.0);
        PhiIC->Initialize(lev, phi_mf);
        Util::Message(INFO, "Regridding on level ", lev);
    } 
} // namespace Integrator
