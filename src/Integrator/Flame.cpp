#include "Flame.H"
#include "IO/ParmParse.H"
#include "BC/Constant.H"
#include "Numeric/Stencil.H"
#include "IC/Laminate.H"
#include "IC/Constant.H"
#include "IC/PSRead.H"
#include "IC/Expression.H"
#include "Base/Mechanics.H"

#include <cmath>

namespace Integrator
{

    Flame::Flame() : Base::Mechanics<Model::Solid::Affine::Isotropic>() {}

    Flame::Flame(IO::ParmParse &pp) : Flame()
    {pp.queryclass(*this);}

    void 
    Flame::Parse(Flame &value, IO::ParmParse &pp)
    {
        BL_PROFILE("Integrator::Flame::Flame()");
        {
            pp.query("pf.eps", value.pf.eps); // Burn width thickness
            pp.query("pf.kappa", value.pf.kappa); // Interface energy param
            pp.query("pf.gamma",value.pf.gamma); // Scaling factor for mobility
            pp.query("pf.lambda",value.pf.lambda); // Chemical potential multiplier
            pp.query("pf.w1", value.pf.w1); // Unburned rest energy
            pp.query("pf.w12", value.pf.w12);  // Barrier energy
            pp.query("pf.w0", value.pf.w0);    // Burned rest energy
            pp.query("pf.P",      value.pf.P);             // Pressure [UNITS?]
            pp.query("pf.r_ap",   value.pf.r_ap);       // AP Power law multiplier
            pp.query("pf.n_ap",   value.pf.n_ap);       // AP Power law exponent
            pp.query("pf.r_htpb", value.pf.r_htpb);   // HTPB Power law multiplier
            pp.query("pf.n_htpb", value.pf.n_htpb);   // HTPB Power law exponent
            pp.query("pf.r_comb", value.pf.r_comb);   // Combination power law multiplier
            pp.query("pf.n_comb", value.pf.n_comb);   // Combination power law exponent

            value.bc_eta = new BC::Constant(1);
            pp.queryclass("pf.eta.bc", *static_cast<BC::Constant *>(value.bc_eta)); // See :ref:`BC::Constant`

            value.RegisterNewFab(value.eta_mf, value.bc_eta, 1, 2, "eta", true);
            value.RegisterNewFab(value.eta_old_mf, value.bc_eta, 1, 2, "eta_old", false);
        }

        {
            // Read in parameters to determine the IC for eta
            std::string type = "constant";
            pp.query("eta.ic.type", type); // IC type - [packedspheres,laminate] - see classes for more information
            if (type == "constant") value.ic_eta = new IC::Constant(value.geom,pp,"eta.ic.constant");
            else if (type == "expression") value.ic_eta  = new IC::Expression(value.geom,pp,"eta.ic.expression");
            else Util::Abort(INFO,"Invalid eta.ic ", type);
        }

        {
            // These parameters are for the **Thermal transport model**
            pp.query("thermal.on",   value.thermal.on);       // Whether to use the thermal model
            pp.query("thermal.rho1", value.thermal.rho1); // Density (before)
            pp.query("thermal.rho0", value.thermal.rho0); // Density (after)
            pp.query("thermal.ka",   value.thermal.ka); // Thermal conductivity (before and after)
            pp.query("thermal.kh",   value.thermal.kh); // Thermal conductivity (before and after)
            pp.query("thermal.k0",   value.thermal.k0); // Thermal conductivity (before and after)
            pp.query("thermal.cp1",  value.thermal.cp1); // Specific heat (before and after)
            pp.query("thermal.cp0",  value.thermal.cp0); // Specific heat (before and after)
            pp.query("thermal.delA", value.thermal.delA); // Thermal flux of each material
            pp.query("thermal.delH", value.thermal.delH); // Thermal flux of each material
            if (value.thermal.on)
            {
                value.bc_temp = new BC::Constant(1,pp,"thermal.bc");
                value.RegisterNewFab(value.temp_mf, value.bc_temp, 1, 2, "temp", true);
                value.RegisterNewFab(value.temp_old_mf, value.bc_temp, 1, 2, "temp_old", false);
            }
        }


        // Refinement criterion for eta field
        pp.query("amr.refinement_criterion", value.m_refinement_criterion); 
        // Refinement criterion for temperature field
        pp.query("amr.refinement_criterion_temp", value.t_refinement_criterion); 
        // Eta value to restrict the refinament for the temperature field
        pp.query("amr.refinament_restriction", value.t_refinement_restriction);



        {
            // The material field is referred to as :math:`\phi(\mathbf{x})` and is 
            // specified using these parameters. 
            //IO::ParmParse pp("phi.ic");
            std::string type = "packedspheres";
            pp.query("phi.ic.type", type); // IC type (psread, laminate, constant)
            if      (type == "psread") value.ic_phi = new IC::PSRead(value.geom,pp,"phi.ic.psread");
            else if (type == "laminate") value.ic_phi = new IC::Laminate(value.geom,pp,"phi.ic.laminate");
            else if (type == "constant") value.ic_phi = new IC::Constant(value.geom,pp,"phi.ic.constant");
            else Util::Abort(INFO,"Invalid IC type ",type);
            
            value.RegisterNewFab(value.phi_mf, value.bc_eta, 1, 2, "phi", true);
        }

        value.m_type = Base::Mechanics<model_type>::Type::Disable;
        pp.queryclass<Base::Mechanics<Model::Solid::Affine::Isotropic>>("elastic",value);
        if (value.m_type  != Type::Disable)
        {
            pp.queryclass("model_ap",value.elastic.model_ap);
            pp.queryclass("model_htpb",value.elastic.model_htpb);
        }

        value.RegisterIntegratedVariable(&value.volume,"volume");
        value.RegisterIntegratedVariable(&value.area,"area");
    }

    void Flame::Initialize(int lev)
    {
        BL_PROFILE("Integrator::Flame::Initialize");
        Base::Mechanics<Model::Solid::Affine::Isotropic>::Initialize(lev);

        if (thermal.on)
        {
            temp_mf[lev]->setVal(0.0);
            temp_old_mf[lev]->setVal(0.0);
        } 

        ic_eta->Initialize(lev,eta_mf);
        ic_eta->Initialize(lev,eta_old_mf);

        ic_phi->Initialize(lev, phi_mf);
        
    }
    
    void Flame::UpdateModel(int /*a_step*/)
    {
        if (m_type == Base::Mechanics<Model::Solid::Affine::Isotropic>::Type::Disable) return;

        for (int lev = 0; lev <= finest_level; ++lev)
        {
            phi_mf[lev]->FillBoundary();
            eta_mf[lev]->FillBoundary();
            temp_mf[lev]->FillBoundary();
            
            for (MFIter mfi(*model_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.nodaltilebox();
                amrex::Array4<model_type>        const &model = model_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const &phi = phi_mf[lev]->array(mfi);
                if (thermal.on)
                {
                    amrex::Array4<const Set::Scalar> const &temp = temp_mf[lev]->array(mfi);
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        Set::Scalar phi_avg = Numeric::Interpolate::CellToNodeAverage(phi,i,j,k,0);
                        Set::Scalar temp_avg = Numeric::Interpolate::CellToNodeAverage(temp,i,j,k,0);
                        model_type model_ap = elastic.model_ap;
                        model_ap.F0 *= temp_avg;
                        model_type model_htpb = elastic.model_htpb;
                        model_htpb.F0 *= temp_avg;
                        model(i,j,k) = model_ap*phi_avg + model_htpb*(1.-phi_avg);
                    });
                }
                else
                {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        Set::Scalar phi_avg = Numeric::Interpolate::CellToNodeAverage(phi,i,j,k,0);
                        model_type model_ap = elastic.model_ap;
                        model_ap.F0 *= Set::Matrix::Zero();
                        model_type model_htpb = elastic.model_htpb;
                        model_htpb.F0 *= Set::Matrix::Zero();
                        model(i,j,k) = model_ap*phi_avg + model_htpb*(1.-phi_avg);
                    });
                }
            }

            Util::RealFillBoundary(*model_mf[lev], geom[lev]);
            
            amrex::MultiFab::Copy(*psi_mf[lev],*eta_mf[lev],0,0,1,psi_mf[lev]->nGrow());
        }
    }
    
    void Flame::TimeStepBegin(Set::Scalar a_time, int a_iter)
    {
        BL_PROFILE("Integrator::Flame::TimeStepBegin");
        Base::Mechanics<Model::Solid::Affine::Isotropic>::TimeStepBegin(a_time,a_iter);
    }


    void Flame::Advance(int lev, Set::Scalar time, Set::Scalar dt)
    {
        const amrex::Real *DX = geom[lev].CellSize();


        //
        // Phase field evolution
        //

        if (lev == finest_level)
        {
            std::swap(eta_old_mf[lev], eta_mf[lev]);
            
            Set::Scalar 
                a0 = pf.w0, 
                a1 = 0.0, 
                a2 = -5.0 * pf.w1 + 16.0 * pf.w12 - 11.0 * a0, 
                a3 = 14.0 * pf.w1 - 32.0 * pf.w12 + 18.0 * a0, 
                a4 = -8.0 * pf.w1 + 16.0 * pf.w12 -  8.0 * a0;

            for (amrex::MFIter mfi(*eta_mf[lev], true); mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.tilebox();

                amrex::Array4<Set::Scalar> const &eta = (*eta_mf[lev]).array(mfi);
                amrex::Array4<const Set::Scalar> const &eta_old = (*eta_old_mf[lev]).array(mfi);
                amrex::Array4<const Set::Scalar> const &phi = (*phi_mf[lev]).array(mfi);

                Set::Scalar fmod_ap   = pf.r_ap * pow(pf.P, pf.n_ap);
                Set::Scalar fmod_htpb = pf.r_htpb * pow(pf.P, pf.n_htpb);
                Set::Scalar fmod_comb = pf.r_comb * pow(pf.P, pf.n_comb);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {

                    Set:: Scalar fs_actual;

                    fs_actual = 
                            fmod_ap * phi(i, j, k) 
                            + fmod_htpb * (1.0 - phi(i, j, k))
                            + 4.0*fmod_comb*phi(i,j,k)*(1.0-phi(i,j,k));
                            
                    Set::Scalar L = fs_actual / pf.gamma / (pf.w1 - pf.w0);

                    Set::Scalar eta_lap = Numeric::Laplacian(eta_old, i, j, k, 0, DX);

                    eta(i, j, k) = 
                        eta_old(i, j, k) 
                        - L * dt * (
                            (pf.lambda/pf.eps)*(a1 + 2.0 * a2 * eta_old(i, j, k) + 3.0 * a3 * eta_old(i, j, k) * eta_old(i, j, k) + 4 * a4 * eta_old(i, j, k) * eta_old(i, j, k) * eta_old(i, j, k))
                            - pf.eps * pf.kappa * eta_lap);
                    
                    if (eta(i,j,k) > eta_old(i,j,k)) eta(i,j,k) = eta_old(i,j,k);
                });
            }
        }

        //
        // Temperature evolution
        //

        Set::Scalar temperature_delay = 0.01; // hard coded for now, need to make input
        if (thermal.on && time >= temperature_delay)
        {
            std::swap(temp_old_mf[lev], temp_mf[lev]);
            for (amrex::MFIter mfi(*temp_mf[lev], true); mfi.isValid(); ++mfi)
            {
                const amrex::Box &bx = mfi.tilebox();

                amrex::Array4<const Set::Scalar> const &eta_old = (*eta_old_mf[lev]).array(mfi);
                amrex::Array4<Set::Scalar> const &temp = (*temp_mf[lev]).array(mfi);
                amrex::Array4<const Set::Scalar> const &temp_old = (*temp_old_mf[lev]).array(mfi);
                amrex::Array4<const Set::Scalar> const &phi = (*phi_mf[lev]).array(mfi);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Vector eta_grad = Numeric::Gradient(eta_old, i, j, k, 0, DX);
                    Set::Vector temp_grad = Numeric::Gradient(temp_old, i, j, k, 0, DX);
                    Set::Scalar temp_lap = Numeric::Laplacian(temp_old, i, j, k, 0, DX);
                    Set::Scalar eta_grad_mag = eta_grad.lpNorm<2>();
                    Set::Vector normvec = eta_grad/eta_old(i,j,k);

                    amrex::Real rho = (thermal.rho1 - thermal.rho0) * eta_old(i,j,k) + thermal.rho0;
                    amrex::Real Ka = (thermal.ka - thermal.k0) * eta_old(i,j,k) + thermal.k0;
                    amrex::Real Kh = (thermal.kh -thermal.k0) * eta_old(i,j,k) + thermal.k0;
                    Set:: Scalar K = Ka*phi(i,j,k) + Kh*(1-phi(i,j,k));

                    amrex::Real cp = (thermal.cp1 - thermal.cp0) * eta_old(i,j,k) + thermal.cp0;

                    Set::Scalar test = normvec.dot(temp_grad);
                    Set:: Scalar neumbound = thermal.delA*phi(i,j,k) + thermal.delH*(1-phi(i,j,k));

                    if (eta_old(i,j,k) > 0.001 && eta_old(i,j,k)<1)
                    {
                        temp(i,j,k) = temp_old(i,j,k) + dt*(K/cp/rho) * (test + temp_lap + eta_grad_mag/eta_old(i,j,k) * neumbound);
                    }
                    else
                    {
                        if(eta_old(i,j,k)<=0.001)
                        {
                            temp(i,j,k)= 0;
                        }
                    else
                    {
                        temp(i,j,k) = temp_old(i,j,k)+ dt*(K/cp/rho) * temp_lap;
                    }
                } });
            }
        }
    }


    void Flame::TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, Set::Scalar time, int ngrow)
    {
        BL_PROFILE("Integrator::Flame::TagCellsForRefinement");
        Base::Mechanics<Model::Solid::Affine::Isotropic>::TagCellsForRefinement(lev,a_tags,time,ngrow);
        
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
                amrex::Array4<const Set::Scalar> const &eta  = (*eta_mf[lev]).array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Vector tempgrad = Numeric::Gradient(temp, i, j, k, 0, DX);
                    if (tempgrad.lpNorm<2>() * dr  > t_refinement_criterion && eta(i,j,k) > t_refinement_restriction)

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

    void Flame::Integrate(int amrlev, Set::Scalar /*time*/, int /*step*/,
                                            const amrex::MFIter &mfi, const amrex::Box &box)
    {
        BL_PROFILE("Flame::Integrate");
        const Set::Scalar *DX = geom[amrlev].CellSize();
        Set::Scalar dv = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);
        amrex::Array4<amrex::Real> const &eta = (*eta_mf[amrlev]).array(mfi);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
        {
            volume += eta(i, j, k, 0) * dv;
            Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX);
            Set::Scalar normgrad = grad.lpNorm<2>();
            Set::Scalar da = normgrad * dv;
            area += da;
        });
    }

} // namespace Integrator
