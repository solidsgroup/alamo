#include "Flame.H"
#include "IO/ParmParse.H"
#include "BC/Constant.H"
#include "Numeric/Stencil.H"
#include "IC/Laminate.H"
#include "IC/Constant.H"
#include "IC/PSRead.H"
#include "IC/Expression.H"

namespace Integrator
{

    Flame::Flame() : Integrator()
    {
    }

    void Flame::Initialize(int lev)
    {
        if (thermal.on) Temp_mf[lev]->setVal(0.0);
        if (thermal.on) Temp_old_mf[lev]->setVal(0.0);

        EtaIC->Initialize(lev,Eta_mf);
        EtaIC->Initialize(lev,Eta_old_mf);

        PhiIC->Initialize(lev, phi_mf);
    }

    void Flame::UpdateModel(int /*a_step*/) 
    {
        if (m_type == Base::Mechanics<Model::Solid::Affine::Isotropic>::Type::Disable) return;

        for (int lev = 0; lev <= finest_level; ++lev)
        {
            Util::RealFillBoundary(*Eta_mf[lev], geom[lev]);
            Util::RealFillBoundary(*phi_mf[lev], geom[lev]);
            if (thermal.on) Util::RealFillBoundary(*Temp_mf[lev], geom[lev]);
            rhs_mf[lev]->setVal(Set::Vector::Zero());
            rhs_mf[lev]->setVal(Set::Vector::Zero());
            disp_mf[lev]->setVal(Set::Vector::Zero());
            model_mf[lev]->setVal(mechanics.model_ap);
            Set::Vector DX(geom[lev].CellSize());
            for (MFIter mfi(*model_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.nodaltilebox();
                bx.grow(1);
                amrex::Array4<model_type> const &model = model_mf[lev]->array(mfi);
                //amrex::Array4<const Set::Scalar> const &eta = Eta_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const &phi = phi_mf[lev]->array(mfi);
                if (thermal.on)
                {
                    amrex::Array4<const Set::Scalar> const &temp = Temp_mf[lev]->array(mfi);
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        Set::Scalar phi_avg = Numeric::Interpolate::CellToNodeAverage(phi,i,j,k,0);
                        Set::Scalar temp_avg = Numeric::Interpolate::CellToNodeAverage(temp,i,j,k,0);
                        model_type model_ap = mechanics.model_ap;
                        model_ap.F0 *= temp_avg;
                        model_type model_htpb = mechanics.model_htpb;
                        model_htpb.F0 *= temp_avg;
                        model(i,j,k) = model_ap*phi_avg + model_htpb*(1.-phi_avg);
                    });
                }
                else
                {
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        Set::Scalar phi_avg = Numeric::Interpolate::CellToNodeAverage(phi,i,j,k,0);
                        model_type model_ap = mechanics.model_ap;
                        model_ap.F0 *= Set::Matrix::Zero();
                        model_type model_htpb = mechanics.model_htpb;
                        model_htpb.F0 *= Set::Matrix::Zero();
                        model(i,j,k) = model_ap*phi_avg + model_htpb*(1.-phi_avg);
                    });
                }
            }

            Util::RealFillBoundary(*model_mf[lev], geom[lev]);
            
            amrex::MultiFab::Copy(*psi_mf[lev],*Eta_mf[lev],0,0,1,psi_mf[lev]->nGrow());
        }
    }

    void Flame::TimeStepBegin(Set::Scalar a_time, int a_iter)
    {
        Base::Mechanics<Model::Solid::Affine::Isotropic>::TimeStepBegin(a_time,a_iter);
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

                    Set::Scalar eta_lap = Numeric::Laplacian(Eta_old, i, j, k, 0, DX);

                    Eta(i, j, k) = 
                        Eta_old(i, j, k) 
                        - L * dt * (
                            (pf.lambda/pf.eps)*(a1 + 2.0 * a2 * Eta_old(i, j, k) + 3.0 * a3 * Eta_old(i, j, k) * Eta_old(i, j, k) + 4 * a4 * Eta_old(i, j, k) * Eta_old(i, j, k) * Eta_old(i, j, k))
                            - pf.eps * pf.kappa * eta_lap);
                    
                    if (Eta(i,j,k) > Eta_old(i,j,k)) Eta(i,j,k) = Eta_old(i,j,k);
                });
            }
        }

        //
        // Temperature evolution
        //

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

                    amrex::Real rho = (thermal.rho1 - thermal.rho0) * Eta_old(i,j,k) + thermal.rho0;
                    amrex::Real Ka = (thermal.ka - thermal.k0) * Eta_old(i,j,k) + thermal.k0;
                    amrex::Real Kh = (thermal.kh -thermal.k0) * Eta_old(i,j,k) + thermal.k0;
                    Set:: Scalar K = Ka*phi(i,j,k) + Kh*(1-phi(i,j,k));

                    amrex::Real cp = (thermal.cp1 - thermal.cp0) * Eta_old(i,j,k) + thermal.cp0;

                    Set::Scalar test = normvec.dot(temp_grad);
                    Set:: Scalar neumbound = thermal.delA*phi(i,j,k) + thermal.delH*(1-phi(i,j,k));

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

    void Flame::Integrate(int amrlev, Set::Scalar /*time*/, int /*step*/,
                                            const amrex::MFIter &mfi, const amrex::Box &box)
    {
        BL_PROFILE("Flame::Integrate");
        const Set::Scalar *DX = geom[amrlev].CellSize();
        Set::Scalar dv = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);
        amrex::Array4<amrex::Real> const &eta = (*Eta_mf[amrlev]).array(mfi);
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
