
#include <eigen3/Eigen/Eigenvalues>

#include <omp.h>
#include <cmath>

#include <AMReX_SPACE.H>

#include "PhaseFieldMicrostructure.H"
#include "Integrator/Base/Mechanics.H"
#include "BC/Constant.H"
#include "BC/Step.H"
#include "Set/Set.H"
#include "Util/Util.H"
#include "IC/Random.H"
#include "IC/Trig.H"
#include "IC/Sphere.H"
#include "IC/Expression.H"
#include "Model/Interface/GB/SH.H"
#include "Numeric/Stencil.H"
#include "Solver/Nonlocal/Linear.H"
#include "Solver/Nonlocal/Newton.H"
#include "IC/Trig.H"
namespace Integrator
{
void PhaseFieldMicrostructure::Advance(int lev, Set::Scalar time, Set::Scalar dt)
{
    BL_PROFILE("PhaseFieldMicrostructure::Advance");
    Base::Mechanics<model_type>::Advance(lev,time,dt);
    /// TODO Make this optional
    //if (lev != max_level) return;
    std::swap(eta_old_mf[lev], eta_new_mf[lev]);
    const Set::Scalar *DX = geom[lev].CellSize();


    Model::Interface::GB::SH gbmodel(0.0, 0.0, anisotropy.sigma0, anisotropy.sigma1);

    for (amrex::MFIter mfi(*eta_new_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        amrex::Box bx = mfi.tilebox();
        //if (m_type == MechanicsBase<model_type>::Type::Static)
        //bx.grow(number_of_ghost_cells-1);
        //bx = bx & domain;
        amrex::Array4<const amrex::Real> const &eta = (*eta_old_mf[lev]).array(mfi);
        amrex::Array4<amrex::Real> const &etanew = (*eta_new_mf[lev]).array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            for (int m = 0; m < number_of_grains; m++)
            {
                Set::Scalar driving_force = 0.0;
                Set::Scalar kappa = NAN, mu = NAN;

                //
                // BOUNDARY TERM and SECOND ORDER REGULARIZATION
                //

                Set::Vector Deta = Numeric::Gradient(eta, i, j, k, m, DX);
                Set::Scalar normgrad = Deta.lpNorm<2>();
                if (normgrad < 1E-4)
                    continue; // This ought to speed things up.

                Set::Matrix DDeta = Numeric::Hessian(eta, i, j, k, m, DX);
                Set::Scalar laplacian = DDeta.trace();

                if (!anisotropy.on || time < anisotropy.tstart)
                {
                    kappa = pf.l_gb * 0.75 * pf.sigma0;
                    mu = 0.75 * (1.0 / 0.23) * pf.sigma0 / pf.l_gb;
                    driving_force += -kappa * laplacian;
                }
                else
                {
                    Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Full> DDDDEta = Numeric::DoubleHessian<AMREX_SPACEDIM>(eta, i, j, k, m, DX);
                    auto anisotropic_df = boundary->DrivingForce(Deta, DDeta, DDDDEta);
                    driving_force += pf.l_gb * 0.75 * std::get<0>(anisotropic_df);
                    if (std::isnan(std::get<0>(anisotropic_df))) Util::Abort(INFO);
                    driving_force += anisotropy.beta * std::get<1>(anisotropic_df);
                    if (std::isnan(std::get<1>(anisotropic_df))) Util::Abort(INFO);
                    mu = 0.75 * (1.0/0.23) * boundary->W(Deta) / pf.l_gb;
                }

                //
                // CHEMICAL POTENTIAL
                //

                Set::Scalar sum_of_squares = 0.;
                for (int n = 0; n < number_of_grains; n++)
                {
                    if (m == n)
                        continue;
                    sum_of_squares += eta(i, j, k, n) * eta(i, j, k, n);
                }
                driving_force += mu * (eta(i, j, k, m) * eta(i, j, k, m) - 1.0 + 2.0 * pf.gamma * sum_of_squares) * eta(i, j, k, m);

                //
                // SYNTHETIC DRIVING FORCE
                //
                if (lagrange.on && m == 0 && time > lagrange.tstart)
                {
                    driving_force += lagrange.lambda * (volume - lagrange.vol0);
                }

                //
                // EVOLVE ETA
                //
                etanew(i, j, k, m) = eta(i, j, k, m) - pf.L * dt * driving_force;
                if (std::isnan(driving_force)) Util::Abort(INFO, "Eta is nan at lev=",lev,", (", i, " ", j, " ", k, ")[",m,"]");
            } 
        });

        //
        // ELASTIC DRIVING FORCE
        //
        if (false) // TODO need to replace
        {
            const amrex::Box &bx = mfi.tilebox();
            amrex::Array4<const Set::Scalar> const &eta = (*eta_old_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const &etanew = (*eta_new_mf[lev]).array(mfi);
            amrex::Array4<const Set::Matrix> const &sigma = (*stress_mf[lev]).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
            {
                for (int m = 0; m < number_of_grains; m++)
                {
                    Set::Scalar driving_force = 0.0;
                    Set::Scalar etasum = 0.0;
                    Set::Matrix F0avg = Set::Matrix::Zero();
                    
                    for (int n = 0; n < number_of_grains; n++)
                    {
                        etasum += eta(i,j,k,n);
                        F0avg  += eta(i,j,k,n) * mechanics.model[n].F0;
                    } 
                    
                    Set::Matrix dF0deta = mechanics.model[m].F0;//(etasum * elastic.model[m].F0 - F0avg) / (etasum * etasum);


                    Set::Matrix sig = Numeric::Interpolate::NodeToCellAverage(sigma,i,j,k,0);

                    Set::Scalar tmpdf = (dF0deta.transpose() * sig).trace();

                    if (tmpdf > pf.elastic_threshold)
                    {
                        driving_force -= pf.elastic_mult * (tmpdf-pf.elastic_threshold);
                    }
                    else if (tmpdf < -pf.elastic_threshold)
                    {
                        driving_force -= pf.elastic_mult * (tmpdf+pf.elastic_threshold);
                    }

                    etanew(i, j, k, m) -= pf.L * dt * driving_force;
                }
            });

        }
    }
}

void PhaseFieldMicrostructure::Initialize(int lev)
{
    BL_PROFILE("PhaseFieldMicrostructure::Initialize");
    Base::Mechanics<model_type>::Initialize(lev);
    ic->Initialize(lev, eta_new_mf);
    ic->Initialize(lev, eta_old_mf);
}

void PhaseFieldMicrostructure::TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, Set::Scalar time, int ngrow)
{
    BL_PROFILE("PhaseFieldMicrostructure::TagCellsForRefinement");
    Base::Mechanics<model_type>::TagCellsForRefinement(lev,a_tags,time,ngrow);
    const amrex::Real *DX = geom[lev].CellSize();
    const Set::Vector dx(DX);
    const Set::Scalar dxnorm = dx.lpNorm<2>();

    for (amrex::MFIter mfi(*eta_new_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.tilebox();
        amrex::Array4<const amrex::Real> const &etanew = (*eta_new_mf[lev]).array(mfi);
        amrex::Array4<char> const &tags = a_tags.array(mfi);

        for (int n = 0; n < number_of_grains; n++)
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                                        Set::Vector grad = Numeric::Gradient(etanew, i, j, k, n, DX);

                                        if (dxnorm * grad.lpNorm<2>() > ref_threshold)
                                            tags(i, j, k) = amrex::TagBox::SET;
                                    });
    }
}

void PhaseFieldMicrostructure::TimeStepComplete(Set::Scalar /*time*/, int /*iter*/)
{
    // TODO: remove this function, it is no longer needed.
}

void PhaseFieldMicrostructure::UpdateModel(int a_step)
{
    BL_PROFILE("PhaseFieldMicrostructure::UpdateModel");
    if (a_step % m_interval) return;

    for (int lev = 0; lev <= finest_level; ++lev)
    {
        amrex::Box domain = geom[lev].Domain();
        domain.convert(amrex::IntVect::TheNodeVector());

        eta_new_mf[lev]->FillBoundary();

        for (MFIter mfi(*model_mf[lev], false); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.grownnodaltilebox() & domain;
            
            amrex::Array4<model_type> const &model = model_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const &eta = eta_new_mf[lev]->array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
            {
                std::vector<Set::Scalar> etas(number_of_grains);
                for (int n = 0; n < number_of_grains; n++) 
                    etas[n] = Numeric::Interpolate::CellToNodeAverage(eta,i,j,k,n);
                model(i, j, k) = model_type::Combine(mechanics.model,etas);
            });
        }

        Util::RealFillBoundary(*model_mf[lev],geom[lev]);
    }

}

void PhaseFieldMicrostructure::TimeStepBegin(Set::Scalar time, int iter)
{
    BL_PROFILE("PhaseFieldMicrostructure::TimeStepBegin");
    Base::Mechanics<model_type>::TimeStepBegin(time,iter);

    if (anisotropy.on && time >= anisotropy.tstart)
    {
        SetTimestep(anisotropy.timestep);
        if (anisotropy.elastic_int > 0) 
            if (iter % anisotropy.elastic_int) return;
    }    
}

void PhaseFieldMicrostructure::Integrate(int amrlev, Set::Scalar time, int step,
                                        const amrex::MFIter &mfi, const amrex::Box &box)
{
    BL_PROFILE("PhaseFieldMicrostructure::Integrate");
    Base::Mechanics<model_type>::Integrate(amrlev,time,step,mfi,box);

    Model::Interface::GB::SH gbmodel(0.0, 0.0, anisotropy.sigma0, anisotropy.sigma1);
    const amrex::Real *DX = geom[amrlev].CellSize();
    Set::Scalar dv = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);

    amrex::Array4<amrex::Real> const &eta = (*eta_new_mf[amrlev]).array(mfi);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        #if AMREX_SPACEDIM == 2
        auto sten = Numeric::GetStencil(i,j,k,box);
        #endif

        volume += eta(i, j, k, 0) * dv;

        Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX);
        Set::Scalar normgrad = grad.lpNorm<2>();

        if (normgrad > 1E-8)
        {
            Set::Vector normal = grad / normgrad;

            Set::Scalar da = normgrad * dv;
            area += da;

            if (!anisotropy.on || time < anisotropy.tstart)
            {
                gbenergy += pf.sigma0 * da;
                Set::Scalar k = 0.75 * pf.sigma0 * pf.l_gb;
                realgbenergy += 0.5 * k * normgrad * normgrad * dv;
                regenergy = 0.0;
            }
            else
            {
            #if AMREX_SPACEDIM == 2
                Set::Scalar theta = atan2(grad(1), grad(0));
                Set::Scalar sigma = boundary->W(theta);
                gbenergy += sigma * da;

                Set::Scalar k = 0.75 * sigma * pf.l_gb;
                realgbenergy += 0.5 * k * normgrad * normgrad * dv;

                Set::Matrix DDeta = Numeric::Hessian(eta, i, j, k, 0, DX, sten);
                Set::Vector tangent(normal[1], -normal[0]);
                Set::Scalar k2 = (DDeta * tangent).dot(tangent);
                regenergy += 0.5 * anisotropy.beta * k2 * k2;
            #elif AMREX_SPACEDIM == 3
                gbenergy += gbmodel.W(normal) * da;
            #endif
            }
        }
    });
}

} // namespace Integrator
