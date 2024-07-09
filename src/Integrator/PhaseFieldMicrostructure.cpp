
#include <eigen3/Eigen/Eigenvalues>

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

#include "Model/Solid/Affine/Cubic.H"
#include "Model/Solid/Affine/Hexagonal.H"

namespace Integrator
{
template<class model_type>
void PhaseFieldMicrostructure<model_type>::Advance(int lev, Set::Scalar time, Set::Scalar dt)
{
    BL_PROFILE("PhaseFieldMicrostructure::Advance");
    Base::Mechanics<model_type>::Advance(lev, time, dt);
    /// TODO Make this optional
    //if (lev != max_level) return;
    //std::swap(eta_old_mf[lev], eta_new_mf[lev]);
    const Set::Scalar* DX = this->geom[lev].CellSize();


    Model::Interface::GB::SH gbmodel(0.0, 0.0, anisotropy.sigma0, anisotropy.sigma1);

    for (amrex::MFIter mfi(*eta_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        amrex::Box bx = mfi.tilebox();
        //if (m_type == MechanicsBase<model_type>::Type::Static)
        //bx.grow(number_of_ghost_cells-1);
        //bx = bx & domain;
        Set::Patch<Set::Scalar> eta = eta_mf.Patch(lev, mfi);
        Set::Patch<Set::Scalar> driving_force = driving_force_mf.Patch(lev, mfi);
        Set::Patch<Set::Scalar> driving_force_threshold = driving_force_threshold_mf.Patch(lev, mfi);// = (*driving_force_mf[lev]).array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            for (int m = 0; m < number_of_grains; m++)
            {
                driving_force(i, j, k, m) = 0.0;
                if (pf.threshold.on) driving_force_threshold(i, j, k, m) = 0.0;
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
                    if (pf.threshold.boundary)  driving_force_threshold(i, j, k, m) += -kappa * laplacian;
                    else                        driving_force(i, j, k, m) += -kappa * laplacian;
                }
                else
                {
                    Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Full> DDDDEta = Numeric::DoubleHessian<AMREX_SPACEDIM>(eta, i, j, k, m, DX);
                    auto anisotropic_df = boundary->DrivingForce(Deta, DDeta, DDDDEta);
                    if (pf.threshold.boundary) driving_force_threshold(i, j, k, m) += pf.l_gb * 0.75 * std::get<0>(anisotropic_df);
                    else                       driving_force(i, j, k, m) += pf.l_gb * 0.75 * std::get<0>(anisotropic_df);
                    if (std::isnan(std::get<0>(anisotropic_df))) Util::Abort(INFO);
                    if (pf.threshold.boundary) driving_force_threshold(i, j, k, m) += anisotropy.beta * std::get<1>(anisotropic_df);
                    else                       driving_force(i, j, k, m) += anisotropy.beta * std::get<1>(anisotropic_df);
                    if (std::isnan(std::get<1>(anisotropic_df))) Util::Abort(INFO);
                    mu = 0.75 * (1.0 / 0.23) * boundary->W(Deta) / pf.l_gb;
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
                if (pf.threshold.chempot)
                    driving_force_threshold(i, j, k, m) += mu * (eta(i, j, k, m) * eta(i, j, k, m) - 1.0 + 2.0 * pf.gamma * sum_of_squares) * eta(i, j, k, m);
                else
                    driving_force(i, j, k, m) += mu * (eta(i, j, k, m) * eta(i, j, k, m) - 1.0 + 2.0 * pf.gamma * sum_of_squares) * eta(i, j, k, m);

                //
                // SYNTHETIC DRIVING FORCE
                //
                if (lagrange.on && m == 0 && time > lagrange.tstart)
                {
                    if (pf.threshold.lagrange)
                        driving_force_threshold(i, j, k, m) += lagrange.lambda * (volume - lagrange.vol0);
                    else
                        driving_force(i, j, k, m) += lagrange.lambda * (volume - lagrange.vol0);
                }
            }
        });

        //
        // ELASTIC DRIVING FORCE
        //
        if (pf.elastic_df)
        {
            amrex::Array4<const Set::Matrix> const& sigma = (*this->stress_mf[lev]).array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                for (int m = 0; m < number_of_grains; m++)
                {
                    Set::Scalar etasum = 0.0;
                    Set::Matrix F0avg = Set::Matrix::Zero();

                    for (int n = 0; n < number_of_grains; n++)
                    {
                        etasum += eta(i, j, k, n);
                        F0avg += eta(i, j, k, n) * mechanics.model[n].F0;
                    }

                    Set::Matrix sig = Numeric::Interpolate::NodeToCellAverage(sigma, i, j, k, 0);

                    //Set::Matrix dF0deta = mechanics.model[m].F0;//(etasum * elastic.model[m].F0 - F0avg) / (etasum * etasum);
                    Set::Matrix dF0deta = Set::Matrix::Zero();

                    for (int n = 0; n < number_of_grains; n++)
                    {
                        if (n == m) continue;
                        Set::Scalar normsq = eta(i, j, k, m) * eta(i, j, k, m) + eta(i, j, k, n) * eta(i, j, k, n);
                        dF0deta += (2.0 * eta(i, j, k, m) * eta(i, j, k, n) * eta(i, j, k, n) * (mechanics.model[m].F0 - mechanics.model[n].F0))
                            / normsq / normsq;
                    }

                    Set::Scalar tmpdf = (dF0deta.transpose() * sig).trace();

                    if (pf.threshold.mechanics)
                        driving_force_threshold(i, j, k, m) -= pf.elastic_mult * tmpdf;
                    else
                        driving_force(i, j, k, m) -= pf.elastic_mult * tmpdf;
                }
            });
        }


        //
        // Update eta
        // (if NOT using anisotropic kinetics)
        //

        if (!anisotropic_kinetics.on || time < anisotropic_kinetics.tstart)
        {
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                for (int m = 0; m < number_of_grains; m++)
                {
                    if (pf.threshold.on)
                    {
                        if (driving_force_threshold(i, j, k, m) > pf.threshold.value)
                        {
                            if (pf.threshold.type == ThresholdType::Continuous)
                                eta(i, j, k, m) -= pf.L * dt * (driving_force_threshold(i, j, k, m) - pf.threshold.value);
                            else if (pf.threshold.type == ThresholdType::Chop)
                                eta(i, j, k, m) -= pf.L * dt * (driving_force_threshold(i, j, k, m));
                        }
                        else if (driving_force_threshold(i, j, k, m) < -pf.threshold.value)
                        {
                            if (pf.threshold.type == ThresholdType::Continuous)
                                eta(i, j, k, m) -= pf.L * dt * (driving_force_threshold(i, j, k, m) + pf.threshold.value);
                            else if (pf.threshold.type == ThresholdType::Chop)
                                eta(i, j, k, m) -= pf.L * dt * (driving_force_threshold(i, j, k, m));
                        }
                    }

                    eta(i, j, k, m) -= pf.L * dt * driving_force(i, j, k, m);
                }
            });
        }
    }

    //
    // Update eta
    // (if we ARE using anisotropic kinetics)
    //

    if (anisotropic_kinetics.on && time >= anisotropic_kinetics.tstart)
    {
        for (amrex::MFIter mfi(*eta_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            amrex::Array4<const Set::Scalar> const& eta = (*eta_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& L = (*anisotropic_kinetics.L_mf[lev]).array(mfi);
            amrex::Array4<Set::Scalar> const& threshold = (*anisotropic_kinetics.threshold_mf[lev]).array(mfi);

            for (int m = 0; m < number_of_grains; m++)
            {
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Vector Deta = Numeric::Gradient(eta, i, j, k, m, DX);
                    Set::Scalar theta = atan2(Deta(1), Deta(0));
                    L(i, j, k, m) = (4. / 3.) * anisotropic_kinetics.mobility(theta) / pf.l_gb;
                    threshold(i, j, k, m) = anisotropic_kinetics.threshold(theta);
                });
            }
        }

        for (amrex::MFIter mfi(*eta_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.tilebox();
            Set::Patch<Set::Scalar>       eta = eta_mf.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> driving_force = driving_force_mf.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> driving_force_threshold = driving_force_threshold_mf.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> L = anisotropic_kinetics.L_mf.Patch(lev, mfi);
            Set::Patch<const Set::Scalar> threshold = anisotropic_kinetics.threshold_mf.Patch(lev, mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                for (int m = 0; m < number_of_grains; m++)
                {
                    if (pf.threshold.on)
                    {
                        if (driving_force_threshold(i, j, k, m) > threshold(i, j, k, m))
                        {
                            if (pf.threshold.type == ThresholdType::Continuous)
                                eta(i, j, k, m) -= L(i, j, k, m) * dt * (driving_force_threshold(i, j, k, m) - threshold(i, j, k, m));
                            else if (pf.threshold.type == ThresholdType::Chop)
                                eta(i, j, k, m) -= L(i, j, k, m) * dt * (driving_force_threshold(i, j, k, m));
                        }
                        else if (driving_force_threshold(i, j, k, m) < -pf.threshold.value)
                        {
                            if (pf.threshold.type == ThresholdType::Continuous)
                                eta(i, j, k, m) -= L(i, j, k, m) * dt * (driving_force_threshold(i, j, k, m) + threshold(i, j, k, m));
                            else if (pf.threshold.type == ThresholdType::Chop)
                                eta(i, j, k, m) -= L(i, j, k, m) * dt * (driving_force_threshold(i, j, k, m));
                        }
                    }

                    eta(i, j, k, m) -= L(i, j, k, m) * dt * driving_force(i, j, k, m);
                }
            });
        }
    }
}

template<class model_type>
void PhaseFieldMicrostructure<model_type>::Initialize(int lev)
{
    BL_PROFILE("PhaseFieldMicrostructure::Initialize");
    Base::Mechanics<model_type>::Initialize(lev);
    ic->Initialize(lev, eta_mf);
    ic->Initialize(lev, eta_mf);
}

template<class model_type>
void PhaseFieldMicrostructure<model_type>::TagCellsForRefinement(int lev, amrex::TagBoxArray& a_tags, Set::Scalar time, int ngrow)
{
    BL_PROFILE("PhaseFieldMicrostructure::TagCellsForRefinement");
    Base::Mechanics<model_type>::TagCellsForRefinement(lev, a_tags, time, ngrow);
    const amrex::Real* DX = this->geom[lev].CellSize();
    const Set::Vector dx(DX);
    const Set::Scalar dxnorm = dx.lpNorm<2>();

    for (amrex::MFIter mfi(*eta_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.tilebox();
        amrex::Array4<const amrex::Real> const& etanew = (*eta_mf[lev]).array(mfi);
        amrex::Array4<char> const& tags = a_tags.array(mfi);

        for (int n = 0; n < number_of_grains; n++)
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Set::Vector grad = Numeric::Gradient(etanew, i, j, k, n, DX);

            if (dxnorm * grad.lpNorm<2>() > ref_threshold)
                tags(i, j, k) = amrex::TagBox::SET;
        });
    }
}

template<class model_type>
void PhaseFieldMicrostructure<model_type>::TimeStepComplete(Set::Scalar /*time*/, int /*iter*/)
{
    // TODO: remove this function, it is no longer needed.
}

template<class model_type>
void PhaseFieldMicrostructure<model_type>::UpdateModel(int a_step, Set::Scalar /*a_time*/)
{
    BL_PROFILE("PhaseFieldMicrostructure::UpdateModel");
    if (a_step % this->m_interval) return;

    for (int lev = 0; lev <= this->finest_level; ++lev)
    {
        amrex::Box domain = this->geom[lev].Domain();
        domain.convert(amrex::IntVect::TheNodeVector());

        eta_mf[lev]->FillBoundary();

        for (MFIter mfi(*this->model_mf[lev], false); mfi.isValid(); ++mfi)
        {
            amrex::Box bx = mfi.grownnodaltilebox() & domain;

            amrex::Array4<model_type> const& model = this->model_mf[lev]->array(mfi);
            amrex::Array4<const Set::Scalar> const& eta = eta_mf[lev]->array(mfi);

            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                std::vector<Set::Scalar> etas(number_of_grains);
                for (int n = 0; n < number_of_grains; n++)
                    etas[n] = Numeric::Interpolate::CellToNodeAverage(eta, i, j, k, n);
                model(i, j, k) = model_type::Combine(mechanics.model, etas);
            });
        }

        Util::RealFillBoundary(*this->model_mf[lev], this->geom[lev]);
    }

}

template<class model_type>
void PhaseFieldMicrostructure<model_type>::TimeStepBegin(Set::Scalar time, int iter)
{
    BL_PROFILE("PhaseFieldMicrostructure::TimeStepBegin");
    Base::Mechanics<model_type>::TimeStepBegin(time, iter);

    if (anisotropy.on && time >= anisotropy.tstart)
    {
        this->SetTimestep(anisotropy.timestep);
        if (anisotropy.elastic_int > 0)
            if (iter % anisotropy.elastic_int) return;
    }
}

template<class model_type>
void PhaseFieldMicrostructure<model_type>::Integrate(int amrlev, Set::Scalar time, int step,
    const amrex::MFIter& mfi, const amrex::Box& box)
{
    BL_PROFILE("PhaseFieldMicrostructure::Integrate");
    Base::Mechanics<model_type>::Integrate(amrlev, time, step, mfi, box);

    Model::Interface::GB::SH gbmodel(0.0, 0.0, anisotropy.sigma0, anisotropy.sigma1);
    const amrex::Real* DX = this->geom[amrlev].CellSize();
    Set::Scalar dv = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);

    amrex::Array4<amrex::Real> const& eta = (*eta_mf[amrlev]).array(mfi);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
#if AMREX_SPACEDIM == 2
        auto sten = Numeric::GetStencil(i, j, k, box);
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

template class PhaseFieldMicrostructure<Model::Solid::Affine::Cubic>;
template class PhaseFieldMicrostructure<Model::Solid::Affine::Hexagonal>;


} // namespace Integrator
