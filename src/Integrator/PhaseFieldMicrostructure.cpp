
#include <eigen3/Eigen/Eigenvalues>

#include <cmath>

#include <AMReX_SPACE.H>

#include "PhaseFieldMicrostructure.H"
#include "Integrator/Base/Mechanics.H"
#include "BC/Constant.H"
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
#include "Model/Solid/Finite/PseudoAffine/Cubic.H"

#include "Util/MPI.H"

namespace Integrator
{
template<class model_type>
void PhaseFieldMicrostructure<model_type>::Advance(int lev, Set::Scalar time, Set::Scalar dt)
{
    BL_PROFILE("PhaseFieldMicrostructure::Advance");
    Base::Mechanics<model_type>::Advance(lev, time, dt);
    const Set::Scalar* DX = this->geom[lev].CellSize();

    std::swap(eta_old_mf[lev], eta_mf[lev]);
    

    Set::Scalar df_max = std::numeric_limits<Set::Scalar>::min();

    Model::Interface::GB::SH gbmodel(0.0, 0.0, anisotropy.sigma0, anisotropy.sigma1);

    for (amrex::MFIter mfi(*eta_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        amrex::Box bx = mfi.tilebox();
        Set::Patch<Set::Scalar> etanew = eta_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> eta    = eta_old_mf.Patch(lev,mfi);

        Set::Scalar *df_max_handle = &df_max;

        Set::Patch<const Set::Matrix> sigma = stress_mf.Patch(lev,mfi); 
        Set::Patch<const Set::Vector> disp  = this->disp_mf.Patch(lev,mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Matrix sig = Set::Matrix::Zero();
            if (pf.elastic_df) sig = Numeric::Interpolate::NodeToCellAverage(sigma, i, j, k, 0);

            for (int m = 0; m < number_of_grains; m++)
            {
                //if (shearcouple.on) etaold(i,j,k,m) = eta(i,j,k,m);

                Set::Scalar driving_force = 0.0;
                Set::Scalar driving_force_threshold = NAN;
                if (pf.threshold.on) driving_force_threshold = 0.0;
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
                    if (pf.threshold.boundary)  driving_force_threshold += -kappa * laplacian;
                    else                        driving_force += -kappa * laplacian;
                }
                else
                {
                    Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Full> DDDDEta = Numeric::DoubleHessian<AMREX_SPACEDIM>(eta, i, j, k, m, DX);
                    auto anisotropic_df = boundary->DrivingForce(Deta, DDeta, DDDDEta);
                    if (pf.threshold.boundary) driving_force_threshold += pf.l_gb * 0.75 * std::get<0>(anisotropic_df);
                    else                       driving_force += pf.l_gb * 0.75 * std::get<0>(anisotropic_df);
                    if (std::isnan(std::get<0>(anisotropic_df))) Util::Abort(INFO);
                    if (pf.threshold.boundary) driving_force_threshold += anisotropy.beta * std::get<1>(anisotropic_df);
                    else                       driving_force += anisotropy.beta * std::get<1>(anisotropic_df);
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
                    driving_force_threshold += mu * (eta(i, j, k, m) * eta(i, j, k, m) - 1.0 + 2.0 * pf.gamma * sum_of_squares) * eta(i, j, k, m);
                else
                    driving_force += mu * (eta(i, j, k, m) * eta(i, j, k, m) - 1.0 + 2.0 * pf.gamma * sum_of_squares) * eta(i, j, k, m);


                //
                // LAGRANGE MULTIPLIER
                //
                if (lagrange.on && m == 0 && time > lagrange.tstart)
                {
                    if (pf.threshold.lagrange)
                        driving_force_threshold += lagrange.lambda * (volume - lagrange.vol0);
                    else
                        driving_force += lagrange.lambda * (volume - lagrange.vol0);
                }

                //
                // SYNTHETIC DRIVING FORCE
                //
                if (sdf.on && time > sdf.tstart)
                {
                    if (pf.threshold.sdf)
                        driving_force_threshold += sdf.val[m](time);
                    else
                        driving_force += sdf.val[m](time);
                }

              
                //
                // ELASTIC DRIVING FORCE
                //
                if (pf.elastic_df)
                {
                    Set::Scalar elastic_df_m = 0.0;
                    
                    if (shearcouple.on)
                    {
                        for (int n = 0; n < number_of_grains; n++)
                        {
                            if (n==m) continue;
                            Set::Scalar etam = eta(i,j,k,m);
                            Set::Scalar etan = eta(i,j,k,n);

                            Set::Scalar sumsq = (etam*etam + etan*etan);
                            sumsq = sumsq * sumsq;

                            //Set::Scalar dgm = 2.0*etan*etan*etam / sumsq;
                            Set::Scalar dgn = 2.0*etan*etam*etam / sumsq;


                            Set::Matrix Fgbn = shearcouple.Fgb[n*number_of_grains + m];

                            if (model_type::kinvar == Model::Solid::KinematicVariable::F)
                            {
                                Set::Matrix F =
                                    Set::Matrix::Identity() +
                                    Numeric::NodeGradientOnCell(disp,i,j,k,DX);
                                Fgbn = F*Fgbn;
                            }

                            elastic_df_m += (sig.transpose() * Fgbn).trace() * dgn;
                        }
                    }
                    else
                    {
                        Set::Matrix F0avg = Set::Matrix::Zero();

                        for (int n = 0; n < number_of_grains; n++)
                            F0avg += eta(i, j, k, n) * mechanics.model[n].F0;

                        Set::Matrix dF0deta = Set::Matrix::Zero();

                        for (int n = 0; n < number_of_grains; n++)
                        {
                            if (n == m) continue;
                            Set::Scalar normsq = eta(i, j, k, m) * eta(i, j, k, m) + eta(i, j, k, n) * eta(i, j, k, n);
                            dF0deta += (2.0 * eta(i, j, k, m) * eta(i, j, k, n) * eta(i, j, k, n) *
                                        (mechanics.model[m].F0 - mechanics.model[n].F0))
                                / normsq / normsq;
                        }

                        elastic_df_m += (dF0deta.transpose() * sig).trace();
                    }


                    if (pf.threshold.mechanics)
                        driving_force_threshold -= pf.elastic_mult * elastic_df_m;
                    else
                        driving_force -= pf.elastic_mult * elastic_df_m;
                }

                //
                // Update eta
                //
                // Final mobility and threshold values:
                Set::Scalar L = NAN, threshold = NAN;
                // (if we are NOT using anisotropic kinetics)
                if (!anisotropic_kinetics.on || time < anisotropic_kinetics.tstart)
                {
                    L = pf.L;
                    threshold = pf.threshold.value;
                }
                // (if we ARE using anisotropic kinetics)
                else
                {
                    Set::Scalar theta = atan2(Deta(1), Deta(0));
                    L = (4. / 3.) * anisotropic_kinetics.mobility(theta) / pf.l_gb;
                    threshold = anisotropic_kinetics.threshold(theta);
                }

                Set::Scalar totaldf = 0.0;
                if (pf.threshold.on)
                {
                    if (driving_force_threshold > threshold)
                    {
                        if (pf.threshold.type == ThresholdType::Continuous)
                            totaldf -= L * (driving_force_threshold - threshold);
                        else if (pf.threshold.type == ThresholdType::Chop)
                            totaldf -= L * (driving_force_threshold);
                    }
                    else if (driving_force_threshold < -threshold)
                    {
                        if (pf.threshold.type == ThresholdType::Continuous)
                            totaldf -= L * (driving_force_threshold + threshold);
                        else if (pf.threshold.type == ThresholdType::Chop)
                            totaldf -= L  * (driving_force_threshold);
                    }
                }
                totaldf -= L * driving_force;

                // Final Eta Update
                etanew(i, j, k, m) = eta(i,j,k,m) +  dt * totaldf;
                *df_max_handle = std::max(df_max, std::fabs(totaldf));
            }
        });
    }

    if (shearcouple.on && time >= mechanics.tstart) UpdateEigenstrain(lev);
    //this->DynamicTimestep_SyncDrivingForce(lev,df_max);
}

template <class model_type>
void PhaseFieldMicrostructure<model_type>::UpdateEigenstrain(int lev)
{
    if (this->m_type == Mechanics<model_type>::Disable) return;
    eta_mf[lev]->FillBoundary();
    eta_old_mf[lev]->FillBoundary();

    amrex::Box domain = this->geom[lev].Domain();
    domain.convert(amrex::IntVect::TheNodeVector());

    for (amrex::MFIter mfi(*this->model_mf[lev], false); mfi.isValid(); ++mfi)
    {
        amrex::Box bx = mfi.grownnodaltilebox() & domain;
        Set::Patch<const Set::Scalar> etaold = eta_old_mf.Patch(lev,mfi);
        Set::Patch<const Set::Scalar> etanew = eta_mf.Patch(lev,mfi);
        Set::Patch<model_type>        model  = model_mf.Patch(lev,mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            for (int m = 0; m < number_of_grains; m++)
                for (int n = 0; n < number_of_grains; n++)
                {
                    if (m==n) continue;
                    Set::Scalar emnew = Numeric::Interpolate::CellToNodeAverage(etanew, i, j, k, m);//, sten);
                    Set::Scalar emold = Numeric::Interpolate::CellToNodeAverage(etaold, i, j, k, m);//, sten);
                    Set::Scalar ennew = Numeric::Interpolate::CellToNodeAverage(etanew, i, j, k, n);//, sten);
                    Set::Scalar enold = Numeric::Interpolate::CellToNodeAverage(etaold, i, j, k, n);//, sten);
                    Set::Scalar gmnew = (emnew * emnew) / (emnew * emnew + ennew * ennew);
                    Set::Scalar gmold = (emold * emold) / (emold * emold + enold * enold);
                    model(i, j, k).F0 -= 0.5 * (gmnew - gmold) * shearcouple.Fgb[m*number_of_grains + n]; 
                }
        });
    }
    
    Util::RealFillBoundary(*this->model_mf[lev],this->geom[lev]);
}

template<class model_type>
void PhaseFieldMicrostructure<model_type>::Initialize(int lev)
{
    BL_PROFILE("PhaseFieldMicrostructure::Initialize");
    Base::Mechanics<model_type>::Initialize(lev);
    ic->Initialize(lev, eta_mf);
    ic->Initialize(lev, eta_old_mf);
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
    this->DynamicTimestep_Update();
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
                model_type mix = model_type::Combine(mechanics.model, etas, mechanics.model_mix_order);
                if (a_step == 0 || !shearcouple.on) // Do a complete model initialization
                    model(i,j,k) = mix ;
                else // Otherwise, we will set the modulus only and leave the eigenstrain alone
                {
                    Set::Matrix F0 = model(i,j,k).F0;
                    model(i,j,k) = mix;
                    model(i,j,k).F0 = F0;
                }
            });
        }

        if (mechanics.model_neumann_boundary)
        {
            Util::AssertException(INFO,TEST(AMREX_SPACEDIM==2),"neumann boundary works in 2d only");
            for (MFIter mfi(*this->model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.grownnodaltilebox() & domain;
                amrex::Array4<model_type> const& model = this->model_mf[lev]->array(mfi);
                const Dim3 lo = amrex::lbound(domain), hi = amrex::ubound(domain);

                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    if      (i==lo.x && j==lo.y)
                        model(i,j,k) = 0.5*(model(i+1,j,k)+model(i,j+1,k));
                    else if (i==lo.x && j==hi.y)
                        model(i,j,k) = 0.5*(model(i+1,j,k)+model(i,j-1,k));
                    else if (i==hi.x && j==lo.y)
                        model(i,j,k) = 0.5*(model(i-1,j,k)+model(i,j+1,k));
                    else if (i==hi.x && j==hi.y)
                        model(i,j,k) = 0.5*(model(i-1,j,k)+model(i,j-1,k));

                    else if (i==lo.x)
                        model(i,j,k) = model(i+1,j,k);
                    else if (i==hi.x)
                        model(i,j,k) = model(i-1,j,k);
                    else if (j==lo.y)
                        model(i,j,k) = model(i,j+1,k);
                    else if (j==hi.y)
                        model(i,j,k) = model(i,j-1,k);
                });
            }
        }

        Util::RealFillBoundary(*this->model_mf[lev], this->geom[lev]);
    }

}

template<class model_type>
void PhaseFieldMicrostructure<model_type>::TimeStepBegin(Set::Scalar time, int iter)
{
    BL_PROFILE("PhaseFieldMicrostructure::TimeStepBegin");
    for (int lev = 0; lev < totaldf_mf.size(); lev++) totaldf_mf[lev]->setVal(0.0);

    // df_limit_max.resize(this->max_level+1);
    // df_limit_min.resize(this->max_level+1);

    //
    // Manual Disconnection Nucleation
    //
    // currently supportedfor 2d only
    //
    #if AMREX_SPACEDIM == 2 
    if (disconnection.on && time > disconnection.tstart && !(iter % disconnection.interval))
    {
        disconnection.sitex.clear();
        disconnection.sitey.clear();
        disconnection.phases.clear();

        int lev = this->max_level;
        const Set::Scalar *DX = this->geom[lev].CellSize();
        Set::Scalar exponent = DX[0] * DX[0] * (this->timestep * disconnection.interval) / disconnection.tau_vol;

        // Determine the nucleation sites in my portion of the mesh
        for (amrex::MFIter mfi(*eta_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box &bx = mfi.tilebox();
            amrex::Array4<amrex::Real> const &eta = (*eta_mf[lev]).array(mfi);
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Set::Scalar E0 = 2.0*disconnection.nucleation_energy;
                E0 /= disconnection.epsilon + 256.0*eta(i,j,k,0)*eta(i,j,k,0)*eta(i,j,k,0)*eta(i,j,k,0)*eta(i,j,k,1)*eta(i,j,k,1)*eta(i,j,k,1)*eta(i,j,k,1);
                Set::Scalar p = std::exp(-E0/(disconnection.K_b*disconnection.temp));
                Set::Scalar P = 1.0 - std::pow(1.0 - p,exponent);
                if (eta(i,j,k,0) < 0 || eta(i,j,k,0) > 1.0 || eta(i,j,k,1) < 0 || eta(i,j,k,1) > 1.0) P = 0.0;
                Set::Scalar q = 0.0;
                q = disconnection.unif_dist(disconnection.rand_num_gen);
                if (q < P)
                {
                    disconnection.sitex.push_back(this->geom[lev].ProbLo()[0] + ((amrex::Real)(i)) * DX[0]);
                    disconnection.sitey.push_back(this->geom[lev].ProbLo()[1] + ((amrex::Real)(j)) * DX[1]);
                    int phase = disconnection.int_dist(disconnection.rand_num_gen);
                    disconnection.phases.push_back(phase);
                } });
        }
        // Sync up all the nucleation sites among processors
        Util::MPI::Allgather(disconnection.sitex);
        Util::MPI::Allgather(disconnection.sitey);
        Util::MPI::Allgather(disconnection.phases);
        Util::Message(INFO, "Nucleating ", disconnection.phases.size(), " disconnections");
        if (disconnection.sitex.size() > 0)
        {
            // Now that we all know the nucleation locations, perform the nucleation
            for (int lev = 0; lev <= this->max_level; lev++)
            {
                amrex::Box domain = this->geom[lev].Domain();
                domain.convert(amrex::IntVect::TheNodeVector());
                const amrex::Real *DX = this->geom[lev].CellSize();
                for (amrex::MFIter mfi(*this->model_mf[lev], TilingIfNotGPU()); mfi.isValid(); ++mfi)
                {
                    const amrex::Box bx = mfi.grownnodaltilebox() & domain;
                    //amrex::Array4<const Set::Scalar> const &etaold = (*eta_old_mf[lev]).array(mfi);
                    amrex::Array4<Set::Scalar> const &eta = (*eta_mf[lev]).array(mfi);
                    amrex::Array4<Set::Scalar> const &disc = (*disc_mf[lev]).array(mfi);
                    amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                    {
                        Set::Vector x;
                        AMREX_D_TERM(
                            x(0) = this->geom[lev].ProbLo()[0] + ((amrex::Real)(i)) * DX[0];,
                            x(1) = this->geom[lev].ProbLo()[1] + ((amrex::Real)(j)) * DX[1];,
                            x(2) = this->geom[lev].ProbLo()[2] + ((amrex::Real)(k)) * DX[2];);
                        disc(i, j, k, 0) = 0.0;
                        for (unsigned int m = 0; m < disconnection.phases.size(); m++)
                        {
                            amrex::Real r_squared = 0;
                            Set::Vector nucleation_site(disconnection.sitex[m], disconnection.sitey[m]);
                            for (int n = 0; n < AMREX_SPACEDIM; n++)
                            {
                                amrex::Real dist = nucleation_site(n) - x(n);
                                r_squared += dist * dist;
                            }
                            amrex::Real bump = exp(-r_squared / disconnection.box_size);
                            disc(i, j, k, 0) += bump;
                            eta(i, j, k, disconnection.phases[m]) = bump * (1 - eta(i, j, k, disconnection.phases[m])) + eta(i, j, k, disconnection.phases[m]);
                            eta(i, j, k, 1 - disconnection.phases[m]) = (1. - bump) * eta(i, j, k, 1 - disconnection.phases[m]);
                        }
                    });
                }
                UpdateEigenstrain(lev);
                Util::RealFillBoundary(*this->model_mf[lev], this->geom[lev]);
            }
        }
    }
    #endif
    
    Base::Mechanics<model_type>::TimeStepBegin(time, iter);

    if (anisotropy.on && time >= anisotropy.tstart)
    {
        Util::AssertException(INFO,TEST(this->dynamictimestep.on == false),
                              "Cannot use anisotropy with dynamic timestep yet");
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
template class PhaseFieldMicrostructure<Model::Solid::Finite::PseudoAffine::Cubic>;


} // namespace Integrator
