
#include <eigen3/Eigen/Eigenvalues>

#include <omp.h>
#include <cmath>

#include <AMReX_SPACE.H>

#include "PhaseFieldMicrostructure.H"
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
    MechanicsBase<model_type>::Advance(lev,time,dt);
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
                    Set::Vector normal = Deta / normgrad;
                    Set::Matrix4<AMREX_SPACEDIM, Set::Sym::Full> DDDDEta = Numeric::DoubleHessian<AMREX_SPACEDIM>(eta, i, j, k, m, DX);

#if AMREX_SPACEDIM == 1
                    Util::Abort(INFO, "Anisotropy is enabled but works in 2D/3D ONLY");
#elif AMREX_SPACEDIM == 2
                    Set::Vector tangent(normal[1],-normal[0]);
                    Set::Scalar Theta = atan2(Deta(1),Deta(0));
                    Set::Scalar kappa = pf.l_gb*0.75*boundary->W(Theta);
                    Set::Scalar Dkappa = pf.l_gb*0.75*boundary->DW(Theta);
                    Set::Scalar DDkappa = pf.l_gb*0.75*boundary->DDW(Theta);
                    mu = 0.75 * (1.0/0.23) * boundary->W(Theta) / pf.l_gb;
                    Set::Scalar sinTheta = sin(Theta);
                    Set::Scalar cosTheta = cos(Theta);
            
                    Set::Scalar Curvature_term =
                        DDDDEta(0,0,0,0)*(    sinTheta*sinTheta*sinTheta*sinTheta) +
                        DDDDEta(0,0,0,1)*(4.0*sinTheta*sinTheta*sinTheta*cosTheta) +
                        DDDDEta(0,0,1,1)*(6.0*sinTheta*sinTheta*cosTheta*cosTheta) +
                        DDDDEta(0,1,1,1)*(4.0*sinTheta*cosTheta*cosTheta*cosTheta) +
                        DDDDEta(1,1,1,1)*(    cosTheta*cosTheta*cosTheta*cosTheta);

                    Set::Scalar Boundary_term =
                        kappa*laplacian +
                        Dkappa*(cos(2.0*Theta)*DDeta(0,1) + 0.5*sin(2.0*Theta)*(DDeta(1,1) - DDeta(0,0)))
                        + 0.5*DDkappa*(sinTheta*sinTheta*DDeta(0,0) - 2.*sinTheta*cosTheta*DDeta(0,1) + cosTheta*cosTheta*DDeta(1,1));
                    if (std::isnan(Boundary_term)) Util::Abort(INFO,"nan at m=",i,",",j,",",k);
            
                    driving_force += - (Boundary_term) + anisotropy.beta*(Curvature_term);
                    if (std::isnan(driving_force)) Util::Abort(INFO,"nan at m=",i,",",j,",",k);

                    #elif AMREX_SPACEDIM == 3
                    // GRAHM-SCHMIDT PROCESS 
                    const Set::Vector e1(1,0,0), e2(0,1,0), e3(0,0,1);
                    Set::Vector _t2, _t3;
                    if      (fabs(normal(0)) > fabs(normal(1)) && fabs(normal(0)) > fabs(normal(2)))
                    {
                        _t2 = e2 - normal.dot(e2)*normal; _t2 /= _t2.lpNorm<2>();
                        _t3 = e3 - normal.dot(e3)*normal - _t2.dot(e3)*_t2; _t3 /= _t3.lpNorm<2>();
                    }
                    else if (fabs(normal(1)) > fabs(normal(0)) && fabs(normal(1)) > fabs(normal(2)))
                    {
                        _t2 = e1 - normal.dot(e1)*normal; _t2 /= _t2.lpNorm<2>();
                        _t3 = e3 - normal.dot(e3)*normal - _t2.dot(e3)*_t2; _t3 /= _t3.lpNorm<2>();
                    }
                    else
                    {
                        _t2 = e1 - normal.dot(e1)*normal; _t2 /= _t2.lpNorm<2>();
                        _t3 = e2 - normal.dot(e2)*normal - _t2.dot(e2)*_t2; _t3 /= _t3.lpNorm<2>();
                    }                            
                    // Compute Hessian projected into tangent space (spanned by _t1,_t2)
                    Eigen::Matrix2d DDeta2D;
                    DDeta2D <<
                        _t2.dot(DDeta*_t2) , _t2.dot(DDeta*_t3),
                        _t3.dot(DDeta*_t2) , _t3.dot(DDeta*_t3);
                    Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eigensolver(2);
                    eigensolver.computeDirect(DDeta2D);
                    Eigen::Matrix2d eigenvecs = eigensolver.eigenvectors();

                    // Compute tangent vectors embedded in R^3
                    Set::Vector t2 = _t2*eigenvecs(0,0) + _t3*eigenvecs(0,1),
                        t3 = _t2*eigenvecs(1,0) + _t3*eigenvecs(1,1);

                    // Compute components of second Hessian in t2,t3 directions
                    Set::Scalar DH2 = 0.0, DH3 = 0.0;
                    Set::Scalar DH23 = 0.0;
                    for (int p = 0; p < 3; p++)
                        for (int q = 0; q < 3; q++)
                            for (int r = 0; r < 3; r++)
                                for (int s = 0; s < 3; s++)
                                {
                                    DH2 += DDDDEta(p,q,r,s)*t2(p)*t2(q)*t2(r)*t2(s);
                                    DH3 += DDDDEta(p,q,r,s)*t3(p)*t3(q)*t3(r)*t3(s);
                                    DH23 += DDDDEta(p,q,r,s)*t2(p)*t2(q)*t3(r)*t3(s);
                                }

                    Set::Scalar gbe = gbmodel.W(normal);
                    //Set::Scalar kappa = l_gb*0.75*gbe;
                    kappa = pf.l_gb*0.75*gbe;
                    mu = 0.75 * (1.0/0.23) * gbe / pf.l_gb;
                    Set::Scalar DDK2 = gbmodel.DDW(normal,_t2) * pf.l_gb * 0.75;
                    Set::Scalar DDK3 = gbmodel.DDW(normal,_t3) * pf.l_gb * 0.75;

                    // GB energy anisotropy term
                    Set::Scalar gbenergy_df = - kappa*laplacian - DDK2*DDeta2D(0,0) - DDK3*DDeta2D(1,1);
                    driving_force += gbenergy_df;
                                  
                    // Second order curvature term
                    Set::Scalar reg_df = NAN;
                    switch(regularization)
                    {
                    case Wilmore:
                        reg_df = anisotropy.beta*(DH2 + DH3 + 2.0*DH23);
                        break;
                    case K12:
                        reg_df = anisotropy.beta*(DH2+DH3);
                        break;
                    }
                    driving_force += reg_df;

                    if (std::isnan(driving_force) || std::isinf(driving_force))
                    {
                        for (int p = 0; p < 3; p++)
                            for (int q = 0; q < 3; q++)
                                for (int r = 0; r < 3; r++)
                                    for (int s = 0; s < 3; s++)
                                    {
                                        Util::Message(INFO,p,q,r,s," ",DDDDEta(p,q,r,s));
                                    }
                        Util::Abort(INFO,"nan/inf detected at amrlev = ", lev," i=",i," j=",j," k=",k);
                    }
                    #endif
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
                if (std::isnan(driving_force))
                    Util::Abort(INFO, i, " ", j, " ", k, " ", m);
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
    MechanicsBase<model_type>::Initialize(lev);
    ic->Initialize(lev, eta_new_mf);
    ic->Initialize(lev, eta_old_mf);
}

void PhaseFieldMicrostructure::TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, Set::Scalar time, int ngrow)
{
    BL_PROFILE("PhaseFieldMicrostructure::TagCellsForRefinement");
    MechanicsBase<model_type>::TagCellsForRefinement(lev,a_tags,time,ngrow);
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
    MechanicsBase<model_type>::TimeStepBegin(time,iter);

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
    MechanicsBase<model_type>::Integrate(amrlev,time,step,mfi,box);

    Model::Interface::GB::SH gbmodel(0.0, 0.0, anisotropy.sigma0, anisotropy.sigma1);
    const amrex::Real *DX = geom[amrlev].CellSize();
    Set::Scalar dv = AMREX_D_TERM(DX[0], *DX[1], *DX[2]);

    amrex::Array4<amrex::Real> const &eta = (*eta_new_mf[amrlev]).array(mfi);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {

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

                                        Set::Matrix DDeta = Numeric::Hessian(eta, i, j, k, 0, DX);
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
