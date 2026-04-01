#include <AMReX_MLPoisson.H>

#include "BC/Constant.H"
#include "DoublyDegenerateCahnHilliard.H"
#include "IC/Ellipse.H"
#include "IC/Random.H"
#include "IO/ParmParse.H"
#include "Numeric/Stencil.H"
#include "Set/Set.H"

namespace Integrator
{
DoublyDegenerateCahnHilliard::DoublyDegenerateCahnHilliard() : Integrator()
{
}
DoublyDegenerateCahnHilliard::~DoublyDegenerateCahnHilliard()
{
    delete ic;
    delete bc;
}

void DoublyDegenerateCahnHilliard::Parse(DoublyDegenerateCahnHilliard &value, IO::ParmParse &pp)
{
    // Interface energy
    pp.query_default("epsilon", value.epsilon, 0.0005);
    // Chemical potential parameter
    pp.query_default("omega", value.omega, 72.0);
    // Energy restriction parameter
    pp.query_default("gamma", value.gamma, 6.0);
    // Energy restriction exponent
    // pp.query_default("p", value.p, 1.0);
    // Energy restriction normalization parameter
    pp.query_default("alpha", value.alpha, 1e-4);
    // Mobility parameter
    pp.query_default("mu", value.mu, 36.0);
    // Regridding criterion
    pp.query_default("refinement_threshold", value.refinement_threshold, 1e100);

    // initial condition for :math:`\eta`
    pp.select_default<IC::Ellipse, IC::Random>("eta.ic", value.ic, value.geom);
    // boundary condition for :math:`\eta`
    pp.select_default<BC::Constant>("eta.bc", value.bc, 1);

    value.RegisterNewFab(value.etanew_mf, value.bc, 1, 1, "eta", true);
    value.RegisterNewFab(value.etaold_mf, value.bc, 1, 1, "eta_old", false);
    value.RegisterNewFab(value.driving_force_mf, value.bc, 1, 1, "driving_Force", false, false);
}

void
DoublyDegenerateCahnHilliard::Initialize(int lev)
{
    driving_force_mf[lev]->setVal(0.0);
    ic->Initialize(lev, etanew_mf);
    ic->Initialize(lev, etaold_mf);
}

void
DoublyDegenerateCahnHilliard::Advance(int lev, Set::Scalar /*time*/, Set::Scalar dt)
{
    std::swap(etaold_mf[lev], etanew_mf[lev]);
    const Set::Scalar *DX = geom[lev].CellSize();
    for (amrex::MFIter mfi(*etanew_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.tilebox();
        Set::Patch<const Set::Scalar> &etaold = etaold_mf[lev]->array(mfi);
        Set::Patch<Set::Scalar> &driving_force = driving_force_mf[lev]->array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Set::Scalar eta = etaold(i, j, k);
            Set::Scalar g_alpha = 1.0 / sqrt(gamma * gamma * (eta * eta * (1 - eta) * (1 - eta)) + alpha * alpha * epsilon * epsilon);
            Set::Scalar g_alphaprime_denominator = sqrt(eta * eta * (1 - eta) * (1 - eta) + alpha * alpha * epsilon * epsilon);
            Set::Scalar g_alphaprime = eta * (1 - eta) * (2 * eta - 1) / (gamma * g_alphaprime_denominator * g_alphaprime_denominator * g_alphaprime_denominator);
            Set::Scalar f = omega / 4 * eta * eta * (1 - eta) * (1 - eta);
            Set::Scalar fprime = omega / 2 * eta * (1 - eta) * (1 - 2 * eta);
            Set::Vector grad_eta = Numeric::Gradient(etaold, i, j, k, 0, DX);
            Set::Scalar norm_squared_grad_eta = grad_eta.squaredNorm();
            Set::Scalar lap_eta = Numeric::Laplacian(etaold, i, j, k, 0, DX);

            driving_force(i, j, k) = g_alpha * fprime / epsilon - epsilon * (g_alphaprime * norm_squared_grad_eta + g_alpha * lap_eta) + g_alphaprime * (1 / epsilon * f + epsilon / 2 * norm_squared_grad_eta);
        });
    }

    driving_force_mf[lev]->FillBoundary(geom[lev].periodicity());

    for (amrex::MFIter mfi(*etanew_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.tilebox();
        Set::Patch<Set::Scalar> &etanew = etanew_mf[lev]->array(mfi);
        Set::Patch<const Set::Scalar> &etaold = etaold_mf[lev]->array(mfi);
        Set::Patch<const Set::Scalar> &driving_force = driving_force_mf[lev]->array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Set::Scalar eta = etaold(i, j, k);
            Set::Vector grad_eta = Numeric::Gradient(etaold, i, j, k, 0, DX);
            Set::Scalar M_alpha = mu * eta * eta * (1 - eta) * (1 - eta) + alpha * epsilon;
            Set::Scalar M_alphaprime = 2.0 * mu * eta * (eta - 1) * (2 * eta - 1);
            Set::Vector grad_driving_force = Numeric::Gradient(driving_force, i, j, k, 0, DX);
            Set::Scalar lap_driving_force = Numeric::Laplacian(driving_force, i, j, k, 0, DX);

            etanew(i, j, k) = 1.0 / epsilon * (M_alphaprime * grad_eta.dot(grad_driving_force) + M_alpha * lap_driving_force) * dt + etaold(i, j, k);
        });
    }
}

void
DoublyDegenerateCahnHilliard::TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, Set::Scalar /*time*/, int /*ngrow*/)
{
    const Set::Scalar *DX = geom[lev].CellSize();
    Set::Scalar dr = sqrt(AMREX_D_TERM(DX[0] * DX[0], +DX[1] * DX[1], +DX[2] * DX[2]));

    for (amrex::MFIter mfi(*etanew_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.tilebox();
        amrex::Array4<char> const &tags = a_tags.array(mfi);
        Set::Patch<const Set::Scalar> eta = (*etanew_mf[lev]).array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Set::Vector grad = Numeric::Gradient(eta, i, j, k, 0, DX);
            if (grad.lpNorm<2>() * dr > refinement_threshold)
                tags(i, j, k) = amrex::TagBox::SET;
        });
    }
}
}
