#include <AMReX_MLPoisson.H>

#ifdef ALAMO_FFT
#include <AMReX_FFT.H>
#endif

#include "BC/Constant.H"
#include "DoublyDegenerateCahnHilliard.H"
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
    pp.select_default<IC::Random>("eta.ic", value.ic, value.geom);
    // boundary condition for :math:`\eta`
    pp.select_default<BC::Constant>("eta.bc", value.bc, 1);

    value.RegisterNewFab(value.etanew_mf, value.bc, 1, 1, "eta", true);
    value.RegisterNewFab(value.etaold_mf, value.bc, 1, 1, "eta_old", false);
    value.RegisterNewFab(value.grad_etaold_mf, value.bc, AMREX_SPACEDIM, 1, "grad_eta_old", false, true, { AMREX_D_DECL("x", "y", "z") });
    value.RegisterNewFab(value.lap_etaold_mf, value.bc, 1, 1, "lap_eta_old", false);
}

void
DoublyDegenerateCahnHilliard::Initialize(int lev)
{
    grad_etaold_mf[lev]->setVal(0.0);
    lap_etaold_mf[lev]->setVal(0.0);
    ic->Initialize(lev, etanew_mf);
    ic->Initialize(lev, etaold_mf);
}

void
DoublyDegenerateCahnHilliard::Advance(int lev, Set::Scalar time, Set::Scalar dt)
{
    std::swap(etaold_mf[lev], etanew_mf[lev]);
    const Set::Scalar *DX = geom[lev].CellSize();
    for (amrex::MFIter mfi(*etanew_mf[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.tilebox();
        Set::Patch<Set::Scalar> &etanew = etanew_mf[lev]->array(mfi);
        Set::Patch<const Set::Scalar> &etaold = etaold_mf[lev]->array(mfi);
        Set::Patch<Set::Scalar> &grad_etaold = grad_etaold_mf[lev]->array(mfi);
        Set::Patch<Set::Scalar> &lap_etaold = lap_etaold_mf[lev]->array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            // Set::Scalar lap_eta = Numeric::Laplacian(eta, i, j, k, 0, DX);

            // inter(i, j, k) = eta(i, j, k) * eta(i, j, k) * eta(i, j, k)
            //                  - eta(i, j, k)
            //                  - gamma * lap_eta;

            // etanew(i, j, k) = eta(i, j, k) - dt * inter(i, j, k); // Allen Cahn

            Set::Scalar eta = etaold(i, j, k);
            Set::Scalar g_alpha = 1.0 / sqrt(gamma * gamma * (eta * eta * (1 - eta) * (1 - eta)) + alpha * alpha * epsilon * epsilon);
            Set::Vector grad_etaold_matrix = g_alpha * Numeric::Gradient(etaold, i, j, k, DX);
            AMREX_D_TERM(
                grad_etaold(i, j, k, 0) = grad_etaold_matrix[0];,
                grad_etaold(i, j, k, 1) = grad_etaold_matrix[1];,
                grad_etaold(i, j, k, 2) = grad_etaold_matrix[2];
            )
        });

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            lap_etaold(i, j, k) = -epsilon * Numeric::Divergence(grad_etaold, i, j, k, 0, DX);
        });

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            // Set::Scalar lap_inter = Numeric::Laplacian(inter, i, j, k, 0, DX);

            // etanew(i, j, k) = eta(i, j, k) + dt * lap_inter;

            Set::Scalar eta = etaold(i, j, k);
            Set::Scalar g_alpha = 1.0 / sqrt(gamma * gamma * (eta * eta * (1 - eta) * (1 - eta)) + alpha * alpha * epsilon * epsilon);
            Set::Scalar g_alphaprime_denominator = sqrt(eta * eta * (1 - eta) * (1 - eta) + alpha * alpha * epsilon * epsilon);
            Set::Scalar g_alphaprime = eta * (1 - eta) * (2 * eta - 1) / (gamma * g_alphaprime_denominator * g_alphaprime_denominator * g_alphaprime_denominator);
            Set::Scalar f = omega / 4 * eta * eta * (1 - eta) * (1 - eta);
            Set::Scalar fprime = omega / 2 * eta * (1 - eta) * (1 - 2 * eta);
            Set::Scalar M_alpha = mu * eta * eta * (1 - eta) * (1 - eta) + alpha * epsilon;
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
