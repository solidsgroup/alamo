#include <utility>

#include "Agglomeration.H"
#include "Flame.H"

#include "BC/Constant.H"
#include "IC/Constant.H"
#include "IC/Ellipse.H"
#include "IC/Expression.H"
#include "IC/PSRead.H"
#include "IC/Random.H"
#include "IO/ParmParse.H"
#include "Numeric/Stencil.H"
#include "Set/Base.H"
#include "Set/Set.H"
#include "Unit/Unit.H"
#include "Util/Util.H"

#include "AMReX_Array4.H"
#include "AMReX_Box.H"
#include "AMReX_GpuLaunchFunctsC.H"
#include "AMReX_GpuQualifiers.H"
#include "AMReX_MFIter.H"
#include "AMReX_MultiFab.H"
#include "AMReX_MultiFabUtil.H"

namespace Integrator
{
void
Agglomeration::Parse(Agglomeration &value, IO::ParmParse &pp)
{
    BL_PROFILE("Integrator::Agglomeration::Agglomeration()")
    pp.queryclass<Flame>("flame", &value);

    if (value.max_level >= 1)
        Util::Abort(INFO, "Agglomeration does not support mesh refinement; set amr.max_level = 0");

    pp.query_default("p", value.agglom.p, "1.0", Unit::Less());
    pp.query_default("epsilon", value.agglom.epsilon, "1.0_m", Unit::Length());
    pp.query_default("omega", value.agglom.omega, "1.0", Unit::Less());
    pp.query_default("gamma", value.agglom.gamma, AMREX_SPACEDIM == 2 ? "6.0_m/J" : "6.0_m^2/J", (Unit::Length() ^ (AMREX_SPACEDIM - 1)) / Unit::Energy());
    pp.query_default("kappa_g", value.agglom.kappa_g, AMREX_SPACEDIM == 2 ? "1.0e-4_1/J" : "1.0e-4_m/J", (Unit::Length() ^ (AMREX_SPACEDIM - 2)) / Unit::Energy());
    pp.query_default("mu", value.agglom.mu, AMREX_SPACEDIM == 2 ? "1.0_m^5/J/s" : "1.0_m^6/J/s", (Unit::Length() ^ (AMREX_SPACEDIM + 3)) / Unit::Energy() / Unit::Time());
    pp.query_default("kappa_M", value.agglom.kappa_M, AMREX_SPACEDIM == 2 ? "1.0e-4_m^4/J/s" : "1.0e-4_m^5/J/s", (Unit::Length() ^ (AMREX_SPACEDIM + 2)) / Unit::Energy() / Unit::Time());
    pp.query_default("q", value.agglom.q, "3.0", Unit::Less());

    pp.select_default<IC::Constant, IC::Random, IC::Ellipse, IC::Expression, IC::PSRead>("alpha.ic", value.agglom.alpha_ic, value.geom);
    pp.select_default<BC::Constant>("alpha.bc", value.agglom.alpha_bc, 1);

    value.RegisterNewFab(value.agglom.alpha_old, value.agglom.alpha_bc, 1, 1, "agglom.alpha_old", false);
    value.RegisterNewFab(value.agglom.alpha, value.agglom.alpha_bc, 1, 1, "agglom.alpha", true);
    value.RegisterNewFab(value.agglom.dalpha, value.agglom.alpha_bc, 1, 1, "agglom.dalpha", true);
    value.RegisterNewFab(value.agglom.driving_force, value.agglom.alpha_bc, 1, 1, "agglom.driving_force", true, false);
    value.AddField<Set::Vector, Set::Hypercube::Cell>(value.agglom.evolution_driving_force, nullptr, 1, 1, "agglom.evolution_driving_force", false, false);

    value.RegisterIntegratedVariable(&value.agglom.V, "agglom.V");
}

void
Agglomeration::ScaleByComplement(MultiFab &dst, const MultiFab &src, int srccomp, int dstcomp, int numcomp, int nghost)
{
    BL_PROFILE("Integrator::Agglomeration::ScaleByComplement")
    int nCompSrc = src.nComp();
    int nGrowSrc = src.nGrow();

    for (int comp = 0; comp < numcomp; ++comp)
    {
        Util::Assert(INFO, TEST(src.min(srccomp + comp) >= 0.0));
        Util::Assert(INFO, TEST(src.max(srccomp + comp) <= 1.0));
    }

    MultiFab complement(src.boxArray(), src.DistributionMap(), nCompSrc, nGrowSrc);
    complement.setVal(1.0);

    MultiFab::Subtract(complement, src, 0, 0, nCompSrc, nGrowSrc);
    MultiFab::Multiply(dst, complement, srccomp, dstcomp, numcomp, nghost);
}

void
Agglomeration::ScaleByPhiComplement(const Set::Field<Set::Scalar> &mf, int lev)
{
    BL_PROFILE("Integrator::Agglomeration::ScaleByPhiComplement")
    int nComp = mf[lev]->nComp();
    int nGrow = mf[lev]->nGrow();
    MultiFab cell_based_phi(mf[lev]->boxArray(), mf[lev]->DistributionMap(), nComp, nGrow);
    average_node_to_cellcenter(cell_based_phi, 0, *phi_mf[lev], 0, nComp, nGrow);
    ScaleByComplement(*mf[lev], cell_based_phi, 0, 0, nComp, nGrow);
}

void
Agglomeration::Initialize(int lev)
{
    BL_PROFILE("Integrator::Agglomeration::Initialize")
    Flame::Initialize(lev);

    agglom.alpha_old[lev]->setVal(0.0);
    agglom.driving_force[lev]->setVal(0.0);
    agglom.alpha_ic->Initialize(lev, agglom.alpha);
    ScaleByPhiComplement(agglom.alpha, lev);
}

void
Agglomeration::Advance(int lev, Set::Scalar time, Set::Scalar dt)
{
    BL_PROFILE("Integrator::Agglomeration::Advance")
    Flame::Advance(lev, time, dt);

    std::swap(agglom.alpha_old[lev], agglom.alpha[lev]);
    const Set::Scalar *dx = geom[lev].CellSize();

    // compute the driving force
    for (amrex::MFIter mfi(*agglom.alpha[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.tilebox();
        Set::Patch<const Set::Scalar> alpha_old = agglom.alpha_old.Patch(lev, mfi);
        Set::Patch<Set::Scalar> driving_force = agglom.driving_force.Patch(lev, mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Set::Scalar alpha = alpha_old(i, j, k);
            Set::Scalar g_kappa_denom = std::sqrt(agglom.gamma * agglom.gamma * std::pow(alpha * alpha * (1 - alpha) * (1 - alpha), agglom.p) + agglom.kappa_g * agglom.kappa_g * agglom.epsilon * agglom.epsilon);
            Set::Scalar g_kappa = 1.0 / g_kappa_denom;
            Set::Scalar g_kappa_prime = -(2 * agglom.p * agglom.gamma * agglom.gamma * std::pow(1 - alpha, 2 * agglom.p) * std::pow(alpha, 2 * agglom.p - 1) - 2 * agglom.p * agglom.gamma * agglom.gamma * std::pow(1 - alpha, 2 * agglom.p - 1) * std::pow(alpha, 2 * agglom.p)) / (2 * g_kappa_denom * g_kappa_denom * g_kappa_denom);
            Set::Scalar w = agglom.omega * alpha * alpha * (1 - alpha) * (1 - alpha);
            Set::Scalar w_prime = 2 * agglom.omega * alpha * (1 - alpha) * (1 - 2 * alpha);
            Set::Vector grad_alpha = Numeric::Gradient(alpha_old, i, j, k, 0, dx);
            Set::Scalar norm_sq_grad_alpha = grad_alpha.squaredNorm();
            Set::Scalar lap_alpha = Numeric::Laplacian(alpha_old, i, j, k, 0, dx);

            driving_force(i, j, k) = (g_kappa_prime * w + g_kappa * w_prime) / agglom.epsilon - agglom.epsilon * g_kappa * lap_alpha - agglom.epsilon / 2 * g_kappa_prime * norm_sq_grad_alpha;
        });
    }

    agglom.driving_force[lev]->FillBoundary(geom[lev].periodicity());

    // compute the evolution of alpha
    for (amrex::MFIter mfi(*agglom.alpha[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.tilebox();
        Set::Patch<const Set::Scalar> alpha_old = agglom.alpha_old.Patch(lev, mfi);
        Set::Patch<Set::Scalar> alpha_new = agglom.alpha.Patch(lev, mfi);
        Set::Patch<Set::Scalar> dalpha = agglom.dalpha.Patch(lev, mfi);
        Set::Patch<const Set::Scalar> driving_force = agglom.driving_force.Patch(lev, mfi);
        Set::Patch<const Set::Scalar> eta_old = eta_mf.Patch(lev, mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Set::Scalar alpha = alpha_old(i, j, k);
            Set::Scalar eta = eta_old(i, j, k);
            Set::Scalar tildeM_kappa = (agglom.mu * alpha * alpha * (1 - alpha) * (1 - alpha) + agglom.kappa_M * agglom.epsilon) * std::pow(1 - eta, agglom.q);
            Set::Scalar dtildeM_kappa_dalpha = 2 * agglom.mu * alpha * (1 - 2 * alpha) * (1 - alpha) * std::pow(1 - eta, agglom.q);
            Set::Scalar dtildeM_kappa_deta = -agglom.q * ((1 - alpha) * (1 - alpha) * alpha * alpha * agglom.mu + agglom.epsilon * agglom.kappa_M) * std::pow(1 - eta, agglom.q - 1);

            Set::Vector grad_alpha = Numeric::Gradient(alpha_old, i, j, k, 0, dx);
            Set::Vector grad_eta = Numeric::Gradient(eta_old, i, j, k, 0, dx);
            Set::Vector grad_driving_force = Numeric::Gradient(driving_force, i, j, k, 0, dx);
            Set::Scalar lap_driving_force = Numeric::Laplacian(driving_force, i, j, k, 0, dx);

            dalpha(i, j, k) = dt / agglom.epsilon * (tildeM_kappa * lap_driving_force + dtildeM_kappa_dalpha * grad_alpha.dot(grad_driving_force) + dtildeM_kappa_deta * grad_eta.dot(grad_driving_force));

            alpha_new(i, j, k) = alpha_old(i, j, k) + dalpha(i, j, k);
        });
    }
}

void
Agglomeration::Integrate(int amrlev, Set::Scalar time, int step, const amrex::MFIter &mfi, const amrex::Box &box)
{
    BL_PROFILE("Integrator::Agglomeration::Integrate");
    Flame::Integrate(amrlev, time, step, mfi, box);

    const Set::Scalar *dx = geom[amrlev].CellSize();
    Set::Scalar dv = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
    Set::Patch<const Set::Scalar> alpha = agglom.alpha.Patch(amrlev, mfi);

    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        agglom.V += alpha(i, j, k) * dv;
    });
}
}
