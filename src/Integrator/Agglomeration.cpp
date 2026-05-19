#include <utility>

#include "Agglomeration.H"
#include "Flame.H"

#include "BC/Constant.H"
#include "IC/Constant.H"
#include "IC/Ellipse.H"
#include "IC/Expression.H"
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

    // Diffuse interface length parameter
    pp.query_default("epsilon", value.agglom.epsilon, "1.0_m", Unit::Length());
    // DDCH bulk chemical potential coefficient
    pp.query_default("omega", value.agglom.omega, "1.0", Unit::Less());
    // DDCH energy restriction parameter
    pp.query_default("gamma", value.agglom.gamma, "6.0", Unit::Less());
    // DDCH regularization normalization parameter
    pp.query_default("kappa", value.agglom.kappa, "1.0e-4", Unit::Less());
    // Bulk coefficient in the DDCH agglomeration mobility
    // M_agglom(alpha) = mu * alpha^2 * (1-alpha)^2 + kappa * epsilon
    pp.query_default("mu", value.agglom.mu, "1.0", Unit::Less());
    // Free-flow reaction mobility (eta-restricted reaction sink uses
    // M_reaction * (1-eta)^3)
    pp.query_default("M_reaction", value.agglom.M_reaction, "1.0_m^3/J/s", Unit::Volume() / Unit::Energy() / Unit::Time());

    // Initial condition for the agglomerate order parameter
    pp.select_default<IC::Constant, IC::Random, IC::Ellipse, IC::Expression>("alpha.ic", value.agglom.alpha_ic, value.geom);
    // Boundary condition for agglomerate order parameter
    pp.select_default<BC::Constant>("alpha.bc", value.agglom.alpha_bc, 1);

    value.RegisterNewFab(value.agglom.alpha_old, value.agglom.alpha_bc, 1, 1, "agglom.alpha_old", false);
    value.RegisterNewFab(value.agglom.alpha, value.agglom.alpha_bc, 1, 1, "agglom.alpha", true);
    value.RegisterNewFab(value.agglom.dalpha_agglom, value.agglom.alpha_bc, 1, 1, "agglom.dalpha_agglom", true);
    value.RegisterNewFab(value.agglom.dalpha_reaction, value.agglom.alpha_bc, 1, 1, "agglom.dalpha_reaction", true);
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

    // Pass 1: compute the DDCH driving force (variational derivative of free energy)
    for (amrex::MFIter mfi(*agglom.alpha[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.tilebox();
        Set::Patch<const Set::Scalar> alpha_old = agglom.alpha_old.Patch(lev, mfi);
        Set::Patch<Set::Scalar> driving_force = agglom.driving_force.Patch(lev, mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Set::Scalar alpha = alpha_old(i, j, k);
            Set::Scalar g_denom = sqrt(agglom.gamma * agglom.gamma * (1 - alpha) * (1 - alpha) * alpha * alpha + agglom.kappa * agglom.kappa * agglom.epsilon * agglom.epsilon);
            Set::Scalar g = 1.0 / g_denom;
            Set::Scalar gprime = -agglom.gamma * agglom.gamma * (alpha - 1) * alpha * (2 * alpha - 1) / (g_denom * g_denom * g_denom);
            Set::Scalar f = agglom.omega / 4 * alpha * alpha * (1 - alpha) * (1 - alpha);
            Set::Scalar fprime = agglom.omega / 2 * alpha * (1 - alpha) * (1 - 2 * alpha);
            Set::Vector grad_alpha = Numeric::Gradient(alpha_old, i, j, k, 0, dx);
            Set::Scalar norm_sq_grad_alpha = grad_alpha.squaredNorm();
            Set::Scalar lap_alpha = Numeric::Laplacian(alpha_old, i, j, k, 0, dx);

            driving_force(i, j, k) = g * fprime / agglom.epsilon
                                     - agglom.epsilon * (gprime * norm_sq_grad_alpha + g * lap_alpha)
                                     + gprime * (f / agglom.epsilon + agglom.epsilon / 2 * norm_sq_grad_alpha);
        });
    }

    agglom.driving_force[lev]->FillBoundary(geom[lev].periodicity());

    // ========================================================================
    // Pass 2: advance alpha. Choose exactly ONE of the two formulations below.
    // The product-rule formulation is active by default. To switch, comment
    // out the entire FORMULATION A block and uncomment the entire
    // FORMULATION B block.
    // ========================================================================

    // ------------------------------------------------------------------------
    // FORMULATION A (active): product-rule expansion of
    //   (1/eps) * div(tildeM_agglom * grad(dF/dalpha))
    //
    // = (1/eps) * [ tildeM_agglom * Delta(dF/dalpha)
    //             + d(tildeM_agglom)/dalpha * grad(alpha) . grad(dF/dalpha)
    //             + d(tildeM_agglom)/deta   * grad(eta)   . grad(dF/dalpha) ]
    //
    // The grad(eta) term is required because tildeM_agglom depends on eta
    // through the (1-eta)^3 binder-restriction factor.
    // ------------------------------------------------------------------------
    for (amrex::MFIter mfi(*agglom.alpha[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.tilebox();
        Set::Patch<const Set::Scalar> alpha_old = agglom.alpha_old.Patch(lev, mfi);
        Set::Patch<Set::Scalar> alpha = agglom.alpha.Patch(lev, mfi);
        Set::Patch<Set::Scalar> dalpha_agglom = agglom.dalpha_agglom.Patch(lev, mfi);
        Set::Patch<Set::Scalar> dalpha_reaction = agglom.dalpha_reaction.Patch(lev, mfi);
        Set::Patch<const Set::Scalar> driving_force = agglom.driving_force.Patch(lev, mfi);
        Set::Patch<const Set::Scalar> eta = eta_mf.Patch(lev, mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Set::Scalar alpha_val = alpha_old(i, j, k);
            Set::Scalar eta_complement = 1 - eta(i, j, k);
            Set::Scalar eta_scale = eta_complement * eta_complement * eta_complement;
            Set::Scalar deta_scale_deta = -3 * eta_complement * eta_complement;

            Set::Scalar M_agglom = agglom.mu * alpha_val * alpha_val * (1 - alpha_val) * (1 - alpha_val) + agglom.kappa * agglom.epsilon;
            Set::Scalar dM_agglom_dalpha = 2.0 * agglom.mu * alpha_val * (alpha_val - 1) * (2 * alpha_val - 1);
            Set::Scalar tildeM_agglom = M_agglom * eta_scale;
            Set::Scalar dtildeM_agglom_dalpha = dM_agglom_dalpha * eta_scale;
            Set::Scalar dtildeM_agglom_deta = M_agglom * deta_scale_deta;
            Set::Scalar tildeM_reaction = agglom.M_reaction * eta_scale;

            Set::Vector grad_alpha = Numeric::Gradient(alpha_old, i, j, k, 0, dx);
            Set::Vector grad_eta = Numeric::Gradient(eta, i, j, k, 0, dx);
            Set::Vector grad_driving_force = Numeric::Gradient(driving_force, i, j, k, 0, dx);
            Set::Scalar lap_driving_force = Numeric::Laplacian(driving_force, i, j, k, 0, dx);

            dalpha_agglom(i, j, k) = dt / agglom.epsilon * (tildeM_agglom * lap_driving_force + dtildeM_agglom_dalpha * grad_alpha.dot(grad_driving_force) + dtildeM_agglom_deta * grad_eta.dot(grad_driving_force));
            dalpha_reaction(i, j, k) = -tildeM_reaction * driving_force(i, j, k) * dt;

            alpha(i, j, k) = alpha_old(i, j, k) + dalpha_agglom(i, j, k) + dalpha_reaction(i, j, k);
        });
    }

    // ------------------------------------------------------------------------
    // FORMULATION B (commented out): divergence form
    //
    //   dalpha/dt|_DDCH = (1/eps) * div(tildeM_agglom * grad(dF/dalpha))
    //
    // The discrete divergence handles the eta-dependence of tildeM_agglom
    // naturally, so no explicit grad(eta) term is needed. Computed in two
    // sub-passes: first build the flux  tildeM_agglom * grad(dF/dalpha)  on
    // cells (with a FillBoundary afterward), then take its discrete
    // divergence in the alpha update.
    // ------------------------------------------------------------------------
    // for (amrex::MFIter mfi(*agglom.alpha[lev], true); mfi.isValid(); ++mfi)
    // {
    //     const amrex::Box &bx = mfi.tilebox();
    //     Set::Patch<const Set::Scalar> alpha_old = agglom.alpha_old.Patch(lev, mfi);
    //     Set::Patch<const Set::Scalar> driving_force = agglom.driving_force.Patch(lev, mfi);
    //     Set::Patch<const Set::Scalar> eta = eta_mf.Patch(lev, mfi);
    //     Set::Patch<Set::Vector> evolution_driving_force = agglom.evolution_driving_force.Patch(lev, mfi);
    //
    //     amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
    //         Set::Scalar alpha_val = alpha_old(i, j, k);
    //         Set::Scalar eta_complement = 1 - eta(i, j, k);
    //         Set::Scalar eta_scale = eta_complement * eta_complement * eta_complement;
    //         Set::Scalar M_agglom = agglom.mu * alpha_val * alpha_val * (1 - alpha_val) * (1 - alpha_val) + agglom.kappa * agglom.epsilon;
    //         Set::Scalar tildeM_agglom = M_agglom * eta_scale;
    //         Set::Vector grad_driving_force = Numeric::Gradient(driving_force, i, j, k, 0, dx);
    //
    //         evolution_driving_force(i, j, k) = tildeM_agglom * grad_driving_force;
    //     });
    // }
    //
    // agglom.evolution_driving_force[lev]->FillBoundary(geom[lev].periodicity());
    //
    // for (amrex::MFIter mfi(*agglom.alpha[lev], true); mfi.isValid(); ++mfi)
    // {
    //     const amrex::Box &bx = mfi.tilebox();
    //     Set::Patch<const Set::Scalar> alpha_old = agglom.alpha_old.Patch(lev, mfi);
    //     Set::Patch<Set::Scalar> alpha = agglom.alpha.Patch(lev, mfi);
    //     Set::Patch<Set::Scalar> dalpha_agglom = agglom.dalpha_agglom.Patch(lev, mfi);
    //     Set::Patch<Set::Scalar> dalpha_reaction = agglom.dalpha_reaction.Patch(lev, mfi);
    //     Set::Patch<const Set::Scalar> driving_force = agglom.driving_force.Patch(lev, mfi);
    //     Set::Patch<const Set::Scalar> eta = eta_mf.Patch(lev, mfi);
    //     Set::Patch<const Set::Vector> evolution_driving_force = agglom.evolution_driving_force.Patch(lev, mfi);
    //
    //     amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
    //         Set::Scalar eta_complement = 1 - eta(i, j, k);
    //         Set::Scalar eta_scale = eta_complement * eta_complement * eta_complement;
    //         Set::Scalar tildeM_reaction = agglom.M_reaction * eta_scale;
    //
    //         dalpha_agglom(i, j, k) = dt / agglom.epsilon * Numeric::Divergence(evolution_driving_force, i, j, k, dx);
    //         dalpha_reaction(i, j, k) = -tildeM_reaction * driving_force(i, j, k) * dt;
    //
    //         alpha(i, j, k) = alpha_old(i, j, k) + dalpha_agglom(i, j, k) + dalpha_reaction(i, j, k);
    //     });
    // }
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
