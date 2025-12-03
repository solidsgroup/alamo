#include <cmath>
#include <utility>

#include "Agglomeration.H"
#include "Flame.H"

#include "BC/Constant.H"
#include "IC/BMP.H"
#include "IC/Constant.H"
#include "IC/Expression.H"
#include "IC/PNG.H"
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
#include "AMReX_TagBox.H"

namespace Integrator
{
void
Agglomeration::Parse(Agglomeration &value, IO::ParmParse &pp)
{
    pp.queryclass<Flame>("flame", &value);

    // Interface energy
    pp.query_default("sigma", value.agglom.sigma, "1.0_J/m^2", Unit::Energy() / Unit::Area());
    // Diffuse interface length of the agglomerate phase
    pp.query_default("epsilon", value.agglom.epsilon, "1.0_m", Unit::Length());

    // Prescribed agglomerate mobility parameter
    pp.query_default("L_0", value.agglom.L_0, "1.0_m^3/J/s", Unit::Volume() / Unit::Energy() / Unit::Time());
    // Agglomeration mobility exponent
    pp.query_default("n", value.agglom.n, "1.0", Unit::Less());

    // Lagrangian multiplier for the agglomerate volume constraint
    pp.query_default("lambda", value.agglom.lambda, "1.0_J/m^3", Unit::Energy() / Unit::Volume());
    // If true, the value of V_0 will be ignored, and the
    // prescribed volume of agglomerate will be calculated just
    // before the first timestep after all initialization has
    // occurred. If false, a value of V_0 must be specified.
    pp.query_default("calculate_initial_volume", value.agglom.calculate_initial_volume, true);
    // Prescribed volume of agglomerate in the domain. This value
    // is only used if calculate_initial_volume is false.
    pp.query_default("V_0", value.agglom.V_0, "1.0", Unit::Volume());

    // If the gradient of alpha is larger than this value, refine/regrid alpha
    pp.query_default("gradient_alpha_refinement_threshold", value.agglom.gradient_alpha_refinement_threshold, "1.0", Unit::Less());
    // If eta is larger than this value, do NOT refine/regrid alpha
    pp.query_default("eta_refinement_threshold", value.agglom.eta_refinement_threshold, "1.0", Unit::Less());

    // Initial conditions for agglomerate order parameter
    pp.select_default<IC::Constant, IC::Expression, IC::BMP, IC::PNG, IC::Random>("alpha.ic", value.agglom.alpha_ic, value.geom);
    // Boundary conditions for agglomerate order parameter
    pp.select_default<BC::Constant>("alpha.bc", value.agglom.alpha_bc, 1);

    value.RegisterNewFab(value.agglom.alpha_old, value.agglom.alpha_bc, 1, 1, "agglom.alpha_old", false);
    value.RegisterNewFab(value.agglom.alpha, value.agglom.alpha_bc, 1, 1, "agglom.alpha", true);

    value.RegisterIntegratedVariable(&value.agglom.V, "agglom.V", true);
}

void
Agglomeration::Initialize(int lev)
{
    Flame::Initialize(lev);

    agglom.alpha_old[lev]->setVal(0.0);
    agglom.alpha_ic->Initialize(lev, agglom.alpha);

    int nComp = agglom.alpha[lev]->nComp();
    int nGrow = agglom.alpha[lev]->nGrow();
    MultiFab cell_based_phi(agglom.alpha[lev]->boxArray(), agglom.alpha[lev]->DistributionMap(), nComp, nGrow);
    average_node_to_cellcenter(cell_based_phi, 0, *phi_mf[lev], 0, nComp, nGrow);

    scaleByComplement(*agglom.alpha[lev], cell_based_phi, 0, 0, nComp, nGrow);
}

void
Agglomeration::scaleByComplement(MultiFab &dst, const MultiFab &src, int srccomp, int dstcomp, int numcomp, int nghost)
{
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
Agglomeration::Advance(int lev, Set::Scalar time, Set::Scalar dt)
{
    Flame::Advance(lev, time, dt);

    std::swap(agglom.alpha_old[lev], agglom.alpha[lev]);
    const Set::Scalar *dx = geom[lev].CellSize();

    if (time == 0.0 && agglom.calculate_initial_volume)
    {
        agglom.V_0 = agglom.V;
        Util::DebugMessage(INFO, "agglom.V_0 dynamically set to ", agglom.V_0);
    }

    for (amrex::MFIter mfi(*agglom.alpha[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.tilebox();
        Set::Patch<const Set::Scalar> alpha_old = agglom.alpha_old.Patch(lev, mfi);
        Set::Patch<Set::Scalar> alpha = agglom.alpha.Patch(lev, mfi);
        Set::Patch<Set::Scalar> eta = eta_mf.Patch(lev, mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            // calculate the Laplacian of alpha
            Set::Scalar laplacian_alpha = Numeric::Laplacian(alpha_old, i, j, k, 0, dx);

            // calculate the variational derivative components
            Set::Scalar free_energy_derivative = 2 * agglom.sigma / agglom.epsilon * alpha_old(i, j, k) * (1 - alpha_old(i, j, k)) * (1 - 2 * alpha_old(i, j, k)) - agglom.sigma * agglom.epsilon * laplacian_alpha;

            // calculate effective mobility L
            Set::Scalar L = agglom.L_0 * std::pow(1 - eta(i, j, k), agglom.n);

            // constrained Allen-Cahn equation
            Set::Scalar lagrange_term = agglom.lambda * (agglom.V / agglom.V_0 - 1);
            Set::Scalar dalpha_dt = -L * (free_energy_derivative + lagrange_term);

            // evolve alpha
            alpha(i, j, k) = alpha_old(i, j, k) + dt * dalpha_dt;
        });
    }
}

void
Agglomeration::TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, Set::Scalar time, int ngrow)
{
    Flame::TagCellsForRefinement(lev, a_tags, time, ngrow);

    const Set::Vector dx(geom[lev].CellSize());
    Set::Scalar dr = dx.lpNorm<2>();

    for (amrex::MFIter mfi(*agglom.alpha[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.tilebox();
        Set::Patch<char> tags = a_tags.array(mfi);
        Set::Patch<const Set::Scalar> alpha = agglom.alpha.Patch(lev, mfi);
        Set::Patch<const Set::Scalar> eta = eta_mf.Patch(lev, mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            if (eta(i, j, k) < agglom.eta_refinement_threshold)
            {
                Set::Vector grad = Numeric::Gradient(alpha, i, j, k, 0, dx.data());
                if (grad.lpNorm<2>() * dr > agglom.gradient_alpha_refinement_threshold)
                    tags(i, j, k) = amrex::TagBox::SET;
            }
        });
    }
}

void
Agglomeration::Integrate(int amrlev, Set::Scalar time, int step, const amrex::MFIter &mfi, const amrex::Box &box)
{
    Flame::Integrate(amrlev, time, step, mfi, box);

    const Set::Scalar *dx = geom[amrlev].CellSize();
    Set::Scalar dv = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
    Set::Patch<const Set::Scalar> alpha = agglom.alpha.Patch(amrlev, mfi);
    amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
        agglom.V += alpha(i, j, k, 0) * dv;
    });
}
}
