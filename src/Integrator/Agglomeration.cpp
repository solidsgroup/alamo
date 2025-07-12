#include <cmath>
#include <utility>

#include "Agglomeration.H"
#include "Flame.H"

#include "BC/Constant.H"
#include "IC/BMP.H"
#include "IC/BetaDistribution.H"
#include "IC/Constant.H"
#include "IC/Expression.H"
#include "IC/PNG.H"
#include "IC/Random.H"
#include "IO/ParmParse.H"
#include "Numeric/Stencil.H"
#include "Set/Base.H"
#include "Set/Set.H"
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

    // Method to use for agglomeration kinetics
    pp.query_enum_case_insensitive("kinetics_method", value.agglom.kinetics_method);

    // Diffuse interface length of the agglomerate phase
    pp.query_default("epsilon", value.agglom.epsilon, 1e-7);
    // Chemical potential of the agglomerating material
    pp.query_default("gamma", value.agglom.gamma, 0.0005);
    // Surface tension of the agglomerating material
    pp.query_default("kappa", value.agglom.kappa, 1.0);
    // Agglomeration mobility
    pp.query_default("L_0", value.agglom.L_0, 1.0);
    // Agglomeration mobility exponent
    pp.query_default("n", value.agglom.n, 1.0);

    if (value.agglom.kinetics_method == AgglomerationKinetics::ConstrainedAllenCahn)
    {
        // Lagrangian multiplier for the agglomerate volume fraction
        // constraint
        pp.query_default("lambda", value.agglom.lambda, 100.0);
        // Prescribed volume fraction of agglomerate in the domain. If
        // unset, it will default to the initial volume fraction of
        // agglomerate in the domain based on the given initial
        // condition.
        pp.query("V_0", value.agglom.V_0);
    }

    // If the gradient of alpha is larger than this value, refine/regrid alpha
    pp.query_default("gradient_refinement_threshold", value.agglom.gradient_refinement_threshold, 1e100);
    // If eta is larger than this value, do not refine/regrid alpha
    pp.query_default("eta_refinement_threshold", value.agglom.eta_refinement_threshold, 1e100);

    // Initial conditions for agglomerate order parameter
    pp.select_default<IC::Constant, IC::Expression, IC::BMP, IC::PNG, IC::Random, IC::BetaDistribution>("alpha.ic", value.agglom.alpha_ic, value.geom);
    // Boundary conditions for agglomerate order parameter
    pp.select_default<BC::Constant>("alpha.bc", value.agglom.alpha_bc, 1);

    value.RegisterNewFab(value.agglom.alpha_old, value.agglom.alpha_bc, 1, 1, "agglom.alpha_old", false);
    value.RegisterNewFab(value.agglom.alpha, value.agglom.alpha_bc, 1, 1, "agglom.alpha", true);
    value.RegisterNewFab(value.agglom.free_energy_derivative, value.agglom.alpha_bc, 1, 1, "agglom.free_energy_derivative", false);

    if (value.agglom.kinetics_method == AgglomerationKinetics::ConstrainedAllenCahn)
        value.RegisterIntegratedVariable(&value.agglom.V, "agglom.V", true);
}

void
Agglomeration::Initialize(int lev)
{
    Flame::Initialize(lev);

    agglom.alpha_old[lev]->setVal(0.0);
    agglom.alpha_ic->Initialize(lev, agglom.alpha);
    agglom.free_energy_derivative[lev]->setVal(0.0);

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

    Set::Scalar agglom_volume_fraction = NAN;
    if (agglom.kinetics_method == AgglomerationKinetics::ConstrainedAllenCahn)
        agglom_volume_fraction = agglom.V / geom[0].ProbSize();

    if (agglom.kinetics_method == AgglomerationKinetics::ConstrainedAllenCahn && time == 0.0 && std::isnan(agglom.V_0))
    {
        agglom.V_0 = agglom_volume_fraction;
        Util::DebugMessage(INFO, "agglom.V_0 dynamically set to ", agglom.V_0);
    }

    for (amrex::MFIter mfi(*agglom.alpha[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.tilebox();
        Set::Patch<const Set::Scalar> alpha_old = agglom.alpha_old.Patch(lev, mfi);
        Set::Patch<Set::Scalar> alpha = agglom.alpha.Patch(lev, mfi);
        Set::Patch<Set::Scalar> free_energy_derivative = agglom.free_energy_derivative.Patch(lev, mfi);
        Set::Patch<Set::Scalar> eta = eta_mf.Patch(lev, mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            // calculate the Laplacian of alpha
            Set::Scalar laplacian_alpha = Numeric::Laplacian(alpha_old, i, j, k, 0, dx);

            // calculate the variational derivative components
            Set::Scalar chem_term = 2 * agglom.gamma / agglom.epsilon * alpha_old(i, j, k) * (1 - alpha_old(i, j, k)) * (1 - 2 * alpha_old(i, j, k));
            Set::Scalar grad_term = -agglom.epsilon * agglom.kappa * laplacian_alpha;

            free_energy_derivative(i, j, k) = chem_term + grad_term;
        });

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            // calculate effective mobility L
            Set::Scalar L = agglom.L_0 * std::pow(1 - eta(i, j, k), agglom.n);

            if (agglom.kinetics_method == AgglomerationKinetics::ConstrainedAllenCahn)
            {
                // constrained Allen-Cahn equation
                Set::Scalar lagrange_term = agglom.lambda * (agglom_volume_fraction - agglom.V_0);
                Set::Scalar dalpha_dt = -L * (free_energy_derivative(i, j, k) + lagrange_term);
                alpha(i, j, k) = alpha_old(i, j, k) + dt * dalpha_dt;
            }
            else if (agglom.kinetics_method == AgglomerationKinetics::CahnHilliard)
            {
                // calculate the Laplacian of the variational derivative
                Set::Scalar free_energy_laplacian = Numeric::Laplacian(free_energy_derivative, i, j, k, 0, dx);

                // Cahn-Hilliard equation
                Set::Scalar dalpha_dt = L * free_energy_laplacian;
                alpha(i, j, k) = alpha_old(i, j, k) + dt * dalpha_dt;
            }
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
                if (grad.lpNorm<2>() * dr > agglom.gradient_refinement_threshold)
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
