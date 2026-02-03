#include <algorithm>
#include <utility>

#include "Agglomeration.H"
#include "Flame.H"

#include "BC/Constant.H"
#include "IC/Constant.H"
#include "IC/Expression.H"
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
#include "AMReX_IntVect.H"
#include "AMReX_MFIter.H"
#include "AMReX_MultiFab.H"
#include "AMReX_MultiFabUtil.H"
#include "AMReX_ParallelDescriptor.H"
#include "AMReX_TagBox.H"

namespace Integrator
{
void
Agglomeration::Parse(Agglomeration &value, IO::ParmParse &pp)
{
    BL_PROFILE("Integrator::Agglomeration::Agglomeration()")
    pp.queryclass<Flame>("flame", &value);

    // Interface energy
    pp.query_default("sigma", value.agglom.sigma, "1.0_J/m^2", Unit::Energy() / Unit::Area());
    // Diffuse interface length parameter
    pp.query_default("epsilon", value.agglom.epsilon, "1.0_m", Unit::Length());

    // Free-flow agglomeration mobility parameter
    pp.query_default("L_0_agglom", value.agglom.L_0_agglom, "1.0_m^3/J/s", Unit::Volume() / Unit::Energy() / Unit::Time());
    // Free-flow reaction mobility parameter
    pp.query_default("L_0_reaction", value.agglom.L_0_reaction, "1.0_m^3/J/s", Unit::Volume() / Unit::Energy() / Unit::Time());

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

    // If eta is larger than this value, do NOT refine/regrid alpha
    pp.query_default("eta_refinement_threshold", value.agglom.eta_refinement_threshold, "1.0", Unit::Less());
    // If the gradient of alpha is larger than this value and eta is smaller than `eta_refinement_threshold`, refine/regrid alpha
    pp.query_default("gradient_alpha_refinement_threshold", value.agglom.gradient_alpha_refinement_threshold, "1.0", Unit::Less());

    // Initial condition for the agglomerate order parameter
    pp.select_default<IC::Constant, IC::Expression>("alpha.ic", value.agglom.alpha_ic, value.geom);
    // Boundary condition for for agglomerate order parameter
    pp.select_default<BC::Constant>("alpha.bc", value.agglom.alpha_bc, 1);

    value.RegisterNewFab(value.agglom.alpha_0, value.agglom.alpha_bc, 1, 1, "agglom.alpha_0", false);
    value.RegisterNewFab(value.agglom.alpha_old, value.agglom.alpha_bc, 1, 1, "agglom.alpha_old", false);
    value.RegisterNewFab(value.agglom.alpha, value.agglom.alpha_bc, 1, 1, "agglom.alpha", true);
    value.RegisterNewFab(value.agglom.dalpha_agglom, value.agglom.alpha_bc, 1, 1, "agglom.dalpha_agglom", true);
    value.RegisterNewFab(value.agglom.dalpha_reaction, value.agglom.alpha_bc, 1, 1, "agglom.dalpha_reaction", true);

    value.RegisterIntegratedVariable(&value.agglom.V, "agglom.V", false);
    value.RegisterIntegratedVariable(&value.agglom.V_0, "agglom.V_0", false);
    value.RegisterIntegratedVariable(&value.agglom.dV, "agglom.dV", false);
    value.RegisterIntegratedVariable(&value.agglom.dV_0, "agglom.dV_0", false);
}

Set::Scalar
Agglomeration::CalculateInitialVolume(int target_resolution_level)
{
    BL_PROFILE("Integrator::Agglomeration::CalculateInitialVolume");

    // calculate refined cell count based on target resolution level
    amrex::IntVect refined_cells = geom[0].Domain().length();
    for (int i = 0; i < target_resolution_level; i++)
        refined_cells *= refRatio(std::min(i, max_level - 1));

    // create fine geometry with same physical domain but higher resolution
    amrex::Box fine_domain(amrex::IntVect(0), refined_cells - amrex::IntVect(1));
    amrex::Geometry fine_geom(fine_domain, geom[0].ProbDomain(), geom[0].Coord(), geom[0].isPeriodic());

    // create BoxArray and DistributionMapping for the fine grid
    amrex::BoxArray fine_ba(fine_domain);
    fine_ba.maxSize(32);
    amrex::DistributionMapping fine_dm(fine_ba);

    // create BoxArray for node-centered data
    amrex::BoxArray fine_nodal_ba = amrex::convert(fine_ba, amrex::IntVect::TheNodeVector());

    // create temporary fields
    Set::Field<Set::Scalar> fine_phi_field(1);
    fine_phi_field.Define(0, fine_nodal_ba, fine_dm, 1, 0);
    Set::Field<Set::Scalar> fine_alpha_field(1);
    fine_alpha_field.Define(0, fine_ba, fine_dm, 1, 0);

    // temporarily replace geom[0] with fine geometry so existing ICs use correct positions
    // this is a bit hacky, but necessary because when the ICs are parsed, they're associated with a geometry
    amrex::Geometry original_geom = geom[0];
    geom[0] = fine_geom;

    // initialize phi and alpha using the existing ICs (which now see fine_geom via geom[0])
    ic_phi->Initialize(0, fine_phi_field);
    agglom.alpha_ic->Initialize(0, fine_alpha_field);

    // restore original geometry
    geom[0] = original_geom;

    // average phi from nodes to cell centers and apply complement scaling
    amrex::MultiFab fine_phi_cc(fine_ba, fine_dm, 1, 0);
    amrex::average_node_to_cellcenter(fine_phi_cc, 0, *fine_phi_field[0], 0, 1, 0);
    ScaleByComplement(*fine_alpha_field[0], fine_phi_cc, 0, 0, 1, 0);

    // compute volume integral
    const Set::Scalar *dx = fine_geom.CellSize();
    Set::Scalar dv = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);
    Set::Scalar volume = fine_alpha_field[0]->sum() * dv;

    Util::Message(INFO, "Calculated initial volume on fine grid (level ", target_resolution_level, ", ", refined_cells[0], "x", refined_cells[1], "x", refined_cells[2], " cells): V_0 = ", volume);

    return volume;
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
    agglom.alpha_ic->Initialize(lev, agglom.alpha);
    ScaleByPhiComplement(agglom.alpha, lev);
}

void
Agglomeration::TimeStepBegin(Set::Scalar a_time, int a_iter)
{
    BL_PROFILE("Integrator::Agglomeration::TimeStepBegin");
    Flame::TimeStepBegin(a_time, a_iter);

    if (a_time == 0.0)
    {
        if (agglom.calculate_initial_volume)
            agglom.V_0 = CalculateInitialVolume(max_level);
        agglom.V = agglom.V_0;
        agglom.prev_finest_ba = grids[finest_level];
        agglom.prev_finest_level = finest_level;
        Util::Message(INFO, "agglom.V_0 = ", agglom.V_0, ", agglom.V = ", agglom.V);
    }

    agglom.dV = 0.0;
    agglom.dV_0 = 0.0;
}

void
Agglomeration::Advance(int lev, Set::Scalar time, Set::Scalar dt)
{
    BL_PROFILE("Integrator::Agglomeration::Advance")
    Flame::Advance(lev, time, dt);

    std::swap(agglom.alpha_old[lev], agglom.alpha[lev]);
    const Set::Scalar *dx = geom[lev].CellSize();
    Set::Scalar dv = AMREX_D_TERM(dx[0], *dx[1], *dx[2]);

    amrex::BoxArray cfba;
    if (lev < finest_level)
        cfba = amrex::coarsen(grids[lev + 1], refRatio(lev));

    for (amrex::MFIter mfi(*agglom.alpha[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.tilebox();
        Set::Patch<const Set::Scalar> alpha_old = agglom.alpha_old.Patch(lev, mfi);
        Set::Patch<Set::Scalar> alpha = agglom.alpha.Patch(lev, mfi);
        Set::Patch<Set::Scalar> dalpha_agglom = agglom.dalpha_agglom.Patch(lev, mfi);
        Set::Patch<Set::Scalar> dalpha_reaction = agglom.dalpha_reaction.Patch(lev, mfi);
        Set::Patch<Set::Scalar> eta = eta_mf.Patch(lev, mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            Set::Scalar laplacian_alpha = Numeric::Laplacian(alpha_old, i, j, k, 0, dx);

            Set::Scalar free_energy_derivative = 2 * agglom.sigma / agglom.epsilon * alpha_old(i, j, k) * (1 - alpha_old(i, j, k)) * (1 - 2 * alpha_old(i, j, k)) - agglom.sigma * agglom.epsilon * laplacian_alpha;

            Set::Scalar eta_complement = 1 - eta(i, j, k);
            Set::Scalar eta_scale = eta_complement * eta_complement * eta_complement;
            Set::Scalar L_agglom = agglom.L_0_agglom * eta_scale;
            Set::Scalar L_reaction = agglom.L_0_reaction * eta_scale;

            Set::Scalar lagrange_term = agglom.lambda * (agglom.V / agglom.V_0 - 1);
            dalpha_agglom(i, j, k) = -L_agglom * (free_energy_derivative + lagrange_term) * dt;
            dalpha_reaction(i, j, k) = -L_reaction * free_energy_derivative * dt;
            Set::Scalar dalpha = dalpha_agglom(i, j, k) + dalpha_reaction(i, j, k);

            alpha(i, j, k) = alpha_old(i, j, k) + dalpha;
        });

        // accumulate volume flux using complement to avoid double-counting with finer levels
        amrex::BoxList boxes_to_integrate;
        if (lev < finest_level)
        {
            amrex::BoxArray comp = amrex::complementIn(bx, cfba);
            boxes_to_integrate = comp.boxList();
        }
        else
            boxes_to_integrate.push_back(bx);

        // we do a manual loop instead of a ParallelFor to avoid
        // threading/race conditions resulting from the += operator on
        // agglom.dV and agglom.dV_0
        for (const amrex::Box &box : boxes_to_integrate)
        {
            const auto lo = amrex::lbound(box);
            const auto hi = amrex::ubound(box);
            for (int k = lo.z; k <= hi.z; ++k)
                for (int j = lo.y; j <= hi.y; ++j)
                    for (int i = lo.x; i <= hi.x; ++i)
                    {
                        agglom.dV += (dalpha_agglom(i, j, k) + dalpha_reaction(i, j, k)) * dv;
                        agglom.dV_0 += dalpha_reaction(i, j, k) * dv;
                    }
        }
    }
}

void
Agglomeration::TimeStepComplete(Set::Scalar a_time, int a_iter)
{
    BL_PROFILE("Integrator::Agglomeration::TimeStepComplete");
    Flame::TimeStepComplete(a_time, a_iter);

    amrex::ParallelDescriptor::ReduceRealSum(agglom.dV);
    amrex::ParallelDescriptor::ReduceRealSum(agglom.dV_0);
    agglom.V += agglom.dV;
    agglom.V_0 += agglom.dV_0;
}

void
Agglomeration::TagCellsForRefinement(int lev, amrex::TagBoxArray &a_tags, Set::Scalar time, int ngrow)
{
    BL_PROFILE("Integrator::Agglomeration::TagCellsForRefinement")
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
Agglomeration::Regrid(int lev, Set::Scalar time)
{
    BL_PROFILE("Integrator::Agglomeration::Regrid")
    Flame::Regrid(lev, time);

    agglom.alpha_ic->Initialize(lev, agglom.alpha_0);
    ScaleByPhiComplement(agglom.alpha_0, lev);

    // The goal here is, on regrid, to set the value of the
    // agglomerate phase to its value based on the IC without
    // disrupting the agglomerate phase near the regression front or
    // in the "burned" region of the simulation. The temp and eta
    // refinement criterion ensure that at the regression front, we'll
    // be at the finest level, so I only want to do reset the values
    // based on the IC when we're a.) above the
    // eta_refinement_threshold and b.) not already on the finest
    // level. There's one small caveat to b.): if a cell has just been
    // refined from finest_level - 1 to the finest_level, we do want
    // to reset it to its value based on the IC.
    for (amrex::MFIter mfi(*agglom.alpha[lev], true); mfi.isValid(); ++mfi)
    {
        const amrex::Box &bx = mfi.tilebox();
        Set::Patch<Set::Scalar> alpha_0 = agglom.alpha_0.Patch(lev, mfi);
        Set::Patch<Set::Scalar> alpha = agglom.alpha.Patch(lev, mfi);
        Set::Patch<Set::Scalar> eta = eta_mf.Patch(lev, mfi);

        amrex::BoxList boxes_to_update;
        if (lev == finest_level && agglom.prev_finest_level == finest_level)
            // only update cells on the finest level that were JUST refined
            boxes_to_update = amrex::complementIn(bx, agglom.prev_finest_ba).boxList();
        else
            // this implies we're on level coarser than the
            // finest_level, or there is just now a new
            // finer_level. in either case, all the boxes should be
            // set based on the IC
            boxes_to_update.push_back(bx);

        for (const amrex::Box &box : boxes_to_update)
            amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
                if (eta(i, j, k) > agglom.eta_refinement_threshold)
                    alpha(i, j, k) = alpha_0(i, j, k);
            });
    }

    if (lev == finest_level)
    {
        agglom.prev_finest_ba = grids[finest_level];
        agglom.prev_finest_level = finest_level;
    }
}
}
