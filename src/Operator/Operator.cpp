
#include <AMReX_MLNodeLinOp.H>
#include <AMReX_MLCellLinOp.H>
#include <AMReX_MLNodeLap_K.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Reduce.H>
#include <cmath>
#include <cstdlib>
#include "Util/Color.H"
#include "Set/Set.H"
#include "Operator.H"

using namespace amrex;

#ifdef ALAMO_GPU
#define ALAMO_OPERATOR_FOR amrex::ParallelFor
#define ALAMO_OPERATOR_DEVICE AMREX_GPU_DEVICE
#else
#define ALAMO_OPERATOR_FOR amrex::LoopConcurrentOnCpu
#define ALAMO_OPERATOR_DEVICE
#endif

namespace
{
int AlamoEnvInt(const char* name)
{
    const char* value = std::getenv(name);
    return value ? std::atoi(value) : 0;
}

int AlamoMLCorrectionDiagLimit()
{
    return AlamoEnvInt("ALAMO_ML_COR_DIAG");
}

int AlamoMLResidualDiagLimit()
{
    return AlamoEnvInt("ALAMO_ML_RESID_DIAG");
}

int AlamoMLSpatialDiagLimit()
{
    return AlamoEnvInt("ALAMO_ML_SPATIAL_DIAG");
}

bool AlamoMLZeroCrseGhosts()
{
    return AlamoEnvInt("ALAMO_ML_ZERO_CRSE_GHOSTS") != 0;
}

void AlamoPrintMFNorms(const char* prefix, const char* label, const amrex::MultiFab& mf)
{
    const amrex::IntVect ng(0);
    const int ncomp = mf.nComp();
    const amrex::Real norminf = mf.norminf(0, ncomp, ng);
    const amrex::Real l2 = std::sqrt(amrex::Dot(mf, 0, ncomp, ng));
    const int valid_contains_nan = mf.contains_nan(0, ncomp, 0);
    const int valid_contains_inf = mf.contains_inf(0, ncomp, 0);
    const int grown_contains_nan = mf.contains_nan(0, ncomp, mf.nGrowVect());
    const int grown_contains_inf = mf.contains_inf(0, ncomp, mf.nGrowVect());

    Util::Message(INFO, prefix, label,
        " ncomp=", ncomp,
        " nboxes=", mf.boxArray().size(),
        " ngrow=", mf.nGrowVect(),
        " norminf=", norminf,
        " l2=", l2,
        " valid_contains_nan=", valid_contains_nan,
        " valid_contains_inf=", valid_contains_inf,
        " grown_contains_nan=", grown_contains_nan,
        " grown_contains_inf=", grown_contains_inf,
        " ghost_contains_nan=", (grown_contains_nan && !valid_contains_nan),
        " ghost_contains_inf=", (grown_contains_inf && !valid_contains_inf));

    for (int n = 0; n < ncomp; ++n)
    {
        const amrex::Real comp_norminf = mf.norminf(n, 1, ng);
        const amrex::Real comp_l2 = std::sqrt(amrex::Dot(mf, n, 1, ng));
        Util::Message(INFO, prefix, label,
            " comp=", n,
            " min=", mf.min(n, 0),
            " max=", mf.max(n, 0),
            " norminf=", comp_norminf,
            " l2=", comp_l2);
    }
}

void AlamoZeroGhostsPreserveValid(amrex::MultiFab& mf)
{
    if (mf.nGrowVect() == amrex::IntVect::TheZeroVector()) return;

    amrex::MultiFab valid(mf.boxArray(), mf.DistributionMap(), mf.nComp(), 0);
    amrex::MultiFab::Copy(valid, mf, 0, 0, mf.nComp(), 0);
    mf.setVal(0.0, mf.nGrowVect());
    amrex::MultiFab::Copy(mf, valid, 0, 0, mf.nComp(), 0);
}

void AlamoPrintCorrectionDiagMF(const char* label, const amrex::MultiFab& mf)
{
    AlamoPrintMFNorms("ALAMO_ML_COR_DIAG ", label, mf);
}

void AlamoPrintSpatialDiagMF(const char* label, int call, int amrlev, int mglev,
    const amrex::MultiFab& mf, const amrex::Geometry& geom)
{
    amrex::Box domain = geom.Domain();
    domain.convert(amrex::IntVect::TheNodeVector());

    Util::Message(INFO, "ALAMO_ML_SPATIAL_DIAG call=", call,
        " amrlev=", amrlev,
        " mglev=", mglev,
        " field=", label,
        " domain=", domain,
        " nboxes=", mf.boxArray().size(),
        " ngrow=", mf.nGrowVect());

    const int domlo0 = domain.smallEnd(0);
    const int domhi0 = domain.bigEnd(0);
#if AMREX_SPACEDIM >= 2
    const int domlo1 = domain.smallEnd(1);
    const int domhi1 = domain.bigEnd(1);
#endif
#if AMREX_SPACEDIM >= 3
    const int domlo2 = domain.smallEnd(2);
    const int domhi2 = domain.bigEnd(2);
#endif

    for (int n = 0; n < mf.nComp(); ++n)
    {
        const amrex::Real min_value = mf.min(n, 0);
        const amrex::Real max_value = mf.max(n, 0);
        const amrex::IntVect min_index = mf.minIndex(n, 0);
        const amrex::IntVect max_index = mf.maxIndex(n, 0);

        amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpMax, amrex::ReduceOpMax,
                        amrex::ReduceOpMax, amrex::ReduceOpSum, amrex::ReduceOpSum,
                        amrex::ReduceOpSum, amrex::ReduceOpSum, amrex::ReduceOpSum,
                        amrex::ReduceOpSum> reduce_op;
        amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Real,
                        amrex::Real, amrex::Real, amrex::Real, amrex::Real,
                        amrex::Real, amrex::Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (amrex::MFIter mfi(mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.tilebox();
            const int bxlo0 = bx.smallEnd(0);
            const int bxhi0 = bx.bigEnd(0);
#if AMREX_SPACEDIM >= 2
            const int bxlo1 = bx.smallEnd(1);
            const int bxhi1 = bx.bigEnd(1);
#endif
#if AMREX_SPACEDIM >= 3
            const int bxlo2 = bx.smallEnd(2);
            const int bxhi2 = bx.bigEnd(2);
#endif
            auto const& fab = mf.const_array(mfi);

            reduce_op.eval(bx, reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                {
                    const amrex::Real value = std::abs(fab(i, j, k, n));
                    const bool on_domain_boundary =
                        (i == domlo0 || i == domhi0)
#if AMREX_SPACEDIM >= 2
                        || (j == domlo1 || j == domhi1)
#endif
#if AMREX_SPACEDIM >= 3
                        || (k == domlo2 || k == domhi2)
#endif
                        ;
                    const bool on_box_boundary =
                        (i == bxlo0 || i == bxhi0)
#if AMREX_SPACEDIM >= 2
                        || (j == bxlo1 || j == bxhi1)
#endif
#if AMREX_SPACEDIM >= 3
                        || (k == bxlo2 || k == bxhi2)
#endif
                        ;

                    amrex::Real neighbor_jump = 0.0;
                    if (i > bxlo0) neighbor_jump = amrex::max(neighbor_jump, std::abs(fab(i, j, k, n) - fab(i - 1, j, k, n)));
                    if (i < bxhi0) neighbor_jump = amrex::max(neighbor_jump, std::abs(fab(i + 1, j, k, n) - fab(i, j, k, n)));
#if AMREX_SPACEDIM >= 2
                    if (j > bxlo1) neighbor_jump = amrex::max(neighbor_jump, std::abs(fab(i, j, k, n) - fab(i, j - 1, k, n)));
                    if (j < bxhi1) neighbor_jump = amrex::max(neighbor_jump, std::abs(fab(i, j + 1, k, n) - fab(i, j, k, n)));
#endif
#if AMREX_SPACEDIM >= 3
                    if (k > bxlo2) neighbor_jump = amrex::max(neighbor_jump, std::abs(fab(i, j, k, n) - fab(i, j, k - 1, n)));
                    if (k < bxhi2) neighbor_jump = amrex::max(neighbor_jump, std::abs(fab(i, j, k + 1, n) - fab(i, j, k, n)));
#endif

                    return {
                        on_domain_boundary ? value : 0.0,
                        on_box_boundary ? value : 0.0,
                        (!on_domain_boundary && !on_box_boundary) ? value : 0.0,
                        neighbor_jump,
                        value,
                        on_domain_boundary ? value : 0.0,
                        on_box_boundary ? value : 0.0,
                        1.0,
                        on_domain_boundary ? 1.0 : 0.0,
                        on_box_boundary ? 1.0 : 0.0
                    };
                });
        }

        ReduceTuple hv = reduce_data.value(reduce_op);
        amrex::Real domain_boundary_max = amrex::get<0>(hv);
        amrex::Real box_boundary_max = amrex::get<1>(hv);
        amrex::Real interior_max = amrex::get<2>(hv);
        amrex::Real neighbor_jump_max = amrex::get<3>(hv);
        amrex::Real sum_abs = amrex::get<4>(hv);
        amrex::Real domain_boundary_sum_abs = amrex::get<5>(hv);
        amrex::Real box_boundary_sum_abs = amrex::get<6>(hv);
        amrex::Real count = amrex::get<7>(hv);
        amrex::Real domain_boundary_count = amrex::get<8>(hv);
        amrex::Real box_boundary_count = amrex::get<9>(hv);

        amrex::ParallelDescriptor::ReduceRealMax(domain_boundary_max);
        amrex::ParallelDescriptor::ReduceRealMax(box_boundary_max);
        amrex::ParallelDescriptor::ReduceRealMax(interior_max);
        amrex::ParallelDescriptor::ReduceRealMax(neighbor_jump_max);
        amrex::ParallelDescriptor::ReduceRealSum(sum_abs);
        amrex::ParallelDescriptor::ReduceRealSum(domain_boundary_sum_abs);
        amrex::ParallelDescriptor::ReduceRealSum(box_boundary_sum_abs);
        amrex::ParallelDescriptor::ReduceRealSum(count);
        amrex::ParallelDescriptor::ReduceRealSum(domain_boundary_count);
        amrex::ParallelDescriptor::ReduceRealSum(box_boundary_count);

        const amrex::Real abs_min = std::abs(min_value);
        const amrex::Real abs_max = std::abs(max_value);
        const amrex::Real max_abs = amrex::max(abs_min, abs_max);
        const amrex::IntVect& max_abs_index = (abs_min > abs_max) ? min_index : max_index;
        const amrex::Real inv_sum_abs = (sum_abs > 0.0) ? (1.0 / sum_abs) : 0.0;
        const amrex::Real inv_max_abs = (max_abs > 0.0) ? (1.0 / max_abs) : 0.0;

        Util::Message(INFO, "ALAMO_ML_SPATIAL_DIAG ", label,
            " comp=", n,
            " min=", min_value,
            " min_index=", min_index,
            " max=", max_value,
            " max_index=", max_index,
            " max_abs=", max_abs,
            " max_abs_index=", max_abs_index,
            " domain_boundary_max_abs=", domain_boundary_max,
            " box_boundary_max_abs=", box_boundary_max,
            " interior_max_abs=", interior_max,
            " neighbor_jump_max=", neighbor_jump_max,
            " neighbor_jump_over_max_abs=", neighbor_jump_max * inv_max_abs,
            " sum_abs=", sum_abs,
            " domain_boundary_sum_abs_frac=", domain_boundary_sum_abs * inv_sum_abs,
            " box_boundary_sum_abs_frac=", box_boundary_sum_abs * inv_sum_abs,
            " count=", count,
            " domain_boundary_count=", domain_boundary_count,
            " box_boundary_count=", box_boundary_count);
    }
}

void AlamoPrintResidualDiagMF(const char* label, int call, int amrlev, int mglev,
    const amrex::MultiFab& mf)
{
    Util::Message(INFO, "ALAMO_ML_RESID_DIAG call=", call,
        " amrlev=", amrlev,
        " mglev=", mglev,
        " field=", label);
    AlamoPrintMFNorms("ALAMO_ML_RESID_DIAG ", label, mf);
}
}

namespace Operator {

// constexpr amrex::IntVect AMREX_D_DECL(Operator<Grid::Node>::dx,Operator<Grid::Node>::dy,Operator<Grid::Node>::dz);
constexpr amrex::IntVect AMREX_D_DECL(Operator<Grid::Cell>::dx, Operator<Grid::Cell>::dy, Operator<Grid::Cell>::dz);

void Operator<Grid::Node>::Diagonal(bool recompute)
{
    BL_PROFILE(Color::FG::Yellow + "Operator::Diagonal()" + Color::Reset);
    if (!recompute && m_diagonal_computed) return;
    m_diagonal_computed = true;

    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            Diagonal(amrlev, mglev, *m_diag[amrlev][mglev]);
        }
    }
}

void Operator<Grid::Node>::Diagonal(int amrlev, int mglev, amrex::MultiFab& diag)
{
    BL_PROFILE("Operator::Diagonal()");
    //Util::Message(INFO);

    int ncomp = diag.nComp();
    int nghost = 0;

    int sep = 2;
    int num = AMREX_D_TERM(sep, *sep, *sep);
    int cntr = 0;

    amrex::MultiFab x(m_diag[amrlev][mglev]->boxArray(), m_diag[amrlev][mglev]->DistributionMap(), ncomp, nghost);
    amrex::MultiFab Ax(m_diag[amrlev][mglev]->boxArray(), m_diag[amrlev][mglev]->DistributionMap(), ncomp, nghost);

    for (MFIter mfi(x, false); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();
        amrex::FArrayBox& diagfab = diag[mfi];
        amrex::FArrayBox& xfab = x[mfi];
        amrex::FArrayBox& Axfab = Ax[mfi];

        diagfab.setVal<amrex::RunOn::Device>(0.0);

        for (int i = 0; i < num; i++)
        {
            for (int n = 0; n < ncomp; n++)
            {
                xfab.setVal<amrex::RunOn::Device>(0.0);
                Axfab.setVal<amrex::RunOn::Device>(0.0);

                //BL_PROFILE_VAR("Operator::Part1", part1); 
                AMREX_D_TERM(for (int m1 = bx.loVect()[0]; m1 <= bx.hiVect()[0]; m1++),
                    for (int m2 = bx.loVect()[1]; m2 <= bx.hiVect()[1]; m2++),
                        for (int m3 = bx.loVect()[2]; m3 <= bx.hiVect()[2]; m3++))
                {
                    amrex::IntVect m(AMREX_D_DECL(m1, m2, m3));

                    if (m1 % sep == i / sep && m2 % sep == i % sep) xfab(m, n) = 1.0;
                    else xfab(m, n) = 0.0;
                }
                //BL_PROFILE_VAR_STOP(part1);

                BL_PROFILE_VAR("Operator::Part2", part2);
                Util::Message(INFO, "Calling fapply...", cntr++);
                Fapply(amrlev, mglev, Ax, x);
                BL_PROFILE_VAR_STOP(part2);

                //BL_PROFILE_VAR("Operator::Part3", part3); 
                Axfab.mult<amrex::RunOn::Device>(xfab, n, n, 1);
                diagfab.plus<amrex::RunOn::Device>(Axfab, n, n, 1);
                //BL_PROFILE_VAR_STOP(part3);
            }
        }
    }
}

void Operator<Grid::Node>::Fsmooth(int amrlev, int mglev, amrex::MultiFab& x, const amrex::MultiFab& b) const
{
    BL_PROFILE("Operator::Fsmooth()");

    amrex::Box domain(m_geom[amrlev][mglev].growPeriodicDomain(1));
    domain.convert(amrex::IntVect::TheNodeVector());

    int ncomp = b.nComp();
    int nghost = 2; //b.nGrow();


    amrex::MultiFab Ax(x.boxArray(), x.DistributionMap(), ncomp, nghost);
    amrex::MultiFab Dx(x.boxArray(), x.DistributionMap(), ncomp, nghost);
    amrex::MultiFab Rx(x.boxArray(), x.DistributionMap(), ncomp, nghost);

    if (!m_diagonal_computed) Util::Abort(INFO, "Operator::Diagonal() must be called before using Fsmooth");

    // This is a JACOBI iteration, not Gauss-Seidel.
    // So we need to do twice the number of iterations to get the same behavior as GS.
    for (int ctr = 0; ctr < 2; ctr++)
    {
        Fapply(amrlev, mglev, Ax, x); // find Ax

        amrex::MultiFab::Copy(Dx, x, 0, 0, ncomp, nghost); // Dx = x
        amrex::MultiFab::Multiply(Dx, *m_diag[amrlev][mglev], 0, 0, ncomp, nghost); // Dx *= diag  (Dx = x*diag)

        amrex::MultiFab::Copy(Rx, Ax, 0, 0, ncomp, nghost); // Rx = Ax
        amrex::MultiFab::Subtract(Rx, Dx, 0, 0, ncomp, nghost); // Rx -= Dx  (Rx = Ax - Dx)

        for (MFIter mfi(x, false); mfi.isValid(); ++mfi)
        {
            Box bx = mfi.grownnodaltilebox();
            
            auto xfab = x.array(mfi);
            auto bfab = b.const_array(mfi);
            auto Rxfab = Rx.const_array(mfi);
            auto diagfab = (*m_diag[amrlev][mglev]).const_array(mfi);


            for (int n = 0; n < ncomp; n++)
            {
                auto m_omega = this->m_omega;
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) 
                {

                    // Skip ghost cells outside problem domain
                    if (!domain.contains(i,j,k))
                    {
                        //continue;
                    }
                    else if ( !bx.strictly_contains(i,j,k))
                    {
                        xfab(i, j, k, n) = 0.0;
                        //continue;
                    }
                    else
                    {
                        xfab(i,j,k,n) = (1. - m_omega) * xfab(i,j,k, n) + m_omega * (bfab(i,j,k, n) - Rxfab(i,j,k, n)) / diagfab(i,j,k,n);
                    }
                });
            }
        }
        amrex::Geometry geom = m_geom[amrlev][mglev];
        x.setMultiGhost(true);
        x.FillBoundary(geom.periodicity());
        nodalSync(amrlev, mglev, x);
    }
}

void Operator<Grid::Node>::normalize(int amrlev, int mglev, MultiFab& a_x) const
{
    if (!m_diagonal_computed)
        Util::Abort(INFO, "Operator::Diagonal() must be called before using normalize");

    a_x.divide(*m_diag[amrlev][mglev],0,getNComp(),2);

    a_x.setMultiGhost(true);
    a_x.FillBoundaryAndSync(Geom(amrlev,mglev).periodicity());
}

Operator<Grid::Node>::Operator(const Vector<Geometry>& a_geom,
    const Vector<BoxArray>& a_grids,
    const Vector<DistributionMapping>& a_dmap,
    const LPInfo& a_info,
    const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    BL_PROFILE("Operator::Operator()");
    Util::Message(INFO);

    if (!(a_grids[0].ixType() == amrex::IndexType::TheNodeType()))
        Util::Abort(INFO, "Operator must be defined using CELL CENTERED boxarrays.");

    define(a_geom, a_grids, a_dmap, a_info, a_factory);
}

Operator<Grid::Node>::~Operator()
{}

void Operator<Grid::Node>::define(const Vector<Geometry>& a_geom,
    const Vector<BoxArray>& a_grids,
    const Vector<DistributionMapping>& a_dmap,
    const LPInfo& a_info,
    const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    BL_PROFILE("Operator::~Operator()");

    // Make sure we're not trying to parallelize in vain.
    if (amrex::ParallelDescriptor::NProcs() > a_grids[0].size())
    {
        Util::Warning(INFO, "There are more processors than there are boxes in the amrlev=0 boxarray!!\n",
            "(NProcs = ", amrex::ParallelDescriptor::NProcs(), ", a_grids[0].size() = ", a_grids[0].size(), ")\n",
            "You should decrease max_grid_size or you will not get proper scaling!");
    }

    // This makes sure grids are node-centered;
    Vector<BoxArray> cc_grids = a_grids;
    for (auto& ba : cc_grids) {
        ba.enclosedCells();
    }

    MLNodeLinOp::define(a_geom, a_grids, a_dmap, a_info, a_factory);

    int nghost = 2;
    // Resize the multifab containing the operator diagonal
    m_diag.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_diag[amrlev].resize(m_num_mg_levels[amrlev]);

        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            m_diag[amrlev][mglev].reset(new MultiFab(amrex::convert(m_grids[amrlev][mglev], amrex::IntVect::TheNodeVector()),
                m_dmap[amrlev][mglev], getNComp(), nghost));
        }
    }

    // We need to instantiate the m_lobc objects.
    // WE DO NOT USE THEM - our BCs are implemented differently.
    // But they need to be the right size or the code will segfault.
    m_lobc.resize(getNComp(), { {AMREX_D_DECL(BCType::bogus,BCType::bogus,BCType::bogus)} });
    m_hibc.resize(getNComp(), { {AMREX_D_DECL(BCType::bogus,BCType::bogus,BCType::bogus)} });
}

void Operator<Grid::Node>::fixUpResidualMask(int amrlev, iMultiFab& resmsk)
{
    BL_PROFILE("Operator::fixUpResidualMask()");

    if (!m_masks_built) buildMasks();

    const iMultiFab& cfmask = *m_nd_fine_mask[amrlev];

#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(resmsk, TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        Array4<int> const& rmsk = resmsk.array(mfi);
        Array4<int const> const& fmsk = cfmask.const_array(mfi);
        AMREX_HOST_DEVICE_PARALLEL_FOR_3D(bx, i, j, k,
            {
                if (fmsk(i,j,k) == amrex::nodelap_detail::crse_fine_node) rmsk(i,j,k) = 1;
            });
    }
}

void Operator<Grid::Node>::prepareForSolve()
{
    BL_PROFILE("Operator::prepareForSolve()");
    MLNodeLinOp::prepareForSolve();
    buildMasks();
    averageDownCoeffs();
    Diagonal(true);
}

void Operator<Grid::Node>::restriction(int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const
{
    BL_PROFILE("Operator::restriction()");

    applyBC(amrlev, cmglev - 1, fine, BCMode::Homogeneous, StateMode::Solution);

    amrex::Box cdomain = m_geom[amrlev][cmglev].growPeriodicDomain(1);
    cdomain = cdomain.convert(amrex::IntVect::TheNodeVector());

    bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
    MultiFab cfine;
    if (need_parallel_copy) {
        const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
        cfine.define(ba, fine.DistributionMap(), fine.nComp(), fine.nGrow());
    }

    MultiFab* pcrse = (need_parallel_copy) ? &cfine : &crse;

    for (MFIter mfi(*pcrse, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.grownnodaltilebox(-1,1) & cdomain;

        amrex::Array4<const amrex::Real> const& fdata = fine.array(mfi);
        amrex::Array4<amrex::Real> const& cdata = pcrse->array(mfi);

        const Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);


        for (int n = 0; n < crse.nComp(); n++)
        {
            // I,J,K == coarse coordinates
            // i,j,k == fine coordinates
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int I, int J, int K) {
                int i = 2 * I, j = 2 * J, k = 2 * K;

#if AMREX_SPACEDIM == 2
                if ((I == lo.x || I == hi.x) && (J == lo.y || J == hi.y)) // Corner
                {
                    cdata(I, J, K, n) = fdata(i, j, k, n);
                }
                else if (J == lo.y || J == hi.y) // Y boundary
                {
                    cdata(I, J, K, n) = 0.25 * fdata(i - 1, j, k, n) + 0.5 * fdata(i, j, k, n) + 0.25 * fdata(i + 1, j, k, n);
                }
                else if (I == lo.x || I == hi.x) // X boundary
                {
                    cdata(I, J, K, n) = 0.25 * fdata(i, j - 1, k, n) + 0.5 * fdata(i, j, k, n) + 0.25 * fdata(i, j + 1, k, n);
                }
                else // Interior
                {
                    cdata(I, J, K, n) =
                        (+fdata(i - 1, j - 1, k, n) + 2.0 * fdata(i, j - 1, k, n) + fdata(i + 1, j - 1, k, n)
                            + 2.0 * fdata(i - 1, j, k, n) + 4.0 * fdata(i, j, k, n) + 2.0 * fdata(i + 1, j, k, n)
                            + fdata(i - 1, j + 1, k, n) + 2.0 * fdata(i, j + 1, k, n) + fdata(i + 1, j + 1, k, n)) / 16.0;
                }
#else
                if ((I == lo.x || I == hi.x) &&
                    (J == lo.y || J == hi.y) &&
                    (K == lo.z || K == hi.z)) // Corner
                {
                    cdata(I, J, K, n) = fdata(i, j, k, n);
                }
                else if ((J == lo.y || J == hi.y) &&
                    (K == lo.z || K == hi.z)) // X edge
                {
                    cdata(I, J, K, n) = 0.25 * fdata(i - 1, j, k, n) + 0.5 * fdata(i, j, k, n) + 0.25 * fdata(i + 1, j, k, n);
                }
                else if ((K == lo.z || K == hi.z) &&
                    (I == lo.x || I == hi.x)) // Y edge
                {
                    cdata(I, J, K, n) = 0.25 * fdata(i, j - 1, k, n) + 0.5 * fdata(i, j, k, n) + 0.25 * fdata(i, j + 1, k, n);
                }
                else if ((I == lo.x || I == hi.x) &&
                    (J == lo.y || J == hi.y)) // Z edge
                {
                    cdata(I, J, K, n) = 0.25 * fdata(i, j, k - 1, n) + 0.5 * fdata(i, j, k, n) + 0.25 * fdata(i, j, k + 1, n);
                }
                else if (I == lo.x || I == hi.x) // X face
                {
                    cdata(I, J, K, n) =
                        (+fdata(i, j - 1, k - 1, n) + 2.0 * fdata(i, j, k - 1, n) + fdata(i, j + 1, k - 1, n)
                            + 2.0 * fdata(i, j - 1, k, n) + 4.0 * fdata(i, j, k, n) + 2.0 * fdata(i, j + 1, k, n)
                            + fdata(i, j - 1, k + 1, n) + 2.0 * fdata(i, j, k + 1, n) + fdata(i, j + 1, k + 1, n)) / 16.0;
                }
                else if (J == lo.y || J == hi.y) // Y face
                {
                    cdata(I, J, K, n) =
                        (+fdata(i - 1, j, k - 1, n) + 2.0 * fdata(i - 1, j, k, n) + fdata(i - 1, j, k + 1, n)
                            + 2.0 * fdata(i, j, k - 1, n) + 4.0 * fdata(i, j, k, n) + 2.0 * fdata(i, j, k + 1, n)
                            + fdata(i + 1, j, k - 1, n) + 2.0 * fdata(i + 1, j, k, n) + fdata(i + 1, j, k + 1, n)) / 16.0;
                }
                else if (K == lo.z || K == hi.z) // Z face
                {
                    cdata(I, J, K, n) =
                        (+fdata(i - 1, j - 1, k, n) + 2.0 * fdata(i, j - 1, k, n) + fdata(i + 1, j - 1, k, n)
                            + 2.0 * fdata(i - 1, j, k, n) + 4.0 * fdata(i, j, k, n) + 2.0 * fdata(i + 1, j, k, n)
                            + fdata(i - 1, j + 1, k, n) + 2.0 * fdata(i, j + 1, k, n) + fdata(i + 1, j + 1, k, n)) / 16.0;
                }
                else // Interior
                    cdata(I, J, K, n) =
                    (fdata(i - 1, j - 1, k - 1, n) + fdata(i - 1, j - 1, k + 1, n) + fdata(i - 1, j + 1, k - 1, n) + fdata(i - 1, j + 1, k + 1, n) +
                        fdata(i + 1, j - 1, k - 1, n) + fdata(i + 1, j - 1, k + 1, n) + fdata(i + 1, j + 1, k - 1, n) + fdata(i + 1, j + 1, k + 1, n)) / 64.0
                    +
                    (fdata(i, j - 1, k - 1, n) + fdata(i, j - 1, k + 1, n) + fdata(i, j + 1, k - 1, n) + fdata(i, j + 1, k + 1, n) +
                        fdata(i - 1, j, k - 1, n) + fdata(i + 1, j, k - 1, n) + fdata(i - 1, j, k + 1, n) + fdata(i + 1, j, k + 1, n) +
                        fdata(i - 1, j - 1, k, n) + fdata(i - 1, j + 1, k, n) + fdata(i + 1, j - 1, k, n) + fdata(i + 1, j + 1, k, n)) / 32.0
                    +
                    (fdata(i - 1, j, k, n) + fdata(i, j - 1, k, n) + fdata(i, j, k - 1, n) +
                        fdata(i + 1, j, k, n) + fdata(i, j + 1, k, n) + fdata(i, j, k + 1, n)) / 16.0
                    +
                    fdata(i, j, k, n) / 8.0;
#endif
            });
        }
    }

    if (need_parallel_copy) {
        crse.ParallelCopy(cfine);
    }

    crse.setMultiGhost(true);
    crse.FillBoundary(Geom(amrlev,cmglev).periodicity());
    nodalSync(amrlev, cmglev, crse);
}

void Operator<Grid::Node>::interpolation(int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const
{
    BL_PROFILE("Operator::interpolation()");
    amrex::Box fdomain = m_geom[amrlev][fmglev].growPeriodicDomain(2);
    fdomain.convert(amrex::IntVect::TheNodeVector());

    static int diag_call = 0;
    const int diag_limit = AlamoMLCorrectionDiagLimit();
    const int spatial_diag_limit = AlamoMLSpatialDiagLimit();
    const bool diag_enabled = diag_call < diag_limit;
    const bool spatial_diag_enabled = diag_call < spatial_diag_limit;
    if (diag_enabled || spatial_diag_enabled)
    {
        Util::Message(INFO, "ALAMO_ML_COR_DIAG interpolation_call=", diag_call,
            " amrlev=", amrlev,
            " fmglev=", fmglev,
            " crse_mglev=", fmglev + 1);
    }
    if (diag_enabled)
    {
        AlamoPrintCorrectionDiagMF("crse_after_bottom", crse);
        AlamoPrintCorrectionDiagMF("fine_before_interp", fine);
    }
    if (spatial_diag_enabled)
    {
        AlamoPrintSpatialDiagMF("crse_after_bottom", diag_call, amrlev, fmglev + 1, crse, m_geom[amrlev][fmglev + 1]);
        AlamoPrintSpatialDiagMF("fine_before_interp", diag_call, amrlev, fmglev, fine, m_geom[amrlev][fmglev]);
    }

    MultiFab sanitized_crse;
    const MultiFab* source_crse = &crse;
    if (AlamoMLZeroCrseGhosts())
    {
        sanitized_crse.define(crse.boxArray(), crse.DistributionMap(),
            crse.nComp(), crse.nGrowVect());
        MultiFab::Copy(sanitized_crse, crse, 0, 0, crse.nComp(), crse.nGrowVect());
        AlamoZeroGhostsPreserveValid(sanitized_crse);
        source_crse = &sanitized_crse;

        Util::Message(INFO, "ALAMO_ML_COR_DIAG zeroed coarse correction ghosts before interpolation",
            " amrlev=", amrlev,
            " fmglev=", fmglev,
            " crse_mglev=", fmglev + 1);
        if (diag_enabled)
            AlamoPrintCorrectionDiagMF("crse_after_bottom_ghosts_zeroed", *source_crse);
        if (spatial_diag_enabled)
            AlamoPrintSpatialDiagMF("crse_after_bottom_ghosts_zeroed", diag_call, amrlev, fmglev + 1, *source_crse, m_geom[amrlev][fmglev + 1]);
    }

    bool need_parallel_copy = !amrex::isMFIterSafe(*source_crse, fine);
    MultiFab cfine;
    const MultiFab* cmf = source_crse;
    if (need_parallel_copy) {
        const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
        cfine.define(ba, fine.DistributionMap(), source_crse->nComp(), source_crse->nGrow());
        cfine.ParallelCopy(*source_crse);
        cmf = &cfine;
    }

    for (MFIter mfi(fine, false); mfi.isValid(); ++mfi)
    {
        Box fine_bx = mfi.validbox() & fdomain;

        const Box& course_bx = amrex::coarsen(fine_bx, 2);
        const Box& tmpbx = amrex::refine(course_bx, 2);
        FArrayBox tmpfab;
        tmpfab.resize(tmpbx, fine.nComp());
        // GPU FIX (root cause of the multi-box MLMG divergence): this per-box
        // temporary lives only within one MFIter iteration, but AMReX's non-tiled
        // GPU MFIter cycles through a pool of CUDA streams across iterations. Without
        // an Elixir, tmpfab's device memory is freed at the end of this iteration and
        // can be reused by a later iteration running on a *different* stream while the
        // interpolation/plus kernels below are still in flight -> a cross-stream race
        // that corrupts the interpolated correction (deterministically wrong on CPU's
        // single stream = never; nondeterministically wrong on GPU). The garbage
        // correction is then amplified by the BC-penalty diagonal into a 1e20 blow-up.
        // The Elixir keeps tmpfab alive until its own stream has finished. (restriction
        // writes its result directly, with no temp, which is why it is unaffected.)
        amrex::Gpu::Elixir tmpfab_eli = tmpfab.elixir();
        tmpfab.setVal<amrex::RunOn::Device>(0.0);
        const amrex::FArrayBox& crsefab = (*cmf)[mfi];

        amrex::Array4<const amrex::Real> const& cdata = crsefab.const_array();
        amrex::Array4<amrex::Real> const& fdata = tmpfab.array();

        for (int n = 0; n < crse.nComp(); n++)
        {
            // I,J,K == coarse coordinates
            // i,j,k == fine coordinates
            ALAMO_OPERATOR_FOR(fine_bx, [=] ALAMO_OPERATOR_DEVICE (int i, int j, int k) {

                int I = i / 2, J = j / 2, K = k / 2;

#if AMREX_SPACEDIM == 2
                if (i % 2 == 0 && j % 2 == 0) // Coincident
                    fdata(i, j, k, n) = cdata(I, J, K, n);
                else if (j % 2 == 0) // X edge
                    fdata(i, j, k, n) = 0.5 * (cdata(I, J, K, n) + cdata(I + 1, J, K, n));
                else if (i % 2 == 0) // Y edge
                    fdata(i, j, k, n) = 0.5 * (cdata(I, J, K, n) + cdata(I, J + 1, K, n));
                else // Center
                    fdata(i, j, k, n) = 0.25 * (cdata(I, J, K, n) + cdata(I + 1, J, K, n) +
                        cdata(I, J + 1, K, n) + cdata(I + 1, J + 1, K, n));
#else
                if (i % 2 == 0 && j % 2 == 0 && k % 2 == 0) // Coincident
                    fdata(i, j, k, n) = cdata(I, J, K, n);
                else if (j % 2 == 0 && k % 2 == 0) // X Edge
                    fdata(i, j, k, n) = 0.5 * (cdata(I, J, K, n) + cdata(I + 1, J, K, n));
                else if (k % 2 == 0 && i % 2 == 0) // Y Edge
                    fdata(i, j, k, n) = 0.5 * (cdata(I, J, K, n) + cdata(I, J + 1, K, n));
                else if (i % 2 == 0 && j % 2 == 0) // Z Edge
                    fdata(i, j, k, n) = 0.5 * (cdata(I, J, K, n) + cdata(I, J, K + 1, n));
                else if (i % 2 == 0) // X Face
                    fdata(i, j, k, n) = 0.25 * (cdata(I, J, K, n) + cdata(I, J + 1, K, n) +
                        cdata(I, J, K + 1, n) + cdata(I, J + 1, K + 1, n));
                else if (j % 2 == 0) // Y Face
                    fdata(i, j, k, n) = 0.25 * (cdata(I, J, K, n) + cdata(I, J, K + 1, n) +
                        cdata(I + 1, J, K, n) + cdata(I + 1, J, K + 1, n));
                else if (k % 2 == 0) // Z Face
                    fdata(i, j, k, n) = 0.25 * (cdata(I, J, K, n) + cdata(I + 1, J, K, n) +
                        cdata(I, J + 1, K, n) + cdata(I + 1, J + 1, K, n));
                else // Center
                    fdata(i, j, k, n) = 0.125 * (cdata(I, J, K, n) +
                        cdata(I + 1, J, K, n) + cdata(I, J + 1, K, n) + cdata(I, J, K + 1, n) +
                        cdata(I, J + 1, K + 1, n) + cdata(I + 1, J, K + 1, n) + cdata(I + 1, J + 1, K, n) +
                        cdata(I + 1, J + 1, K + 1, n));
#endif

            });
        }
        fine[mfi].plus<amrex::RunOn::Device>(tmpfab, fine_bx, fine_bx, 0, 0, fine.nComp());
    }

    fine.setMultiGhost(true);
    fine.FillBoundary(Geom(amrlev,fmglev).periodicity());
    nodalSync(amrlev, fmglev, fine);

    if (diag_enabled)
    {
        AlamoPrintCorrectionDiagMF("fine_after_interp_sync", fine);
    }
    if (spatial_diag_enabled)
    {
        AlamoPrintSpatialDiagMF("fine_after_interp_sync", diag_call, amrlev, fmglev, fine, m_geom[amrlev][fmglev]);
    }
    if (diag_enabled || spatial_diag_enabled)
    {
        ++diag_call;
    }
}

void Operator<Grid::Node>::averageDownSolutionRHS(int camrlev, MultiFab& crse_sol, MultiFab& /*crse_rhs*/,
    const MultiFab& fine_sol, const MultiFab& /*fine_rhs*/)
{
    BL_PROFILE("Operator::averageDownSolutionRHS()");
    const auto& amrrr = AMRRefRatio(camrlev);
    amrex::average_down(fine_sol, crse_sol, 0, crse_sol.nComp(), amrrr);

    if (isSingular(0))
    {
        Util::Abort(INFO, "Singular operators not supported!");
    }

}

void Operator<Grid::Node>::realFillBoundary(MultiFab& phi, const Geometry& geom)
{
    Util::RealFillBoundary(phi, geom);
}

void Operator<Grid::Node>::applyBC(int amrlev, int mglev, MultiFab& phi, BCMode/* bc_mode*/,
    amrex::MLLinOp::StateMode /**/, bool skip_fillboundary) const
{
    BL_PROFILE("Operator::applyBC()");

    const Geometry& geom = m_geom[amrlev][mglev];

    if (!skip_fillboundary) {
        //phi.FillBoundary(geom.periodicity());
        //phi.setMultiGhost(true);
        phi.FillBoundaryAndSync(geom.periodicity());
    }
}

const amrex::FArrayBox&
Operator<Grid::Node>::GetFab(const int num, const int amrlev, const int mglev, const amrex::MFIter& mfi) const
{
    BL_PROFILE("Operator::GetFab()");
    Util::Message(INFO);
    return m_a_coeffs[num][amrlev][mglev][mfi];
}

void Operator<Grid::Node>::RegisterNewFab(amrex::Vector<amrex::MultiFab>& input)
{
    BL_PROFILE("Operator::RegisterNewFab()");
    Util::Message(INFO);
    /// \todo assertions here
    m_a_coeffs.resize(m_a_coeffs.size() + 1);
    m_a_coeffs[m_num_a_fabs].resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_a_coeffs[m_num_a_fabs][amrlev].resize(m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
            m_a_coeffs[m_num_a_fabs][amrlev][mglev].define(m_grids[amrlev][mglev],
                m_dmap[amrlev][mglev],
                input[amrlev].nComp(),
                input[amrlev].nGrow());

        amrex::MultiFab::Copy(m_a_coeffs[m_num_a_fabs][amrlev][0],
            input[amrlev], 0, 0,
            input[amrlev].nComp(),
            input[amrlev].nGrow());
    }
    m_num_a_fabs++;
}


void Operator<Grid::Node>::RegisterNewFab(amrex::Vector<std::unique_ptr<amrex::MultiFab> >& input)
{
    BL_PROFILE("Operator::RegisterNewFab()");
    Util::Message(INFO);
    /// \todo assertions here
    m_a_coeffs.resize(m_a_coeffs.size() + 1);
    m_a_coeffs[m_num_a_fabs].resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_a_coeffs[m_num_a_fabs][amrlev].resize(m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
            m_a_coeffs[m_num_a_fabs][amrlev][mglev].define(m_grids[amrlev][mglev],
                m_dmap[amrlev][mglev],
                input[amrlev]->nComp(),
                input[amrlev]->nGrow());

        amrex::MultiFab::Copy(m_a_coeffs[m_num_a_fabs][amrlev][0],
            *input[amrlev], 0, 0,
            input[amrlev]->nComp(),
            input[amrlev]->nGrow());
    }
    m_num_a_fabs++;
}

void Operator<Grid::Node>::reflux(int crse_amrlev,
    MultiFab& res, const MultiFab& /*crse_sol*/, const MultiFab& /*crse_rhs*/,
    MultiFab& fine_res, MultiFab& /*fine_sol*/, const MultiFab& /*fine_rhs*/) const
{
    BL_PROFILE("Operator::Elastic::reflux()");

    int ncomp = AMREX_SPACEDIM;

    amrex::Box cdomain(m_geom[crse_amrlev][0].growPeriodicDomain(2));
    cdomain.convert(amrex::IntVect::TheNodeVector());

    const Geometry& cgeom = m_geom[crse_amrlev][0];

    const BoxArray& fba = fine_res.boxArray();
    const DistributionMapping& fdm = fine_res.DistributionMap();

    MultiFab fine_res_for_coarse(amrex::coarsen(fba, 2), fdm, ncomp, 2);
    fine_res_for_coarse.ParallelCopy(res, 0, 0, ncomp, 0, 0, cgeom.periodicity());

    applyBC(crse_amrlev + 1, 0, fine_res, BCMode::Inhomogeneous, StateMode::Solution);

    /// \todo Replace with Enum
    // const int coarse_coarse_node = 0;
    const int coarse_fine_node = 1;
    const int fine_fine_node = 2;

    amrex::iMultiFab nodemask(amrex::coarsen(fba, 2), fdm, 1, 2);
    nodemask.ParallelCopy(*m_nd_fine_mask[crse_amrlev], 0, 0, 1, 0, 0, cgeom.periodicity());

    for (MFIter mfi(fine_res_for_coarse, false); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.grownnodaltilebox(-1,1) & cdomain;

        amrex::Array4<const int> const& nmask = nodemask.array(mfi);
        //amrex::Array4<const int> const& cmask = cellmask.array(mfi);

        amrex::Array4<amrex::Real> const& cdata = fine_res_for_coarse.array(mfi);
        amrex::Array4<const amrex::Real> const& fdata = fine_res.array(mfi);

        const Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);

        for (int n = 0; n < fine_res.nComp(); n++)
        {
            // I,J,K == coarse coordinates
            // i,j,k == fine coordinates
            ALAMO_OPERATOR_FOR(bx, [=] ALAMO_OPERATOR_DEVICE (int I, int J, int K) {
                int i = I * 2, j = J * 2, k = K * 2;

                if (nmask(I, J, K) == fine_fine_node || nmask(I, J, K) == coarse_fine_node)
                {
                    if ((I == lo.x || I == hi.x) &&
                        (J == lo.y || J == hi.y) &&
                        (K == lo.z || K == hi.z)) // Corner
                        cdata(I, J, K, n) = fdata(i, j, k, n);
                    else if ((J == lo.y || J == hi.y) &&
                        (K == lo.z || K == hi.z)) // X edge
                        cdata(I, J, K, n) = 0.25 * fdata(i - 1, j, k, n) + 0.5 * fdata(i, j, k, n) + 0.25 * fdata(i + 1, j, k, n);
                    else if ((K == lo.z || K == hi.z) &&
                        (I == lo.x || I == hi.x)) // Y edge
                        cdata(I, J, K, n) = 0.25 * fdata(i, j - 1, k, n) + 0.5 * fdata(i, j, k, n) + 0.25 * fdata(i, j + 1, k, n);
                    else if ((I == lo.x || I == hi.x) &&
                        (J == lo.y || J == hi.y)) // Z edge
                        cdata(I, J, K, n) = 0.25 * fdata(i, j, k - 1, n) + 0.5 * fdata(i, j, k, n) + 0.25 * fdata(i, j, k + 1, n);
                    else if (I == lo.x || I == hi.x) // X face
                        cdata(I, J, K, n) =
                        (+fdata(i, j - 1, k - 1, n) + 2.0 * fdata(i, j, k - 1, n) + fdata(i, j + 1, k - 1, n)
                            + 2.0 * fdata(i, j - 1, k, n) + 4.0 * fdata(i, j, k, n) + 2.0 * fdata(i, j + 1, k, n)
                            + fdata(i, j - 1, k + 1, n) + 2.0 * fdata(i, j, k + 1, n) + fdata(i, j + 1, k + 1, n)) / 16.0;
                    else if (J == lo.y || J == hi.y) // Y face
                        cdata(I, J, K, n) =
                        (+fdata(i - 1, j, k - 1, n) + 2.0 * fdata(i - 1, j, k, n) + fdata(i - 1, j, k + 1, n)
                            + 2.0 * fdata(i, j, k - 1, n) + 4.0 * fdata(i, j, k, n) + 2.0 * fdata(i, j, k + 1, n)
                            + fdata(i + 1, j, k - 1, n) + 2.0 * fdata(i + 1, j, k, n) + fdata(i + 1, j, k + 1, n)) / 16.0;
                    else if (K == lo.z || K == hi.z) // Z face
                        cdata(I, J, K, n) =
                        (+fdata(i - 1, j - 1, k, n) + 2.0 * fdata(i, j - 1, k, n) + fdata(i + 1, j - 1, k, n)
                            + 2.0 * fdata(i - 1, j, k, n) + 4.0 * fdata(i, j, k, n) + 2.0 * fdata(i + 1, j, k, n)
                            + fdata(i - 1, j + 1, k, n) + 2.0 * fdata(i, j + 1, k, n) + fdata(i + 1, j + 1, k, n)) / 16.0;
                    else // Interior
                        cdata(I, J, K, n) =
                        (fdata(i - 1, j - 1, k - 1, n) + fdata(i - 1, j - 1, k + 1, n) + fdata(i - 1, j + 1, k - 1, n) + fdata(i - 1, j + 1, k + 1, n) +
                            fdata(i + 1, j - 1, k - 1, n) + fdata(i + 1, j - 1, k + 1, n) + fdata(i + 1, j + 1, k - 1, n) + fdata(i + 1, j + 1, k + 1, n)) / 64.0
                        +
                        (fdata(i, j - 1, k - 1, n) + fdata(i, j - 1, k + 1, n) + fdata(i, j + 1, k - 1, n) + fdata(i, j + 1, k + 1, n) +
                            fdata(i - 1, j, k - 1, n) + fdata(i + 1, j, k - 1, n) + fdata(i - 1, j, k + 1, n) + fdata(i + 1, j, k + 1, n) +
                            fdata(i - 1, j - 1, k, n) + fdata(i - 1, j + 1, k, n) + fdata(i + 1, j - 1, k, n) + fdata(i + 1, j + 1, k, n)) / 32.0
                        +
                        (fdata(i - 1, j, k, n) + fdata(i, j - 1, k, n) + fdata(i, j, k - 1, n) +
                            fdata(i + 1, j, k, n) + fdata(i, j + 1, k, n) + fdata(i, j, k + 1, n)) / 16.0
                        +
                        fdata(i, j, k, n) / 8.0;
                }

            });
        }
    }

    // Copy the fine residual restricted onto the coarse grid
    // into the final residual.
    res.ParallelCopy(fine_res_for_coarse, 0, 0, ncomp, 0, 0, cgeom.periodicity());

    // Sync up ghost nodes
    res.setMultiGhost(true);
    res.FillBoundaryAndSync(Geom(crse_amrlev).periodicity());
    nodalSync(crse_amrlev, 0, res);

    return;
}

void
Operator<Grid::Node>::solutionResidual(int amrlev, MultiFab& resid, MultiFab& x, const MultiFab& b,
    const MultiFab* /*crse_bcdata*/)
{
    const int mglev = 0;
    const int ncomp = b.nComp();
    apply(amrlev, mglev, resid, x, BCMode::Inhomogeneous, StateMode::Solution);
    MultiFab::Xpay(resid, -1.0, b, 0, 0, ncomp, 2);
    resid.setMultiGhost(true);
    resid.FillBoundaryAndSync(Geom(amrlev).periodicity());
}

void
Operator<Grid::Node>::correctionResidual(int amrlev, int mglev, MultiFab& resid, MultiFab& x, const MultiFab& b,
    BCMode /*bc_mode*/, const MultiFab* /*crse_bcdata*/)
{
    static int resid_diag_call = 0;
    const int resid_diag_limit = AlamoMLResidualDiagLimit();
    const int spatial_diag_limit = AlamoMLSpatialDiagLimit();
    const bool resid_diag_enabled = resid_diag_call < resid_diag_limit;
    const bool spatial_diag_enabled = resid_diag_call < spatial_diag_limit;
    const int this_resid_diag_call = resid_diag_call;

    if (resid_diag_enabled)
    {
        Util::Message(INFO, "ALAMO_ML_RESID_DIAG correctionResidual_begin call=", this_resid_diag_call,
            " amrlev=", amrlev,
            " mglev=", mglev,
            " resid_ngrow=", resid.nGrow(),
            " x_ngrow=", x.nGrow(),
            " b_ngrow=", b.nGrow());
        AlamoPrintResidualDiagMF("correction_x_before_apply", this_resid_diag_call, amrlev, mglev, x);
        AlamoPrintResidualDiagMF("correction_b_rhs", this_resid_diag_call, amrlev, mglev, b);
    }
    if (spatial_diag_enabled)
    {
        AlamoPrintSpatialDiagMF("correction_x_before_apply", this_resid_diag_call, amrlev, mglev, x, m_geom[amrlev][mglev]);
    }

    resid.setVal(0.0);
    apply(amrlev, mglev, resid, x, BCMode::Homogeneous, StateMode::Correction);

    if (resid_diag_enabled)
    {
        AlamoPrintResidualDiagMF("correction_Ax_after_apply", this_resid_diag_call, amrlev, mglev, resid);
    }
    if (spatial_diag_enabled)
    {
        AlamoPrintSpatialDiagMF("correction_Ax_after_apply", this_resid_diag_call, amrlev, mglev, resid, m_geom[amrlev][mglev]);
    }

    int ncomp = b.nComp();
    MultiFab::Xpay(resid, -1.0, b, 0, 0, ncomp, resid.nGrow());

    if (resid_diag_enabled)
    {
        AlamoPrintResidualDiagMF("correction_resid_after_xpay", this_resid_diag_call, amrlev, mglev, resid);
    }

    resid.setMultiGhost(true);
    resid.FillBoundaryAndSync(Geom(amrlev).periodicity());

    if (resid_diag_enabled)
    {
        AlamoPrintResidualDiagMF("correction_resid_after_sync", this_resid_diag_call, amrlev, mglev, resid);
    }
    if (resid_diag_enabled || spatial_diag_enabled)
    {
        ++resid_diag_call;
    }
}




Operator<Grid::Cell>::Operator()
{
    m_ixtype = amrex::IntVect::TheCellVector();
}

void
Operator<Grid::Cell>::define(amrex::Vector<amrex::Geometry> a_geom,
    const amrex::Vector<amrex::BoxArray>& a_grids,
    const amrex::Vector<amrex::DistributionMapping>& a_dmap,
    BC::BC<Set::Scalar>& a_bc,
    const amrex::LPInfo& a_info,
    const amrex::Vector<amrex::FabFactory<amrex::FArrayBox> const*>& a_factory)
{
    m_bc = &a_bc;

    std::array<int, AMREX_SPACEDIM> is_periodic = m_bc->IsPeriodic();

    MLCellLinOp::define(a_geom, a_grids, a_dmap, a_info, a_factory);

    Util::Warning(INFO, "This section of code has not been tested.");
    for (int n = 0; n < getNComp(); n++)
    {
        m_lobc.push_back({ AMREX_D_DECL(is_periodic[0] ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Dirichlet,
                        is_periodic[1] ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Dirichlet,
                        is_periodic[2] ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Dirichlet) });
        m_hibc.push_back({ AMREX_D_DECL(is_periodic[0] ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Dirichlet,
                        is_periodic[1] ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Dirichlet,
                        is_periodic[2] ? amrex::LinOpBCType::Periodic : amrex::LinOpBCType::Dirichlet) });
    }

    for (int ilev = 0; ilev < a_geom.size(); ++ilev)
        setLevelBC(ilev, nullptr);

}


void
Operator<Grid::Cell>::prepareForSolve()
{
    MLCellLinOp::prepareForSolve();
}

Operator<Grid::Cell>::BndryCondLoc::BndryCondLoc(const amrex::BoxArray& ba, const amrex::DistributionMapping& dm)
    : bcond(ba, dm),
    bcloc(ba, dm)
{
}

void
Operator<Grid::Cell>::BndryCondLoc::setLOBndryConds(const amrex::Geometry& /*geom*/, const amrex::Real* /*dx*/,
    const amrex::Array<BCType, AMREX_SPACEDIM>& /*lobc*/,
    const amrex::Array<BCType, AMREX_SPACEDIM>& /*hibc*/,
    int /*ratio*/, const amrex::RealVect& /*a_loc*/)
{
    Util::Warning(INFO, "This code has not been properlyt tested");
}


void
Operator<Grid::Cell>::averageDownCoeffs()
{
    for (int i = 0; i < m_num_a_fabs; i++)
    {
        for (int amrlev = m_num_amr_levels - 1; amrlev > 0; --amrlev)
        {
            auto& fine_a_coeffs = m_a_coeffs[i][amrlev];
            averageDownCoeffsSameAmrLevel(fine_a_coeffs);
        }
        averageDownCoeffsSameAmrLevel(m_a_coeffs[i][0]);
    }
}

void
Operator<Grid::Cell>::averageDownCoeffsSameAmrLevel(amrex::Vector<amrex::MultiFab>& a)
{
    int nmglevs = a.size();
    for (int mglev = 1; mglev < nmglevs; ++mglev)
    {
        amrex::average_down(a[mglev - 1], a[mglev], 0, a[0].nComp(), mg_coarsen_ratio);
    }
}



const amrex::FArrayBox&
Operator<Grid::Cell>::GetFab(const int num, const int amrlev, const int mglev, const amrex::MFIter& mfi) const
{
    return m_a_coeffs[num][amrlev][mglev][mfi];
}


void
Operator<Grid::Cell>::RegisterNewFab(amrex::Vector<amrex::MultiFab>& input)
{
    /// \todo assertions here
    m_a_coeffs.resize(m_a_coeffs.size() + 1);
    m_a_coeffs[m_num_a_fabs].resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_a_coeffs[m_num_a_fabs][amrlev].resize(m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
            m_a_coeffs[m_num_a_fabs][amrlev][mglev].define(m_grids[amrlev][mglev],
                m_dmap[amrlev][mglev],
                input[amrlev].nComp(),
                input[amrlev].nGrow());

        amrex::MultiFab::Copy(m_a_coeffs[m_num_a_fabs][amrlev][0],
            input[amrlev], 0, 0,
            input[amrlev].nComp(),
            input[amrlev].nGrow());
    }
    m_num_a_fabs++;
}


void
Operator<Grid::Cell>::RegisterNewFab(amrex::Vector<std::unique_ptr<amrex::MultiFab> >& input)
{
    /// \todo assertions here
    m_a_coeffs.resize(m_a_coeffs.size() + 1);
    m_a_coeffs[m_num_a_fabs].resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_a_coeffs[m_num_a_fabs][amrlev].resize(m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
            m_a_coeffs[m_num_a_fabs][amrlev][mglev].define(m_grids[amrlev][mglev],
                m_dmap[amrlev][mglev],
                input[amrlev]->nComp(),
                input[amrlev]->nGrow());

        amrex::MultiFab::Copy(m_a_coeffs[m_num_a_fabs][amrlev][0],
            *input[amrlev], 0, 0,
            input[amrlev]->nComp(),
            input[amrlev]->nGrow());
    }
    m_num_a_fabs++;
}


}
