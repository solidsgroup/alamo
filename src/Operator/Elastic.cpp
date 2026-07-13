// TODO: Remove these 

#include "Elastic.H"
#include "AMReX_Loop.H"
#include "AMReX_Reduce.H"
#include "Set/Set.H"

#include "Numeric/Stencil.H"
#include <cstdlib>
#include <cstring>

namespace
{
int AlamoCoeffDiagLimit()
{
    const char* value = std::getenv("ALAMO_ML_COEFF_DIAG");
    return value ? std::atoi(value) : 0;
}

// Dumps a Frobenius-norm-squared and max-abs-entry summary of the elasticity
// coefficient field at one (amrlev, mglev) so it can be diffed CPU-vs-GPU.
template<int SYM>
void AlamoPrintCoeffDiag(int amrlev, int mglev, amrex::FabArray<amrex::BaseFab<Set::Matrix4<AMREX_SPACEDIM, SYM>>>& ddw)
{
    using MATRIX4 = Set::Matrix4<AMREX_SPACEDIM, SYM>;

    amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpMax, amrex::ReduceOpSum> reduce_op;
    amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;

    for (amrex::MFIter mfi(ddw, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const amrex::Box& bx = mfi.validbox();
        amrex::Array4<const MATRIX4> const& C = ddw.const_array(mfi);

        // Read the tensor's raw storage rather than calling MATRIX4::operator(),
        // which is host-only for some symmetry specializations (e.g. Major) and
        // would fail to compile in a device lambda. MATRIX4 is a trivially
        // copyable struct of doubles for every specialization in use here.
        reduce_op.eval(bx, reduce_data,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
            {
                const MATRIX4& m = C(i, j, k);
                const double* raw = reinterpret_cast<const double*>(&m);
                constexpr int nraw = sizeof(MATRIX4) / sizeof(double);
                amrex::Real sumsq = 0.0, maxabs = 0.0;
                for (int t = 0; t < nraw; ++t)
                {
                    const amrex::Real v = raw[t];
                    sumsq += v * v;
                    maxabs = amrex::max(maxabs, std::abs(v));
                }
                return { sumsq, maxabs, 1.0 };
            });
    }

    ReduceTuple hv = reduce_data.value(reduce_op);
    amrex::Real sumsq = amrex::get<0>(hv);
    amrex::Real maxabs = amrex::get<1>(hv);
    amrex::Real count = amrex::get<2>(hv);

    amrex::ParallelDescriptor::ReduceRealSum(sumsq);
    amrex::ParallelDescriptor::ReduceRealMax(maxabs);
    amrex::ParallelDescriptor::ReduceRealSum(count);

    Util::Message(INFO, "ALAMO_ML_COEFF_DIAG amrlev=", amrlev,
        " mglev=", mglev,
        " nboxes=", ddw.boxArray().size(),
        " nnodes=", count,
        " frob_norm=", std::sqrt(sumsq),
        " max_abs_entry=", maxabs);
}

int AlamoDiagProbeLimit()
{
    const char* value = std::getenv("ALAMO_ML_DIAG_PROBE");
    return value ? std::atoi(value) : 0;
}

// Dumps min/max/frobenius-norm/nan-inf-zero counts of m_diag at one
// (amrlev, mglev), separately over the valid region and over the full grown
// (ghost-inclusive) region, so the two can be diffed to see whether ghost
// values differ from interior ones (the latter is independently verified
// correct here; the former is unverified prior to this probe).
void AlamoPrintDiagProbe(int amrlev, int mglev, const char* label, const amrex::MultiFab& diag)
{
    const int ncomp = diag.nComp();

    auto summarize = [&](const amrex::IntVect& ng, const char* region)
    {
        amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpMin, amrex::ReduceOpMax, amrex::ReduceOpSum> reduce_op;
        amrex::ReduceData<amrex::Real, amrex::Real, amrex::Real, amrex::Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        amrex::Long nancount = 0, infcount = 0, zerocount = 0;
        for (amrex::MFIter mfi(diag, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box bx = mfi.growntilebox(ng);
            amrex::Array4<const amrex::Real> const& d = diag.const_array(mfi);

            amrex::Gpu::DeviceScalar<amrex::Long> nan_d(0), inf_d(0), zero_d(0);
            amrex::Long* nan_p = nan_d.dataPtr();
            amrex::Long* inf_p = inf_d.dataPtr();
            amrex::Long* zero_p = zero_d.dataPtr();

            reduce_op.eval(bx, reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                {
                    amrex::Real sumsq = 0.0, vmin = 1e300, vmax = -1e300, sum = 0.0;
                    for (int n = 0; n < ncomp; ++n)
                    {
                        const amrex::Real v = d(i, j, k, n);
                        if (std::isnan(v)) amrex::Gpu::Atomic::Add(nan_p, (amrex::Long)1);
                        else if (std::isinf(v)) amrex::Gpu::Atomic::Add(inf_p, (amrex::Long)1);
                        else if (v == 0.0) amrex::Gpu::Atomic::Add(zero_p, (amrex::Long)1);
                        if (!std::isnan(v) && !std::isinf(v))
                        {
                            sumsq += v * v;
                            sum += v;
                            vmin = amrex::min(vmin, v);
                            vmax = amrex::max(vmax, v);
                        }
                    }
                    return { sumsq, vmin, vmax, sum };
                });
            nancount += nan_d.dataValue();
            infcount += inf_d.dataValue();
            zerocount += zero_d.dataValue();
        }

        ReduceTuple hv = reduce_data.value(reduce_op);
        amrex::Real sumsq = amrex::get<0>(hv);
        amrex::Real vmin = amrex::get<1>(hv);
        amrex::Real vmax = amrex::get<2>(hv);
        amrex::Real sum = amrex::get<3>(hv);

        amrex::ParallelDescriptor::ReduceRealSum(sumsq);
        amrex::ParallelDescriptor::ReduceRealMin(vmin);
        amrex::ParallelDescriptor::ReduceRealMax(vmax);
        amrex::ParallelDescriptor::ReduceRealSum(sum);
        amrex::ParallelDescriptor::ReduceLongSum(nancount);
        amrex::ParallelDescriptor::ReduceLongSum(infcount);
        amrex::ParallelDescriptor::ReduceLongSum(zerocount);

        Util::Message(INFO, "ALAMO_ML_DIAG_PROBE ", label,
            " amrlev=", amrlev, " mglev=", mglev, " region=", region,
            " sum=", sum, " min=", vmin, " max=", vmax,
            " nan=", nancount, " inf=", infcount, " zero=", zerocount,
            " frob_norm=", std::sqrt(sumsq));
    };

    summarize(amrex::IntVect(0), "valid");
    summarize(amrex::IntVect(1), "one_ghost");
    summarize(diag.nGrowVect(), "grown");
}
}

#ifdef ALAMO_GPU
#define ALAMO_ELASTIC_OP_FOR amrex::ParallelFor
#define ALAMO_ELASTIC_OP_CAPTURE [=]
#define ALAMO_ELASTIC_OP_DEVICE AMREX_GPU_DEVICE
#define ALAMO_ELASTIC_OP_BC_EVAL(bc, bc_type, u, gradu, sigma, i, j, k, bx) \
    ::BC::Operator::Elastic::Elastic::eval(bc_type, u, gradu, sigma, i, j, k, bx)
#else
#define ALAMO_ELASTIC_OP_FOR amrex::LoopConcurrentOnCpu
#define ALAMO_ELASTIC_OP_CAPTURE [=]
#define ALAMO_ELASTIC_OP_DEVICE
#define ALAMO_ELASTIC_OP_BC_EVAL(bc, bc_type, u, gradu, sigma, i, j, k, bx) \
    (*(bc))(u, gradu, sigma, i, j, k, bx)
#endif

namespace Operator
{
template<int SYM>
Elastic<SYM>::Elastic(const Vector<Geometry>& a_geom,
    const Vector<BoxArray>& a_grids,
    const Vector<DistributionMapping>& a_dmap,
    const LPInfo& a_info)
{
    BL_PROFILE("Operator::Elastic::Elastic()");

    define(a_geom, a_grids, a_dmap, a_info);
}

template<int SYM>
Elastic<SYM>::~Elastic()
{}

template<int SYM>
void
Elastic<SYM>::define(const Vector<Geometry>& a_geom,
    const Vector<BoxArray>& a_grids,
    const Vector<DistributionMapping>& a_dmap,
    const LPInfo& a_info,
    const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
    BL_PROFILE("Operator::Elastic::define()");

    Operator::define(a_geom, a_grids, a_dmap, a_info, a_factory);

    // D- G+ has the ordinary one-cell Laplacian high-frequency spectrum.
    // The generic nodal default (2/3) is too aggressive for 3-D elastic
    // cross-coupling; use a stable damped-Jacobi default.  Users can still
    // override it through elastic.solver.omega.
    SetOmega(0.5);

    int model_nghost = 2;

    m_ddw_mf.resize(m_num_amr_levels);
    m_psi_mf.resize(m_num_amr_levels);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        m_ddw_mf[amrlev].resize(m_num_mg_levels[amrlev]);
        m_psi_mf[amrlev].resize(m_num_mg_levels[amrlev]);
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            m_ddw_mf[amrlev][mglev].reset(new MultiTab(amrex::convert(m_grids[amrlev][mglev],
                amrex::IntVect::TheNodeVector()),
                m_dmap[amrlev][mglev], 1, model_nghost));
            m_psi_mf[amrlev][mglev].reset(new MultiFab(m_grids[amrlev][mglev],
                m_dmap[amrlev][mglev], 1, model_nghost));

            if (!m_psi_set) m_psi_mf[amrlev][mglev]->setVal(1.0);
        }
    }
}

template <int SYM>
void
Elastic<SYM>::SetModel(MATRIX4& a_model)
{
    for (int amrlev = 0; amrlev < m_num_amr_levels; amrlev++)
    {
        amrex::Box domain(m_geom[amrlev][0].Domain());
        domain.convert(amrex::IntVect::TheNodeVector());

#ifdef AMREX_DEBUG
#ifdef ALAMO_GPU
        Util::DeviceErrorFlag setmodel_error;
        int* setmodel_error_flag = setmodel_error.dataPtr();
#endif
#endif
        for (MFIter mfi(*m_ddw_mf[amrlev][0], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box bx = mfi.grownnodaltilebox();

            amrex::Array4<MATRIX4> const& ddw = (*(m_ddw_mf[amrlev][0])).array(mfi);

            ALAMO_ELASTIC_OP_FOR(bx, ALAMO_ELASTIC_OP_CAPTURE ALAMO_ELASTIC_OP_DEVICE (int i, int j, int k) {
                ddw(i, j, k) = a_model;

#ifdef AMREX_DEBUG
#ifdef ALAMO_GPU
                if (ddw(i, j, k).contains_nan()) Util::SetDeviceError(setmodel_error_flag);
#else
                if (ddw(i, j, k).contains_nan()) Util::Abort(INFO, "model is nan at (", i, ",", j, ",", k, "), amrlev=", amrlev);
#endif
#endif
            });
        }
#ifdef AMREX_DEBUG
#ifdef ALAMO_GPU
        Util::AbortIfDeviceError(setmodel_error, INFO, "Operator::Elastic::SetModel() detected a NaN in the device coefficient field");
#endif
#endif
        (*(m_ddw_mf[amrlev][0])).setMultiGhost(true);
        (*(m_ddw_mf[amrlev][0])).FillBoundary( Geom(amrlev,0).periodicity());
    }
    m_model_set = true;
}

template <int SYM>
void
Elastic<SYM>::SetModel(int amrlev, const amrex::FabArray<amrex::BaseFab<MATRIX4> >& a_model)
{
    BL_PROFILE("Operator::Elastic::SetModel()");

    amrex::Box domain(m_geom[amrlev][0].Domain());
    domain.convert(amrex::IntVect::TheNodeVector());

    if (a_model.boxArray() != m_ddw_mf[amrlev][0]->boxArray()) Util::Abort(INFO, "Inconsistent box arrays\n", "a_model.boxArray()=\n", a_model.boxArray(), "\n but the current box array is \n", m_ddw_mf[amrlev][0]->boxArray());
    if (a_model.DistributionMap() != m_ddw_mf[amrlev][0]->DistributionMap()) Util::Abort(INFO, "Inconsistent distribution maps");
    if (a_model.nComp() != m_ddw_mf[amrlev][0]->nComp()) Util::Abort(INFO, "Inconsistent # of components - should be ", m_ddw_mf[amrlev][0]->nComp());
    if (a_model.nGrow() != m_ddw_mf[amrlev][0]->nGrow()) Util::Abort(INFO, "Inconsistent # of ghost nodes, should be ", m_ddw_mf[amrlev][0]->nGrow());


    for (MFIter mfi(a_model, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.grownnodaltilebox();

        amrex::Array4<MATRIX4> const& C = (*(m_ddw_mf[amrlev][0])).array(mfi);
        amrex::Array4<const MATRIX4> const& a_C = a_model.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            C(i, j, k) = a_C(i, j, k);
        });
    }
    m_ddw_mf[amrlev][0]->setMultiGhost(true);
    m_ddw_mf[amrlev][0]->FillBoundaryAndSync(Geom(amrlev,0).periodicity());
    m_model_set = true;
}

template <int SYM>
void
Elastic<SYM>::SetPsi(int amrlev, const amrex::MultiFab& a_psi_mf)
{
    BL_PROFILE("Operator::Elastic::SetPsi()");
    amrex::Box domain(m_geom[amrlev][0].Domain());

    for (MFIter mfi(a_psi_mf, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.growntilebox() & domain;

        amrex::Array4<Set::Scalar> const& m_psi = (*(m_psi_mf[amrlev][0])).array(mfi);
        amrex::Array4<const Set::Scalar> const& a_psi = a_psi_mf.array(mfi);

        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
            m_psi(i, j, k) = a_psi(i, j, k);
        });
    }
    m_psi_set = true;
}

template<int SYM>
void
Elastic<SYM>::Fapply(int amrlev, int mglev, MultiFab& a_f, const MultiFab& a_u) const
{
    BL_PROFILE("Operator::Elastic::Fapply()");

    amrex::Box domain(m_geom[amrlev][mglev].growPeriodicDomain(1));
    domain.convert(amrex::IntVect::TheNodeVector());

    amrex::Box stencilbox(m_geom[amrlev][mglev].growPeriodicDomain(2));
    stencilbox.convert(amrex::IntVect::TheNodeVector());

    const amrex::Box psi_cell_domain(m_geom[amrlev][mglev].Domain());
    const amrex::Dim3 psi_cell_lo = amrex::lbound(psi_cell_domain);
    const amrex::Dim3 psi_cell_hi = amrex::ubound(psi_cell_domain);
    const std::array<bool, AMREX_SPACEDIM> psi_periodic = {
        AMREX_D_DECL(m_geom[amrlev][mglev].isPeriodic(0),
            m_geom[amrlev][mglev].isPeriodic(1),
            m_geom[amrlev][mglev].isPeriodic(2))};

    Set::Vector DX(m_geom[amrlev][mglev].CellSize());
#ifdef AMREX_DEBUG
#ifdef ALAMO_GPU
    Util::DeviceErrorFlag fapply_error;
    int* fapply_error_flag = fapply_error.dataPtr();
#endif
#endif

    for (MFIter mfi(a_f, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.validbox().grow(1) & domain;
        amrex::Box tilebox = mfi.grownnodaltilebox() & bx;

        amrex::Array4<MATRIX4> const& DDW = (*(m_ddw_mf[amrlev][mglev])).array(mfi);
        amrex::Array4<const amrex::Real> const& U = a_u.array(mfi);
        amrex::Array4<amrex::Real> const& F = a_f.array(mfi);
        amrex::Array4<const Set::Scalar> const& psi = m_psi_mf[amrlev][mglev]->const_array(mfi);

        const Dim3 lo = amrex::lbound(stencilbox), hi = amrex::ubound(stencilbox);
        const bool m_psi_set = this->m_psi_set;
        const Set::Scalar psi_regularization = this->m_psi_small;
#ifdef ALAMO_GPU
        auto m_bc_type = this->m_bc->GetBcTypeArray();
#else
        auto m_bc = this->m_bc;
#endif

        ALAMO_ELASTIC_OP_FOR(tilebox, ALAMO_ELASTIC_OP_CAPTURE ALAMO_ELASTIC_OP_DEVICE (int i, int j, int k)
        {
            const auto sten = Numeric::GetStencil(i, j, k, stencilbox);
            const bool on_boundary = AMREX_D_TERM(
                !psi_periodic[0] && (i == lo.x || i == hi.x),
                || !psi_periodic[1] && (j == lo.y || j == hi.y),
                || !psi_periodic[2] && (k == lo.z || k == hi.z));

            // The bulk is the conservative D- G+ pair.  At physical
            // traction rows PairedGradient keeps the one-sided normal but
            // restores centered tangential derivatives; that same closure
            // supplies the boundary flux consumed by the adjacent D- row.
            const Set::Matrix gradu = Numeric::PairedGradient(U, i, j, k, DX.data(), sten, psi_periodic);

            Set::Scalar theta = 1.0;
            if (m_psi_set)
                theta = (1.0 - psi_regularization)
                    * Numeric::Interpolate::CellToNodeAverageBounded(psi, i, j, k, 0,
                        psi_cell_lo, psi_cell_hi, psi_periodic)
                    + psi_regularization;

            const Set::Matrix flux = (DDW(i, j, k) * gradu) * theta;
            Set::Vector u;
            for (int p = 0; p < AMREX_SPACEDIM; ++p) u(p) = U(i, j, k, p);

            Set::Vector f = Set::Vector::Zero();
            if (on_boundary)
            {
                f = ALAMO_ELASTIC_OP_BC_EVAL(m_bc, m_bc_type, u, gradu, flux, i, j, k, stencilbox);
            }
            else
            {
#if AMREX_SPACEDIM > 0
                const auto sten_xlo = Numeric::GetStencil(i - 1, j, k, stencilbox);
                const Set::Matrix gradu_xlo = Numeric::PairedGradient(U, i - 1, j, k, DX.data(), sten_xlo, psi_periodic);
                Set::Scalar theta_xlo = 1.0;
                if (m_psi_set)
                    theta_xlo = (1.0 - psi_regularization)
                        * Numeric::Interpolate::CellToNodeAverageBounded(psi, i - 1, j, k, 0,
                            psi_cell_lo, psi_cell_hi, psi_periodic)
                        + psi_regularization;
                const Set::Matrix flux_xlo = (DDW(i - 1, j, k) * gradu_xlo) * theta_xlo;
                f += (flux.col(0) - flux_xlo.col(0)) / DX[0];
#endif
#if AMREX_SPACEDIM > 1
                const auto sten_ylo = Numeric::GetStencil(i, j - 1, k, stencilbox);
                const Set::Matrix gradu_ylo = Numeric::PairedGradient(U, i, j - 1, k, DX.data(), sten_ylo, psi_periodic);
                Set::Scalar theta_ylo = 1.0;
                if (m_psi_set)
                    theta_ylo = (1.0 - psi_regularization)
                        * Numeric::Interpolate::CellToNodeAverageBounded(psi, i, j - 1, k, 0,
                            psi_cell_lo, psi_cell_hi, psi_periodic)
                        + psi_regularization;
                const Set::Matrix flux_ylo = (DDW(i, j - 1, k) * gradu_ylo) * theta_ylo;
                f += (flux.col(1) - flux_ylo.col(1)) / DX[1];
#endif
#if AMREX_SPACEDIM > 2
                const auto sten_zlo = Numeric::GetStencil(i, j, k - 1, stencilbox);
                const Set::Matrix gradu_zlo = Numeric::PairedGradient(U, i, j, k - 1, DX.data(), sten_zlo, psi_periodic);
                Set::Scalar theta_zlo = 1.0;
                if (m_psi_set)
                    theta_zlo = (1.0 - psi_regularization)
                        * Numeric::Interpolate::CellToNodeAverageBounded(psi, i, j, k - 1, 0,
                            psi_cell_lo, psi_cell_hi, psi_periodic)
                        + psi_regularization;
                const Set::Matrix flux_zlo = (DDW(i, j, k - 1) * gradu_zlo) * theta_zlo;
                f += (flux.col(2) - flux_zlo.col(2)) / DX[2];
#endif
            }

#ifdef AMREX_DEBUG
#ifdef ALAMO_GPU
            if (std::isnan(f(0)) || std::isnan(f(1))
#if AMREX_SPACEDIM == 3
                || std::isnan(f(2))
#endif
            )
            {
                Util::SetDeviceError(fapply_error_flag);
            }
#else
            if (std::isnan(f(0)) || std::isnan(f(1))
#if AMREX_SPACEDIM == 3
                || std::isnan(f(2))
#endif
            )
                Util::Abort(INFO, "Elastic::Fapply produced NaN at (", i, ",", j, ",", k,
                    "), amrlev=", amrlev, ", mglev=", mglev);
#endif
#endif

            AMREX_D_TERM(F(i, j, k, 0) = f[0];,
                F(i, j, k, 1) = f[1];,
                F(i, j, k, 2) = f[2];);
        });
#ifdef AMREX_DEBUG
#ifdef ALAMO_GPU
        Util::AbortIfDeviceError(fapply_error, INFO,
            "Operator::Elastic::Fapply() detected an invalid value on device");
#endif
#endif
    }
}



template<int SYM>
void
Elastic<SYM>::Diagonal(int amrlev, int mglev, MultiFab& a_diag)
{
    BL_PROFILE("Operator::Elastic::Diagonal()");

    amrex::Box domain(m_geom[amrlev][mglev].growPeriodicDomain(1));
    domain.convert(amrex::IntVect::TheNodeVector());

    amrex::Box stencilbox(m_geom[amrlev][mglev].growPeriodicDomain(2));
    stencilbox.convert(amrex::IntVect::TheNodeVector());

    const amrex::Box psi_cell_domain(m_geom[amrlev][mglev].Domain());
    const amrex::Dim3 psi_cell_lo = amrex::lbound(psi_cell_domain);
    const amrex::Dim3 psi_cell_hi = amrex::ubound(psi_cell_domain);
    const std::array<bool, AMREX_SPACEDIM> psi_periodic = {
        AMREX_D_DECL(m_geom[amrlev][mglev].isPeriodic(0),
            m_geom[amrlev][mglev].isPeriodic(1),
            m_geom[amrlev][mglev].isPeriodic(2))};

    Set::Vector DX(m_geom[amrlev][mglev].CellSize());
#ifdef AMREX_DEBUG
#ifdef ALAMO_GPU
    Util::DeviceErrorFlag diagonal_error;
    int* diagonal_error_flag = diagonal_error.dataPtr();
#endif
#endif

    for (MFIter mfi(a_diag, false); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.validbox().grow(1) & domain;
        amrex::Box tilebox = mfi.grownnodaltilebox() & bx;

        amrex::Array4<MATRIX4> const& DDW = (*(m_ddw_mf[amrlev][mglev])).array(mfi);
        amrex::Array4<Set::Scalar> const& diag = a_diag.array(mfi);
        amrex::Array4<const Set::Scalar> const& psi = m_psi_mf[amrlev][mglev]->const_array(mfi);

        const Dim3 lo = amrex::lbound(stencilbox), hi = amrex::ubound(stencilbox);
        const bool m_psi_set = this->m_psi_set;
        const Set::Scalar psi_regularization = this->m_psi_small;
#ifdef ALAMO_GPU
        auto m_bc_type = this->m_bc->GetBcTypeArray();
#else
        auto m_bc = this->m_bc;
#endif

        amrex::ParallelFor(tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            const auto sten = Numeric::GetStencil(i, j, k, stencilbox);
            auto physical_sten = sten;
            AMREX_D_TERM(
                if (psi_periodic[0]) physical_sten[0] = Numeric::StencilType::Central;,
                if (psi_periodic[1]) physical_sten[1] = Numeric::StencilType::Central;,
                if (psi_periodic[2]) physical_sten[2] = Numeric::StencilType::Central;);
            const bool on_boundary = AMREX_D_TERM(
                physical_sten[0] != Numeric::StencilType::Central,
                || physical_sten[1] != Numeric::StencilType::Central,
                || physical_sten[2] != Numeric::StencilType::Central);

            Set::Scalar theta = 1.0;
            if (m_psi_set)
                theta = (1.0 - psi_regularization)
                    * Numeric::Interpolate::CellToNodeAverageBounded(psi, i, j, k, 0,
                        psi_cell_lo, psi_cell_hi, psi_periodic)
                    + psi_regularization;

            for (int p = 0; p < AMREX_SPACEDIM; ++p)
            {
                if (on_boundary)
                {
                    // Derivative of PairedGradient(e_p) at this node.  The
                    // physical-boundary closure has a centered tangential
                    // derivative (zero self coefficient), a forward low-side
                    // normal derivative (-1/dx), and a backward high-side
                    // normal derivative (+1/dx).
                    Set::Matrix gradu = Set::Matrix::Zero();
                    for (int d = 0; d < AMREX_SPACEDIM; ++d)
                        gradu(p, d) = physical_sten[d] == Numeric::StencilType::Lo
                            ? 1.0 / DX[d]
                            : (physical_sten[d] == Numeric::StencilType::Hi
                                ? -1.0 / DX[d] : 0.0);

                    const Set::Matrix flux = (DDW(i, j, k) * gradu) * theta;
                    Set::Vector u = Set::Vector::Zero();
                    u(p) = 1.0;
                    const Set::Vector f = ALAMO_ELASTIC_OP_BC_EVAL(
                        m_bc, m_bc_type, u, gradu, flux, i, j, k, stencilbox);
                    diag(i, j, k, p) = f(p);
                }
                else
                {
                    // Exact local impulse diagonal of D- [theta DDW Gpair].
                    // At i the impulse's forward difference is -1/dx; at
                    // each preceding flux anchor it is +1/dx in that
                    // divergence direction.
                    Set::Scalar value = 0.0;
                    Set::Matrix e0 = Set::Matrix::Zero();
                    for (int d = 0; d < AMREX_SPACEDIM; ++d)
                        e0(p, d) = -1.0 / DX[d];
                    const Set::Matrix flux0 = (DDW(i, j, k) * e0) * theta;

#if AMREX_SPACEDIM > 0
                    Set::Matrix exlo = Set::Matrix::Zero();
                    exlo(p, 0) = 1.0 / DX[0];
                    const Set::Scalar theta_xlo = m_psi_set
                        ? (1.0 - psi_regularization)
                            * Numeric::Interpolate::CellToNodeAverageBounded(psi, i - 1, j, k, 0,
                                psi_cell_lo, psi_cell_hi, psi_periodic)
                            + psi_regularization
                        : 1.0;
                    const Set::Matrix flux_xlo = (DDW(i - 1, j, k) * exlo) * theta_xlo;
                    value += (flux0(p, 0) - flux_xlo(p, 0)) / DX[0];
#endif
#if AMREX_SPACEDIM > 1
                    Set::Matrix eylo = Set::Matrix::Zero();
                    eylo(p, 1) = 1.0 / DX[1];
                    const Set::Scalar theta_ylo = m_psi_set
                        ? (1.0 - psi_regularization)
                            * Numeric::Interpolate::CellToNodeAverageBounded(psi, i, j - 1, k, 0,
                                psi_cell_lo, psi_cell_hi, psi_periodic)
                            + psi_regularization
                        : 1.0;
                    const Set::Matrix flux_ylo = (DDW(i, j - 1, k) * eylo) * theta_ylo;
                    value += (flux0(p, 1) - flux_ylo(p, 1)) / DX[1];
#endif
#if AMREX_SPACEDIM > 2
                    Set::Matrix ezlo = Set::Matrix::Zero();
                    ezlo(p, 2) = 1.0 / DX[2];
                    const Set::Scalar theta_zlo = m_psi_set
                        ? (1.0 - psi_regularization)
                            * Numeric::Interpolate::CellToNodeAverageBounded(psi, i, j, k - 1, 0,
                                psi_cell_lo, psi_cell_hi, psi_periodic)
                            + psi_regularization
                        : 1.0;
                    const Set::Matrix flux_zlo = (DDW(i, j, k - 1) * ezlo) * theta_zlo;
                    value += (flux0(p, 2) - flux_zlo(p, 2)) / DX[2];
#endif

                    diag(i, j, k, p) = value;
                }

#ifdef AMREX_DEBUG
#ifdef ALAMO_GPU
                if (std::isnan(diag(i, j, k, p)) || std::isinf(diag(i, j, k, p)))
                    Util::SetDeviceError(diagonal_error_flag);
#else
                if (std::isnan(diag(i, j, k, p))
                    || std::isinf(diag(i, j, k, p)))
                    Util::Abort(INFO, "Elastic::Diagonal produced an invalid value at (",
                        i, ",", j, ",", k, "), amrlev=", amrlev, ", mglev=", mglev);
#endif
#endif
            }
        });
#ifdef AMREX_DEBUG
#ifdef ALAMO_GPU
        Util::AbortIfDeviceError(diagonal_error, INFO,
            "Operator::Elastic::Diagonal() detected an invalid value on device");
#endif
#endif
    }

    a_diag.FillBoundaryAndSync(Geom(amrlev, mglev).periodicity());
    nodalSync(amrlev, mglev, a_diag);

    // A raw exact-zero mask has genuinely unconstrained rows; Fsmooth and
    // normalize divide by this diagonal.  Reject that unsupported formulation
    // before MLMG can silently turn it into NaN/Inf.  An inactive-DOF/nullspace
    // policy is deliberately separate from this paired-stencil correction.
    if (m_psi_set && m_psi_small == 0.0)
    {
        amrex::ReduceOps<amrex::ReduceOpSum> zero_op;
        amrex::ReduceData<amrex::Real> zero_data(zero_op);
        using ZeroTuple = typename decltype(zero_data)::Type;
        for (MFIter mfi(a_diag, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box bx = mfi.validbox() & domain;
            const auto diag = a_diag.const_array(mfi);
            const int ncomp = a_diag.nComp();
            zero_op.eval(bx, zero_data,
                [=] AMREX_GPU_DEVICE(int i, int j, int k) -> ZeroTuple
                {
                    amrex::Real count = 0.0;
                    for (int p = 0; p < ncomp; ++p)
                        if (diag(i, j, k, p) == 0.0) count += 1.0;
                    return {count};
                });
        }
        amrex::Real zero_rows = amrex::get<0>(zero_data.value(zero_op));
        amrex::ParallelDescriptor::ReduceRealSum(zero_rows);
        if (zero_rows > 0.0)
            Util::Abort(INFO, "Raw psi produced ", zero_rows,
                " zero-support elastic diagonal entries. Exact-zero masked elasticity requires an explicit inactive-DOF/nullspace policy; use elastic.use_psi=0 with a finite model_void, or explicitly opt into elasticop.small.");
    }

    if (AlamoDiagProbeLimit()) AlamoPrintDiagProbe(amrlev, mglev, "paired", a_diag);
}


template<int SYM>
void
Elastic<SYM>::Error0x(int amrlev, int mglev, MultiFab& R0x, const MultiFab& x) const
{
    BL_PROFILE("Operator::Elastic::Error0x()");
    Util::Message(INFO);

    int ncomp = x.nComp();//getNComp();
    int nghost = x.nGrow();

    if (!m_diagonal_computed)
        Util::Abort(INFO, "Operator::Diagonal() must be called before using normalize");

    amrex::MultiFab D0x(x.boxArray(), x.DistributionMap(), ncomp, nghost);
    amrex::MultiFab AD0x(x.boxArray(), x.DistributionMap(), ncomp, nghost);

    amrex::MultiFab::Copy(D0x, x, 0, 0, ncomp, nghost); // D0x = x
    amrex::MultiFab::Divide(D0x, *m_diag[amrlev][mglev], 0, 0, ncomp, 0); // D0x = x/diag
    amrex::MultiFab::Copy(AD0x, D0x, 0, 0, ncomp, nghost); // AD0x = D0x

    Fapply(amrlev, mglev, AD0x, D0x);    // AD0x = A * D0 * x

    amrex::MultiFab::Copy(R0x, x, 0, 0, ncomp, nghost); // R0x = x
    amrex::MultiFab::Subtract(R0x, AD0x, 0, 0, ncomp, nghost); // R0x = x - AD0x
}


template<int SYM>
void
Elastic<SYM>::FFlux(int /*amrlev*/, const MFIter& /*mfi*/,
    const std::array<FArrayBox*, AMREX_SPACEDIM>& sigmafab,
    const FArrayBox& /*ufab*/, const int /*face_only*/) const
{
    BL_PROFILE("Operator::Elastic::FFlux()");
    Util::Message(INFO);
    amrex::BaseFab<amrex::Real> AMREX_D_DECL(&fxfab = *sigmafab[0],
        &fyfab = *sigmafab[1],
        &fzfab = *sigmafab[2]);
    AMREX_D_TERM(fxfab.setVal<amrex::RunOn::Device>(0.0);,
        fyfab.setVal<amrex::RunOn::Device>(0.0);,
        fzfab.setVal<amrex::RunOn::Device>(0.0););

}

template<int SYM>
void
Elastic<SYM>::Strain(int amrlev,
    amrex::MultiFab& a_eps,
    const amrex::MultiFab& a_u,
    bool voigt) const
{
    BL_PROFILE("Operator::Elastic::Strain()");

    Set::Vector DX(m_geom[amrlev][0].CellSize()); // device-safe cell size (copied into device closure)
    amrex::Box domain(m_geom[amrlev][0].Domain());
    domain.convert(amrex::IntVect::TheNodeVector());


    for (MFIter mfi(a_u, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        amrex::Array4<amrex::Real> const& epsilon = a_eps.array(mfi);
        amrex::Array4<const amrex::Real> const& u = a_u.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Matrix gradu;

            std::array<Numeric::StencilType, AMREX_SPACEDIM> sten
                = Numeric::GetStencil(i, j, k, domain);

            // Fill gradu
            for (int p = 0; p < AMREX_SPACEDIM; p++)
            {
                AMREX_D_TERM(gradu(p, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(u, i, j, k, p, DX.data(), sten));,
                    gradu(p, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(u, i, j, k, p, DX.data(), sten));,
                    gradu(p, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(u, i, j, k, p, DX.data(), sten)););
            }

            Set::Matrix eps = 0.5 * (gradu + gradu.transpose());

            if (voigt)
            {
                AMREX_D_PICK(epsilon(i, j, k, 0) = eps(0, 0);
                ,
                    epsilon(i, j, k, 0) = eps(0, 0); epsilon(i, j, k, 1) = eps(1, 1); epsilon(i, j, k, 2) = eps(0, 1);
                ,
                    epsilon(i, j, k, 0) = eps(0, 0); epsilon(i, j, k, 1) = eps(1, 1); epsilon(i, j, k, 2) = eps(2, 2);
                epsilon(i, j, k, 3) = eps(1, 2); epsilon(i, j, k, 4) = eps(2, 0); epsilon(i, j, k, 5) = eps(0, 1););
            }
            else
            {
                AMREX_D_PICK(epsilon(i, j, k, 0) = eps(0, 0);
                ,
                    epsilon(i, j, k, 0) = eps(0, 0); epsilon(i, j, k, 1) = eps(0, 1);
                epsilon(i, j, k, 2) = eps(1, 0); epsilon(i, j, k, 3) = eps(1, 1);
                ,
                    epsilon(i, j, k, 0) = eps(0, 0); epsilon(i, j, k, 1) = eps(0, 1); epsilon(i, j, k, 2) = eps(0, 2);
                epsilon(i, j, k, 3) = eps(1, 0); epsilon(i, j, k, 4) = eps(1, 1); epsilon(i, j, k, 5) = eps(1, 2);
                epsilon(i, j, k, 6) = eps(2, 0); epsilon(i, j, k, 7) = eps(2, 1); epsilon(i, j, k, 8) = eps(2, 2););
            }
        });
    }
}


template<int SYM>
void
Elastic<SYM>::Stress(int amrlev,
    amrex::MultiFab& a_sigma,
    const amrex::MultiFab& a_u,
    bool voigt, bool a_homogeneous)
{
    BL_PROFILE("Operator::Elastic::Stress()");
    SetHomogeneous(a_homogeneous);

    Set::Vector DX(m_geom[amrlev][0].CellSize()); // device-safe cell size (copied into device closure)
    amrex::Box domain(m_geom[amrlev][0].Domain());
    domain.convert(amrex::IntVect::TheNodeVector());

    for (MFIter mfi(a_u, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        amrex::Array4<Set::Matrix4<AMREX_SPACEDIM, SYM>> const& DDW = (*(m_ddw_mf[amrlev][0])).array(mfi);
        amrex::Array4<amrex::Real> const& sigma = a_sigma.array(mfi);
        amrex::Array4<Set::Scalar> const& psi = m_psi_mf[amrlev][0]->array(mfi);
        amrex::Array4<const amrex::Real> const& u = a_u.array(mfi);
        auto m_psi_set = this->m_psi_set;
        auto m_psi_small = this->m_psi_small;
        ALAMO_ELASTIC_OP_FOR(bx, ALAMO_ELASTIC_OP_CAPTURE ALAMO_ELASTIC_OP_DEVICE (int i, int j, int k)
        {
            Set::Matrix gradu;

            std::array<Numeric::StencilType, AMREX_SPACEDIM> sten
                = Numeric::GetStencil(i, j, k, domain);

            // Fill gradu
            for (int p = 0; p < AMREX_SPACEDIM; p++)
            {
                AMREX_D_TERM(gradu(p, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(u, i, j, k, p, DX.data(), sten));,
                    gradu(p, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(u, i, j, k, p, DX.data(), sten));,
                    gradu(p, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(u, i, j, k, p, DX.data(), sten)););
            }

            Set::Scalar psi_avg = 1.0;
            if (m_psi_set) psi_avg = (1.0 - m_psi_small) * Numeric::Interpolate::CellToNodeAverage(psi, i, j, k, 0) + m_psi_small;
            Set::Matrix sig = (DDW(i, j, k) * gradu) * psi_avg;

            if (voigt)
            {
                AMREX_D_PICK(sigma(i, j, k, 0) = sig(0, 0);
                ,
                    sigma(i, j, k, 0) = sig(0, 0); sigma(i, j, k, 1) = sig(1, 1); sigma(i, j, k, 2) = sig(0, 1);
                ,
                    sigma(i, j, k, 0) = sig(0, 0); sigma(i, j, k, 1) = sig(1, 1); sigma(i, j, k, 2) = sig(2, 2);
                sigma(i, j, k, 3) = sig(1, 2); sigma(i, j, k, 4) = sig(2, 0); sigma(i, j, k, 5) = sig(0, 1););
            }
            else
            {
                AMREX_D_PICK(sigma(i, j, k, 0) = sig(0, 0);
                ,
                    sigma(i, j, k, 0) = sig(0, 0); sigma(i, j, k, 1) = sig(0, 1);
                sigma(i, j, k, 2) = sig(1, 0); sigma(i, j, k, 3) = sig(1, 1);
                ,
                    sigma(i, j, k, 0) = sig(0, 0); sigma(i, j, k, 1) = sig(0, 1); sigma(i, j, k, 2) = sig(0, 2);
                sigma(i, j, k, 3) = sig(1, 0); sigma(i, j, k, 4) = sig(1, 1); sigma(i, j, k, 5) = sig(1, 2);
                sigma(i, j, k, 6) = sig(2, 0); sigma(i, j, k, 7) = sig(2, 1); sigma(i, j, k, 8) = sig(2, 2););
            }
        });
    }
}


template<int SYM>
void
Elastic<SYM>::Energy(int amrlev,
    amrex::MultiFab& a_energy,
    const amrex::MultiFab& a_u, bool a_homogeneous)
{
    BL_PROFILE("Operator::Elastic::Energy()");
    SetHomogeneous(a_homogeneous);

    amrex::Box domain(m_geom[amrlev][0].Domain());
    domain.convert(amrex::IntVect::TheNodeVector());

    Set::Vector DX(m_geom[amrlev][0].CellSize()); // device-safe cell size (copied into device closure)

    for (MFIter mfi(a_u, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.tilebox();
        amrex::Array4<Set::Matrix4<AMREX_SPACEDIM, SYM>> const& DDW = (*(m_ddw_mf[amrlev][0])).array(mfi);
        amrex::Array4<amrex::Real> const& energy = a_energy.array(mfi);
        amrex::Array4<const amrex::Real> const& u = a_u.array(mfi);
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {
            Set::Matrix gradu;

            std::array<Numeric::StencilType, AMREX_SPACEDIM> sten
                = Numeric::GetStencil(i, j, k, domain);

            // Fill gradu
            for (int p = 0; p < AMREX_SPACEDIM; p++)
            {
                AMREX_D_TERM(gradu(p, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(u, i, j, k, p, DX.data(), sten));,
                    gradu(p, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(u, i, j, k, p, DX.data(), sten));,
                    gradu(p, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(u, i, j, k, p, DX.data(), sten)););
            }

            Set::Matrix eps = .5 * (gradu + gradu.transpose());
            Set::Matrix sig = DDW(i, j, k) * gradu;

            // energy(i,j,k) = (gradu.transpose() * sig).trace();

            //Util::Abort(INFO,"Fix this"); //
            //energy(i,j,k) = C(i,j,k).W(gradu);
            for (int m = 0; m < AMREX_SPACEDIM; m++)
            {
                for (int n = 0; n < AMREX_SPACEDIM; n++)
                {
                    energy(i, j, k) += .5 * sig(m, n) * eps(m, n);
                }
            }
        });
    }
}

template<int SYM>
void
Elastic<SYM>::averageDownCoeffs()
{
    BL_PROFILE("Elastic::averageDownCoeffs()");

    if (m_average_down_coeffs)
        for (int amrlev = m_num_amr_levels - 1; amrlev > 0; --amrlev)
            averageDownCoeffsDifferentAmrLevels(amrlev);

    averageDownCoeffsSameAmrLevel(0);
    for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
    {
        for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
        {
            if (m_ddw_mf[amrlev][mglev]) {
                FillBoundaryCoeff(*m_ddw_mf[amrlev][mglev], Geom(amrlev,mglev).periodicity());
                FillBoundaryCoeff(*m_psi_mf[amrlev][mglev], Geom(amrlev,mglev).periodicity());
            }
        }
    }

    const int coeff_diag_limit = AlamoCoeffDiagLimit();
    if (coeff_diag_limit > 0)
    {
        for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
            for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
                if (m_ddw_mf[amrlev][mglev])
                    AlamoPrintCoeffDiag<SYM>(amrlev, mglev, *m_ddw_mf[amrlev][mglev]);
    }
}

template<int SYM>
void
Elastic<SYM>::averageDownCoeffsDifferentAmrLevels(int fine_amrlev)
{
    BL_PROFILE("Operator::Elastic::averageDownCoeffsDifferentAmrLevels()");
    Util::Assert(INFO, TEST(fine_amrlev > 0));

    const int crse_amrlev = fine_amrlev - 1;
    const int ncomp = 1;

    MultiTab& crse_ddw = *m_ddw_mf[crse_amrlev][0];
    MultiTab& fine_ddw = *m_ddw_mf[fine_amrlev][0];

    amrex::Box cdomain(m_geom[crse_amrlev][0].Domain());
    cdomain.convert(amrex::IntVect::TheNodeVector());

    const Geometry& cgeom = m_geom[crse_amrlev][0];

    const BoxArray& fba = fine_ddw.boxArray();
    const DistributionMapping& fdm = fine_ddw.DistributionMap();

    MultiTab fine_ddw_for_coarse(amrex::coarsen(fba, 2), fdm, ncomp, 2);
    fine_ddw_for_coarse.ParallelCopy(crse_ddw, 0, 0, ncomp, 0, 0, cgeom.periodicity());

    const int coarse_fine_node = 1;
    const int fine_fine_node = 2;

    amrex::iMultiFab nodemask(amrex::coarsen(fba, 2), fdm, 1, 2);
    nodemask.ParallelCopy(*m_nd_fine_mask[crse_amrlev], 0, 0, 1, 0, 0, cgeom.periodicity());

    amrex::iMultiFab cellmask(amrex::convert(amrex::coarsen(fba, 2), amrex::IntVect::TheCellVector()), fdm, 1, 2);
    cellmask.ParallelCopy(*m_cc_fine_mask[crse_amrlev], 0, 0, 1, 1, 1, cgeom.periodicity());

    for (MFIter mfi(fine_ddw_for_coarse, false); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.validbox();

        amrex::Array4<const int> const& nmask = nodemask.array(mfi);
        //amrex::Array4<const int> const& cmask = cellmask.array(mfi);

        amrex::Array4<MATRIX4> const& cdata = fine_ddw_for_coarse.array(mfi);
        amrex::Array4<const MATRIX4> const& fdata = fine_ddw.array(mfi);

        const Dim3 lo = amrex::lbound(cdomain), hi = amrex::ubound(cdomain);

        for (int n = 0; n < fine_ddw.nComp(); n++)
        {
            // I,J,K == coarse coordinates
            // i,j,k == fine coordinates
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int I, int J, int K) {
                int i = I * 2, j = J * 2, k = K * 2;

                if (nmask(I, J, K) == fine_fine_node || nmask(I, J, K) == coarse_fine_node)
                {
                    if ((I == lo.x || I == hi.x) &&
                        (J == lo.y || J == hi.y) &&
                        (K == lo.z || K == hi.z)) // Corner
                        cdata(I, J, K, n) = fdata(i, j, k, n);
                    else if ((J == lo.y || J == hi.y) &&
                        (K == lo.z || K == hi.z)) // X edge
                        cdata(I, J, K, n) = fdata(i - 1, j, k, n) * 0.25 + fdata(i, j, k, n) * 0.5 + fdata(i + 1, j, k, n) * 0.25;
                    else if ((K == lo.z || K == hi.z) &&
                        (I == lo.x || I == hi.x)) // Y edge
                        cdata(I, J, K, n) = fdata(i, j - 1, k, n) * 0.25 + fdata(i, j, k, n) * 0.5 + fdata(i, j + 1, k, n) * 0.25;
                    else if ((I == lo.x || I == hi.x) &&
                        (J == lo.y || J == hi.y)) // Z edge
                        cdata(I, J, K, n) = fdata(i, j, k - 1, n) * 0.25 + fdata(i, j, k, n) * 0.5 + fdata(i, j, k + 1, n) * 0.25;
                    else if (I == lo.x || I == hi.x) // X face
                        cdata(I, J, K, n) =
                        (fdata(i, j - 1, k - 1, n) + fdata(i, j, k - 1, n) * 2.0 + fdata(i, j + 1, k - 1, n)
                            + fdata(i, j - 1, k, n) * 2.0 + fdata(i, j, k, n) * 4.0 + fdata(i, j + 1, k, n) * 2.0
                            + fdata(i, j - 1, k + 1, n) + fdata(i, j, k + 1, n) * 2.0 + fdata(i, j + 1, k + 1, n)) / 16.0;
                    else if (J == lo.y || J == hi.y) // Y face
                        cdata(I, J, K, n) =
                        (fdata(i - 1, j, k - 1, n) + fdata(i - 1, j, k, n) * 2.0 + fdata(i - 1, j, k + 1, n)
                            + fdata(i, j, k - 1, n) * 2.0 + fdata(i, j, k, n) * 4.0 + fdata(i, j, k + 1, n) * 2.0
                            + fdata(i + 1, j, k - 1, n) + fdata(i + 1, j, k, n) * 2.0 + fdata(i + 1, j, k + 1, n)) / 16.0;
                    else if (K == lo.z || K == hi.z) // Z face
                        cdata(I, J, K, n) =
                        (fdata(i - 1, j - 1, k, n) + fdata(i, j - 1, k, n) * 2.0 + fdata(i + 1, j - 1, k, n)
                            + fdata(i - 1, j, k, n) * 2.0 + fdata(i, j, k, n) * 4.0 + fdata(i + 1, j, k, n) * 2.0
                            + fdata(i - 1, j + 1, k, n) + fdata(i, j + 1, k, n) * 2.0 + fdata(i + 1, j + 1, k, n)) / 16.0;
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

#ifdef AMREX_DEBUG
#ifndef ALAMO_GPU
                    if (cdata(I, J, K).contains_nan()) Util::Abort(INFO, "restricted model is nan at (", i, ",", j, ",", k, "), fine_amrlev=", fine_amrlev);
#endif
#endif
                }

            });
        }
    }

    // Copy the fine residual restricted onto the coarse grid
    // into the final residual.

    crse_ddw.ParallelCopy(fine_ddw_for_coarse, 0, 0, ncomp, 0, 0, cgeom.periodicity());
    //const int mglev = 0;
    //Util::RealFillBoundary(crse_ddw, m_geom[crse_amrlev][mglev]);

    FillBoundaryCoeff(crse_ddw,Geom(fine_amrlev,0).periodicity());

}



template<int SYM>
void
Elastic<SYM>::averageDownCoeffsSameAmrLevel(int amrlev)
{
    BL_PROFILE("Elastic::averageDownCoeffsSameAmrLevel()");

    for (int mglev = 1; mglev < m_num_mg_levels[amrlev]; ++mglev)
    {
        amrex::Box cdomain(m_geom[amrlev][mglev].growPeriodicDomain(2));
        cdomain.convert(amrex::IntVect::TheNodeVector());
        amrex::Box fdomain(m_geom[amrlev][mglev - 1].Domain());
        fdomain.convert(amrex::IntVect::TheNodeVector());

        MultiTab& crse = *m_ddw_mf[amrlev][mglev];
        MultiTab& fine = *m_ddw_mf[amrlev][mglev - 1];

        amrex::BoxArray crseba = crse.boxArray();
        amrex::BoxArray fineba = fine.boxArray();

        BoxArray newba = crseba;
        newba.refine(2);
        MultiTab fine_on_crseba;
        fine_on_crseba.define(newba, crse.DistributionMap(), 1, 4);
        fine_on_crseba.ParallelCopy(fine, 0, 0, 1, 2, 4, m_geom[amrlev][mglev-1].periodicity());
        /* ine_on_crseba.FillBoundaryAndSync(m_geom[amrlev][mglev-1].periodicity()); */

        for (MFIter mfi(crse, false); mfi.isValid(); ++mfi)
        {

            Box bx = mfi.grownnodaltilebox() & cdomain;
            /*Box bx = mfi.grownnodaltilebox(-1,1) & cdomain;*/

            amrex::Array4<const Set::Matrix4<AMREX_SPACEDIM, SYM>> const& fdata = fine_on_crseba.array(mfi);
            amrex::Array4<Set::Matrix4<AMREX_SPACEDIM, SYM>> const& cdata = crse.array(mfi);

            const Dim3 lo = amrex::lbound(bx), hi = amrex::ubound(bx);
            /*const Dim3 lo = amrex::lbound(cdomain), hi = amrex::ubound(cdomain);*/

            // I,J,K == coarse coordinates
            // i,j,k == fine coordinates
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int I, int J, int K) {
                int i = 2 * I, j = 2 * J, k = 2 * K;

                if ((I == lo.x || I == hi.x) &&
                    (J == lo.y || J == hi.y) &&
                    (K == lo.z || K == hi.z)) // Corner
                    cdata(I, J, K) = fdata(i, j, k);
                else if ((J == lo.y || J == hi.y) &&
                    (K == lo.z || K == hi.z)) // X edge
                    cdata(I, J, K) = fdata(i - 1, j, k) * 0.25 + fdata(i, j, k) * 0.5 + fdata(i + 1, j, k) * 0.25;
                else if ((K == lo.z || K == hi.z) &&
                    (I == lo.x || I == hi.x)) // Y edge
                    cdata(I, J, K) = fdata(i, j - 1, k) * 0.25 + fdata(i, j, k) * 0.5 + fdata(i, j + 1, k) * 0.25;
                else if ((I == lo.x || I == hi.x) &&
                    (J == lo.y || J == hi.y)) // Z edge
                    cdata(I, J, K) = fdata(i, j, k - 1) * 0.25 + fdata(i, j, k) * 0.5 + fdata(i, j, k + 1) * 0.25;
                else if (I == lo.x || I == hi.x) // X face
                    cdata(I, J, K) =
                    (fdata(i, j - 1, k - 1) + fdata(i, j, k - 1) * 2.0 + fdata(i, j + 1, k - 1)
                        + fdata(i, j - 1, k) * 2.0 + fdata(i, j, k) * 4.0 + fdata(i, j + 1, k) * 2.0
                        + fdata(i, j - 1, k + 1) + fdata(i, j, k + 1) * 2.0 + fdata(i, j + 1, k + 1)) / 16.0;
                else if (J == lo.y || J == hi.y) // Y face
                    cdata(I, J, K) =
                    (fdata(i - 1, j, k - 1) + fdata(i - 1, j, k) * 2.0 + fdata(i - 1, j, k + 1)
                        + fdata(i, j, k - 1) * 2.0 + fdata(i, j, k) * 4.0 + fdata(i, j, k + 1) * 2.0
                        + fdata(i + 1, j, k - 1) + fdata(i + 1, j, k) * 2.0 + fdata(i + 1, j, k + 1)) / 16.0;
                else if (K == lo.z || K == hi.z) // Z face
                    cdata(I, J, K) =
                    (fdata(i - 1, j - 1, k) + fdata(i, j - 1, k) * 2.0 + fdata(i + 1, j - 1, k)
                        + fdata(i - 1, j, k) * 2.0 + fdata(i, j, k) * 4.0 + fdata(i + 1, j, k) * 2.0
                        + fdata(i - 1, j + 1, k) + fdata(i, j + 1, k) * 2.0 + fdata(i + 1, j + 1, k)) / 16.0;
                else // Interior
                    cdata(I, J, K) =
                    (fdata(i - 1, j - 1, k - 1) + fdata(i - 1, j - 1, k + 1) + fdata(i - 1, j + 1, k - 1) + fdata(i - 1, j + 1, k + 1) +
                        fdata(i + 1, j - 1, k - 1) + fdata(i + 1, j - 1, k + 1) + fdata(i + 1, j + 1, k - 1) + fdata(i + 1, j + 1, k + 1)) / 64.0
                    +
                    (fdata(i, j - 1, k - 1) + fdata(i, j - 1, k + 1) + fdata(i, j + 1, k - 1) + fdata(i, j + 1, k + 1) +
                        fdata(i - 1, j, k - 1) + fdata(i + 1, j, k - 1) + fdata(i - 1, j, k + 1) + fdata(i + 1, j, k + 1) +
                        fdata(i - 1, j - 1, k) + fdata(i - 1, j + 1, k) + fdata(i + 1, j - 1, k) + fdata(i + 1, j + 1, k)) / 32.0
                    +
                    (fdata(i - 1, j, k) + fdata(i, j - 1, k) + fdata(i, j, k - 1) +
                        fdata(i + 1, j, k) + fdata(i, j + 1, k) + fdata(i, j, k + 1)) / 16.0
                    +
                    fdata(i, j, k) / 8.0;

#ifdef AMREX_DEBUG
#ifndef ALAMO_GPU
                if (cdata(I, J, K).contains_nan()) Util::Abort(INFO, "restricted model is nan at crse coordinates (I=", I, ",J=", J, ",K=", k, "), amrlev=", amrlev, " interpolating from mglev", mglev - 1, " to ", mglev);
#endif
#endif
            });
        }
        FillBoundaryCoeff(crse, Geom(amrlev,mglev).periodicity());


        if (!m_psi_set) continue;

        amrex::Box cdomain_cell(m_geom[amrlev][mglev].Domain());
        amrex::Box fdomain_cell(m_geom[amrlev][mglev - 1].Domain());
        MultiFab& crse_psi = *m_psi_mf[amrlev][mglev];
        MultiFab& fine_psi = *m_psi_mf[amrlev][mglev - 1];
        MultiFab fine_psi_on_crseba;
        fine_psi_on_crseba.define(newba.convert(amrex::IntVect::TheCellVector()), crse_psi.DistributionMap(), 1, 1);
        fine_psi_on_crseba.ParallelCopy(fine_psi, 0, 0, 1, 1, 1, m_geom[amrlev][mglev].periodicity());

        for (MFIter mfi(crse_psi, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            Box bx = mfi.tilebox();
            bx = bx & cdomain_cell;

            amrex::Array4<const Set::Scalar> const& fdata = fine_psi_on_crseba.array(mfi);
            amrex::Array4<Set::Scalar> const& cdata = crse_psi.array(mfi);

            const Dim3 lo = amrex::lbound(cdomain), hi = amrex::ubound(cdomain);

            // I,J,K == coarse coordinates
            // i,j,k == fine coordinates
            amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int I, int J, int K) {
                int i = 2 * I, j = 2 * J, k = 2 * K;

                if ((I == lo.x || I == hi.x) &&
                    (J == lo.y || J == hi.y) &&
                    (K == lo.z || K == hi.z)) // Corner
                    cdata(I, J, K) = fdata(i, j, k);
                else if ((J == lo.y || J == hi.y) &&
                    (K == lo.z || K == hi.z)) // X edge
                    cdata(I, J, K) = fdata(i - 1, j, k) * 0.25 + fdata(i, j, k) * 0.5 + fdata(i + 1, j, k) * 0.25;
                else if ((K == lo.z || K == hi.z) &&
                    (I == lo.x || I == hi.x)) // Y edge
                    cdata(I, J, K) = fdata(i, j - 1, k) * 0.25 + fdata(i, j, k) * 0.5 + fdata(i, j + 1, k) * 0.25;
                else if ((I == lo.x || I == hi.x) &&
                    (J == lo.y || J == hi.y)) // Z edge
                    cdata(I, J, K) = fdata(i, j, k - 1) * 0.25 + fdata(i, j, k) * 0.5 + fdata(i, j, k + 1) * 0.25;
                else if (I == lo.x || I == hi.x) // X face
                    cdata(I, J, K) =
                    (fdata(i, j - 1, k - 1) + fdata(i, j, k - 1) * 2.0 + fdata(i, j + 1, k - 1)
                        + fdata(i, j - 1, k) * 2.0 + fdata(i, j, k) * 4.0 + fdata(i, j + 1, k) * 2.0
                        + fdata(i, j - 1, k + 1) + fdata(i, j, k + 1) * 2.0 + fdata(i, j + 1, k + 1)) / 16.0;
                else if (J == lo.y || J == hi.y) // Y face
                    cdata(I, J, K) =
                    (fdata(i - 1, j, k - 1) + fdata(i - 1, j, k) * 2.0 + fdata(i - 1, j, k + 1)
                        + fdata(i, j, k - 1) * 2.0 + fdata(i, j, k) * 4.0 + fdata(i, j, k + 1) * 2.0
                        + fdata(i + 1, j, k - 1) + fdata(i + 1, j, k) * 2.0 + fdata(i + 1, j, k + 1)) / 16.0;
                else if (K == lo.z || K == hi.z) // Z face
                    cdata(I, J, K) =
                    (fdata(i - 1, j - 1, k) + fdata(i, j - 1, k) * 2.0 + fdata(i + 1, j - 1, k)
                        + fdata(i - 1, j, k) * 2.0 + fdata(i, j, k) * 4.0 + fdata(i + 1, j, k) * 2.0
                        + fdata(i - 1, j + 1, k) + fdata(i, j + 1, k) * 2.0 + fdata(i + 1, j + 1, k)) / 16.0;
                else // Interior
                    cdata(I, J, K) =
                    (fdata(i - 1, j - 1, k - 1) + fdata(i - 1, j - 1, k + 1) + fdata(i - 1, j + 1, k - 1) + fdata(i - 1, j + 1, k + 1) +
                        fdata(i + 1, j - 1, k - 1) + fdata(i + 1, j - 1, k + 1) + fdata(i + 1, j + 1, k - 1) + fdata(i + 1, j + 1, k + 1)) / 64.0
                    +
                    (fdata(i, j - 1, k - 1) + fdata(i, j - 1, k + 1) + fdata(i, j + 1, k - 1) + fdata(i, j + 1, k + 1) +
                        fdata(i - 1, j, k - 1) + fdata(i + 1, j, k - 1) + fdata(i - 1, j, k + 1) + fdata(i + 1, j, k + 1) +
                        fdata(i - 1, j - 1, k) + fdata(i - 1, j + 1, k) + fdata(i + 1, j - 1, k) + fdata(i + 1, j + 1, k)) / 32.0
                    +
                    (fdata(i - 1, j, k) + fdata(i, j - 1, k) + fdata(i, j, k - 1) +
                        fdata(i + 1, j, k) + fdata(i, j + 1, k) + fdata(i, j, k + 1)) / 16.0
                    +
                    fdata(i, j, k) / 8.0;
            });
        }
        FillBoundaryCoeff(crse_psi, Geom(amrlev,mglev).periodicity());

    }
}

template<int SYM>
void
Elastic<SYM>::FillBoundaryCoeff(MultiTab& sigma, const amrex::Periodicity& p)
{
    sigma.setMultiGhost(true);
    sigma.FillBoundaryAndSync(p);
}

template<int SYM>
void
Elastic<SYM>::FillBoundaryCoeff(MultiFab& psi, const amrex::Periodicity& p)
{
    psi.setMultiGhost(true);
    psi.FillBoundaryAndSync(p);
}

template class Elastic<Set::Sym::Major>;
template class Elastic<Set::Sym::Isotropic>;
template class Elastic<Set::Sym::MajorMinor>;
template class Elastic<Set::Sym::Diagonal>;

}
