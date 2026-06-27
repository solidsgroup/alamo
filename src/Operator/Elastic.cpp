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

    Set::Vector DX(m_geom[amrlev][mglev].CellSize()); // device-safe cell size (copied into device closure)
#ifdef AMREX_DEBUG
#ifdef ALAMO_GPU
    Util::DeviceErrorFlag fapply_error;
    int* fapply_error_flag = fapply_error.dataPtr();
#endif
#endif

    // ALAMO_ML_FAPPLY_PROBE=<mglev> (temporary diagnostic): on the first call
    // at the requested mglev, override the input field with an analytic
    // displacement field so the expected output can be computed in closed
    // form, then dump a boundary-vs-interior norm summary and stop. Lets us
    // check a single Fapply evaluation at a coarsened MG level without
    // running a full multi-hundred-iteration solve.
    //
    // ALAMO_ML_FAPPLY_PROBE_MODE selects the field:
    //   "linear" (default) - component 0 = global x-coordinate, ramped to a
    //     physically-representative displacement scale (~1e-6 m over the
    //     whole domain). gradgradu is exactly zero analytically, so this
    //     exercises the derivative-stencil/one-sided-difference selection
    //     logic but lies in the operator's near-null-space - it CANNOT see a
    //     defect that only manifests on rough/high-frequency inputs.
    //   "sine" - component 0 = amp*sin(ax*i)*sin(ay*j), a high-frequency
    //     field whose discrete second-difference/mixed-derivative stencil
    //     output has an exact closed form (sampled sinusoids are eigenfunctions
    //     of the standard central-difference stencils used here), so Fapply's
    //     actual output can be diffed against a true expected value rather
    //     than just "should be near zero". This is the discriminating test:
    //     it is the only mode that can see a defect invisible to constant/
    //     linear inputs. Valid as an exact ground truth only when the operator
    //     is uniform and psi is disabled (matches the documented reproduction
    //     baseline in benchmark/elastic_sensitivity_20260621/); a warning is
    //     printed if that does not hold.
    static bool fapply_probe_fired = false;
    const char* fapply_probe_env = std::getenv("ALAMO_ML_FAPPLY_PROBE");
    const char* fapply_probe_mode_env = std::getenv("ALAMO_ML_FAPPLY_PROBE_MODE");
    const bool fapply_probe_sine_mode =
        fapply_probe_mode_env && std::strcmp(fapply_probe_mode_env, "sine") == 0;
    const bool fapply_probe_this_call =
        fapply_probe_env && !fapply_probe_fired && mglev == std::atoi(fapply_probe_env);
    amrex::Real probe_sine_ax = 0.0, probe_sine_ay = 0.0, probe_sine_amp = 0.0;
    if (fapply_probe_this_call)
    {
        fapply_probe_fired = true;
        MultiFab& u_mut = const_cast<MultiFab&>(a_u);
        if (u_mut.nComp() > 1) u_mut.setVal(0.0, 1, u_mut.nComp() - 1, u_mut.nGrowVect());

        if (fapply_probe_sine_mode)
        {
            // Wavelength ~4 cells in each direction: high-frequency relative
            // to the linear probe, but not exactly Nyquist (which would zero
            // out the field on integer indices).
            const amrex::Box probe_cell_domain = m_geom[amrlev][mglev].Domain();
            const int Nx = probe_cell_domain.length(0);
            const int Ny = probe_cell_domain.length(1);
            const int mx = std::max(1, Nx / 4);
            const int my = std::max(1, Ny / 4);
            probe_sine_ax = 2.0 * 3.14159265358979323846 * mx / (amrex::Real)Nx;
            probe_sine_ay = 2.0 * 3.14159265358979323846 * my / (amrex::Real)Ny;
            probe_sine_amp = 1e-6;
            const amrex::Real ax = probe_sine_ax, ay = probe_sine_ay, amp = probe_sine_amp;

            for (amrex::MFIter mfi(u_mut, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box bx = mfi.growntilebox(u_mut.nGrowVect());
                auto const& ufab = u_mut.array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    ufab(i, j, k, 0) = amp * std::sin(ax * i) * std::sin(ay * j);
                });
            }

            amrex::ReduceOps<amrex::ReduceOpMax> verify_op;
            amrex::ReduceData<amrex::Real> verify_data(verify_op);
            using VerifyTuple = typename decltype(verify_data)::Type;
            for (amrex::MFIter mfi(u_mut, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box bx = mfi.growntilebox(u_mut.nGrowVect());
                auto const& ufab = u_mut.const_array(mfi);
                verify_op.eval(bx, verify_data,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) -> VerifyTuple
                    {
                        return { std::abs(ufab(i, j, k, 0) - amp * std::sin(ax * i) * std::sin(ay * j)) };
                    });
            }
            amrex::Real verify_max_diff = amrex::get<0>(verify_data.value(verify_op));
            amrex::ParallelDescriptor::ReduceRealMax(verify_max_diff);
            Util::Message(INFO, "ALAMO_ML_FAPPLY_PROBE(sine) override self-consistency check: max|u-expected|=",
                verify_max_diff, " ax=", ax, " ay=", ay, " amp=", amp, " (expect exactly 0)");
            if (!m_uniform || m_psi_set)
                Util::Message(INFO, "ALAMO_ML_FAPPLY_PROBE(sine) WARNING: m_uniform=", m_uniform,
                    " m_psi_set=", m_psi_set, " - the closed-form comparison below is only an exact ground",
                    " truth for a uniform operator with psi disabled; treat results outside that regime as",
                    " indicative only.");
        }
        else
        {
            // Linear field (component 0 = global x-coordinate): gradgradu is exactly
            // zero analytically, but unlike a constant field this exercises the
            // actual derivative-stencil/one-sided-difference selection logic.
            // Scale the linear ramp to a physically-representative displacement
            // magnitude (~1e-6 m over the whole domain) instead of O(domain size)
            // coordinate values, to avoid 1/dx^2-amplified round-off in the
            // second-derivative stencil swamping the analytic zero signal.
            const amrex::Real probe_x0 = 0.0;
            const amrex::Real probe_dx0 = 1e-6 * m_geom[amrlev][mglev].CellSize(0) / m_geom[amrlev][mglev].ProbLength(0);
            for (amrex::MFIter mfi(u_mut, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box bx = mfi.growntilebox(u_mut.nGrowVect());
                auto const& ufab = u_mut.array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    ufab(i, j, k, 0) = probe_x0 + i * probe_dx0;
                });
            }

            // Verify the override itself is globally self-consistent (every box's
            // own ghost storage matches the same global linear formula) before
            // trusting any Fapply discrepancy as a real bug.
            amrex::ReduceOps<amrex::ReduceOpMax> verify_op;
            amrex::ReduceData<amrex::Real> verify_data(verify_op);
            using VerifyTuple = typename decltype(verify_data)::Type;
            for (amrex::MFIter mfi(u_mut, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box bx = mfi.growntilebox(u_mut.nGrowVect());
                auto const& ufab = u_mut.const_array(mfi);
                verify_op.eval(bx, verify_data,
                    [=] AMREX_GPU_DEVICE(int i, int j, int k) -> VerifyTuple
                    {
                        return { std::abs(ufab(i, j, k, 0) - (probe_x0 + i * probe_dx0)) };
                    });
            }
            amrex::Real verify_max_diff = amrex::get<0>(verify_data.value(verify_op));
            amrex::ParallelDescriptor::ReduceRealMax(verify_max_diff);
            Util::Message(INFO, "ALAMO_ML_FAPPLY_PROBE override self-consistency check: max|u-expected|=", verify_max_diff, " (expect exactly 0)");
        }
    }

    amrex::MultiFab fapply_probe_gg;
    if (fapply_probe_this_call)
        fapply_probe_gg.define(a_f.boxArray(), a_f.DistributionMap(), 1, a_f.nGrowVect());

    amrex::MultiFab fapply_probe_expected;
    if (fapply_probe_this_call && fapply_probe_sine_mode)
        fapply_probe_expected.define(a_f.boxArray(), a_f.DistributionMap(), AMREX_SPACEDIM, a_f.nGrowVect());

    for (MFIter mfi(a_f, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
    {
        Box bx = mfi.validbox().grow(1) & domain;
        amrex::Box tilebox = mfi.grownnodaltilebox() & bx;

        amrex::Array4<MATRIX4> const& DDW = (*(m_ddw_mf[amrlev][mglev])).array(mfi);
        amrex::Array4<const amrex::Real> const& U = a_u.array(mfi);
        amrex::Array4<amrex::Real> const& F = a_f.array(mfi);
        amrex::Array4<Set::Scalar> const& psi = m_psi_mf[amrlev][mglev]->array(mfi);
        amrex::Array4<amrex::Real> const& probe_gg =
            fapply_probe_this_call ? fapply_probe_gg.array(mfi) : amrex::Array4<amrex::Real>();
        amrex::Array4<amrex::Real> const& probe_expected =
            (fapply_probe_this_call && fapply_probe_sine_mode) ? fapply_probe_expected.array(mfi) : amrex::Array4<amrex::Real>();

        const Dim3 lo = amrex::lbound(stencilbox), hi = amrex::ubound(stencilbox);
        auto m_psi_set = this->m_psi_set;
        auto m_psi_small = this->m_psi_small;
        auto m_uniform = this->m_uniform;
        const bool probe_capture_gg = fapply_probe_this_call;
        const bool probe_sine_this_call = fapply_probe_this_call && fapply_probe_sine_mode;
        const amrex::Real probe_ax = probe_sine_ax, probe_ay = probe_sine_ay, probe_amp = probe_sine_amp;
#ifdef ALAMO_GPU
        auto m_bc_type = this->m_bc->GetBcTypeArray();
#else
        auto m_bc = this->m_bc;
#endif

        ALAMO_ELASTIC_OP_FOR(tilebox, ALAMO_ELASTIC_OP_CAPTURE ALAMO_ELASTIC_OP_DEVICE (int i, int j, int k) {

            Set::Vector f = Set::Vector::Zero();

            Set::Vector u;
            for (int p = 0; p < AMREX_SPACEDIM; p++) u(p) = U(i, j, k, p);


            bool AMREX_D_DECL(xmin = (i == lo.x), ymin = (j == lo.y), zmin = (k == lo.z)),
                AMREX_D_DECL(xmax = (i == hi.x), ymax = (j == hi.y), zmax = (k == hi.z));

            // Determine if a special stencil will be necessary for first derivatives
            std::array<Numeric::StencilType, AMREX_SPACEDIM>
                sten = Numeric::GetStencil(i, j, k, stencilbox);

            // The displacement gradient tensor
            Set::Matrix gradu; // gradu(i,j) = u_{i,j)

            // Fill gradu
            for (int p = 0; p < AMREX_SPACEDIM; p++)
            {
                AMREX_D_TERM(gradu(p, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(U, i, j, k, p, DX.data(), sten));,
                    gradu(p, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(U, i, j, k, p, DX.data(), sten));,
                    gradu(p, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(U, i, j, k, p, DX.data(), sten)););
            }

            Set::Scalar psi_avg = 1.0;
            if (m_psi_set) psi_avg = (1.0 - m_psi_small) * Numeric::Interpolate::CellToNodeAverage(psi, i, j, k, 0) + m_psi_small;

            // Stress tensor computed using the model fab
            Set::Matrix sig = (DDW(i, j, k) * gradu) * psi_avg;

            // Boundary conditions
            /// \todo Important: we need a way to handle corners and edges.
            amrex::IntVect m(AMREX_D_DECL(i, j, k));
            if (AMREX_D_TERM(xmax || xmin, || ymax || ymin, || zmax || zmin))
            {
                f = ALAMO_ELASTIC_OP_BC_EVAL(m_bc, m_bc_type, u, gradu, sig, i, j, k, stencilbox);
            }
            else
            {


                // The gradient of the displacement gradient tensor
                // TODO - replace with this call. But not for this PR
                //Set::Matrix3 gradgradu = Numeric::Hessian(U,i,j,k,DX,sten); // gradgradu[k](l,j) = u_{k,lj}
                Set::Matrix3 gradgradu; // gradgradu[k](l,j) = u_{k,lj}

                // Fill gradu and gradgradu
                for (int p = 0; p < AMREX_SPACEDIM; p++)
                {
                    // Diagonal terms:
                    AMREX_D_TERM(
                        gradgradu(p, 0, 0) = (Numeric::Stencil<Set::Scalar, 2, 0, 0>::D(U, i, j, k, p, DX.data()));,
                        gradgradu(p, 1, 1) = (Numeric::Stencil<Set::Scalar, 0, 2, 0>::D(U, i, j, k, p, DX.data()));,
                        gradgradu(p, 2, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 2>::D(U, i, j, k, p, DX.data())););

                    // Off-diagonal terms:
                    AMREX_D_TERM(
                        ,// 2D
                        gradgradu(p, 0, 1) = (Numeric::Stencil<Set::Scalar, 1, 1, 0>::D(U, i, j, k, p, DX.data()));
                        gradgradu(p, 1, 0) = gradgradu(p, 0, 1);
                        ,// 3D
                        gradgradu(p, 0, 2) = (Numeric::Stencil<Set::Scalar, 1, 0, 1>::D(U, i, j, k, p, DX.data()));
                        gradgradu(p, 1, 2) = (Numeric::Stencil<Set::Scalar, 0, 1, 1>::D(U, i, j, k, p, DX.data()));
                        gradgradu(p, 2, 0) = gradgradu(p, 0, 2);
                        gradgradu(p, 2, 1) = gradgradu(p, 1, 2););
                }

                //
                // Operator
                //
                // The return value is
                //    f = C(grad grad u) + grad(C)*grad(u)
                // In index notation
                //    f_i = C_{ijkl,j} u_{k,l}  +  C_{ijkl}u_{k,lj}
                //

                if (probe_capture_gg) probe_gg(i, j, k) = gradgradu.norm();

                if (probe_sine_this_call)
                {
                    // Closed-form gradgradu for u_0 = amp*sin(ax*i)*sin(ay*j), u_{p>0} = 0
                    // (exact for these standard central-difference stencils, not an
                    // approximation: a uniformly-sampled sinusoid is an eigenfunction
                    // of both the central second-difference and the central mixed-
                    // derivative stencil). Component 0 is the only nonzero row since
                    // the override leaves every other displacement component at zero
                    // everywhere, including ghost nodes, so all of their derivatives
                    // are exactly zero too.
                    const amrex::Real cx = std::cos(probe_ax * i), sx = std::sin(probe_ax * i);
                    const amrex::Real cy = std::cos(probe_ay * j), sy = std::sin(probe_ay * j);
                    const amrex::Real eigx = 2.0 * (std::cos(probe_ax) - 1.0) / (DX[0] * DX[0]);
                    const amrex::Real eigy = 2.0 * (std::cos(probe_ay) - 1.0) / (DX[1] * DX[1]);
                    const amrex::Real eigxy = std::sin(probe_ax) * std::sin(probe_ay) / (DX[0] * DX[1]);

                    Set::Matrix3 gradgradu_expected = Set::Matrix3::Zero();
                    gradgradu_expected(0, 0, 0) = probe_amp * eigx * sx * sy;
                    gradgradu_expected(0, 1, 1) = probe_amp * eigy * sx * sy;
                    gradgradu_expected(0, 0, 1) = probe_amp * eigxy * cx * cy;
                    gradgradu_expected(0, 1, 0) = gradgradu_expected(0, 0, 1);

                    // Only valid as an exact ground truth when DDW is uniform and psi
                    // is disabled (no grad(C) or grad(psi) correction terms below) -
                    // matches this probe's documented reproduction baseline.
                    Set::Vector f_expected = (DDW(i, j, k) * gradgradu_expected) * psi_avg;
                    AMREX_D_TERM(probe_expected(i, j, k, 0) = f_expected[0];,
                        probe_expected(i, j, k, 1) = f_expected[1];,
                        probe_expected(i, j, k, 2) = f_expected[2];);
                }

                f = (DDW(i, j, k) * gradgradu) * psi_avg;

                if (!m_uniform)
                {
                    MATRIX4
                        AMREX_D_DECL(Cgrad1 = (Numeric::Stencil<MATRIX4, 1, 0, 0>::D(DDW, i, j, k, 0, DX.data(), sten)),
                            Cgrad2 = (Numeric::Stencil<MATRIX4, 0, 1, 0>::D(DDW, i, j, k, 0, DX.data(), sten)),
                            Cgrad3 = (Numeric::Stencil<MATRIX4, 0, 0, 1>::D(DDW, i, j, k, 0, DX.data(), sten)));
                    f += (AMREX_D_TERM((Cgrad1 * gradu).col(0),
                        +(Cgrad2 * gradu).col(1),
                        +(Cgrad3 * gradu).col(2))) * (psi_avg);
                }
                if (m_psi_set)
                {
                    Set::Vector gradpsi = Numeric::CellGradientOnNode(psi, i, j, k, 0, DX.data());
                    gradpsi *= (1.0 - m_psi_small);
                    f += (DDW(i, j, k) * gradu) * gradpsi;
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
                if (std::isnan(f(0)) || std::isnan(f(1)))
                {
                    Util::Message(INFO,"  =================  ");
                    Util::Message(INFO,"amrlev=",amrlev);
                    Util::Message(INFO,"mglev=",mglev);
                    Util::Message(INFO,"i=",i," j=",j);
                    Util::Message(INFO,"f:            ",f.transpose());
                    Util::Message(INFO,"U(i,j,k):     ",U(i,j,k,0)," ",U(i,j,k,1));
                    Util::Message(INFO,"U(i-1,j,k):   ",U(i-1,j,k,0)," ",U(i-1,j,k,1));
                    Util::Message(INFO,"U(i+1,j,k):   ",U(i+1,j,k,0)," ",U(i-1,j,k,1));
                    Util::Message(INFO,"U(i,j-1,k):     ",U(i,j-1,k,0)," ",U(i,j-1,k,1));
                    Util::Message(INFO,"U(i,j+1,k):     ",U(i,j+1,k,0)," ",U(i,j+1,k,1));
                    Util::Message(INFO,"gradu:        ",gradu);
                    Util::Message(INFO,"gradgradu[0]: ",gradgradu[0]);
                    Util::Message(INFO,"gradgradu[1]: ",gradgradu[1]);
                    Util::Message(INFO,"DDW (i  ,j  ): ",DDW(i,j,k));
                    Util::Message(INFO,"DDW (i-1,j  ): ",DDW(i-1,j,k));
                    Util::Message(INFO,"DDW (i+1,j  ): ",DDW(i+1,j,k));
                    Util::Message(INFO,"DDW (i  ,j+1): ",DDW(i,j+1,k));
                    Util::Message(INFO,"DDW (i  ,j-1): ",DDW(i,j-1,k));
                    Util::Message(INFO,"psi_av: ",psi_avg);
                    Util::Message(INFO,"psi_set: ",m_psi_set);
                    Util::Message(INFO,"  =================  ");
                    Util::Abort(INFO);
                }
#endif
#endif
            }
            AMREX_D_TERM(F(i, j, k, 0) = f[0];, F(i, j, k, 1) = f[1];, F(i, j, k, 2) = f[2];);
        });
#ifdef AMREX_DEBUG
#ifdef ALAMO_GPU
        Util::AbortIfDeviceError(fapply_error, INFO, "Operator::Elastic::Fapply() detected an invalid value on device");
#endif
#endif
    }

    if (fapply_probe_this_call)
    {
        amrex::Box probe_domain(m_geom[amrlev][mglev].growPeriodicDomain(1));
        probe_domain.convert(amrex::IntVect::TheNodeVector());
        const int dlo0 = stencilbox.smallEnd(0), dhi0 = stencilbox.bigEnd(0);
        const int dlo1 = stencilbox.smallEnd(1), dhi1 = stencilbox.bigEnd(1);

        amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpMax> reduce_op;
        amrex::ReduceData<amrex::Real, amrex::Real> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (amrex::MFIter mfi(a_f, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const amrex::Box& bx = mfi.validbox() & probe_domain;
            auto const& fab = a_f.const_array(mfi);
            const int ncomp = a_f.nComp();

            reduce_op.eval(bx, reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                {
                    amrex::Real maxabs = 0.0;
                    for (int n = 0; n < ncomp; ++n) maxabs = amrex::max(maxabs, std::abs(fab(i, j, k, n)));
                    const bool on_boundary = (i == dlo0 || i == dhi0 || j == dlo1 || j == dhi1);
                    return { on_boundary ? maxabs : 0.0, (!on_boundary) ? maxabs : 0.0 };
                });
        }

        ReduceTuple hv = reduce_data.value(reduce_op);
        amrex::Real boundary_max_abs = amrex::get<0>(hv);
        amrex::Real interior_max_abs = amrex::get<1>(hv);
        amrex::ParallelDescriptor::ReduceRealMax(boundary_max_abs);
        amrex::ParallelDescriptor::ReduceRealMax(interior_max_abs);

        if (!fapply_probe_sine_mode)
        {
            amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpMax> gg_reduce_op;
            amrex::ReduceData<amrex::Real, amrex::Real> gg_reduce_data(gg_reduce_op);
            using GgReduceTuple = typename decltype(gg_reduce_data)::Type;
            for (amrex::MFIter mfi(fapply_probe_gg, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox() & probe_domain;
                auto const& ggfab = fapply_probe_gg.const_array(mfi);
                gg_reduce_op.eval(bx, gg_reduce_data,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) -> GgReduceTuple
                    {
                        const bool on_boundary = (i == dlo0 || i == dhi0 || j == dlo1 || j == dhi1);
                        return { on_boundary ? ggfab(i, j, k) : 0.0, (!on_boundary) ? ggfab(i, j, k) : 0.0 };
                    });
            }
            GgReduceTuple gghv = gg_reduce_data.value(gg_reduce_op);
            amrex::Real gg_boundary_max = amrex::get<0>(gghv);
            amrex::Real gg_interior_max = amrex::get<1>(gghv);
            amrex::ParallelDescriptor::ReduceRealMax(gg_boundary_max);
            amrex::ParallelDescriptor::ReduceRealMax(gg_interior_max);
            Util::Message(INFO, "ALAMO_ML_FAPPLY_PROBE gradgradu.norm() interior_max=", gg_interior_max,
                " (expect ~0 for a linear field) boundary_max=", gg_boundary_max);

            Util::Message(INFO, "ALAMO_ML_FAPPLY_PROBE amrlev=", amrlev, " mglev=", mglev,
                " constant-input-field test: interior_max_abs=", interior_max_abs,
                " (expect ~0) boundary_max_abs=", boundary_max_abs, " (BC nodes, ignore)");
            amrex::Abort("ALAMO_ML_FAPPLY_PROBE: probe complete, stopping run");
        }
        else
        {
            // THE discriminating test (priority 1, GPU_ELASTIC_DEBUG_PLAN.md): diff
            // Fapply's actual output against the closed-form expected output computed
            // from the sinusoidal field above. Unlike the constant/linear probes,
            // this field is not in the operator's near-null-space, so a coarse-level
            // operator-application defect that only manifests on rough inputs will
            // show up here as a large interior_max_abs_diff.
            amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpMax> diff_reduce_op;
            amrex::ReduceData<amrex::Real, amrex::Real> diff_reduce_data(diff_reduce_op);
            using DiffReduceTuple = typename decltype(diff_reduce_data)::Type;
            for (amrex::MFIter mfi(a_f, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox() & probe_domain;
                auto const& fab = a_f.const_array(mfi);
                auto const& expfab = fapply_probe_expected.const_array(mfi);
                const int ncomp = a_f.nComp();
                diff_reduce_op.eval(bx, diff_reduce_data,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) -> DiffReduceTuple
                    {
                        amrex::Real maxdiff = 0.0;
                        for (int n = 0; n < ncomp; ++n)
                            maxdiff = amrex::max(maxdiff, std::abs(fab(i, j, k, n) - expfab(i, j, k, n)));
                        const bool on_boundary = (i == dlo0 || i == dhi0 || j == dlo1 || j == dhi1);
                        return { on_boundary ? maxdiff : 0.0, (!on_boundary) ? maxdiff : 0.0 };
                    });
            }
            DiffReduceTuple dhv = diff_reduce_data.value(diff_reduce_op);
            amrex::Real boundary_max_diff = amrex::get<0>(dhv);
            amrex::Real interior_max_diff = amrex::get<1>(dhv);
            amrex::ParallelDescriptor::ReduceRealMax(boundary_max_diff);
            amrex::ParallelDescriptor::ReduceRealMax(interior_max_diff);

            Util::Message(INFO, "ALAMO_ML_FAPPLY_PROBE(sine) amrlev=", amrlev, " mglev=", mglev,
                " ax=", probe_sine_ax, " ay=", probe_sine_ay, " amp=", probe_sine_amp,
                " high-frequency closed-form test: interior_max_abs=", interior_max_abs,
                " interior_max_abs_diff(F-expected)=", interior_max_diff,
                " relative=", interior_max_abs > 0.0 ? interior_max_diff / interior_max_abs : 0.0,
                " (expect relative ~1e-12, roundoff only) boundary_max_abs_diff=", boundary_max_diff, " (BC nodes, ignore)");
            amrex::Abort("ALAMO_ML_FAPPLY_PROBE: sine probe complete, stopping run");
        }
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

    Set::Vector DX(m_geom[amrlev][mglev].CellSize()); // device-safe cell size (copied into device closure)
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
        amrex::Array4<Set::Scalar> const& psi = m_psi_mf[amrlev][mglev]->array(mfi);

        const Dim3 lo = amrex::lbound(stencilbox), hi = amrex::ubound(stencilbox);
        auto m_psi_set = this->m_psi_set;
        auto m_psi_small = this->m_psi_small;
#ifdef ALAMO_GPU
        auto m_bc_type = this->m_bc->GetBcTypeArray();
#else
        auto m_bc = this->m_bc;
#endif

        amrex::ParallelFor(tilebox, [=] AMREX_GPU_DEVICE(int i, int j, int k)
        {

            // Determine if a special stencil will be necessary for first derivatives
            std::array<Numeric::StencilType, AMREX_SPACEDIM>
                sten = Numeric::GetStencil(i, j, k, stencilbox);

            // gradu(i,j) = u_{i,j)
            std::array<Set::Matrix,AMREX_SPACEDIM> gradu = Numeric::Gradient_Diagonal<Set::Matrix>(DX.data(), sten);  

            // gradgradu[k](l,j) = u_{k,lj}
            std::array<Set::Matrix3,AMREX_SPACEDIM>  gradgradu = Numeric::Gradient_Diagonal<Set::Matrix3>(DX.data()); 


            Set::Vector f = Set::Vector::Zero();

            bool    
                AMREX_D_DECL(xmin = (i == lo.x), ymin = (j == lo.y), zmin = (k == lo.z)),
                AMREX_D_DECL(xmax = (i == hi.x), ymax = (j == hi.y), zmax = (k == hi.z));

            Set::Scalar psi_avg = 1.0;
            if (m_psi_set) psi_avg = (1.0 - m_psi_small) * Numeric::Interpolate::CellToNodeAverage(psi, i, j, k, 0) + m_psi_small;



            for (int p = 0; p < AMREX_SPACEDIM; p++)
            {

                diag(i, j, k, p) = 0.0;


                amrex::IntVect m(AMREX_D_DECL(i, j, k));
                if (AMREX_D_TERM(xmax || xmin, || ymax || ymin, || zmax || zmin))
                {
                    Set::Matrix sig = DDW(i, j, k) * gradu[p] * psi_avg;
                    Set::Vector u = Set::Vector::Zero();
                    u(p) = 1.0;
                    f = ALAMO_ELASTIC_OP_BC_EVAL(m_bc, m_bc_type, u, gradu[p], sig, i, j, k, stencilbox);
                    diag(i, j, k, p) = f(p);
                }
                else
                {
                    Set::Vector f = (DDW(i, j, k) * gradgradu[p]) * psi_avg;
                    diag(i, j, k, p) += f(p);
                }

#ifdef AMREX_DEBUG
#ifdef ALAMO_GPU
                if (std::isnan(diag(i, j, k, p)) || std::isinf(diag(i, j, k, p)) || diag(i, j, k, p) == 0)
                {
                    Util::SetDeviceError(diagonal_error_flag);
                }
#else
                if (std::isnan(diag(i, j, k, p))) Util::Abort(INFO, "diagonal is nan at (", i, ",", j, ",", k, "), amrlev=", amrlev, ", mglev=", mglev);
                if (std::isinf(diag(i, j, k, p))) Util::Abort(INFO, "diagonal is inf at (", i, ",", j, ",", k, "), amrlev=", amrlev, ", mglev=", mglev);
                if (diag(i, j, k, p) == 0) Util::Abort(INFO, "diagonal is zero at (", i, ",", j, ",", k, "), amrlev=", amrlev, ", mglev=", mglev);
#endif
#endif

            }
        });
#ifdef AMREX_DEBUG
#ifdef ALAMO_GPU
        Util::AbortIfDeviceError(diagonal_error, INFO, "Operator::Elastic::Diagonal() detected an invalid value on device");
#endif
#endif
    }

    if (AlamoDiagProbeLimit()) AlamoPrintDiagProbe(amrlev, mglev, "presync", a_diag);

    a_diag.FillBoundaryAndSync(Geom(amrlev,mglev).periodicity());
    nodalSync(amrlev,mglev,a_diag);

    if (AlamoDiagProbeLimit()) AlamoPrintDiagProbe(amrlev, mglev, "postsync", a_diag);
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
