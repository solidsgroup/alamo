#ifndef SOLVER_NONLOCAL_NEWTON
#define SOLVER_NONLOCAL_NEWTON

#include <cmath>
#include <limits>
#include <string>
#include <AMReX_Math.H>
#include <AMReX_ParallelContext.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Reduce.H>
#include "Set/Set.H"
#include "Operator/Elastic.H"
#include "Solver/Nonlocal/Linear.H"
#include "IO/ParmParse.H"
#include "Numeric/Stencil.H"

#ifdef ALAMO_GPU
#define ALAMO_NEWTON_FOR amrex::ParallelFor
#define ALAMO_NEWTON_DEVICE AMREX_GPU_DEVICE
#define ALAMO_NEWTON_BC_EVAL(elastic, bc_type, u, gradu, sigma, i, j, k, bx) \
    ::BC::Operator::Elastic::Elastic::eval(bc_type, u, gradu, sigma, i, j, k, bx)
#else
#define ALAMO_NEWTON_FOR amrex::LoopConcurrentOnCpu
#define ALAMO_NEWTON_DEVICE
#define ALAMO_NEWTON_BC_EVAL(elastic, bc_type, u, gradu, sigma, i, j, k, bx) \
    (elastic)->GetBC()(u, gradu, sigma, i, j, k, bx)
#endif

namespace Solver
{
namespace Nonlocal
{
template <typename T>
class Newton: public Linear
{
public:
    Newton() {};

    Newton(Operator::Elastic<T::sym>& a_op)
    {
        this->Define(a_op);
    };

    ~Newton()
    {
        if (m_defined) Clear();
    }

    void Define(Operator::Elastic<T::sym>& a_op)
    {
        Linear::Define(a_op);
        m_elastic = dynamic_cast<Operator::Elastic<T::sym> *>(linop);
        //m_bc = &m_elastic->GetBC();
    }
    void Clear()
    {
        Linear::Clear();
        m_elastic = nullptr;
        //m_bc = nullptr;
    }


    void setNRIters(int a_nriters) { m_nriters = a_nriters; }

    void setPsi(Set::Field<Set::Scalar>& a_psi)
    {
        m_psi = &a_psi;
    }

    void prepareForSolve(const Set::Field<Set::Scalar>& a_u_mf,
        const Set::Field<Set::Scalar>& a_b_mf,
        Set::Field<Set::Scalar>& a_rhs_mf,
        Set::Field<Set::Matrix>& a_dw_mf,
        Set::Field<Set::Matrix4<AMREX_SPACEDIM, T::sym>>& a_ddw_mf,
        Set::Field<T>& a_model_mf)
    {
        for (int lev = 0; lev <= a_b_mf.finest_level; ++lev)
        {
            amrex::Box domain(linop->Geom(lev).growPeriodicDomain(2));
            domain.convert(amrex::IntVect::TheNodeVector());
            // Device-safe cell size (see note on the Set::Vector overload below).
            Set::Vector DX(linop->Geom(lev).CellSize());
            const std::array<bool, AMREX_SPACEDIM> periodic = {
                AMREX_D_DECL(linop->Geom(lev).isPeriodic(0),
                    linop->Geom(lev).isPeriodic(1),
                    linop->Geom(lev).isPeriodic(2))};
#ifdef ALAMO_GPU
            Util::DeviceErrorFlag prepare_error;
            int* prepare_error_flag = prepare_error.dataPtr();
#endif
            const amrex::IntVect required_ghost(2);
            if (!a_u_mf[lev]->nGrowVect().allGE(required_ghost)
                || !a_model_mf[lev]->nGrowVect().allGE(required_ghost)
                || !a_dw_mf[lev]->nGrowVect().allGE(required_ghost)
                || !a_ddw_mf[lev]->nGrowVect().allGE(required_ghost))
                Util::Abort(INFO, "The paired elastic stencil requires two coefficient, model, and displacement ghosts.");

            for (MFIter mfi(*a_model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                // D- G+ needs flux anchors through two low and one high
                // ghosts.  A fine C/F model ghost is not guaranteed to be
                // initialized, so clamp only the material selection to the
                // nearest valid node.  Crucially, derive DW/DDW from the
                // *actual* paired gradient at that ghost anchor; copying an
                // already-derived DW would make the residual's derivative
                // disagree with Fapply.
                const amrex::Box valid = mfi.validbox();
                const amrex::Dim3 valid_lo = amrex::lbound(valid);
                const amrex::Dim3 valid_hi = amrex::ubound(valid);
                amrex::Box bx = mfi.nodaltilebox();
                AMREX_D_TERM(bx.growLo(0, 2);, bx.growLo(1, 2);, bx.growLo(2, 2););
                AMREX_D_TERM(bx.growHi(0, 1);, bx.growHi(1, 1);, bx.growHi(2, 1););
                bx &= domain;

                amrex::Array4<const T>           const& model = a_model_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const& u = a_u_mf[lev]->array(mfi);
                amrex::Array4<Set::Matrix>       const& dw = a_dw_mf[lev]->array(mfi);
                amrex::Array4<Set::Matrix4<AMREX_SPACEDIM, T::sym>>  const& ddw = a_ddw_mf[lev]->array(mfi);

                // Set model internal dw and ddw.
                ALAMO_NEWTON_FOR(bx, [=] ALAMO_NEWTON_DEVICE (int i, int j, int k) {
                    auto sten = Numeric::GetStencil(i, j, k, domain);
                    const int ii = i < valid_lo.x ? valid_lo.x : (i > valid_hi.x ? valid_hi.x : i);
#if AMREX_SPACEDIM > 1
                    const int jj = j < valid_lo.y ? valid_lo.y : (j > valid_hi.y ? valid_hi.y : j);
#else
                    const int jj = j;
#endif
#if AMREX_SPACEDIM > 2
                    const int kk = k < valid_lo.z ? valid_lo.z : (k > valid_hi.z ? valid_hi.z : k);
#else
                    const int kk = k;
#endif
                    const T& material = model(ii, jj, kk);

                    Set::Matrix gradu = Numeric::PairedGradient(u, i, j, k, DX.data(), sten, periodic);

                    if (material.kinvar == Model::Solid::KinematicVariable::gradu)
                    {
                        dw(i, j, k) = material.DW(gradu);
                        ddw(i, j, k) = material.DDW(gradu);
                    }
                    else if (material.kinvar == Model::Solid::KinematicVariable::epsilon)
                    {
                        Set::Matrix eps = 0.5 * (gradu + gradu.transpose());
                        dw(i, j, k) = material.DW(eps);
                        ddw(i, j, k) = material.DDW(eps);
                    }
                    else if (material.kinvar == Model::Solid::KinematicVariable::F)
                    {
                        Set::Matrix F = gradu + Set::Matrix::Identity();
                        dw(i, j, k) = material.DW(F);
                        ddw(i, j, k) = material.DDW(F);
                    }
#ifdef ALAMO_GPU
                    else
                    {
                        Util::SetDeviceError(prepare_error_flag);
                    }
#endif
                });
            }
#ifdef ALAMO_GPU
            Util::AbortIfDeviceError(prepare_error, INFO,
                "Solver::Nonlocal::Newton::prepareForSolve() detected an invalid scalar model state on device");
#endif

            a_dw_mf[lev]->setMultiGhost(true);
            a_dw_mf[lev]->FillBoundaryAndSync(m_elastic->Geom(lev).periodicity());
            a_ddw_mf[lev]->setMultiGhost(true);
            a_ddw_mf[lev]->FillBoundaryAndSync(m_elastic->Geom(lev).periodicity());
        }

        m_elastic->SetModel(a_ddw_mf);

        for (int lev = 0; lev <= a_b_mf.finest_level; ++lev)
        {
            amrex::Box domain(linop->Geom(lev).growPeriodicDomain(2));
            domain.convert(amrex::IntVect::TheNodeVector());
            // Device-safe cell size (see note on the Set::Vector overload below).
            Set::Vector DX(linop->Geom(lev).CellSize());
            const std::array<bool, AMREX_SPACEDIM> periodic = {
                AMREX_D_DECL(linop->Geom(lev).isPeriodic(0),
                    linop->Geom(lev).isPeriodic(1),
                    linop->Geom(lev).isPeriodic(2))};
            const amrex::Dim3 lo = amrex::lbound(domain), hi = amrex::ubound(domain);
            for (MFIter mfi(*a_model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.nodaltilebox();
                bx = bx & domain;
                amrex::Array4<const Set::Scalar>  const& u = a_u_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar>  const& b = a_b_mf[lev]->array(mfi);
                amrex::Array4<const Set::Matrix>  const& dw = a_dw_mf[lev]->array(mfi);
                amrex::Array4<Set::Scalar>        const& rhs = a_rhs_mf[lev]->array(mfi);
#ifdef ALAMO_GPU
                auto m_bc_type = m_elastic->GetBC().GetBcTypeArray();
#endif
                ALAMO_NEWTON_FOR(bx, [=] ALAMO_NEWTON_DEVICE (int i, int j, int k)
                {
                    auto sten = Numeric::GetStencil(i, j, k, domain);
                    // Do this if on the domain boundary
                    if (AMREX_D_TERM(i == lo.x || i == hi.x, || j == lo.y || j == hi.y, || k == lo.z || k == hi.z))
                    {
                        Set::Matrix gradu = Numeric::PairedGradient(u, i, j, k, DX.data(), sten, periodic);
                        Set::Vector U(AMREX_D_DECL(u(i, j, k, 0), u(i, j, k, 1), u(i, j, k, 2)));
                        Set::Vector ret = ALAMO_NEWTON_BC_EVAL(m_elastic, m_bc_type, U, gradu, dw(i, j, k), i, j, k, domain);
                        for (int d = 0; d < AMREX_SPACEDIM; d++) rhs(i, j, k, d) = b(i, j, k, d) - ret(d);
                    }
                    else
                    {
                        Set::Vector divdw = Numeric::BackwardDivergence(dw, i, j, k, DX.data());
                        for (int p = 0; p < AMREX_SPACEDIM; p++) rhs(i, j, k, p) = b(i, j, k, p) - divdw(p);
                    }
                });
            }
            Util::RealFillBoundary(*a_ddw_mf[lev], m_elastic->Geom(lev));
            Util::RealFillBoundary(*a_rhs_mf[lev], m_elastic->Geom(lev));
        }
    }

    void prepareForSolve(const Set::Field<Set::Vector>& a_u_mf,
        const Set::Field<Set::Vector>& a_b_mf,
        Set::Field<Set::Scalar>& a_rhs_mf,
        Set::Field<Set::Matrix>& a_dw_mf,
        Set::Field<Set::Matrix4<AMREX_SPACEDIM, T::sym>>& a_ddw_mf,
        Set::Field<T>& a_model_mf)
    {
        for (int lev = 0; lev <= a_b_mf.finest_level; ++lev)
        {
            amrex::Box domain(linop->Geom(lev).growPeriodicDomain(2));
            domain.convert(amrex::IntVect::TheNodeVector());
            // NB: this kernel captures DX (a Set::Vector copied into the device
            // closure) and uses DX.data(); a raw Geom().CellSize() host pointer
            // would be invalid on device. The old `dx` local is therefore gone.
            Set::Vector DX(linop->Geom(lev).CellSize());
            const std::array<bool, AMREX_SPACEDIM> periodic = {
                AMREX_D_DECL(linop->Geom(lev).isPeriodic(0),
                    linop->Geom(lev).isPeriodic(1),
                    linop->Geom(lev).isPeriodic(2))};
#ifdef ALAMO_GPU
            Util::DeviceErrorFlag prepare_error;
            int* prepare_error_flag = prepare_error.dataPtr();
#endif
            const amrex::IntVect required_ghost(2);
            if (!a_u_mf[lev]->nGrowVect().allGE(required_ghost)
                || !a_model_mf[lev]->nGrowVect().allGE(required_ghost)
                || !a_dw_mf[lev]->nGrowVect().allGE(required_ghost)
                || !a_ddw_mf[lev]->nGrowVect().allGE(required_ghost))
                Util::Abort(INFO, "The paired elastic stencil requires two coefficient, model, and displacement ghosts.");

            for (MFIter mfi(*a_model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                // See the scalar overload: C/F model ghosts are not reliable.
                // Clamp the material choice, but evaluate DW/DDW with the
                // actual paired gradient at every flux anchor so the residual
                // and Fapply retain the same discrete derivative.
                const amrex::Box valid = mfi.validbox();
                const amrex::Dim3 valid_lo = amrex::lbound(valid);
                const amrex::Dim3 valid_hi = amrex::ubound(valid);
                amrex::Box bx = mfi.nodaltilebox();
                AMREX_D_TERM(bx.growLo(0, 2);, bx.growLo(1, 2);, bx.growLo(2, 2););
                AMREX_D_TERM(bx.growHi(0, 1);, bx.growHi(1, 1);, bx.growHi(2, 1););
                bx &= domain;

                amrex::Array4<const T>           const& model = a_model_mf[lev]->array(mfi);
                amrex::Array4<const Set::Vector> const& u = a_u_mf[lev]->array(mfi);
                amrex::Array4<Set::Matrix>       const& dw = a_dw_mf[lev]->array(mfi);
                amrex::Array4<Set::Matrix4<AMREX_SPACEDIM, T::sym>>  const& ddw = a_ddw_mf[lev]->array(mfi);

                // Set model internal dw and ddw.
                amrex::ParallelFor (bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    auto sten = Numeric::GetStencil(i, j, k, domain);
                    const int ii = i < valid_lo.x ? valid_lo.x : (i > valid_hi.x ? valid_hi.x : i);
#if AMREX_SPACEDIM > 1
                    const int jj = j < valid_lo.y ? valid_lo.y : (j > valid_hi.y ? valid_hi.y : j);
#else
                    const int jj = j;
#endif
#if AMREX_SPACEDIM > 2
                    const int kk = k < valid_lo.z ? valid_lo.z : (k > valid_hi.z ? valid_hi.z : k);
#else
                    const int kk = k;
#endif
                    const T& material = model(ii, jj, kk);

                    Set::Matrix gradu = Numeric::PairedGradient(u, i, j, k, DX.data(), sten, periodic);
                    Set::Matrix kinvar;
                    if (material.kinvar == Model::Solid::KinematicVariable::gradu)
                        kinvar = gradu; // gradu
                    else if (material.kinvar == Model::Solid::KinematicVariable::epsilon)
                        kinvar = 0.5 * (gradu + gradu.transpose()); // epsilon
                    else if (material.kinvar == Model::Solid::KinematicVariable::F)
                        kinvar = gradu + Set::Matrix::Identity(); // F
#ifdef ALAMO_GPU
                    else
                    {
                        kinvar = Set::Matrix::Zero();
                        Util::SetDeviceError(prepare_error_flag);
                    }
#endif

                    dw(i, j, k) = material.DW(kinvar);
                    ddw(i, j, k) = material.DDW(kinvar);

                });
            }
#ifdef ALAMO_GPU
            Util::AbortIfDeviceError(prepare_error, INFO,
                "Solver::Nonlocal::Newton::prepareForSolve() detected an invalid vector model state on device");
#endif

            a_dw_mf[lev]->setMultiGhost(true);
            a_dw_mf[lev]->FillBoundaryAndSync(m_elastic->Geom(lev).periodicity());
            a_ddw_mf[lev]->setMultiGhost(true);
            a_ddw_mf[lev]->FillBoundaryAndSync(m_elastic->Geom(lev).periodicity());
        }

        m_elastic->SetModel(a_ddw_mf);
        if (m_psi)
            for (int i = 0; i <= a_b_mf.finest_level; i++)
                m_elastic->SetPsi(i, *(*m_psi)[i]);


        for (int lev = 0; lev <= a_b_mf.finest_level; ++lev)
        {
            amrex::Box domain(linop->Geom(lev).growPeriodicDomain(2));
            domain.convert(amrex::IntVect::TheNodeVector());
            // Device-safe cell size: DX is a Set::Vector copied into the device
            // closure; DX.data() is valid on device. A raw Geom().CellSize()
            // host pointer would fault when dereferenced in these kernels.
            Set::Vector DX(linop->Geom(lev).CellSize());
            const amrex::Dim3 lo = amrex::lbound(domain), hi = amrex::ubound(domain);
            const amrex::Box psi_cell_domain = linop->Geom(lev).Domain();
            const amrex::Dim3 psi_cell_lo = amrex::lbound(psi_cell_domain);
            const amrex::Dim3 psi_cell_hi = amrex::ubound(psi_cell_domain);
            const std::array<bool, AMREX_SPACEDIM> psi_periodic = {
                AMREX_D_DECL(linop->Geom(lev).isPeriodic(0),
                    linop->Geom(lev).isPeriodic(1),
                    linop->Geom(lev).isPeriodic(2))};
            for (MFIter mfi(*a_model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.nodaltilebox() & domain;
                amrex::Array4<const Set::Vector>  const& u = a_u_mf[lev]->array(mfi);
                amrex::Array4<const Set::Vector>  const& b = a_b_mf[lev]->array(mfi);
                amrex::Array4<const Set::Matrix>  const& dw = a_dw_mf[lev]->array(mfi);
                amrex::Array4<Set::Scalar>        const& rhs = a_rhs_mf[lev]->array(mfi);
#ifdef ALAMO_GPU
                auto m_bc_type = m_elastic->GetBC().GetBcTypeArray();
#endif

                // The paired residual is D- [theta * DW(G+ u)].  The same
                // global-domain stencil is used when preparing DW/DDW and
                // here, so an MFIter tile edge never changes the equation.
                if (m_psi == nullptr)
                {
                    ALAMO_NEWTON_FOR (bx, [=] ALAMO_NEWTON_DEVICE (int i, int j, int k)
                    {
                        auto sten = Numeric::GetStencil(i, j, k, domain);
                        // Do this if on the domain boundary
                        if (AMREX_D_TERM(i == lo.x || i == hi.x, || j == lo.y || j == hi.y, || k == lo.z || k == hi.z))
                        {
                            Set::Matrix gradu = Numeric::PairedGradient(u, i, j, k, DX.data(), sten, psi_periodic);
                            Set::Vector ret = ALAMO_NEWTON_BC_EVAL(m_elastic, m_bc_type, u(i, j, k), gradu, dw(i, j, k), i, j, k, domain);
                            for (int d = 0; d < AMREX_SPACEDIM; d++) rhs(i, j, k, d) = b(i, j, k)(d) - ret(d);
                        }
                        else
                        {
                            Set::Vector divdw = Numeric::BackwardDivergence(dw, i, j, k, DX.data());
                            for (int d = 0; d < AMREX_SPACEDIM; d++) rhs(i, j, k, d) = b(i, j, k)(d) - divdw(d);
                        }
                    });
                }
                else
                {
                    amrex::Array4<const Set::Scalar>  const& psi = (*m_psi)[lev]->array(mfi);
                    const Set::Scalar psi_regularization = m_elastic->PsiRegularization();
                    ALAMO_NEWTON_FOR(bx, [=] ALAMO_NEWTON_DEVICE (int i, int j, int k)
                    {
                        auto sten = Numeric::GetStencil(i, j, k, domain);
                        const Set::Scalar theta = (1.0 - psi_regularization)
                            * Numeric::Interpolate::CellToNodeAverageBounded(psi, i, j, k, 0,
                                psi_cell_lo, psi_cell_hi, psi_periodic)
                            + psi_regularization;

                        // Do this if on the domain boundary
                        if (AMREX_D_TERM(i == lo.x || i == hi.x, || j == lo.y || j == hi.y, || k == lo.z || k == hi.z))
                        {
                            Set::Matrix gradu = Numeric::PairedGradient(u, i, j, k, DX.data(), sten, psi_periodic);
                            Set::Vector ret = ALAMO_NEWTON_BC_EVAL(m_elastic, m_bc_type, u(i, j, k), gradu, dw(i, j, k) * theta, i, j, k, domain);
                            for (int d = 0; d < AMREX_SPACEDIM; d++) rhs(i, j, k, d) = b(i, j, k)(d) - ret(d);
                        }
                        else
                        {
                            Set::Matrix flux = dw(i, j, k) * theta;
                            Set::Vector divdw = Set::Vector::Zero();

#if AMREX_SPACEDIM > 0
                            const Set::Scalar theta_xlo = (1.0 - psi_regularization)
                                * Numeric::Interpolate::CellToNodeAverageBounded(psi, i - 1, j, k, 0,
                                    psi_cell_lo, psi_cell_hi, psi_periodic)
                                + psi_regularization;
                            divdw += (flux.col(0) - (dw(i - 1, j, k) * theta_xlo).col(0)) / DX[0];
#endif
#if AMREX_SPACEDIM > 1
                            const Set::Scalar theta_ylo = (1.0 - psi_regularization)
                                * Numeric::Interpolate::CellToNodeAverageBounded(psi, i, j - 1, k, 0,
                                    psi_cell_lo, psi_cell_hi, psi_periodic)
                                + psi_regularization;
                            divdw += (flux.col(1) - (dw(i, j - 1, k) * theta_ylo).col(1)) / DX[1];
#endif
#if AMREX_SPACEDIM > 2
                            const Set::Scalar theta_zlo = (1.0 - psi_regularization)
                                * Numeric::Interpolate::CellToNodeAverageBounded(psi, i, j, k - 1, 0,
                                    psi_cell_lo, psi_cell_hi, psi_periodic)
                                + psi_regularization;
                            divdw += (flux.col(2) - (dw(i, j, k - 1) * theta_zlo).col(2)) / DX[2];
#endif

                            for (int d = 0; d < AMREX_SPACEDIM; d++) rhs(i, j, k, d) = b(i, j, k)(d) - divdw(d);
                        }
                    });
                }
            }

            a_rhs_mf[lev]->setMultiGhost(true);
            a_rhs_mf[lev]->FillBoundaryAndSync(m_elastic->Geom(lev).periodicity());
            a_ddw_mf[lev]->setMultiGhost(true);
            a_ddw_mf[lev]->FillBoundaryAndSync(m_elastic->Geom(lev).periodicity());
        }
    }


    Set::Scalar solve(const Set::Field<Set::Vector>& a_u_mf,
        const Set::Field<Set::Vector>& a_b_mf,
        Set::Field<T>& a_model_mf,
        Real a_tol_rel, Real a_tol_abs, const char* checkpoint_file = nullptr)
    {
        Set::Field<Set::Scalar> dsol_mf, rhs_mf, u0_mf, utrial_mf;
        Set::Field<Set::Matrix> dw_mf;
        Set::Field<Set::Matrix4<AMREX_SPACEDIM, T::sym>> ddw_mf;

        dsol_mf.resize(a_u_mf.finest_level + 1); dsol_mf.finest_level = a_u_mf.finest_level;
        dw_mf.resize(a_u_mf.finest_level + 1); dw_mf.finest_level = a_u_mf.finest_level;
        ddw_mf.resize(a_u_mf.finest_level + 1); ddw_mf.finest_level = a_u_mf.finest_level;
        rhs_mf.resize(a_u_mf.finest_level + 1); rhs_mf.finest_level = a_u_mf.finest_level;
        u0_mf.resize(a_u_mf.finest_level + 1); u0_mf.finest_level = a_u_mf.finest_level;
        utrial_mf.resize(a_u_mf.finest_level + 1); utrial_mf.finest_level = a_u_mf.finest_level;
        for (int lev = 0; lev <= a_u_mf.finest_level; lev++)
        {
            dsol_mf.Define(lev, a_u_mf[lev]->boxArray(),
                a_u_mf[lev]->DistributionMap(),
                a_u_mf.NComp(),
                a_u_mf[lev]->nGrow());
            dw_mf.Define(lev, a_b_mf[lev]->boxArray(),
                a_b_mf[lev]->DistributionMap(),
                1,
                a_b_mf[lev]->nGrow());
            ddw_mf.Define(lev, a_b_mf[lev]->boxArray(),
                a_b_mf[lev]->DistributionMap(),
                1,
                a_b_mf[lev]->nGrow());
            rhs_mf.Define(lev, a_b_mf[lev]->boxArray(),
                a_b_mf[lev]->DistributionMap(),
                a_b_mf.NComp(),
                a_b_mf[lev]->nGrow());
            u0_mf.Define(lev, a_u_mf[lev]->boxArray(),
                a_u_mf[lev]->DistributionMap(),
                a_u_mf.NComp(),
                a_u_mf[lev]->nGrow());
            utrial_mf.Define(lev, a_u_mf[lev]->boxArray(),
                a_u_mf[lev]->DistributionMap(),
                a_u_mf.NComp(),
                a_u_mf[lev]->nGrow());

            dsol_mf[lev]->setVal(0.0);
            dw_mf[lev]->setVal(Set::Matrix::Zero());
            ddw_mf[lev]->setVal(Set::Matrix4<AMREX_SPACEDIM, T::sym>::Zero());

            a_b_mf.Copy(lev, *rhs_mf[lev], 0, 2);
            //amrex::MultiFab::Copy(*rhs_mf[lev], *a_b_mf[lev], 0, 0, AMREX_SPACEDIM, 2);
        }

        Set::Scalar initial_resnorm = -1.0;
        Set::Scalar final_resnorm = -1.0;
        Set::Scalar final_resid_rel = std::numeric_limits<Set::Scalar>::infinity();
        for (int nriter = 0; nriter < m_nriters; nriter++)
        {
            if (verbose > 0 && nriter < m_nriters) Util::Message(INFO, "Newton Iteration ", nriter + 1, " of ", m_nriters);

            prepareForSolve(a_u_mf, a_b_mf, rhs_mf, dw_mf, ddw_mf, a_model_mf);

            // The persistent MLMG runs linop.prepareForSolve() only on its
            // first solve. SetModel (inside prepareForSolve above) updates
            // only the mglev=0 coefficients, so on relinearization the coarse
            // MG-level operators and the smoother/normalize diagonal must be
            // re-derived; otherwise the V-cycle smooths the new fine-level
            // operator against a stale hierarchy, which detonates at high
            // modulus contrast (low psi_floor).
            if (m_resync_coeffs && nriter > 0) m_elastic->SyncCoefficients();

            if (nriter == m_nriters) break;

            const bool need_residual_norm = m_line_search || m_nr_diagnostics
                || UsesResidualConvergence();
            const bool need_weighted_update = m_nr_diagnostics || UsesPsiUpdateConvergence();

            Set::Scalar resnorm0 = -1.0;
            if (need_residual_norm)
            {
                if (m_nr_diagnostics) LogFieldNorms(rhs_mf, "raw nonlinear rhs before linear solve");
            }

            // MLMG solves the composite AMR system, not independent equations
            // at coarse nodes covered by finer grids.  Capture the nonlinear
            // right-hand side with that same residual masking/reflux policy
            // after the solve has prepared its AMR/MG hierarchy.  The RHS is
            // passed to MLMG as const, so it is still the pre-update residual.
            Solver::Nonlocal::Linear::solve(dsol_mf, rhs_mf, a_tol_rel, a_tol_abs, checkpoint_file);

            if (need_residual_norm)
            {
                resnorm0 = CompositeResidualNorm(rhs_mf);
                if (nriter == 0) initial_resnorm = resnorm0;
                if (UsesResidualConvergence()
                    && (!std::isfinite(resnorm0) || resnorm0 < 0.0))
                    Util::Abort(INFO, "Newton residual convergence received a non-finite ",
                        "initial composite residual.");
            }

            Set::Scalar cornorm = 0, solnorm = 0;
            for (int lev = 0; lev < dsol_mf.size(); ++lev)
            {
                for (int comp = 0; comp < AMREX_SPACEDIM; comp++)
                {
                    Set::Scalar tmpcornorm = dsol_mf[lev]->norm0(comp, 0);
                    if (tmpcornorm > cornorm) cornorm = tmpcornorm;

                    //Set::Scalar tmpsolnorm = a_u_mf[lev]->norm0(comp,0);
                    //if (tmpsolnorm > solnorm) solnorm = tmpsolnorm;
                }

            }
            NRUpdateNorms update_norms;
            if (need_weighted_update) update_norms = ComputeWeightedUpdateNorms(dsol_mf);

            Set::Scalar alpha = 1.0;
            Set::Scalar resnorm = resnorm0;
            if (m_line_search)
            {
                for (int lev = 0; lev <= a_u_mf.finest_level; lev++)
                    a_u_mf.Copy(lev, *u0_mf[lev], 0, 2);

                constexpr int max_backtrack = 8;
                for (int bt = 0; bt <= max_backtrack; bt++)
                {
                    for (int lev = 0; lev <= a_u_mf.finest_level; lev++)
                    {
                        amrex::MultiFab::Copy(*utrial_mf[lev], *u0_mf[lev], 0, 0, AMREX_SPACEDIM, 2);
                        amrex::MultiFab::Saxpy(*utrial_mf[lev], alpha, *dsol_mf[lev], 0, 0, AMREX_SPACEDIM, 2);
                        a_u_mf.CopyFrom(lev, *utrial_mf[lev], 0, 2);
                    }

                    prepareForSolve(a_u_mf, a_b_mf, rhs_mf, dw_mf, ddw_mf, a_model_mf);

                    resnorm = CompositeResidualNorm(rhs_mf);
                    if (m_nr_diagnostics) LogFieldNorms(rhs_mf, "raw nonlinear rhs after line-search trial");

                    constexpr Set::Scalar bt_accept_tol = 1e-4;
                    if (std::isfinite(resnorm)
                        && resnorm <= resnorm0 * (1.0 + bt_accept_tol)) break;

                    if (bt == max_backtrack)
                    {
                        if (UsesResidualConvergence())
                        {
                            for (int lev = 0; lev <= a_u_mf.finest_level; lev++)
                                a_u_mf.Copy(lev, *u0_mf[lev], 0, 2);
                            prepareForSolve(a_u_mf, a_b_mf, rhs_mf, dw_mf, ddw_mf, a_model_mf);
                            Util::Abort(INFO, "Newton residual line search failed after ",
                                max_backtrack + 1, " trial steps: initial composite residual=",
                                resnorm0, ", final trial residual=", resnorm,
                                ". The displacement was restored and the update was not accepted.");
                        }
                        break;
                    }

                    alpha *= 0.5;
                }
            }
            else
            {
                for (int lev = 0; lev < dsol_mf.size(); ++lev)
                    a_u_mf.AddFrom(lev, *dsol_mf[lev], 0, 2);

                if (need_residual_norm)
                {
                    prepareForSolve(a_u_mf, a_b_mf, rhs_mf, dw_mf, ddw_mf, a_model_mf);
                    resnorm = CompositeResidualNorm(rhs_mf);
                    if (m_nr_diagnostics) LogFieldNorms(rhs_mf, "raw nonlinear rhs after update");
                }
            }

            // dsol_mf (and friends) are about to be either overwritten by the next
            // Linear::solve() call or destroyed if this is the final iteration. The
            // AddFrom/Saxpy above launched an async kernel that reads dsol_mf on
            // whichever stream its MFIter box landed on; without a full-device sync
            // here that read can still be in flight when the buffer is freed/reused
            // on another stream, causing an intermittent CUDA 700 illegal access
            // (same UAF family as the tmpfab.elixir() fix in Operator::interpolation).
            amrex::Gpu::streamSynchronizeAll();

            Set::Scalar relnorm;
            if (solnorm == 0) relnorm = cornorm;
            else relnorm = cornorm / solnorm;
            if (verbose > 0) Util::Message(INFO, "NR iteration ", nriter + 1, ", alpha = ", alpha, ", relative norm(ddisp) = ", relnorm);

            const Set::Scalar accepted_relnorm = alpha * relnorm;
            const Set::Scalar psi_relnorm = alpha * update_norms.psi_update;
            const Set::Scalar solid_relnorm = alpha * update_norms.solid_update;
            Set::Scalar resid_rel = -1.0;
            if (need_residual_norm)
            {
                if (initial_resnorm > 0.0 && std::isfinite(initial_resnorm)
                    && std::isfinite(resnorm))
                    resid_rel = resnorm / initial_resnorm;
                else if (initial_resnorm == 0.0 && resnorm == 0.0)
                    resid_rel = 0.0;
                else
                    resid_rel = std::numeric_limits<Set::Scalar>::infinity();
            }
            final_resnorm = resnorm;
            final_resid_rel = resid_rel;

            const bool converged = CheckNRConvergence(relnorm, psi_relnorm, resid_rel);
            if (m_nr_diagnostics || CustomNRConvergence())
                Util::Message(INFO, "NR convergence metrics: mode=", m_nr_convergence,
                    ", max_update=", cornorm,
                    ", accepted_max_update=", alpha * cornorm,
                    ", accepted_relative_update=", accepted_relnorm,
                    ", psi_update=", psi_relnorm,
                    ", solid_update=", solid_relnorm,
                    ", nonlinear_resid=", resnorm,
                    ", nonlinear_resid_rel=", resid_rel,
                    ", converged=", converged);

            if (converged)
                return CustomNRConvergence() ? SelectedNRMetric(relnorm, psi_relnorm, resid_rel) : relnorm;

        }

        if (UsesResidualConvergence())
            Util::Abort(INFO, "Newton residual convergence failed after ", m_nriters,
                " iterations (solver.nr_convergence=residual, solver.nrtolerance=",
                m_nrtolerance, "): initial composite residual=", initial_resnorm,
                ", final composite residual=", final_resnorm,
                ", relative residual=", final_resid_rel, ".");

        return 0.0;
    }

    Set::Scalar solve(const Set::Field<Set::Scalar>& a_u_mf,
        const Set::Field<Set::Scalar>& a_b_mf,
        Set::Field<T>& a_model_mf,
        Real a_tol_rel, Real a_tol_abs, const char* checkpoint_file = nullptr)
    {
        Set::Field<Set::Scalar> dsol_mf, rhs_mf;
        Set::Field<Set::Matrix> dw_mf;
        Set::Field<Set::Matrix4<AMREX_SPACEDIM, T::sym>> ddw_mf;

        dsol_mf.resize(a_u_mf.finest_level + 1); dsol_mf.finest_level = a_u_mf.finest_level;
        dw_mf.resize(a_u_mf.finest_level + 1); dw_mf.finest_level = a_u_mf.finest_level;
        ddw_mf.resize(a_u_mf.finest_level + 1); ddw_mf.finest_level = a_u_mf.finest_level;
        rhs_mf.resize(a_u_mf.finest_level + 1); rhs_mf.finest_level = a_u_mf.finest_level;
        for (int lev = 0; lev <= a_u_mf.finest_level; lev++)
        {
            dsol_mf.Define(lev, a_u_mf[lev]->boxArray(),
                a_u_mf[lev]->DistributionMap(),
                a_u_mf[lev]->nComp(),
                a_u_mf[lev]->nGrow());
            dw_mf.Define(lev, a_b_mf[lev]->boxArray(),
                a_b_mf[lev]->DistributionMap(),
                1,
                a_b_mf[lev]->nGrow());
            ddw_mf.Define(lev, a_b_mf[lev]->boxArray(),
                a_b_mf[lev]->DistributionMap(),
                1,
                a_b_mf[lev]->nGrow());
            rhs_mf.Define(lev, a_b_mf[lev]->boxArray(),
                a_b_mf[lev]->DistributionMap(),
                a_b_mf[lev]->nComp(),
                a_b_mf[lev]->nGrow());

            dsol_mf[lev]->setVal(0.0);
            dw_mf[lev]->setVal(Set::Matrix::Zero());
            ddw_mf[lev]->setVal(Set::Matrix4<AMREX_SPACEDIM, T::sym>::Zero());

            amrex::MultiFab::Copy(*rhs_mf[lev], *a_b_mf[lev], 0, 0, AMREX_SPACEDIM, 2);
        }

        for (int nriter = 0; nriter < m_nriters; nriter++)
        {
            if (verbose > 0 && nriter < m_nriters) Util::Message(INFO, "Newton Iteration ", nriter + 1, " of ", m_nriters);

            prepareForSolve(a_u_mf, a_b_mf, rhs_mf, dw_mf, ddw_mf, a_model_mf);

            // See the Set::Vector solve overload: re-derive coarse MG-level
            // coefficients and the diagonal after relinearization.
            if (m_resync_coeffs && nriter > 0) m_elastic->SyncCoefficients();

            if (nriter == m_nriters) break;

            Solver::Nonlocal::Linear::solve(dsol_mf, rhs_mf, a_tol_rel, a_tol_abs, checkpoint_file);

            Set::Scalar cornorm = 0, solnorm = 0;
            for (int lev = 0; lev < dsol_mf.size(); ++lev)
            {
                for (int comp = 0; comp < AMREX_SPACEDIM; comp++)
                {
                    Set::Scalar tmpcornorm = dsol_mf[lev]->norm0(comp, 0);
                    if (tmpcornorm > cornorm) cornorm = tmpcornorm;

                    Set::Scalar tmpsolnorm = a_u_mf[lev]->norm0(comp, 0);
                    if (tmpsolnorm > solnorm) solnorm = tmpsolnorm;
                }

            }
            Set::Scalar relnorm;
            if (solnorm == 0) relnorm = cornorm;
            else relnorm = cornorm / solnorm;
            if (verbose > 0) Util::Message(INFO, "NR iteration ", nriter + 1, ", relative norm(ddisp) = ", relnorm);

            for (int lev = 0; lev < dsol_mf.size(); ++lev)
                amrex::MultiFab::Add(*a_u_mf[lev], *dsol_mf[lev], 0, 0, AMREX_SPACEDIM, 2);

            if (relnorm < m_nrtolerance)
                return relnorm;

        }

        return 0.0;
    }
    Set::Scalar solve(const Set::Field<Set::Scalar>& a_u_mf,
        const Set::Field<Set::Scalar>& a_b_mf,
        Set::Field<T>& a_model_mf)
    {
        return solve(a_u_mf, a_b_mf, a_model_mf, tol_rel, tol_abs);
    }

    void compResidual(Set::Field<Set::Scalar>& a_res_mf,
        Set::Field<Set::Scalar>& a_u_mf,
        Set::Field<Set::Scalar>& a_b_mf,
        Set::Field<T>& a_model_mf)
    {
        Set::Field<Set::Matrix> dw_mf;
        Set::Field<Set::Matrix4<AMREX_SPACEDIM, T::sym>> ddw_mf;
        dw_mf.resize(a_u_mf.size());
        ddw_mf.resize(a_u_mf.size());
        for (int lev = 0; lev < a_u_mf.size(); lev++)
        {
            dw_mf.Define(lev, a_b_mf[lev]->boxArray(),
                a_b_mf[lev]->DistributionMap(),
                1, a_b_mf[lev]->nGrow());
            ddw_mf.Define(lev, a_b_mf[lev]->boxArray(),
                a_b_mf[lev]->DistributionMap(),
                1, a_b_mf[lev]->nGrow());
            dw_mf[lev]->setVal(Set::Matrix::Zero());
        }

        prepareForSolve(a_u_mf, a_b_mf, a_res_mf, dw_mf, ddw_mf, a_model_mf);
    }

    void compResidual(Set::Field<Set::Vector>& a_res_mf,
        Set::Field<Set::Vector>& a_u_mf,
        Set::Field<Set::Vector>& a_b_mf,
        Set::Field<T>& a_model_mf)
    {
        Set::Field<Set::Matrix> dw_mf;
        Set::Field<Set::Matrix4<AMREX_SPACEDIM, T::sym>> ddw_mf;
        Set::Field<Set::Scalar> res_mf(a_res_mf.size());
        dw_mf.resize(a_u_mf.size());
        ddw_mf.resize(a_u_mf.size());
        res_mf.resize(a_u_mf.size());
        for (int lev = 0; lev < a_u_mf.size(); lev++)
        {
            dw_mf.Define(lev, a_b_mf[lev]->boxArray(),
                a_b_mf[lev]->DistributionMap(),
                1, a_b_mf[lev]->nGrow());
            ddw_mf.Define(lev, a_b_mf[lev]->boxArray(),
                a_b_mf[lev]->DistributionMap(),
                1, a_b_mf[lev]->nGrow());
            res_mf.Define(lev, a_b_mf[lev]->boxArray(),
                a_b_mf[lev]->DistributionMap(),
                AMREX_SPACEDIM, a_b_mf[lev]->nGrow());
            dw_mf[lev]->setVal(Set::Matrix::Zero());
        }

        prepareForSolve(a_u_mf, a_b_mf, res_mf, dw_mf, ddw_mf, a_model_mf);

        for (int lev = 0; lev < a_res_mf.size(); ++lev)
        {
            Util::RealFillBoundary(*res_mf[lev], m_elastic->Geom(lev));
            a_res_mf.CopyFrom(lev, *res_mf[lev], 0, 2);
        }
    }

    void compLinearSolverResidual(Set::Field<Set::Vector>& a_res_mf,
        Set::Field<Set::Vector>& a_u_mf,
        Set::Field<Set::Vector>& a_b_mf)
    {
        Util::Assert(INFO, TEST(a_res_mf.finest_level == a_u_mf.finest_level));
        Util::Assert(INFO, TEST(a_u_mf.finest_level == a_b_mf.finest_level));
        Set::Field<Set::Scalar> res_mf(a_res_mf.finest_level + 1), sol_mf(a_u_mf.finest_level + 1), rhs_mf(a_b_mf.finest_level + 1);
        for (int lev = 0; lev <= a_u_mf.finest_level; lev++)
        {
            res_mf.Define(lev, a_b_mf[lev]->boxArray(), a_b_mf[lev]->DistributionMap(), AMREX_SPACEDIM, a_b_mf[lev]->nGrow());
            sol_mf.Define(lev, a_b_mf[lev]->boxArray(), a_b_mf[lev]->DistributionMap(), AMREX_SPACEDIM, a_b_mf[lev]->nGrow());
            rhs_mf.Define(lev, a_b_mf[lev]->boxArray(), a_b_mf[lev]->DistributionMap(), AMREX_SPACEDIM, a_b_mf[lev]->nGrow());

            a_u_mf.Copy(lev, *sol_mf[lev], 0, 2);
            a_b_mf.Copy(lev, *rhs_mf[lev], 0, 2);
        }

        mlmg->compResidual(amrex::GetVecOfPtrs(res_mf), amrex::GetVecOfPtrs(sol_mf), amrex::GetVecOfConstPtrs(rhs_mf));

        for (int lev = 0; lev <= a_res_mf.finest_level; ++lev)
        {
            a_res_mf.CopyFrom(lev, *res_mf[lev], 0, 2);
        }
    }


    void W(Set::Field<Set::Scalar>& a_w_mf,
        Set::Field<Set::Scalar>& a_u_mf,
        Set::Field<T>& a_model_mf)
    {
        for (int lev = 0; lev < a_u_mf.size(); lev++)
        {
            BL_PROFILE("Solver::Nonlocal::Newton::DW()");

            // Device-safe cell size: a Set::Vector copied into the device closure.
            Set::Vector DX(linop->Geom(lev).CellSize());
            amrex::Box domain(linop->Geom(lev).Domain());
            domain.convert(amrex::IntVect::TheNodeVector());

            for (MFIter mfi(*a_u_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                amrex::Array4<T> const& C = a_model_mf[lev]->array(mfi);
                amrex::Array4<amrex::Real> const& w = a_w_mf[lev]->array(mfi);
                amrex::Array4<const amrex::Real> const& u = a_u_mf[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Matrix gradu;

                    auto sten = Numeric::GetStencil(i, j, k, bx);

                    // Fill gradu
                    for (int p = 0; p < AMREX_SPACEDIM; p++)
                    {
                        AMREX_D_TERM(gradu(p, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(u, i, j, k, p, DX.data(), sten));,
                            gradu(p, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(u, i, j, k, p, DX.data(), sten));,
                            gradu(p, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(u, i, j, k, p, DX.data(), sten)););
                    }

                    if (C(i, j, k).kinvar == Model::Solid::KinematicVariable::gradu)
                        w(i, j, k) = C(i, j, k).W(gradu);
                    else if (C(i, j, k).kinvar == Model::Solid::KinematicVariable::epsilon)
                        w(i, j, k) = C(i, j, k).W(0.5 * (gradu + gradu.transpose()));
                    else if (C(i, j, k).kinvar == Model::Solid::KinematicVariable::F)
                        w(i, j, k) = C(i, j, k).W(gradu + Set::Matrix::Identity());
                });
            }
        }
    }

    void DW(Set::Field<Set::Scalar>& a_dw_mf,
        Set::Field<Set::Scalar>& a_u_mf,
        Set::Field<T>& a_model_mf)
    {
        for (int lev = 0; lev < a_u_mf.size(); lev++)
        {
            BL_PROFILE("Solver::Nonlocal::Newton::DW()");

            // Device-safe cell size: a Set::Vector copied into the device closure.
            Set::Vector DX(linop->Geom(lev).CellSize());
            amrex::Box domain(linop->Geom(lev).Domain());
            domain.convert(amrex::IntVect::TheNodeVector());

            for (MFIter mfi(*a_u_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.tilebox();
                amrex::Array4<T> const& C = a_model_mf[lev]->array(mfi);
                amrex::Array4<amrex::Real> const& dw = a_dw_mf[lev]->array(mfi);
                amrex::Array4<const amrex::Real> const& u = a_u_mf[lev]->array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
                {
                    Set::Matrix gradu;

                    auto sten = Numeric::GetStencil(i, j, k, bx);

                    // Fill gradu
                    for (int p = 0; p < AMREX_SPACEDIM; p++)
                    {
                        AMREX_D_TERM(gradu(p, 0) = (Numeric::Stencil<Set::Scalar, 1, 0, 0>::D(u, i, j, k, p, DX.data(), sten));,
                            gradu(p, 1) = (Numeric::Stencil<Set::Scalar, 0, 1, 0>::D(u, i, j, k, p, DX.data(), sten));,
                            gradu(p, 2) = (Numeric::Stencil<Set::Scalar, 0, 0, 1>::D(u, i, j, k, p, DX.data(), sten)););
                    }

                    Set::Matrix sig = Set::Matrix::Zero();

                    if (C(i, j, k).kinvar == Model::Solid::KinematicVariable::gradu)
                        sig = C(i, j, k).DW(gradu);
                    else if (C(i, j, k).kinvar == Model::Solid::KinematicVariable::epsilon)
                        sig = C(i, j, k).DW(0.5 * (gradu + gradu.transpose()));
                    else if (C(i, j, k).kinvar == Model::Solid::KinematicVariable::F)
                        sig = C(i, j, k).DW(gradu + Set::Matrix::Identity());

                    // = C(i,j,k)(gradu,m_homogeneous);

                    AMREX_D_PICK(dw(i, j, k, 0) = sig(0, 0);
                    ,
                        dw(i, j, k, 0) = sig(0, 0); dw(i, j, k, 1) = sig(0, 1);
                    dw(i, j, k, 2) = sig(1, 0); dw(i, j, k, 3) = sig(1, 1);
                    ,
                        dw(i, j, k, 0) = sig(0, 0); dw(i, j, k, 1) = sig(0, 1); dw(i, j, k, 2) = sig(0, 2);
                    dw(i, j, k, 3) = sig(1, 0); dw(i, j, k, 4) = sig(1, 1); dw(i, j, k, 5) = sig(1, 2);
                    dw(i, j, k, 6) = sig(2, 0); dw(i, j, k, 7) = sig(2, 1); dw(i, j, k, 8) = sig(2, 2););

                });
            }
        }
    }


    int m_nriters = 1;
    Set::Scalar m_nrtolerance = 0.0;
    bool m_line_search = false;
    bool m_resync_coeffs = false;
    std::string m_nr_convergence = "update";
    Set::Scalar m_nr_psi_weight_floor = 0.0;
    Set::Scalar m_nr_psi_power = 1.0;
    Set::Scalar m_nr_solid_psi = 0.5;
    bool m_nr_diagnostics = false;
    Operator::Elastic<T::sym>* m_elastic;
    //BC::Operator::Elastic::Elastic *m_bc;

    Set::Field<Set::Scalar>* m_psi = nullptr;

    struct NRUpdateNorms
    {
        Set::Scalar psi_update = 0.0;
        Set::Scalar solid_update = 0.0;
    };

    Set::Scalar FieldNorm0(const Set::Field<Set::Scalar>& a_mf) const
    {
        Set::Scalar norm = 0.0;
        for (int lev = 0; lev < a_mf.size(); ++lev)
        {
            const int ncomp = a_mf[lev]->nComp();
            for (int comp = 0; comp < ncomp; ++comp)
            {
                Set::Scalar tmp = a_mf[lev]->norm0(comp, 0);
                if (tmp > norm) norm = tmp;
            }
        }
        return norm;
    }

    Set::Scalar CompositeResidualNorm(const Set::Field<Set::Scalar>& a_rhs_mf)
    {
        Set::Field<Set::Scalar> zero_mf(a_rhs_mf.size());
        Set::Field<Set::Scalar> composite_mf(a_rhs_mf.size());
        zero_mf.finest_level = a_rhs_mf.finest_level;
        composite_mf.finest_level = a_rhs_mf.finest_level;

        for (int lev = 0; lev < a_rhs_mf.size(); ++lev)
        {
            zero_mf.Define(lev, a_rhs_mf[lev]->boxArray(),
                a_rhs_mf[lev]->DistributionMap(),
                a_rhs_mf[lev]->nComp(), a_rhs_mf[lev]->nGrow());
            composite_mf.Define(lev, a_rhs_mf[lev]->boxArray(),
                a_rhs_mf[lev]->DistributionMap(),
                a_rhs_mf[lev]->nComp(), a_rhs_mf[lev]->nGrow());
            zero_mf[lev]->setVal(0.0);
            composite_mf[lev]->setVal(0.0);
        }

        mlmg->compResidual(amrex::GetVecOfPtrs(composite_mf),
            amrex::GetVecOfPtrs(zero_mf),
            amrex::GetVecOfConstPtrs(a_rhs_mf));
        amrex::Gpu::streamSynchronizeAll();

        // MLMG's residual norm ignores coarse nodes covered by finer AMR
        // levels.  Use the operator's masked norm (rather than MultiFab's
        // raw norm0) and the same subcommunicator reduction as MLMG::solve.
        Set::Scalar norm = 0.0;
        for (int lev = 0; lev < composite_mf.size(); ++lev)
            norm = amrex::max(norm, linop->normInf(lev, *composite_mf[lev], true));
        amrex::ParallelAllReduce::Max(norm, amrex::ParallelContext::CommunicatorSub());
        return norm;
    }

    void LogFieldNorms(const Set::Field<Set::Scalar>& a_mf, const char* label) const
    {
        for (int lev = 0; lev < a_mf.size(); ++lev)
        {
            Set::Scalar level_norm = 0.0;
            for (int comp = 0; comp < a_mf[lev]->nComp(); ++comp)
                level_norm = amrex::max(level_norm, a_mf[lev]->norm0(comp, 0));
            Util::Message(INFO, "NR residual level norm: ", label,
                ", lev=", lev, ", norm0=", level_norm);
        }
    }

    NRUpdateNorms ComputeWeightedUpdateNorms(const Set::Field<Set::Scalar>& a_dsol_mf) const
    {
        NRUpdateNorms norms;
        if (!m_psi)
        {
            norms.psi_update = FieldNorm0(a_dsol_mf);
            norms.solid_update = norms.psi_update;
            return norms;
        }

        for (int lev = 0; lev < a_dsol_mf.size(); ++lev)
        {
            amrex::Box stencilbox(linop->Geom(lev).growPeriodicDomain(2));
            stencilbox.convert(amrex::IntVect::TheNodeVector());

            amrex::ReduceOps<amrex::ReduceOpMax, amrex::ReduceOpMax> reduce_op;
            amrex::ReduceData<amrex::Real, amrex::Real> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;

            const Set::Scalar psi_floor = m_nr_psi_weight_floor;
            const Set::Scalar psi_power = m_nr_psi_power;
            const Set::Scalar solid_psi = m_nr_solid_psi;

            for (MFIter mfi(*a_dsol_mf[lev], amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
            {
                const amrex::Box& bx = mfi.validbox();
                amrex::Array4<const Set::Scalar> const& dsol = a_dsol_mf[lev]->const_array(mfi);
                amrex::Array4<const Set::Scalar> const& psi = (*m_psi)[lev]->const_array(mfi);
                const int ncomp = a_dsol_mf[lev]->nComp();

                reduce_op.eval(bx, reduce_data,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                    {
                        const auto sten = Numeric::GetStencil(i, j, k, stencilbox);
                        Set::Scalar psi_avg = Numeric::Interpolate::CellToNodeAverage(psi, i, j, k, 0, sten);
                        psi_avg = amrex::max(Set::Scalar(0.0), psi_avg);

                        Set::Scalar weight = amrex::max(psi_floor, psi_avg);
                        if (psi_power != 1.0) weight = std::pow(weight, psi_power);

                        Set::Scalar weighted_max = 0.0;
                        Set::Scalar solid_max = 0.0;
                        for (int n = 0; n < ncomp; ++n)
                        {
                            const Set::Scalar absdu = amrex::Math::abs(dsol(i, j, k, n));
                            weighted_max = amrex::max(weighted_max, weight * absdu);
                            if (psi_avg >= solid_psi) solid_max = amrex::max(solid_max, absdu);
                        }
                        return {weighted_max, solid_max};
                    });
            }

            ReduceTuple hv = reduce_data.value(reduce_op);
            Set::Scalar weighted = amrex::get<0>(hv);
            Set::Scalar solid = amrex::get<1>(hv);
            amrex::ParallelDescriptor::ReduceRealMax(weighted);
            amrex::ParallelDescriptor::ReduceRealMax(solid);
            if (weighted > norms.psi_update) norms.psi_update = weighted;
            if (solid > norms.solid_update) norms.solid_update = solid;
        }

        return norms;
    }

    bool UsesPsiUpdateConvergence() const
    {
        return m_nr_convergence == "psi_update";
    }

    bool UsesResidualConvergence() const
    {
        return m_nr_convergence == "residual";
    }

    bool CustomNRConvergence() const
    {
        return m_nr_convergence != "update";
    }

    bool CheckNRConvergence(Set::Scalar relnorm, Set::Scalar psi_relnorm,
                            Set::Scalar resid_rel) const
    {
        if (m_nr_convergence == "update") return relnorm < m_nrtolerance;
        if (m_nr_convergence == "psi_update") return psi_relnorm < m_nrtolerance;
        if (m_nr_convergence == "residual")
            return std::isfinite(resid_rel) && resid_rel < m_nrtolerance;
        return false;
    }

    Set::Scalar SelectedNRMetric(Set::Scalar relnorm, Set::Scalar psi_relnorm,
                                 Set::Scalar resid_rel) const
    {
        if (m_nr_convergence == "psi_update") return psi_relnorm;
        if (m_nr_convergence == "residual") return resid_rel;
        return relnorm;
    }

    void ValidateNRConvergenceInputs() const
    {
        if (m_nr_convergence != "update" &&
            m_nr_convergence != "psi_update" &&
            m_nr_convergence != "residual")
            Util::Abort(INFO, "Invalid solver.nr_convergence=", m_nr_convergence,
                ". Valid values: update, psi_update, residual");
        if (!std::isfinite(m_nr_psi_weight_floor) || m_nr_psi_weight_floor < 0.0)
            Util::Abort(INFO, "solver.nr_psi_weight_floor must be finite and nonnegative");
        if (!std::isfinite(m_nr_psi_power) || m_nr_psi_power < 0.0)
            Util::Abort(INFO, "solver.nr_psi_power must be finite and nonnegative");
        if (!std::isfinite(m_nr_solid_psi) || m_nr_solid_psi < 0.0)
            Util::Abort(INFO, "solver.nr_solid_psi must be finite and nonnegative");
    }

    // These paramters control a standard Newton-Raphson solve.
    //
    // **Note**:
    // This class inherits all of the linear solve paramters
    // from its parent class, :ref:`Solver::Nonlocal::Linear`
    static void Parse(Newton<T>& value, amrex::ParmParse& pp)
    {
        Linear::Parse(value, pp);

        // Number of newton-raphson iterations.
        pp_query("nriters", value.m_nriters);

        // Tolerance to use for newton-raphson convergence
        pp_query("nrtolerance", value.m_nrtolerance);

        // Enable backtracking line-search damping of the Newton step.
        pp_query("line_search", value.m_line_search);

        // Re-derive coarse MG-level coefficients and the smoother diagonal
        // after each Newton relinearization (the persistent MLMG otherwise
        // reuses the hierarchy prepared on the first linear solve). Off by
        // default to preserve historical iteration trajectories exactly.
        pp_query("resync_coeffs", value.m_resync_coeffs);

        // Nonlinear convergence mode. Historical default "update" uses the max
        // displacement correction over all DOFs. "psi_update" weights the max
        // update by local psi so low-stiffness void displacement drift does not
        // dominate convergence after the stress-bearing region has settled.
        // "residual" uses the AMR composite nonlinear residual; unlike the
        // compatibility modes it aborts when the requested residual tolerance
        // is not achieved within nriters.
        pp_query("nr_convergence", value.m_nr_convergence);
        pp_query("nr_psi_weight_floor", value.m_nr_psi_weight_floor);
        pp_query("nr_psi_power", value.m_nr_psi_power);
        pp_query("nr_solid_psi", value.m_nr_solid_psi);
        pp_query("nr_diagnostics", value.m_nr_diagnostics);
        value.ValidateNRConvergenceInputs();
    }

};
} // namespace Nonlocal
} // namespace Solver


#endif
