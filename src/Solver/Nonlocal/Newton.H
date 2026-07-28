#ifndef SOLVER_NONLOCAL_NEWTON
#define SOLVER_NONLOCAL_NEWTON

#include <cmath>
#include <sstream>
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
namespace Detail
{
inline bool
AcceptLineSearchResidual(Set::Scalar initial, Set::Scalar trial,
                        Set::Scalar relative_tolerance = 1.0e-4)
{
    if (!std::isfinite(initial) || !std::isfinite(trial) ||
        !std::isfinite(relative_tolerance) || initial < 0.0 || trial < 0.0 ||
        relative_tolerance < 0.0)
        return false;
    if (initial == 0.0) return trial == 0.0;
    return trial / initial <= 1.0 + relative_tolerance;
}

enum class NRTermination
{
    Continue,
    Converged,
    FixedIterationsComplete,
    InvalidMetric,
    Exhausted
};

inline NRTermination
DecideNRTermination(Set::Scalar metric, Set::Scalar tolerance,
                    bool final_iteration)
{
    if (!std::isfinite(metric) || !std::isfinite(tolerance) || metric < 0.0 ||
        tolerance < 0.0)
        return NRTermination::InvalidMetric;
    if (tolerance == 0.0)
        return final_iteration
            ? NRTermination::FixedIterationsComplete
            : NRTermination::Continue;
    if (metric < tolerance) return NRTermination::Converged;
    return final_iteration ? NRTermination::Exhausted : NRTermination::Continue;
}

inline Set::Scalar
AcceptedUpdateMetric(Set::Scalar full_step_metric, Set::Scalar alpha)
{
    return alpha * full_step_metric;
}

template <typename T>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
T
FaceModel(const amrex::Array4<const T>& model,
        const int i, const int j, const int k, const int face)
{
    const int ip = i + (face == 0);
    const int jp = j + (face == 1);
    const int kp = k + (face == 2);
    return 0.5 * (model(i, j, k) + model(ip, jp, kp));
}

template <typename T, typename U>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Matrix
FaceKinematicVariable(const amrex::Array4<const T>& model,
                    const U& u,
                    const int i, const int j, const int k, const int face,
                    const Set::Scalar dx[AMREX_SPACEDIM])
{
    const T face_model = FaceModel(model, i, j, k, face);
    Set::Matrix gradu = Numeric::FaceGradient(u, i, j, k, face, dx);
    if (face_model.kinvar == Model::Solid::KinematicVariable::gradu)
        return gradu;
    if (face_model.kinvar == Model::Solid::KinematicVariable::epsilon)
        return 0.5 * (gradu + gradu.transpose());
    if (face_model.kinvar == Model::Solid::KinematicVariable::F)
        return gradu + Set::Matrix::Identity();
    return Set::Matrix::Zero();
}

template <typename T, typename U>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Matrix
FaceStress(const amrex::Array4<const T>& model,
            const U& u,
            const int i, const int j, const int k, const int face,
            const Set::Scalar dx[AMREX_SPACEDIM])
{
    const T face_model = FaceModel(model, i, j, k, face);
    return face_model.DW(
        FaceKinematicVariable(model, u, i, j, k, face, dx));
}

/// Reconstruct a nodal stress tensor from the conservative face tractions.
/// Each tensor column is the average of the matching low/high face-traction
/// columns, so its divergence uses the same directional fluxes as the solve.
template <typename T, typename U>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Matrix
NodalStressFromFaces(const amrex::Array4<const T>& model,
                    const U& u,
                    const int i, const int j, const int k,
                    const Set::Scalar dx[AMREX_SPACEDIM])
{
    Set::Matrix ret = Set::Matrix::Zero();
    for (int face = 0; face < AMREX_SPACEDIM; ++face)
    {
        const int im = i - (face == 0);
        const int jm = j - (face == 1);
        const int km = k - (face == 2);
        const Set::Matrix hi = FaceStress(model, u, i, j, k, face, dx);
        const Set::Matrix lo = FaceStress(model, u, im, jm, km, face, dx);
        ret.col(face) = 0.5 * (lo.col(face) + hi.col(face));
    }
    return ret;
}

/// Even extension of a material model through explicitly selected symmetry
/// planes. Material parameters do not change sign under reflection.
template <typename T>
struct SymmetryModelExtension
{
    amrex::Array4<const T> model;
    amrex::Dim3 lo, hi;
    amrex::GpuArray<int, AMREX_SPACEDIM> symmetry_lo;
    amrex::GpuArray<int, AMREX_SPACEDIM> symmetry_hi;

    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    T operator()(int i, int j, int k) const
    {
        AMREX_D_TERM(
            if (symmetry_lo[0] && i < lo.x) i = 2 * lo.x - i;
            if (symmetry_hi[0] && i > hi.x) i = 2 * hi.x - i;,
            if (symmetry_lo[1] && j < lo.y) j = 2 * lo.y - j;
            if (symmetry_hi[1] && j > hi.y) j = 2 * hi.y - j;,
            if (symmetry_lo[2] && k < lo.z) k = 2 * lo.z - k;
            if (symmetry_hi[2] && k > hi.z) k = 2 * hi.z - k;);
        return model(i, j, k);
    }
};

/// Symmetry extension of displacement: the component normal to each reflected
/// plane is odd, while tangential components are even.
struct SymmetryDisplacementExtension
{
    amrex::Array4<const Set::Vector> displacement;
    amrex::Dim3 lo, hi;
    amrex::GpuArray<int, AMREX_SPACEDIM> symmetry_lo;
    amrex::GpuArray<int, AMREX_SPACEDIM> symmetry_hi;

    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    Set::Vector operator()(int i, int j, int k) const
    {
        amrex::GpuArray<int, AMREX_SPACEDIM> reflected{};
        AMREX_D_TERM(
            if (symmetry_lo[0] && i < lo.x) { i = 2 * lo.x - i; reflected[0] = 1; }
            if (symmetry_hi[0] && i > hi.x) { i = 2 * hi.x - i; reflected[0] = 1; },
            if (symmetry_lo[1] && j < lo.y) { j = 2 * lo.y - j; reflected[1] = 1; }
            if (symmetry_hi[1] && j > hi.y) { j = 2 * hi.y - j; reflected[1] = 1; },
            if (symmetry_lo[2] && k < lo.z) { k = 2 * lo.z - k; reflected[2] = 1; }
            if (symmetry_hi[2] && k > hi.z) { k = 2 * hi.z - k; reflected[2] = 1; });
        Set::Vector ret = displacement(i, j, k);
        for (int d = 0; d < AMREX_SPACEDIM; ++d)
            if (reflected[d]) ret(d) = -ret(d);
        return ret;
    }
};

[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Matrix
SymmetryFaceGradient(const SymmetryDisplacementExtension& displacement,
                    const int i, const int j, const int k, const int face,
                    const Set::Scalar dx[AMREX_SPACEDIM])
{
    Set::Matrix ret = Set::Matrix::Zero();
    const int ip = i + (face == 0);
    const int jp = j + (face == 1);
    const int kp = k + (face == 2);
    ret.col(face) =
        (displacement(ip, jp, kp) - displacement(i, j, k)) / dx[face];
#if AMREX_SPACEDIM > 1
    if (face != 0)
        ret.col(0) = (displacement(i + 1, j, k) - displacement(i - 1, j, k)
                    + displacement(ip + 1, jp, kp) - displacement(ip - 1, jp, kp))
                    / (4.0 * dx[0]);
    if (face != 1)
        ret.col(1) = (displacement(i, j + 1, k) - displacement(i, j - 1, k)
                    + displacement(ip, jp + 1, kp) - displacement(ip, jp - 1, kp))
                    / (4.0 * dx[1]);
#endif
#if AMREX_SPACEDIM > 2
    if (face != 2)
        ret.col(2) = (displacement(i, j, k + 1) - displacement(i, j, k - 1)
                    + displacement(ip, jp, kp + 1) - displacement(ip, jp, kp - 1))
                    / (4.0 * dx[2]);
#endif
    return ret;
}

template <typename T>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Matrix
SymmetryFaceStress(const SymmetryModelExtension<T>& model,
                    const SymmetryDisplacementExtension& displacement,
                    const int i, const int j, const int k, const int face,
                    const Set::Scalar dx[AMREX_SPACEDIM])
{
    const int ip = i + (face == 0);
    const int jp = j + (face == 1);
    const int kp = k + (face == 2);
    const T face_model = 0.5 * (model(i, j, k) + model(ip, jp, kp));
    const Set::Matrix gradu =
        SymmetryFaceGradient(displacement, i, j, k, face, dx);
    if (face_model.kinvar == Model::Solid::KinematicVariable::gradu)
        return face_model.DW(gradu);
    if (face_model.kinvar == Model::Solid::KinematicVariable::epsilon)
        return face_model.DW(0.5 * (gradu + gradu.transpose()));
    if (face_model.kinvar == Model::Solid::KinematicVariable::F)
        return face_model.DW(gradu + Set::Matrix::Identity());
    return Set::Matrix::Zero();
}

/// Reconstruct stress on a physical symmetry plane by reflecting model and
/// displacement, then applying the identical centered face reconstruction
/// used in the domain interior.
template <typename T>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Matrix
NodalStressFromSymmetryFaces(
    const amrex::Array4<const T>& model,
    const amrex::Array4<const Set::Vector>& displacement,
    const int i, const int j, const int k,
    const Set::Scalar dx[AMREX_SPACEDIM],
    const amrex::Dim3 lo, const amrex::Dim3 hi,
    const amrex::GpuArray<int, AMREX_SPACEDIM>& symmetry_lo,
    const amrex::GpuArray<int, AMREX_SPACEDIM>& symmetry_hi)
{
    const SymmetryModelExtension<T> extended_model{
        model, lo, hi, symmetry_lo, symmetry_hi};
    const SymmetryDisplacementExtension extended_displacement{
        displacement, lo, hi, symmetry_lo, symmetry_hi};
    Set::Matrix ret = Set::Matrix::Zero();
    for (int face = 0; face < AMREX_SPACEDIM; ++face)
    {
        const int im = i - (face == 0);
        const int jm = j - (face == 1);
        const int km = k - (face == 2);
        const Set::Matrix hi_stress = SymmetryFaceStress(
            extended_model, extended_displacement, i, j, k, face, dx);
        const Set::Matrix lo_stress = SymmetryFaceStress(
            extended_model, extended_displacement, im, jm, km, face, dx);
        ret.col(face) = 0.5 * (lo_stress.col(face) + hi_stress.col(face));
    }
    return ret;
}

template <typename T, typename U>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Matrix4<AMREX_SPACEDIM, T::sym>
FaceTangent(const amrex::Array4<const T>& model,
            const U& u,
            const int i, const int j, const int k, const int face,
            const Set::Scalar dx[AMREX_SPACEDIM])
{
    const T face_model = FaceModel(model, i, j, k, face);
    return face_model.DDW(
        FaceKinematicVariable(model, u, i, j, k, face, dx));
}

template <typename T, typename U>
[[nodiscard]]
AMREX_GPU_HOST_DEVICE
AMREX_FORCE_INLINE
Set::Vector
FaceStressDivergence(const amrex::Array4<const T>& model,
                    const U& u,
                    const int i, const int j, const int k,
                    const Set::Scalar dx[AMREX_SPACEDIM])
{
    Set::Vector ret = Set::Vector::Zero();
    for (int face = 0; face < AMREX_SPACEDIM; ++face)
    {
        const int im = i - (face == 0);
        const int jm = j - (face == 1);
        const int km = k - (face == 2);
        const Set::Matrix hi = FaceStress(model, u, i, j, k, face, dx);
        const Set::Matrix lo = FaceStress(model, u, im, jm, km, face, dx);
        ret += (hi.col(face) - lo.col(face)) / dx[face];
    }
    return ret;
}
} // namespace Detail

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
        m_elastic->SetConservativeFaceFlux(m_conservative_face_flux);
    }
    void Clear()
    {
        Linear::Clear();
        m_elastic = nullptr;
    }


    void setNRIters(int a_nriters) { m_nriters = a_nriters; }

    void setPsi(Set::Field<Set::Scalar>& a_psi)
    {
        m_psi = &a_psi;
    }

    void setConservativeFaceFlux(bool a_enabled)
    {
        m_conservative_face_flux = a_enabled;
        if (m_elastic) m_elastic->SetConservativeFaceFlux(a_enabled);
    }

    [[nodiscard]] bool usesConservativeFaceFlux() const
    {
        return m_conservative_face_flux && m_psi == nullptr;
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
            amrex::Box domain(linop->Geom(lev).growPeriodicDomain(1));
            domain.convert(amrex::IntVect::TheNodeVector());
            // Device-safe cell size (see note on the Set::Vector overload below).
            Set::Vector DX(linop->Geom(lev).CellSize());
#ifdef ALAMO_GPU
            Util::DeviceErrorFlag prepare_error;
            int* prepare_error_flag = prepare_error.dataPtr();
#endif
            for (MFIter mfi(*a_model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = mfi.grownnodaltilebox();
                bx = bx & domain;

                amrex::Array4<const T>           const& model = a_model_mf[lev]->array(mfi);
                amrex::Array4<const Set::Scalar> const& u = a_u_mf[lev]->array(mfi);
                amrex::Array4<Set::Matrix>       const& dw = a_dw_mf[lev]->array(mfi);
                amrex::Array4<Set::Matrix4<AMREX_SPACEDIM, T::sym>>  const& ddw = a_ddw_mf[lev]->array(mfi);

                // Set model internal dw and ddw.
                ALAMO_NEWTON_FOR(bx, [=] ALAMO_NEWTON_DEVICE (int i, int j, int k) {
                    auto sten = Numeric::GetStencil(i, j, k, bx);

                    Set::Matrix gradu = Numeric::Gradient(u, i, j, k, DX.data(), sten);

                    if (model(i, j, k).kinvar == Model::Solid::KinematicVariable::gradu)
                    {
                        dw(i, j, k) = model(i, j, k).DW(gradu);
                        ddw(i, j, k) = model(i, j, k).DDW(gradu);
                    }
                    else if (model(i, j, k).kinvar == Model::Solid::KinematicVariable::epsilon)
                    {
                        Set::Matrix eps = 0.5 * (gradu + gradu.transpose());
                        dw(i, j, k) = model(i, j, k).DW(eps);
                        ddw(i, j, k) = model(i, j, k).DDW(eps);
                    }
                    else if (model(i, j, k).kinvar == Model::Solid::KinematicVariable::F)
                    {
                        Set::Matrix F = gradu + Set::Matrix::Identity();
                        dw(i, j, k) = model(i, j, k).DW(F);
                        ddw(i, j, k) = model(i, j, k).DDW(F);
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

            Util::RealFillBoundary(*a_dw_mf[lev], m_elastic->Geom(lev));
            Util::RealFillBoundary(*a_ddw_mf[lev], m_elastic->Geom(lev));
        }

        m_elastic->SetModel(a_ddw_mf);

        for (int lev = 0; lev <= a_b_mf.finest_level; ++lev)
        {
            amrex::Box domain(linop->Geom(lev).growPeriodicDomain(1));
            domain.convert(amrex::IntVect::TheNodeVector());
            // Device-safe cell size (see note on the Set::Vector overload below).
            Set::Vector DX(linop->Geom(lev).CellSize());
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
                    auto sten = Numeric::GetStencil(i, j, k, bx);
                    // Do this if on the domain boundary
                    if (AMREX_D_TERM(i == lo.x || i == hi.x, || j == lo.y || j == hi.y, || k == lo.z || k == hi.z))
                    {
                        Set::Matrix gradu = Numeric::Gradient(u, i, j, k, DX.data(), sten);
                        Set::Vector U(AMREX_D_DECL(u(i, j, k, 0), u(i, j, k, 1), u(i, j, k, 2)));
                        Set::Vector ret = ALAMO_NEWTON_BC_EVAL(m_elastic, m_bc_type, U, gradu, dw(i, j, k), i, j, k, bx);
                        for (int d = 0; d < AMREX_SPACEDIM; d++) rhs(i, j, k, d) = b(i, j, k, d) - ret(d);
                    }
                    else
                    {
                        Set::Vector divdw = Numeric::Divergence(dw, i, j, k, DX.data());
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
        const bool conservative_solve = m_conservative_face_flux && m_psi == nullptr;
        for (int lev = 0; lev <= a_b_mf.finest_level; ++lev)
        {
            amrex::Box domain(linop->Geom(lev).growPeriodicDomain(2));
            domain.convert(amrex::IntVect::TheNodeVector());
            amrex::Box row_domain(linop->Geom(lev).growPeriodicDomain(1));
            row_domain.convert(amrex::IntVect::TheNodeVector());
            Set::Vector DX(linop->Geom(lev).CellSize()); // device-safe: captured by value, DX.data() valid on device
            const amrex::Dim3 hi = amrex::ubound(domain);
#ifdef ALAMO_GPU
            Util::DeviceErrorFlag prepare_error;
            int* prepare_error_flag = prepare_error.dataPtr();
#endif
            for (MFIter mfi(*a_model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = conservative_solve
                    ? (mfi.grownnodaltilebox(-1, 1) & row_domain)
                    : (mfi.grownnodaltilebox() & domain);

                amrex::Array4<const T>           const& model = a_model_mf[lev]->array(mfi);
                amrex::Array4<const Set::Vector> const& u = a_u_mf[lev]->array(mfi);
                amrex::Array4<Set::Matrix>       const& dw = a_dw_mf[lev]->array(mfi);
                amrex::Array4<Set::Matrix4<AMREX_SPACEDIM, T::sym>>  const& ddw = a_ddw_mf[lev]->array(mfi);

                // Set model internal dw and ddw.
                amrex::ParallelFor (bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    auto sten = Numeric::GetStencil(i, j, k, bx);

                    Set::Matrix gradu = Numeric::Gradient(u, i, j, k, DX.data(), sten);
                    Set::Matrix kinvar;
                    if (model(i, j, k).kinvar == Model::Solid::KinematicVariable::gradu)
                        kinvar = gradu; // gradu
                    else if (model(i, j, k).kinvar == Model::Solid::KinematicVariable::epsilon)
                        kinvar = 0.5 * (gradu + gradu.transpose()); // epsilon
                    else if (model(i, j, k).kinvar == Model::Solid::KinematicVariable::F)
                        kinvar = gradu + Set::Matrix::Identity(); // F
#ifdef ALAMO_GPU
                    else
                    {
                        kinvar = Set::Matrix::Zero();
                        Util::SetDeviceError(prepare_error_flag);
                    }
#endif

                    dw(i, j, k) = model(i, j, k).DW(kinvar);
                    ddw(i, j, k, 0) = model(i, j, k).DDW(kinvar);

                    if (!conservative_solve)
                    {
                        // Extra components unused by legacy masked operator; init deterministically without reading past support ghost.
                        for (int face = 0; face < AMREX_SPACEDIM; ++face)
                            ddw(i, j, k, face + 1) = ddw(i, j, k, 0);
                    }
                    else
                    {
                        for (int face = 0; face < AMREX_SPACEDIM; ++face)
                        {
                            const int ip = i + (face == 0);
                            const int jp = j + (face == 1);
                            [[maybe_unused]] const int kp = k + (face == 2);
                            if (AMREX_D_TERM(ip <= hi.x, && jp <= hi.y, && kp <= hi.z))
                                ddw(i, j, k, face + 1) = Detail::FaceTangent(
                                    model, u, i, j, k, face, DX.data());
                            else
                                ddw(i, j, k, face + 1) = ddw(i, j, k, 0);
                        }
                    }

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
            amrex::Box row_domain(linop->Geom(lev).growPeriodicDomain(1));
            row_domain.convert(amrex::IntVect::TheNodeVector());
            Set::Vector DX(linop->Geom(lev).CellSize()); // device-safe cell size
            const amrex::Dim3 lo = amrex::lbound(domain), hi = amrex::ubound(domain);
            for (MFIter mfi(*a_model_mf[lev], false); mfi.isValid(); ++mfi)
            {
                amrex::Box bx = conservative_solve
                    ? (mfi.grownnodaltilebox(-1, 1) & row_domain)
                    : (mfi.grownnodaltilebox() & domain);
                amrex::Array4<const T>            const& model = a_model_mf[lev]->array(mfi);
                amrex::Array4<const Set::Vector>  const& u = a_u_mf[lev]->array(mfi);
                amrex::Array4<const Set::Vector>  const& b = a_b_mf[lev]->array(mfi);
                amrex::Array4<const Set::Matrix>  const& dw = a_dw_mf[lev]->array(mfi);
                amrex::Array4<Set::Scalar>        const& rhs = a_rhs_mf[lev]->array(mfi);
#ifdef ALAMO_GPU
                auto m_bc_type = m_elastic->GetBC().GetBcTypeArray();
#endif

                if (conservative_solve)
                {
                    ALAMO_NEWTON_FOR (bx, [=] ALAMO_NEWTON_DEVICE (int i, int j, int k)
                    {
                        auto sten = Numeric::GetStencil(i, j, k, bx);
                        // Do this if on the domain boundary
                        if (AMREX_D_TERM(i == lo.x || i == hi.x, || j == lo.y || j == hi.y, || k == lo.z || k == hi.z))
                        {
                            Set::Matrix gradu = Numeric::Gradient(u, i, j, k, DX.data(), sten);
                            Set::Vector ret = ALAMO_NEWTON_BC_EVAL(m_elastic, m_bc_type, u(i, j, k), gradu, dw(i, j, k), i, j, k, bx);
                            for (int d = 0; d < AMREX_SPACEDIM; d++) rhs(i, j, k, d) = b(i, j, k)(d) - ret(d);
                        }
                        else
                        {
                            Set::Vector divdw = Detail::FaceStressDivergence(
                                model, u, i, j, k, DX.data());
                            for (int d = 0; d < AMREX_SPACEDIM; d++) rhs(i, j, k, d) = b(i, j, k)(d) - divdw(d);
                        }
                    });
                }
                else
                {
                    const amrex::Dim3 boxlo = amrex::lbound(bx), boxhi = amrex::ubound(bx);
                    amrex::Array4<const Set::Scalar> const psi = m_psi
                        ? (*m_psi)[lev]->array(mfi)
                        : amrex::Array4<const Set::Scalar>();
                    const bool use_psi = m_psi != nullptr;
                    ALAMO_NEWTON_FOR(bx, [=] ALAMO_NEWTON_DEVICE (int i, int j, int k)
                    {
                        auto sten = Numeric::GetStencil(i, j, k, bx);
                        // Do this if on the domain boundary
                        if (AMREX_D_TERM(i == lo.x || i == hi.x, || j == lo.y || j == hi.y, || k == lo.z || k == hi.z))
                        {
                            Set::Matrix gradu = Numeric::Gradient(u, i, j, k, DX.data(), sten);
                            Set::Scalar psiavg = 1.0;
                            if (use_psi && AMREX_D_TERM(i > boxlo.x && i<boxhi.x, && j>boxlo.y && j<boxhi.y, && k>boxlo.z && k < boxhi.z))
                            {
                                psiavg = Numeric::Interpolate::CellToNodeAverage(psi, i, j, k, 0);
                            }
                            Set::Vector ret = ALAMO_NEWTON_BC_EVAL(m_elastic, m_bc_type, u(i, j, k), gradu, dw(i, j, k) * psiavg, i, j, k, bx);
                            for (int d = 0; d < AMREX_SPACEDIM; d++) rhs(i, j, k, d) = b(i, j, k)(d) - ret(d);
                        }
                        else
                        {
                            Set::Vector divdw = Numeric::Divergence(dw, i, j, k, DX.data(), sten);
                            if (use_psi && AMREX_D_TERM(i > boxlo.x && i<boxhi.x, && j>boxlo.y && j<boxhi.y, && k>boxlo.z && k < boxhi.z))
                            {
                                divdw *= Numeric::Interpolate::CellToNodeAverage(psi, i, j, k, 0);
                                Set::Vector gradpsi = Numeric::CellGradientOnNode(psi, i, j, k, 0, DX.data());
                                divdw += dw(i, j, k) * gradpsi;
                            }
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
                AMREX_SPACEDIM + 1,
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
        }

        Set::Scalar initial_resnorm = -1.0;
        Set::Scalar last_nr_metric = -1.0;
        for (int nriter = 0; nriter < m_nriters; nriter++)
        {
            if (verbose > 0 && nriter < m_nriters) Util::Message(INFO, "Newton Iteration ", nriter + 1, " of ", m_nriters);

            prepareForSolve(a_u_mf, a_b_mf, rhs_mf, dw_mf, ddw_mf, a_model_mf);
            AbortIfNonFinite(rhs_mf, "Newton residual");

            // Re-derive coarse MG coefficients on relinearization; a stale hierarchy can detonate at high modulus contrast.
            if (m_resync_coeffs && nriter > 0) m_elastic->SyncCoefficients();

            if (nriter == m_nriters) break;

            const bool need_residual_norm = m_line_search || m_nr_diagnostics;
            const bool need_weighted_update = m_nr_diagnostics || UsesPsiUpdateConvergence();

            Set::Scalar resnorm0 = -1.0;
            if (need_residual_norm)
            {
                resnorm0 = FieldNorm0(rhs_mf);
                if (nriter == 0) initial_resnorm = resnorm0;
            }

            // Zero the correction each iteration; reusing the prior one skews MLMG's resid/resid0 denominator.
            for (int lev = 0; lev < dsol_mf.size(); ++lev)
                dsol_mf[lev]->setVal(0.0);
            Solver::Nonlocal::Linear::solve(dsol_mf, rhs_mf, a_tol_rel, a_tol_abs, checkpoint_file);
            ProjectInvariantPeriodicDirections(dsol_mf);
            AbortIfNonFinite(dsol_mf, "Newton correction", true);

            Set::Scalar cornorm = 0, solnorm = 0;
            for (int lev = 0; lev < dsol_mf.size(); ++lev)
            {
                for (int comp = 0; comp < AMREX_SPACEDIM; comp++)
                {
                    Set::Scalar tmpcornorm = dsol_mf[lev]->norm0(comp, 0);
                    if (tmpcornorm > cornorm) cornorm = tmpcornorm;
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
                AbortIfNonFinite(u0_mf, "Newton line-search baseline", true);

                auto restore_baseline = [&]()
                {
                    for (int lev = 0; lev <= a_u_mf.finest_level; lev++)
                        a_u_mf.CopyFrom(lev, *u0_mf[lev], 0, 2);
                    amrex::Gpu::streamSynchronizeAll();
                };
                auto restore_and_abort = [&](const std::string& reason)
                {
                    restore_baseline();
                    prepareForSolve(a_u_mf, a_b_mf, rhs_mf, dw_mf, ddw_mf, a_model_mf);
                    AbortIfNonFinite(rhs_mf, "restored Newton residual");
                    Util::Abort(INFO, reason);
                };
                auto restore_baseline_state = [&]()
                {
                    restore_baseline();
                    prepareForSolve(a_u_mf, a_b_mf, rhs_mf, dw_mf, ddw_mf, a_model_mf);
                    AbortIfNonFinite(rhs_mf, "restored Newton residual");
                };

                constexpr int max_backtrack = 8;
                bool accepted = false;
                for (int bt = 0; bt <= max_backtrack; bt++)
                {
                    try
                    {
                        for (int lev = 0; lev <= a_u_mf.finest_level; lev++)
                        {
                            amrex::MultiFab::Copy(*utrial_mf[lev], *u0_mf[lev],
                                0, 0, AMREX_SPACEDIM, 2);
                            amrex::MultiFab::Saxpy(*utrial_mf[lev], alpha,
                                *dsol_mf[lev], 0, 0, AMREX_SPACEDIM, 2);
                        }
                        if (HasNonFinite(utrial_mf, true))
                            restore_and_abort("Newton line-search trial contains a non-finite value");
                        for (int lev = 0; lev <= a_u_mf.finest_level; lev++)
                            a_u_mf.CopyFrom(lev, *utrial_mf[lev], 0, 2);

                        prepareForSolve(a_u_mf, a_b_mf, rhs_mf, dw_mf, ddw_mf, a_model_mf);
                    }
                    catch (...)
                    {
                        restore_baseline_state();
                        throw;
                    }
                    if (HasNonFinite(rhs_mf))
                        restore_and_abort("Newton line-search residual contains a non-finite value");

                    resnorm = 0;
                    resnorm = FieldNorm0(rhs_mf);

                    if (Detail::AcceptLineSearchResidual(resnorm0, resnorm))
                    {
                        accepted = true;
                        break;
                    }

                    alpha *= 0.5;
                }
                if (!accepted)
                {
                    std::ostringstream reason;
                    reason << "Newton line search failed after " << max_backtrack
                            << " backtracks: initial residual=" << resnorm0
                            << ", final trial residual=" << resnorm;
                    restore_and_abort(reason.str());
                }
            }
            else
            {
                for (int lev = 0; lev < dsol_mf.size(); ++lev)
                    a_u_mf.AddFrom(lev, *dsol_mf[lev], 0, 2);
                for (int lev = 0; lev <= a_u_mf.finest_level; ++lev)
                {
                    a_u_mf.Copy(lev, *u0_mf[lev], 0, 2);
                }
                AbortIfNonFinite(u0_mf, "accepted Newton solution", true);

                if (need_residual_norm)
                {
                    prepareForSolve(a_u_mf, a_b_mf, rhs_mf, dw_mf, ddw_mf, a_model_mf);
                    AbortIfNonFinite(rhs_mf, "accepted Newton residual");
                    resnorm = FieldNorm0(rhs_mf);
                }
            }

            // Sync before dsol_mf is overwritten/destroyed: an async kernel above may still be reading it (UAF).
            amrex::Gpu::streamSynchronizeAll();

            Set::Scalar relnorm;
            if (solnorm == 0) relnorm = cornorm;
            else relnorm = cornorm / solnorm;
            if (!std::isfinite(cornorm) || !std::isfinite(relnorm))
                Util::Abort(INFO, "Newton update norm is non-finite");
            if (verbose > 0) Util::Message(INFO, "NR iteration ", nriter + 1,
                ", alpha = ", alpha, ", max norm(ddisp) = ", relnorm);

            const Set::Scalar accepted_relnorm =
                Detail::AcceptedUpdateMetric(relnorm, alpha);
            const Set::Scalar psi_relnorm =
                Detail::AcceptedUpdateMetric(update_norms.psi_update, alpha);
            const Set::Scalar solid_relnorm =
                Detail::AcceptedUpdateMetric(update_norms.solid_update, alpha);
            if (!std::isfinite(psi_relnorm) || !std::isfinite(solid_relnorm))
                Util::Abort(INFO, "Newton weighted update norm is non-finite");
            Set::Scalar resid_rel = -1.0;
            if (need_residual_norm)
            {
                if (initial_resnorm > 0.0) resid_rel = resnorm / initial_resnorm;
                else resid_rel = resnorm;
            }

            last_nr_metric = SelectedNRMetric(accepted_relnorm, psi_relnorm);
            const Detail::NRTermination termination = Detail::DecideNRTermination(
                last_nr_metric, m_nrtolerance, nriter + 1 == m_nriters);
            if (termination == Detail::NRTermination::InvalidMetric)
                Util::Abort(INFO, "Newton convergence metric is invalid");
            const bool converged = termination == Detail::NRTermination::Converged;
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
                return CustomNRConvergence()
                    ? SelectedNRMetric(accepted_relnorm, psi_relnorm)
                    : accepted_relnorm;
            if (termination == Detail::NRTermination::Exhausted)
                Util::Abort(INFO, "Newton solver failed to converge in ", m_nriters,
                    " iterations: mode=", m_nr_convergence,
                    ", final metric=", last_nr_metric,
                    ", tolerance=", m_nrtolerance);

        }

        if (m_nrtolerance > 0.0)
            Util::Abort(INFO, "Newton solver failed to converge in ", m_nriters,
                " iterations: mode=", m_nr_convergence,
                ", final metric=", last_nr_metric,
                ", tolerance=", m_nrtolerance);
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

        Set::Scalar last_relnorm = -1.0;
        for (int nriter = 0; nriter < m_nriters; nriter++)
        {
            if (verbose > 0 && nriter < m_nriters) Util::Message(INFO, "Newton Iteration ", nriter + 1, " of ", m_nriters);

            prepareForSolve(a_u_mf, a_b_mf, rhs_mf, dw_mf, ddw_mf, a_model_mf);
            AbortIfNonFinite(rhs_mf, "Newton residual");

            // Re-derive coarse MG coefficients/diagonal after relinearization.
            if (m_resync_coeffs && nriter > 0) m_elastic->SyncCoefficients();

            if (nriter == m_nriters) break;

            for (int lev = 0; lev < dsol_mf.size(); ++lev)
                dsol_mf[lev]->setVal(0.0);
            Solver::Nonlocal::Linear::solve(dsol_mf, rhs_mf, a_tol_rel, a_tol_abs, checkpoint_file);
            ProjectInvariantPeriodicDirections(dsol_mf);
            AbortIfNonFinite(dsol_mf, "Newton correction", true);

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
            if (!std::isfinite(cornorm) || !std::isfinite(relnorm))
                Util::Abort(INFO, "Newton update norm is non-finite");
            last_relnorm = relnorm;
            const Detail::NRTermination termination = Detail::DecideNRTermination(
                relnorm, m_nrtolerance, nriter + 1 == m_nriters);
            if (termination == Detail::NRTermination::InvalidMetric)
                Util::Abort(INFO, "Newton convergence metric is invalid");
            if (verbose > 0) Util::Message(INFO, "NR iteration ", nriter + 1, ", relative norm(ddisp) = ", relnorm);

            for (int lev = 0; lev < dsol_mf.size(); ++lev)
                amrex::MultiFab::Add(*a_u_mf[lev], *dsol_mf[lev], 0, 0, AMREX_SPACEDIM, 2);
            amrex::Gpu::streamSynchronizeAll();
            AbortIfNonFinite(a_u_mf, "accepted Newton solution", true);

            if (termination == Detail::NRTermination::Converged)
                return relnorm;
            if (termination == Detail::NRTermination::Exhausted)
                Util::Abort(INFO, "Newton solver failed to converge in ", m_nriters,
                    " iterations: final relative update=", last_relnorm,
                    ", tolerance=", m_nrtolerance);

        }

        if (m_nrtolerance > 0.0)
            Util::Abort(INFO, "Newton solver failed to converge in ", m_nriters,
                " iterations: final relative update=", last_relnorm,
                ", tolerance=", m_nrtolerance);
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
        const int nlev = a_u_mf.finest_level + 1;
        dw_mf.resize(nlev);
        ddw_mf.resize(nlev);
        dw_mf.finest_level = a_u_mf.finest_level;
        ddw_mf.finest_level = a_u_mf.finest_level;
        for (int lev = 0; lev < nlev; lev++)
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
        const int nlev = a_u_mf.finest_level + 1;
        Set::Field<Set::Scalar> res_mf(nlev);
        dw_mf.resize(nlev);
        ddw_mf.resize(nlev);
        res_mf.resize(nlev);
        dw_mf.finest_level = a_u_mf.finest_level;
        ddw_mf.finest_level = a_u_mf.finest_level;
        res_mf.finest_level = a_u_mf.finest_level;
        for (int lev = 0; lev < nlev; lev++)
        {
            dw_mf.Define(lev, a_b_mf[lev]->boxArray(),
                a_b_mf[lev]->DistributionMap(),
                1, a_b_mf[lev]->nGrow());
            ddw_mf.Define(lev, a_b_mf[lev]->boxArray(),
                a_b_mf[lev]->DistributionMap(),
                AMREX_SPACEDIM + 1, a_b_mf[lev]->nGrow());
            res_mf.Define(lev, a_b_mf[lev]->boxArray(),
                a_b_mf[lev]->DistributionMap(),
                AMREX_SPACEDIM, a_b_mf[lev]->nGrow());
            dw_mf[lev]->setVal(Set::Matrix::Zero());
        }

        prepareForSolve(a_u_mf, a_b_mf, res_mf, dw_mf, ddw_mf, a_model_mf);

        // Reflux to match the composite residual MLMG solves and norms.
        for (int lev = a_u_mf.finest_level - 1; lev >= 0; --lev)
            m_elastic->Reflux(lev, *res_mf[lev], *res_mf[lev], *res_mf[lev],
                *res_mf[lev + 1], *res_mf[lev + 1], *res_mf[lev + 1]);

        for (int lev = 0; lev < nlev; ++lev)
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
    std::vector<int> m_invariant_periodic;
    Operator::Elastic<T::sym>* m_elastic = nullptr;

    Set::Field<Set::Scalar>* m_psi = nullptr;
    bool m_conservative_face_flux = false;

    struct NRUpdateNorms
    {
        Set::Scalar psi_update = 0.0;
        Set::Scalar solid_update = 0.0;
    };

    void ProjectInvariantPeriodicDirections(Set::Field<Set::Scalar>& field) const
    {
        if (m_invariant_periodic.empty()) return;

        for (int dir = 0; dir < AMREX_SPACEDIM; ++dir)
        {
            if (!m_invariant_periodic[dir]) continue;
            for (int lev = 0; lev < field.size(); ++lev)
            {
                const amrex::Geometry& geom = m_elastic->Geom(lev);
                if (!geom.isPeriodic(dir))
                    Util::Abort(INFO, "solver.invariant_periodic direction ", dir,
                        " is not periodic");

                amrex::Box domain = geom.Domain();
                domain.convert(amrex::IntVect::TheNodeVector());
                const int source_index = domain.smallEnd(dir);
                for (amrex::MFIter mfi(*field[lev], false); mfi.isValid(); ++mfi)
                {
                    // Include all ghosts or machine-scale noise can re-enter the residual via face gradients.
                    const amrex::Box bx = mfi.grownnodaltilebox();
                    const amrex::Box valid = mfi.validbox();
                    if (valid.smallEnd(dir) != domain.smallEnd(dir) ||
                        valid.bigEnd(dir) != domain.bigEnd(dir))
                        Util::Abort(INFO,
                            "solver.invariant_periodic requires every grid to span ",
                            "the periodic direction; direction=", dir,
                            ", level=", lev, ", grid=", valid,
                            ", domain=", domain);

                    amrex::Array4<amrex::Real> const& fab = field[lev]->array(mfi);
                    const int ncomp = field[lev]->nComp();
                    amrex::ParallelFor(bx, ncomp,
                        [=] AMREX_GPU_DEVICE(int i, int j, int k, int n)
                        {
                            int is = i, js = j, ks = k;
                            if (dir == 0) is = source_index;
#if AMREX_SPACEDIM > 1
                            if (dir == 1) js = source_index;
#endif
#if AMREX_SPACEDIM > 2
                            if (dir == 2) ks = source_index;
#endif
                            fab(i, j, k, n) = fab(is, js, ks, n);
                        });
                }
                field[lev]->setMultiGhost(true);
                field[lev]->FillBoundaryAndSync(geom.periodicity());
            }
        }
        amrex::Gpu::streamSynchronizeAll();
    }

    bool HasNonFinite(const Set::Field<Set::Scalar>& a_mf,
                    bool include_ghosts = false) const
    {
        amrex::Gpu::streamSynchronizeAll();
        if (!include_ghosts)
        {
            for (int lev = 0; lev < a_mf.size(); ++lev)
                if (!a_mf[lev]->is_finite(0, a_mf[lev]->nComp(), 0)) return true;
            return false;
        }

        bool nonfinite = false;
        for (int lev = 0; lev < a_mf.size(); ++lev)
        {
            const int ncomp = a_mf[lev]->nComp();
            amrex::Box active_domain = m_elastic->Geom(lev).growPeriodicDomain(
                a_mf[lev]->nGrowVect().max());
            active_domain.convert(amrex::IntVect::TheNodeVector());
            for (amrex::MFIter mfi(*a_mf[lev], false); mfi.isValid(); ++mfi)
            {
                const amrex::Box bx = mfi.growntilebox(a_mf[lev]->nGrowVect()) &
                    active_domain;
                if (!bx.ok()) continue;

                amrex::ReduceOps<amrex::ReduceOpLogicalOr> reduce_op;
                amrex::ReduceData<int> reduce_data(reduce_op);
                using ReduceTuple = typename decltype(reduce_data)::Type;
                const auto field = a_mf[lev]->const_array(mfi);
                reduce_op.eval(bx, reduce_data,
                    [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                    {
                        bool bad = false;
                        for (int n = 0; n < ncomp; ++n)
                            bad = bad || !amrex::Math::isfinite(field(i, j, k, n));
                        return {static_cast<int>(bad)};
                    });
                nonfinite = nonfinite || (amrex::get<0>(reduce_data.value(reduce_op)) != 0);
            }
        }
        amrex::ParallelAllReduce::Or(nonfinite,
            amrex::ParallelContext::CommunicatorSub());
        return nonfinite;
    }

    void AbortIfNonFinite(const Set::Field<Set::Scalar>& a_mf,
                        const char* description,
                        bool include_ghosts = false) const
    {
        if (HasNonFinite(a_mf, include_ghosts))
            Util::Abort(INFO, description, " contains a non-finite value");
    }

    Set::Scalar FieldNorm0(const Set::Field<Set::Scalar>& a_mf) const
    {
        amrex::Vector<std::unique_ptr<amrex::MultiFab>> composite(a_mf.size());
        for (int lev = 0; lev < a_mf.size(); ++lev)
        {
            composite[lev] = std::make_unique<amrex::MultiFab>(
                a_mf[lev]->boxArray(), a_mf[lev]->DistributionMap(),
                a_mf[lev]->nComp(), a_mf[lev]->nGrowVect());
            amrex::MultiFab::Copy(*composite[lev], *a_mf[lev], 0, 0,
                a_mf[lev]->nComp(), a_mf[lev]->nGrowVect());
        }

        // Reflux to match MLMG's C/F row replacement, or line search can reject a good Newton step.
        for (int lev = static_cast<int>(composite.size()) - 2; lev >= 0; --lev)
            m_elastic->Reflux(lev, *composite[lev], *composite[lev],
                *composite[lev], *composite[lev + 1], *composite[lev + 1],
                *composite[lev + 1]);

        Set::Scalar norm = 0.0;
        for (int lev = 0; lev < composite.size(); ++lev)
        {
            const Set::Scalar tmp = m_elastic->CompositeNormInf(
                lev, *composite[lev]);
            if (tmp > norm) norm = tmp;
        }
        return norm;
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

    bool CustomNRConvergence() const
    {
        return m_nr_convergence != "update";
    }

    Set::Scalar SelectedNRMetric(Set::Scalar relnorm, Set::Scalar psi_relnorm) const
    {
        if (m_nr_convergence == "psi_update") return psi_relnorm;
        return relnorm;
    }

    void ValidateNRConvergenceInputs() const
    {
        if (!std::isfinite(m_nrtolerance) || m_nrtolerance < 0.0)
            Util::Abort(INFO, "solver.nrtolerance must be finite and nonnegative");
        if (m_nr_convergence != "update" &&
            m_nr_convergence != "psi_update")
            Util::Abort(INFO, "Invalid solver.nr_convergence=", m_nr_convergence,
                ". Valid values: update, psi_update");
        if (!std::isfinite(m_nr_psi_weight_floor) || m_nr_psi_weight_floor < 0.0)
            Util::Abort(INFO, "solver.nr_psi_weight_floor must be finite and nonnegative");
        if (!std::isfinite(m_nr_psi_power) || m_nr_psi_power < 0.0)
            Util::Abort(INFO, "solver.nr_psi_power must be finite and nonnegative");
        if (!std::isfinite(m_nr_solid_psi) || m_nr_solid_psi < 0.0)
            Util::Abort(INFO, "solver.nr_solid_psi must be finite and nonnegative");
    }

    // Also inherits linear solve parameters from Solver::Nonlocal::Linear.
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
        pp_query("nr_convergence", value.m_nr_convergence);
        // Minimum psi used when weighting displacement updates.
        pp_query("nr_psi_weight_floor", value.m_nr_psi_weight_floor);
        // Exponent applied to psi in the weighted Newton convergence metric.
        pp_query("nr_psi_power", value.m_nr_psi_power);
        // Psi threshold defining the stress-bearing (solid) region.
        pp_query("nr_solid_psi", value.m_nr_solid_psi);
        // Enable diagnostic output for nonlinear convergence metrics.
        pp_query("nr_diagnostics", value.m_nr_diagnostics);
        // Per-coordinate flags enabling periodic invariant-preserving solves.
        pp_queryarr("invariant_periodic", value.m_invariant_periodic);
        if (!value.m_invariant_periodic.empty() &&
            value.m_invariant_periodic.size() != AMREX_SPACEDIM)
            Util::Abort(INFO, "solver.invariant_periodic must contain exactly ",
                AMREX_SPACEDIM, " entries");
        value.ValidateNRConvergenceInputs();
    }

};
} // namespace Nonlocal
} // namespace Solver


#endif
