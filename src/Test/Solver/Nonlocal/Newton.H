#ifndef TEST_SOLVER_NONLOCAL_NEWTON_H
#define TEST_SOLVER_NONLOCAL_NEWTON_H

#include <limits>

#include "Solver/Nonlocal/Newton.H"

namespace Test
{
namespace Solver
{
namespace Nonlocal
{
struct FaceOutputTestModel
{
    static constexpr auto kinvar = ::Model::Solid::KinematicVariable::gradu;

    ::Set::Scalar scale = 0.0;

    AMREX_GPU_HOST_DEVICE FaceOutputTestModel operator+(const FaceOutputTestModel& rhs) const
    {
        FaceOutputTestModel ret;
        ret.scale = scale + rhs.scale;
        return ret;
    }

    AMREX_GPU_HOST_DEVICE friend FaceOutputTestModel operator*(
        ::Set::Scalar lhs, const FaceOutputTestModel& rhs)
    {
        FaceOutputTestModel ret;
        ret.scale = lhs * rhs.scale;
        return ret;
    }

    AMREX_GPU_HOST_DEVICE ::Set::Matrix DW(const ::Set::Matrix& gradu) const
    {
        return scale * gradu;
    }
};

inline int
LineSearchResidualAcceptance()
{
    using ::Solver::Nonlocal::Detail::AcceptLineSearchResidual;

    const ::Set::Scalar nan = std::numeric_limits<::Set::Scalar>::quiet_NaN();
    const ::Set::Scalar inf = std::numeric_limits<::Set::Scalar>::infinity();
    int failed = 0;

    failed += !AcceptLineSearchResidual(10.0, 9.0);
    failed += !AcceptLineSearchResidual(10.0, 10.0);
    failed += !AcceptLineSearchResidual(10.0, 10.0005);
    failed += !AcceptLineSearchResidual(1.0, 1.0001);
    failed += !AcceptLineSearchResidual(0.0, 0.0);
    failed += !AcceptLineSearchResidual(
        std::numeric_limits<::Set::Scalar>::max(),
        std::numeric_limits<::Set::Scalar>::max());

    failed += AcceptLineSearchResidual(10.0, 10.01);
    failed += AcceptLineSearchResidual(0.0, 1.0e-12);
    failed += AcceptLineSearchResidual(-1.0, 0.0);
    failed += AcceptLineSearchResidual(1.0, -1.0);
    failed += AcceptLineSearchResidual(nan, 1.0);
    failed += AcceptLineSearchResidual(1.0, nan);
    failed += AcceptLineSearchResidual(inf, 1.0);
    failed += AcceptLineSearchResidual(1.0, inf);
    failed += AcceptLineSearchResidual(1.0, 1.0, -1.0);
    failed += AcceptLineSearchResidual(1.0, 1.0, nan);
    failed += AcceptLineSearchResidual(1.0, 1.0, inf);

    return failed;
}

inline int
NewtonTerminationDecision()
{
    using ::Solver::Nonlocal::Detail::DecideNRTermination;
    using ::Solver::Nonlocal::Detail::NRTermination;

    const ::Set::Scalar nan = std::numeric_limits<::Set::Scalar>::quiet_NaN();
    int failed = 0;

    failed += DecideNRTermination(1.0, 0.0, false) != NRTermination::Continue;
    failed += DecideNRTermination(1.0, 0.0, true) != NRTermination::FixedIterationsComplete;
    failed += DecideNRTermination(1.0e-6, 1.0e-5, false) != NRTermination::Converged;
    failed += DecideNRTermination(1.0e-5, 1.0e-5, false) != NRTermination::Continue;
    failed += DecideNRTermination(1.0e-5, 1.0e-5, true) != NRTermination::Exhausted;
    failed += DecideNRTermination(-1.0, 1.0e-5, false) != NRTermination::InvalidMetric;
    failed += DecideNRTermination(nan, 1.0e-5, false) != NRTermination::InvalidMetric;
    failed += DecideNRTermination(1.0, nan, false) != NRTermination::InvalidMetric;
    failed += DecideNRTermination(1.0, -1.0, false) != NRTermination::InvalidMetric;

    return failed;
}

inline int
ConservativeNodalStressReconstruction()
{
    const amrex::IntVect small(AMREX_D_DECL(-3, -3, -3));
    const amrex::IntVect big(AMREX_D_DECL(5, 5, 5));
    const amrex::Box box(small, big);
    amrex::BaseFab<FaceOutputTestModel> model(box, 1);
    amrex::FArrayBox displacement(box, AMREX_SPACEDIM);
    auto model_array = model.array();
    auto u = displacement.array();
    constexpr ::Set::Scalar scale = 2.5;
    amrex::LoopOnCpu(box, [=] (int i, int j, int k)
    {
        model_array(i, j, k).scale = scale;
        u(i, j, k, 0) = i * i + 2.0 * j;
        u(i, j, k, 1) = 3.0 * i - j * j;
#if AMREX_SPACEDIM == 3
        u(i, j, k, 2) = i + 2.0 * j + k * k;
#endif
    });

    const ::Set::Scalar dx[AMREX_SPACEDIM] = {AMREX_D_DECL(1.0, 1.0, 1.0)};
    const int i = 2, j = 2, k = AMREX_D_PICK(0, 0, 2);
    const ::Set::Matrix actual =
        ::Solver::Nonlocal::Detail::NodalStressFromFaces(
            model.const_array(), displacement.const_array(), i, j, k, dx);
    ::Set::Matrix expected = ::Set::Matrix::Zero();
    expected(0, 0) = 2.0 * i;
    expected(0, 1) = 2.0;
    expected(1, 0) = 3.0;
    expected(1, 1) = -2.0 * j;
#if AMREX_SPACEDIM == 3
    expected(2, 0) = 1.0;
    expected(2, 1) = 2.0;
    expected(2, 2) = 2.0 * k;
#endif
    expected *= scale;
    return (actual - expected).lpNorm<Eigen::Infinity>() < 1.0e-12 ? 0 : 1;
}

inline int
SymmetryNodalStressReconstruction()
{
    const amrex::IntVect small(AMREX_D_DECL(0, 0, 0));
    const amrex::IntVect big(AMREX_D_DECL(4, 4, 4));
    const amrex::Box box(small, big);
    amrex::BaseFab<FaceOutputTestModel> model(box, 1);
    amrex::BaseFab<::Set::Vector> displacement(box, 1);
    auto model_array = model.array();
    auto u = displacement.array();
    constexpr ::Set::Scalar scale = 2.5;
    amrex::LoopOnCpu(box, [=] (int i, int j, int k)
    {
        model_array(i, j, k).scale = scale;
        u(i, j, k) = ::Set::Vector::Zero();
        // About y=0: u_x is even and u_y is odd.
        u(i, j, k)(0) = i * i + j * j;
        u(i, j, k)(1) = (3.0 * i + 1.0) * j;
    });

    const ::Set::Scalar dx[AMREX_SPACEDIM] = {AMREX_D_DECL(1.0, 1.0, 1.0)};
    const amrex::Dim3 lo = amrex::lbound(box), hi = amrex::ubound(box);
    amrex::GpuArray<int, AMREX_SPACEDIM> symmetry_lo{};
    amrex::GpuArray<int, AMREX_SPACEDIM> symmetry_hi{};
    symmetry_lo[1] = 1;
    const int k = AMREX_D_PICK(0, 0, 2);
    ::Set::Matrix actual = ::Solver::Nonlocal::Detail::NodalStressFromSymmetryFaces(
        model.const_array(), displacement.const_array(), 2, 0, k, dx,
        lo, hi, symmetry_lo, symmetry_hi);
    ::Set::Matrix expected = ::Set::Matrix::Zero();
    expected(0, 0) = 4.0;
    expected(1, 1) = 7.0;
    expected *= scale;
    int failed = (actual - expected).lpNorm<Eigen::Infinity>() < 1.0e-12 ? 0 : 1;

    amrex::LoopOnCpu(box, [=] (int i, int j, int k2)
    {
        const ::Set::Scalar q = i - hi.x;
        u(i, j, k2) = ::Set::Vector::Zero();
        // About x=hi.x: u_x is odd and u_y is even.
        u(i, j, k2)(0) = (2.0 * j + 1.0) * q;
        u(i, j, k2)(1) = q * q + j * j;
    });
    symmetry_lo[1] = 0;
    symmetry_hi[0] = 1;
    actual = ::Solver::Nonlocal::Detail::NodalStressFromSymmetryFaces(
        model.const_array(), displacement.const_array(), hi.x, 2, k, dx,
        lo, hi, symmetry_lo, symmetry_hi);
    expected = ::Set::Matrix::Zero();
    expected(0, 0) = 5.0;
    expected(1, 1) = 4.0;
    expected *= scale;
    failed += (actual - expected).lpNorm<Eigen::Infinity>() < 1.0e-12 ? 0 : 1;
    return failed;
}
} // namespace Nonlocal
} // namespace Solver
} // namespace Test

#endif
