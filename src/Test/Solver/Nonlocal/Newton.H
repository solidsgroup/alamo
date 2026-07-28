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
inline int
LineSearchResidualAcceptance()
{
    using ::Solver::Nonlocal::Detail::AcceptLineSearchResidual;

    const ::Set::Scalar nan = std::numeric_limits<::Set::Scalar>::quiet_NaN();
    const ::Set::Scalar inf = std::numeric_limits<::Set::Scalar>::infinity();
    int failed = 0;

    failed += !AcceptLineSearchResidual(10.0, 9.0);
    failed += !AcceptLineSearchResidual(10.0, 10.0);
    failed += !AcceptLineSearchResidual(1.0, 1.0001);
    failed += !AcceptLineSearchResidual(0.0, 0.0);
    failed += !AcceptLineSearchResidual(
        std::numeric_limits<::Set::Scalar>::max(),
        std::numeric_limits<::Set::Scalar>::max());

    failed += AcceptLineSearchResidual(10.0, 10.01);
    failed += AcceptLineSearchResidual(0.0, 1.0e-12);
    failed += AcceptLineSearchResidual(nan, 1.0);
    failed += AcceptLineSearchResidual(1.0, inf);
    failed += AcceptLineSearchResidual(1.0, 1.0, -1.0);

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
    failed += DecideNRTermination(1.0e-5, 1.0e-5, true) != NRTermination::Exhausted;
    failed += DecideNRTermination(nan, 1.0e-5, false) != NRTermination::InvalidMetric;

    failed += DecideNRTermination(5.0e-6, 1.0e-5, false) !=
        NRTermination::Converged;

    return failed;
}
} // namespace Nonlocal
} // namespace Solver
} // namespace Test

#endif
