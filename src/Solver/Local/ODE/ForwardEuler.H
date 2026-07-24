#ifndef SOLVER_LOCAL_ODE_FORWARDEULER_H
#define SOLVER_LOCAL_ODE_FORWARDEULER_H

#include <cmath>

#include "Solver/Local/ODE/ODE.H"

namespace Solver
{
namespace Local
{
namespace ODE
{

struct ForwardEuler
{
    static constexpr const char* name = "forward_euler";
    static constexpr bool implicit = false;
    static constexpr bool requires_jacobian = false;
    static constexpr bool supports_jacobian = false;

    int nsubsteps = 1;

    template<int MAX_SIZE, typename RHS, typename Jacobian>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Solver::Local::ODE::Advance::Result
    Advance(
        Set::Scalar* state,
        int size,
        Set::Scalar dt,
        const RHS& rhs,
        const Jacobian& /*jacobian*/,
        bool /*has_jacobian*/) const
    {
        Solver::Local::ODE::Advance::Result result;
        if (size <= 0 || size > MAX_SIZE || !(dt > 0.0)) return result;

        Set::Scalar rate[MAX_SIZE]{};
        Set::Scalar trial[MAX_SIZE]{};
        if (!rhs(state, rate)) return result;

        for (int n = 0; n < size; ++n)
        {
            trial[n] = state[n] + dt * rate[n];
            if (!std::isfinite(trial[n])) return result;
        }

        // A second RHS evaluation validates model-specific state constraints.
        if (!rhs(trial, rate)) return result;
        for (int n = 0; n < size; ++n) state[n] = trial[n];

        result.converged = true;
        result.residual_norm = 0.0;
        return result;
    }

    static void Parse(ForwardEuler& value, IO::ParmParse& pp)
    {
        Util::Message(INFO, pp.getPrefix());
        pp.query_default("nsubsteps", value.nsubsteps, 1);
        Util::Assert(INFO, TEST(value.nsubsteps > 0),
            "Forward Euler nsubsteps must be positive.");
    }
};

} // namespace ODE
} // namespace Local
} // namespace Solver

#endif
