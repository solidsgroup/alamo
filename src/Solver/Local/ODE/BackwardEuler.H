#ifndef SOLVER_LOCAL_ODE_BACKWARDEULER_H
#define SOLVER_LOCAL_ODE_BACKWARDEULER_H

#include <cmath>
#include <limits>

#include "Solver/Local/ODE/ODE.H"

namespace Solver
{
namespace Local
{
namespace ODE
{

struct BackwardEuler
{
    static constexpr const char* name = "backward_euler";
    static constexpr bool implicit = true;
    static constexpr bool requires_jacobian = false;
    static constexpr bool supports_jacobian = true;

    int nsubsteps = 1;
    int max_iterations = 100;
    Set::Scalar relative_tolerance = 1.0e-8;
    Set::Scalar absolute_tolerance = 1.0e-12;

    void Configure(
        int a_max_iterations,
        Set::Scalar a_relative_tolerance,
        Set::Scalar a_absolute_tolerance)
    {
        max_iterations = a_max_iterations;
        relative_tolerance = a_relative_tolerance;
        absolute_tolerance = a_absolute_tolerance;
    }

    template<int MAX_SIZE, typename RHS, typename Jacobian>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Solver::Local::ODE::Advance::Result
    Advance(
        Set::Scalar* state,
        int size,
        Set::Scalar dt,
        const RHS& rhs,
        const Jacobian& rhs_jacobian,
        bool has_jacobian) const
    {
        Solver::Local::ODE::Advance::Result result;
        if (size <= 0 || size > MAX_SIZE || !(dt > 0.0) ||
            max_iterations <= 0) return result;

        Set::Scalar old_state[MAX_SIZE]{};
        Set::Scalar residual[MAX_SIZE]{};
        for (int n = 0; n < size; ++n) old_state[n] = state[n];

        auto evaluate_residual = [=] AMREX_GPU_HOST_DEVICE(
            const Set::Scalar* trial, Set::Scalar* values)
        {
            Set::Scalar rate[MAX_SIZE]{};
            if (!rhs(trial, rate)) return false;
            for (int n = 0; n < size; ++n)
            {
                values[n] = trial[n] - old_state[n] - dt * rate[n];
                if (!std::isfinite(values[n])) return false;
            }
            return true;
        };

        if (!evaluate_residual(state, residual)) return result;

        auto norm = [=] AMREX_GPU_HOST_DEVICE(const Set::Scalar* values)
        {
            Set::Scalar value = 0.0;
            for (int n = 0; n < size; ++n)
                value = value > std::abs(values[n]) ? value : std::abs(values[n]);
            return value;
        };

        result.residual_norm = norm(residual);
        const Set::Scalar tolerance = absolute_tolerance +
            relative_tolerance * result.residual_norm;
        result.converged = result.residual_norm <= tolerance;

        for (int iter = 0; iter < max_iterations && !result.converged; ++iter)
        {
            Set::Scalar jacobian[MAX_SIZE][MAX_SIZE]{};
            Set::Scalar update[MAX_SIZE]{};
            bool failed = false;

            if (has_jacobian)
            {
                Set::Scalar derivative[MAX_SIZE * MAX_SIZE]{};
                if (!rhs_jacobian(state, derivative)) return result;
                for (int row = 0; row < size; ++row)
                    for (int column = 0; column < size; ++column)
                    {
                        const Set::Scalar value = derivative[row * size + column];
                        if (!std::isfinite(value)) return result;
                        jacobian[row][column] =
                            (row == column ? 1.0 : 0.0) - dt * value;
                    }
            }
            else
            {
                for (int column = 0; column < size; ++column)
                {
                    Set::Scalar perturbed_state[MAX_SIZE]{};
                    Set::Scalar perturbed_residual[MAX_SIZE]{};
                    for (int n = 0; n < size; ++n)
                        perturbed_state[n] = state[n];

                    Set::Scalar delta =
                        std::sqrt(std::numeric_limits<Set::Scalar>::epsilon()) *
                        (std::abs(state[column]) > 1.0 ?
                            std::abs(state[column]) : 1.0);
                    perturbed_state[column] += delta;
                    if (!evaluate_residual(perturbed_state, perturbed_residual))
                    {
                        delta = -delta;
                        perturbed_state[column] = state[column] + delta;
                        if (!evaluate_residual(perturbed_state, perturbed_residual))
                        {
                            failed = true;
                            break;
                        }
                    }

                    for (int row = 0; row < size; ++row)
                        jacobian[row][column] =
                            (perturbed_residual[row] - residual[row]) / delta;
                }
            }
            for (int n = 0; n < size; ++n) update[n] = -residual[n];

            for (int column = 0; column < size && !failed; ++column)
            {
                int pivot = column;
                for (int row = column + 1; row < size; ++row)
                    if (std::abs(jacobian[row][column]) >
                        std::abs(jacobian[pivot][column])) pivot = row;

                for (int n = column; n < size; ++n)
                {
                    const Set::Scalar tmp = jacobian[column][n];
                    jacobian[column][n] = jacobian[pivot][n];
                    jacobian[pivot][n] = tmp;
                }
                {
                    const Set::Scalar tmp = update[column];
                    update[column] = update[pivot];
                    update[pivot] = tmp;
                }

                if (!(std::abs(jacobian[column][column]) > 1.0e-30))
                {
                    failed = true;
                    break;
                }

                for (int row = column + 1; row < size; ++row)
                {
                    const Set::Scalar factor =
                        jacobian[row][column] / jacobian[column][column];
                    for (int n = column + 1; n < size; ++n)
                        jacobian[row][n] -= factor * jacobian[column][n];
                    update[row] -= factor * update[column];
                }
            }
            if (failed) return result;

            for (int row = size - 1; row >= 0; --row)
            {
                for (int n = row + 1; n < size; ++n)
                    update[row] -= jacobian[row][n] * update[n];
                update[row] /= jacobian[row][row];
            }

            bool accepted = false;
            Set::Scalar alpha = 1.0;
            for (int line_search = 0; line_search < 24; ++line_search)
            {
                Set::Scalar trial[MAX_SIZE]{};
                Set::Scalar trial_residual[MAX_SIZE]{};
                for (int n = 0; n < size; ++n)
                    trial[n] = state[n] + alpha * update[n];

                if (evaluate_residual(trial, trial_residual))
                {
                    const Set::Scalar trial_norm = norm(trial_residual);
                    if (trial_norm <= tolerance || trial_norm < result.residual_norm)
                    {
                        for (int n = 0; n < size; ++n)
                        {
                            state[n] = trial[n];
                            residual[n] = trial_residual[n];
                        }
                        result.residual_norm = trial_norm;
                        result.converged = result.residual_norm <= tolerance;
                        accepted = true;
                        break;
                    }
                }
                alpha *= 0.5;
            }
            if (!accepted) return result;
        }

        return result;
    }

    static void Parse(BackwardEuler& value, IO::ParmParse& pp)
    {
        Util::Message(INFO, pp.getPrefix());
        pp.query_default("nsubsteps", value.nsubsteps, 1);
        pp.query_default("max_iter", value.max_iterations, 100);
        pp.query_default("reltol", value.relative_tolerance, 1.0e-8);
        pp.query_default("abstol", value.absolute_tolerance, 1.0e-12);
        Util::Assert(INFO, TEST(value.nsubsteps > 0),
            "Backward Euler nsubsteps must be positive.");
        Util::Assert(INFO, TEST(value.max_iterations > 0),
            "Backward Euler max_iter must be positive.");
        Util::Assert(INFO,
            TEST(value.relative_tolerance >= 0.0 && value.absolute_tolerance >= 0.0),
            "Backward Euler tolerances must be nonnegative.");
        Util::Assert(INFO,
            TEST(value.relative_tolerance > 0.0 || value.absolute_tolerance > 0.0),
            "Backward Euler requires a positive tolerance.");
    }
};

} // namespace ODE
} // namespace Local
} // namespace Solver

#endif
