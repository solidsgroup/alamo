#ifndef SOLVER_LOCAL_ODE_ODE_H
#define SOLVER_LOCAL_ODE_ODE_H

#include <limits>
#include <tuple>
#include <type_traits>

#include "AMReX_Extension.H"
#include "IO/ParmParse.H"
#include "Set/Base.H"

namespace Solver
{
namespace Local
{
namespace ODE
{

struct Advance
{
    struct Result
    {
        bool converged = false;
        Set::Scalar residual_norm =
            std::numeric_limits<Set::Scalar>::infinity();
    };
};

template<typename... OBJ>
class ODE
{
private:
    friend IO::ParmParse;
    std::tuple<OBJ...> obj;
    static constexpr const char* names[] = {OBJ::name...};
    int selected = -1;

    template<typename T, int I>
    void SelectImplementation()
    {
        if constexpr (I < sizeof...(OBJ))
        {
            using Candidate = std::tuple_element_t<I, std::tuple<OBJ...>>;
            if constexpr (std::is_same_v<T, Candidate>) selected = I;
            else SelectImplementation<T, I + 1>();
        }
        else
        {
            static_assert(!std::is_same_v<T, T>, "ODE solver is not in this ODE dispatch set");
        }
    }

public:
    template<typename T>
    void Select()
    {
        SelectImplementation<T, 0>();
    }

    template<typename T>
    T& Get()
    {
        return std::get<T>(obj);
    }

    template<int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    bool Implicit() const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::tuple_element_t<I, std::tuple<OBJ...>>::implicit;
            return Implicit<I + 1>();
        }
        else
        {
            Util::Abort(INFO, selected);
            return false;
        }
    }

    template<int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    bool RequiresJacobian() const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::tuple_element_t<I, std::tuple<OBJ...>>::requires_jacobian;
            return RequiresJacobian<I + 1>();
        }
        else
        {
            Util::Abort(INFO, selected);
            return false;
        }
    }

    template<int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    bool SupportsJacobian() const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::tuple_element_t<I, std::tuple<OBJ...>>::supports_jacobian;
            return SupportsJacobian<I + 1>();
        }
        else
        {
            Util::Abort(INFO, selected);
            return false;
        }
    }

    template<int MAX_SIZE, typename RHS, typename Jacobian, int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Solver::Local::ODE::Advance::Result
    Advance(
        Set::Scalar* state,
        int size,
        Set::Scalar dt,
        const RHS& rhs,
        const Jacobian& jacobian,
        bool has_jacobian) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
            {
                Solver::Local::ODE::Advance::Result result;
                const int nsubsteps = std::get<I>(obj).nsubsteps;
                if (nsubsteps <= 0 || size <= 0 || size > MAX_SIZE ||
                    !(dt > 0.0)) return result;

                Set::Scalar substep_state[MAX_SIZE]{};
                for (int n = 0; n < size; ++n) substep_state[n] = state[n];

                const Set::Scalar substep_dt = dt / nsubsteps;
                for (int substep = 0; substep < nsubsteps; ++substep)
                {
                    const auto substep_result =
                        std::get<I>(obj).template Advance<MAX_SIZE>(
                            substep_state, size, substep_dt, rhs, jacobian,
                            has_jacobian);
                    if (!substep_result.converged) return substep_result;
                    if (substep == 0 ||
                        substep_result.residual_norm > result.residual_norm)
                        result.residual_norm = substep_result.residual_norm;
                }

                for (int n = 0; n < size; ++n) state[n] = substep_state[n];
                result.converged = true;
                return result;
            }
            return Advance<MAX_SIZE, RHS, Jacobian, I + 1>(
                state, size, dt, rhs, jacobian, has_jacobian);
        }
        else
        {
            Util::Abort(INFO, selected);
            return {};
        }
    }
};

} // namespace ODE
} // namespace Local
} // namespace Solver

#endif
