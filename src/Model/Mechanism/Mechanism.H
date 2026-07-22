#ifndef MODEL_MECHANISM_MECHANISM_H
#define MODEL_MECHANISM_MECHANISM_H

#include <tuple>

#include "AMReX_Extension.H"
#include "IO/ParmParse.H"
#include "Set/Set.H"
#include "Util/Util.H"

namespace Model
{
namespace Mechanism
{
// GPU-safe static dispatch for named material interactions. Parsing selects a
// concrete model; kernels use the selected tuple entry without virtual calls.
struct State
{
    amrex::Array4<const Set::Scalar> component_density;
    amrex::Array4<const Set::Scalar> rigid_eta;
    Set::Scalar temperature;
    Set::Scalar pressure;
};

template<typename... OBJ>
class Mechanism
{
public:
    template<int I = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar VolumeSource(const State& state,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).VolumeSource(
                    state, i, j, k, DX);
            return VolumeSource<I + 1>(state, i, j, k, DX);
        }
        Util::Abort(INFO, selected);
        return NAN;
    }

    template<int I = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Vector MomentumSource(const State& state,
                            const Set::Vector& gas_velocity,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).MomentumSource(
                    state, gas_velocity, i, j, k, DX);
            return MomentumSource<I + 1>(
                state, gas_velocity, i, j, k, DX);
        }
        Util::Abort(INFO, selected);
        return Set::Vector::Zero();
    }

    template<int I = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar GasMassSource(const State& state,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).GasMassSource(
                    state, i, j, k, DX);
            return GasMassSource<I + 1>(state, i, j, k, DX);
        }
        Util::Abort(INFO, selected);
        return NAN;
    }

    template<int I = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar GasSpeciesSource(const State& state, const int species,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).GasSpeciesSource(
                    state, species, i, j, k, DX);
            return GasSpeciesSource<I + 1>(
                state, species, i, j, k, DX);
        }
        Util::Abort(INFO, selected);
        return NAN;
    }

    template<int I = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar GasInjectionTemperature(
                            const State& state, const int species,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).GasInjectionTemperature(
                    state, species, i, j, k, DX);
            return GasInjectionTemperature<I + 1>(
                state, species, i, j, k, DX);
        }
        Util::Abort(INFO, selected);
        return NAN;
    }

    template<int I = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Vector InterfaceVelocity(const State& state,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).InterfaceVelocity(
                    state, i, j, k, DX);
            return InterfaceVelocity<I + 1>(state, i, j, k, DX);
        }
        Util::Abort(INFO, selected);
        return Set::Vector::Zero();
    }

    template<int I = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar RegressionSpeed(const State& state,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).RegressionSpeed(
                    state, i, j, k, DX);
            return RegressionSpeed<I + 1>(state, i, j, k, DX);
        }
        Util::Abort(INFO, selected);
        return NAN;
    }

    template<int I = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Apply(const amrex::Array4<Set::Scalar>& component_density_rhs,
                const State& state,
                const int i, const int j, const int k,
                const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
            {
                std::get<I>(obj).Apply(component_density_rhs, state, i, j, k, DX);
                return;
            }
            Apply<I + 1>(component_density_rhs, state, i, j, k, DX);
            return;
        }
        Util::Abort(INFO, selected);
    }

    template<int I = 0>
    Set::Scalar StabilityRate(const Set::Scalar dx,
                            const Set::Scalar temperature) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).StabilityRate(dx, temperature);
            return StabilityRate<I + 1>(dx, temperature);
        }
        Util::Abort(INFO, selected);
        return NAN;
    }

private:
    friend IO::ParmParse;
    std::tuple<OBJ...> obj;
    static constexpr const char* names[] = {OBJ::name...};
    int selected = -1;
};
}
}

#endif
