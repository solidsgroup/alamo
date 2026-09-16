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
    amrex::Array4<const Set::Scalar> rigid_species_eta;
    Set::Scalar temperature;
    Set::Scalar pressure;
    Set::Scalar dt;
};

// Compile-time mechanism interface.  Models inherit these neutral operations
// and define only the physics they actually supply; Mechanism below still
// dispatches to a concrete tuple entry with no virtual calls in GPU kernels.
class MechanismModel
{
public:
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void Apply(const amrex::Array4<Set::Scalar>&, const State&,
               const int, const int, const int,
               const Set::Scalar[AMREX_SPACEDIM]) const {}
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar VolumeSource(const State&, const int, const int, const int,
        const Set::Scalar[AMREX_SPACEDIM]) const { return 0.0; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar HeatSource(const State&, const int, const int, const int,
        const Set::Scalar[AMREX_SPACEDIM]) const { return 0.0; }
    Set::Scalar StabilityRate(const Set::Scalar,
        const Set::Scalar) const { return 0.0; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar LocalRate(const State&, const int, const int,
        const int) const { return 0.0; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar LocalStabilityRate(const State&, const int, const int,
        const int) const { return 0.0; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar GradientCoefficient(const State&, const int, const int,
        const int) const { return 0.0; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar ApplyImplicitChange(
        const amrex::Array4<Set::Scalar>&, const State&, const Set::Scalar,
        Set::Scalar&, const int, const int, const int) const { return 0.0; }
    int RigidComponent() const { return -1; }
    bool Equilibrium() const { return false; }
    bool Kinetic() const { return false; }
    Set::Scalar EquilibriumTemperature() const { return NAN; }
    Set::Scalar LatentHeat() const { return 0.0; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::GpuTuple<Set::Scalar,Set::Scalar,Set::Scalar>
    ApplyEquilibriumChange(
        const amrex::Array4<Set::Scalar>&, const State&, const Set::Scalar,
        const Set::Scalar, const Set::Scalar,
        const int, const int, const int) const { return {0.0, 0.0, 0.0}; }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::GpuTuple<Set::Scalar,Set::Scalar,Set::Scalar,Set::Scalar>
    ApplyKineticChange(
        const amrex::Array4<Set::Scalar>&, const State&, const Set::Scalar,
        const Set::Scalar, Set::Scalar&, Set::Scalar&,
        const int, const int, const int) const
    {
        return {0.0, 0.0, 0.0, 0.0};
    }
    template<class Advection>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::GpuTuple<Set::Scalar,Set::Vector>
    KineticDiffuseSurfaceGeometry(
        const amrex::Array4<const Set::Scalar>&,
        const int, const int, const int,
        const Set::Scalar[AMREX_SPACEDIM], const Advection&) const
    {
        return {0.0, Set::Vector::Zero()};
    }
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar RecoilPressure(const Set::Scalar, const Set::Scalar,
        const Set::Scalar) const { return 0.0; }
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
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar HeatSource(const State& state,
                            const int i, const int j, const int k,
                            const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).HeatSource(
                    state, i, j, k, DX);
            return HeatSource<I + 1>(state, i, j, k, DX);
        }
        Util::Abort(INFO, selected);
        return NAN;
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

    template<int I = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar LocalRate(
        const State& state, const int i, const int j, const int k) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).LocalRate(state, i, j, k);
            return LocalRate<I + 1>(state, i, j, k);
        }
        Util::Abort(INFO, selected);
        return NAN;
    }

    template<int I = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar LocalStabilityRate(
        const State& state, const int i, const int j, const int k) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).LocalStabilityRate(
                    state, i, j, k);
            return LocalStabilityRate<I + 1>(state, i, j, k);
        }
        Util::Abort(INFO, selected);
        return NAN;
    }

    template<int I = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar GradientCoefficient(
        const State& state, const int i, const int j, const int k) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).GradientCoefficient(
                    state, i, j, k);
            return GradientCoefficient<I + 1>(state, i, j, k);
        }
        Util::Abort(INFO, selected);
        return NAN;
    }

    template<int I = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar ApplyImplicitChange(
        const amrex::Array4<Set::Scalar>& component_density,
        const State& state, const Set::Scalar eta_change,
        Set::Scalar& integrated_heat,
        const int i, const int j, const int k) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).ApplyImplicitChange(
                    component_density, state, eta_change, integrated_heat,
                    i, j, k);
            return ApplyImplicitChange<I + 1>(
                component_density, state, eta_change, integrated_heat,
                i, j, k);
        }
        Util::Abort(INFO, selected);
        return NAN;
    }

    template<int I = 0>
    int RigidComponent() const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).RigidComponent();
            return RigidComponent<I + 1>();
        }
        Util::Abort(INFO, selected);
        return -1;
    }

    template<int I = 0>
    bool Equilibrium() const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::get<I>(obj).Equilibrium();
            return Equilibrium<I + 1>();
        }
        Util::Abort(INFO, selected);
        return false;
    }

    template<int I = 0>
    bool Kinetic() const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::get<I>(obj).Kinetic();
            return Kinetic<I + 1>();
        }
        Util::Abort(INFO, selected);
        return false;
    }

    template<int I = 0>
    Set::Scalar EquilibriumTemperature() const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).EquilibriumTemperature();
            return EquilibriumTemperature<I + 1>();
        }
        Util::Abort(INFO, selected);
        return NAN;
    }

    template<int I = 0>
    Set::Scalar LatentHeat() const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::get<I>(obj).LatentHeat();
            return LatentHeat<I + 1>();
        }
        Util::Abort(INFO, selected);
        return NAN;
    }

    template<int I = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::GpuTuple<Set::Scalar,Set::Scalar,Set::Scalar>
    ApplyEquilibriumChange(
        const amrex::Array4<Set::Scalar>& component_density,
        const State& state, const Set::Scalar requested_mass_change,
        const Set::Scalar heat_capacity,
        const Set::Scalar gas_volume_fraction,
        const int i, const int j, const int k) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).ApplyEquilibriumChange(
                    component_density, state, requested_mass_change,
                    heat_capacity, gas_volume_fraction,
                    i, j, k);
            return ApplyEquilibriumChange<I + 1>(
                component_density, state, requested_mass_change,
                heat_capacity, gas_volume_fraction, i, j, k);
        }
        Util::Abort(INFO, selected);
        return {NAN, NAN, NAN};
    }

    template<class Advection, int I = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::GpuTuple<Set::Scalar,Set::Vector>
    KineticDiffuseSurfaceGeometry(
        const amrex::Array4<const Set::Scalar>& component_density,
        const int i, const int j, const int k,
        const Set::Scalar DX[AMREX_SPACEDIM], const Advection& advection) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).KineticDiffuseSurfaceGeometry(
                    component_density, i, j, k, DX, advection);
            return KineticDiffuseSurfaceGeometry<Advection, I + 1>(
                component_density, i, j, k, DX, advection);
        }
        Util::Abort(INFO, selected);
        return {NAN, Set::Vector::Zero()};
    }

    template<int I = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar RecoilPressure(const Set::Scalar mass_flux,
        const Set::Scalar temperature, const Set::Scalar pressure) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).RecoilPressure(
                    mass_flux, temperature, pressure);
            return RecoilPressure<I + 1>(
                mass_flux, temperature, pressure);
        }
        Util::Abort(INFO, selected);
        return NAN;
    }

    template<int I = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    amrex::GpuTuple<Set::Scalar,Set::Scalar,Set::Scalar,Set::Scalar>
    ApplyKineticChange(
        const amrex::Array4<Set::Scalar>& component_density,
        const State& state, const Set::Scalar heat_capacity,
        const Set::Scalar surface_measure,
        Set::Scalar& temperature, Set::Scalar& integrated_heat,
        const int i, const int j, const int k) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).ApplyKineticChange(
                    component_density, state, heat_capacity,
                    surface_measure,
                    temperature, integrated_heat, i, j, k);
            return ApplyKineticChange<I + 1>(
                component_density, state, heat_capacity,
                surface_measure,
                temperature, integrated_heat, i, j, k);
        }
        Util::Abort(INFO, selected);
        return {NAN, NAN, NAN, NAN};
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
