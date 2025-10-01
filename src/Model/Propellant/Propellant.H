//
// This header defines the :code:`Model::Propellant::Propellant` class, which introduces a static
// dispatch mechanism to replicate traditional dynamic polymorphism without incurring the performance
// cost of virtual table lookups.
//
// The :code:`Propellant` class acts as the static equivalent of a pure abstract interface for
// composite solid propellant models. Each function (e.g., :code:`get_K`, :code:`get_cp`,
// :code:`get_L`) statically dispatches to the appropriate underlying implementation, encapsulating
// key thermophysical properties such as  density, thermal conductivity, and specific heat. The
// :code:`set_pressure` and :code:`get_*` functions are resolved at compile-time for GPU
// compatibility and runtime performance.
//
// This system enables polymorphic behavior without relying on runtime virtual functions, making it
// suitable for high-performance computing environments.
//

#ifndef MODEL_PROPELLANT_H
#define MODEL_PROPELLANT_H

#include "AMReX_Extension.H"
#include "AMReX_GpuQualifiers.H"
#include "IO/ParmParse.H"
#include "Set/Base.H"
#include "Util/Util.H"

#include <array>
#include <cmath>
#include <cstddef>
#include <tuple>

namespace Model
{
namespace Propellant
{
template <std::size_t N_SPECIES>
using PhiArray = std::array<Set::Scalar, N_SPECIES - 1>;

template <typename... OBJ>
class Propellant
{
private:
    friend IO::ParmParse;
    std::tuple<OBJ...> obj;
    static constexpr const char *names[] = { OBJ::name... };
    int selected = -1;

public:
    template <int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
        std::size_t
        get_n_species() const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
            {
                using ModelType = std::tuple_element_t<I, std::tuple<OBJ...>>;
                return ModelType::n_species;
            }
            else
                return get_n_species<I + 1>();
        }
        else
            Util::Abort(INFO, selected);
        return 0;
    }

    template <int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE void set_pressure(Set::Scalar P)
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).set_pressure(P);
            else
                return set_pressure<I + 1>(P);
        }
        else
            Util::Abort(INFO, selected);
    }

    template <std::size_t N_SPECIES, int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
        Set::Scalar
        get_K(const PhiArray<N_SPECIES> &phi)
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
            {
                using ModelType = std::tuple_element_t<I, std::tuple<OBJ...>>;
                static_assert(ModelType::n_species == N_SPECIES,
                              "Mismatch between provided N_SPECIES and model's n_species");
                return std::get<I>(obj).get_K(phi);
            }
            else
                return get_K<N_SPECIES, I + 1>(phi);
        }
        else
            Util::Abort(INFO, selected);
        return NAN;
    }

    template <std::size_t N_SPECIES, int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
        Set::Scalar
        get_rho(const PhiArray<N_SPECIES> &phi)
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
            {
                using ModelType = std::tuple_element_t<I, std::tuple<OBJ...>>;
                static_assert(ModelType::n_species == N_SPECIES,
                              "Mismatch between provided N_SPECIES and model's n_species");
                return std::get<I>(obj).get_rho(phi);
            }
            else
                return get_rho<N_SPECIES, I + 1>(phi);
        }
        else
            Util::Abort(INFO, selected);
        return NAN;
    }

    template <std::size_t N_SPECIES, int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
        Set::Scalar
        get_cp(const PhiArray<N_SPECIES> &phi)
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
            {
                using ModelType = std::tuple_element_t<I, std::tuple<OBJ...>>;
                static_assert(ModelType::n_species == N_SPECIES,
                              "Mismatch between provided N_SPECIES and model's n_species");
                return std::get<I>(obj).get_cp(phi);
            }
            else
                return get_cp<N_SPECIES, I + 1>(phi);
        }
        else
            Util::Abort(INFO, selected);
        return NAN;
    }

    template <std::size_t N_SPECIES, int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
        Set::Scalar
        get_qdot(const Set::Scalar mdot, const PhiArray<N_SPECIES> &phi)
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
            {
                using ModelType = std::tuple_element_t<I, std::tuple<OBJ...>>;
                static_assert(ModelType::n_species == N_SPECIES,
                              "Mismatch between provided N_SPECIES and model's n_species");
                return std::get<I>(obj).get_qdot(mdot, phi);
            }
            else
                return get_qdot<N_SPECIES, I + 1>(mdot, phi);
        }
        else
            Util::Abort(INFO, selected);
        return NAN;
    }

    template <std::size_t N_SPECIES, int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
        Set::Scalar
        get_L(const PhiArray<N_SPECIES> &phi, const Set::Scalar T)
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
            {
                using ModelType = std::tuple_element_t<I, std::tuple<OBJ...>>;
                static_assert(ModelType::n_species == N_SPECIES,
                              "Mismatch between provided N_SPECIES and model's n_species");
                return std::get<I>(obj).get_L(phi, T);
            }
            else
                return get_L<N_SPECIES, I + 1>(phi, T);
        }
        else
            Util::Abort(INFO, selected);
        return NAN;
    }
};
}
}

#endif
