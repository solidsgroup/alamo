//
// This header defines the :code:Model::Chemistry::Chemistry class, which introduces a static
// dispatch mechanism to replicate traditional dynamic polymorphism without incurring the performance
// cost of virtual table lookups.
//
// The :code:Chemistry class acts as the static equivalent of a pure abstract interface for
// chemical reaction models. The function :code:ComputeChemistrySources statically dispatches to the
// appropriate underlying implementation, encapsulating key properties for reaction rates
// and heat release. The :code:ComputeChemistrySources function is resolved at compile-time for GPU
// compatibility and runtime performance, and returns the species production rates and heat release.
//
// This system enables polymorphic behavior without relying on runtime virtual functions, making it
// suitable for high-performance computing environments.
//

#ifndef MODEL_CHEMISTRY_H
#define MODEL_CHEMISTRY_H

#include "AMReX_Extension.H"
#include "IO/ParmParse.H"
#include "Model/Gas/Gas.H"

namespace Model
{
namespace Chemistry
{
static constexpr int MAX_SPECIES = 32;
using SpeciesArray = std::array<Set::Scalar, MAX_SPECIES>;
using Source = std::pair<SpeciesArray, Set::Scalar>;

template<typename... OBJ>
class Chemistry
{
private:
    friend IO::ParmParse;
    std::tuple<OBJ...> obj;
    static constexpr const char* names[] = { OBJ::name... };
    int selected = -1;

public:

    template<int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    bool Reactive() const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::tuple_element_t<I, std::tuple<OBJ...>>::reactive;
            return Reactive<I+1>();
        }
        else
        {
            Util::Abort(INFO,selected);
            return false;
        }
    }

    template<int I = 0>
    bool Implicit() const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::get<I>(obj).Implicit();
            return Implicit<I+1>();
        }
        else
        {
            Util::Abort(INFO,selected);
            return false;
        }
    }

    template<int I = 0>
    int MaxIterations() const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::get<I>(obj).MaxIterations();
            return MaxIterations<I+1>();
        }
        else
        {
            Util::Abort(INFO,selected);
            return 0;
        }
    }

    template<int I = 0>
    std::pair<Set::Scalar,Set::Scalar> Tolerances() const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::get<I>(obj).Tolerances();
            return Tolerances<I+1>();
        }
        else
        {
            Util::Abort(INFO,selected);
            return {0.0,0.0};
        }
    }

    template<int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Source
    ComputeChemistrySources (const Set::Scalar P, const Set::Scalar T, const SpeciesArray& rhoY,
                const Set::Scalar dt, const Model::Gas::Gas* gas = nullptr) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::get<I>(obj).ComputeChemistrySources(P,T,rhoY,dt,gas);
            else return ComputeChemistrySources<I+1>(P,T,rhoY,dt,gas);
        }
        else
        {
            Util::Abort(INFO,selected);
            return {};
        }
    }
};
}
}


#endif
