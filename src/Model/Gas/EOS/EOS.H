//
// This header defines the :code:Model::Gas::EOS::EOS:code: static dispatch wrapper
// for gas equation-of-state models.
//
// The :code:EOS:code: class acts as the static equivalent of a polymorphic interface
// for pressure, temperature, and energy conversions. It dispatches to the selected EOS
// implementation without virtual calls, preserving GPU compatibility for per-cell gas
// state reconstruction.
//

#ifndef MODEL_GAS_EOS_H_
#define MODEL_GAS_EOS_H_

#include <vector>
#include <memory>
#include <cmath>

#include "Set/Set.H"
#include "IO/ParmParse.H"

namespace Model {
namespace Gas {
namespace EOS {

template<typename... OBJ>
class EOS {
private:
    friend IO::ParmParse;
    std::tuple<OBJ...> obj;
    static constexpr const char* names[] = { OBJ::name... };
    int selected = -1;

public:
    const char* model_name() const
    {
        return selected >= 0 ? names[selected] : "";
    }

    // Temperature, K
    template<typename GAS, class Composition, int I = 0>
    AMREX_FORCE_INLINE
    Set::Scalar ComputeT(const GAS &gas, Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar E, Set::Scalar Tguess,
                    const Composition& X, int i, int j, int k, Set::Scalar rtol) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).ComputeT(gas,density,momentumx,momentumy,E,Tguess,X,i,j,k,rtol);
            else
                return ComputeT<GAS, Composition, I+1>(gas,density,momentumx,momentumy,E,Tguess,X,i,j,k,rtol);
        }
        else {
            Util::Abort(INFO, selected);
            return NAN;
        }
    }

    // Temperature, K
    template<typename GAS, class Composition, int I = 0>
    AMREX_FORCE_INLINE
    Set::Scalar ComputeT(const GAS &gas, Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar momentumz, Set::Scalar E, Set::Scalar Tguess,
                    const Composition& X, int i, int j, int k, Set::Scalar rtol) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).ComputeT(gas,density,momentumx,momentumy,momentumz,E,Tguess,X,i,j,k,rtol);
            else
                return ComputeT<GAS, Composition, I+1>(gas,density,momentumx,momentumy,momentumz,E,Tguess,X,i,j,k,rtol);
        }
        else {
            Util::Abort(INFO, selected);
            return NAN;
        }
    }

    // Temperature, K
    template<int I = 0>
    AMREX_FORCE_INLINE
    AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeT_from_primitives(Set::Scalar pressure, Set::Scalar density,
                    Set::Scalar R) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).ComputeT_from_primitives(pressure, density, R);
            else
                return ComputeT_from_primitives<I+1>(pressure, density, R);
        }
        else {
            Util::Abort(INFO, selected);
            return NAN;
        }
    }


    // Pressure, Pa
    template<int I = 0>
    AMREX_FORCE_INLINE
    AMREX_GPU_HOST_DEVICE
    Set::Scalar ComputeP(Set::Scalar density, Set::Scalar T, Set::Scalar R) const
    // Energy, J/m^3
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).ComputeP(density, T, R);
            else
                return ComputeP<I+1>(density, T, R);
        }
        else {
            Util::Abort(INFO, selected);
            return NAN;
        }
    }

    // Energy, J
    template<typename GAS, class Composition, int I = 0>
    AMREX_FORCE_INLINE
    Set::Scalar ComputeE(const GAS &gas, Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar T,
                    const Composition& X, int i, int j, int k) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).ComputeE(gas, density, momentumx, momentumy, T, X, i, j, k);
            else
                return ComputeE<GAS, Composition, I+1>(gas, density, momentumx, momentumy, T, X, i, j, k);
        }
        else {
            Util::Abort(INFO, selected);
            return NAN;
        }
    }

    // Energy, J
    template<typename GAS, class Composition, int I = 0>
    AMREX_FORCE_INLINE
    Set::Scalar ComputeE(const GAS &gas, Set::Scalar density, Set::Scalar momentumx, Set::Scalar momentumy, Set::Scalar momentumz, Set::Scalar T,
                    const Composition& X, int i, int j, int k) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).ComputeE(gas, density, momentumx, momentumy, momentumz, T, X, i, j, k);
            else
                return ComputeE<GAS, Composition, I+1>(gas, density, momentumx, momentumy, momentumz, T, X, i, j, k);
        }
        else {
            Util::Abort(INFO, selected);
            return NAN;
        }
    }
}; // class EOS

} // namespace EOS
} // namespace Gas
} // namespace Model

#endif
