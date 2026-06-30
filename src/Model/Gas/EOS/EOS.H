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
    template<typename GAS, int I = 0>
    AMREX_FORCE_INLINE
    AMREX_GPU_HOST_DEVICE
    double ComputeT(const GAS &gas, double density, double momentumx, double momentumy, double E, double Tguess,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k, double rtol) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).ComputeT(gas,density,momentumx,momentumy,E,Tguess,X,i,j,k,rtol);
            else
                return ComputeT<GAS, I+1>(gas,density,momentumx,momentumy,E,Tguess,X,i,j,k,rtol);
        }
        else {
            Util::Abort(INFO, selected);
            return NAN;
        }
    }

    // Temperature, K
    template<typename GAS, int I = 0>
    AMREX_FORCE_INLINE
    AMREX_GPU_HOST_DEVICE
    double ComputeT(const GAS &gas, double density, double momentumx, double momentumy, double momentumz, double E, double Tguess,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k, double rtol) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).ComputeT(gas,density,momentumx,momentumy,momentumz,E,Tguess,X,i,j,k,rtol);
            else
                return ComputeT<GAS, I+1>(gas,density,momentumx,momentumy,momentumz,E,Tguess,X,i,j,k,rtol);
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
    double ComputeT_from_primitives(double pressure, double density, 
                    double R) const
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
    double ComputeP(double density, double T, double R) const
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
    template<typename GAS, int I = 0>
    AMREX_FORCE_INLINE
    AMREX_GPU_HOST_DEVICE
    double ComputeE(const GAS &gas, double density, double momentumx, double momentumy, double T,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).ComputeE(gas, density, momentumx, momentumy, T, X, i, j, k);
            else
                return ComputeE<GAS, I+1>(gas, density, momentumx, momentumy, T, X, i, j, k);
        }
        else {
            Util::Abort(INFO, selected);
            return NAN;
        }
    }

    // Energy, J
    template<typename GAS, int I = 0>
    AMREX_FORCE_INLINE
    AMREX_GPU_HOST_DEVICE
    double ComputeE(const GAS &gas, double density, double momentumx, double momentumy, double momentumz, double T,
                    Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).ComputeE(gas, density, momentumx, momentumy, momentumz, T, X, i, j, k);
            else
                return ComputeE<GAS, I+1>(gas, density, momentumx, momentumy, momentumz, T, X, i, j, k);
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
