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
    // Temperature, K
    template<int I = 0>
    AMREX_FORCE_INLINE
    AMREX_GPU_HOST_DEVICE
    double ComputeT(double density, double momentumx, double momentumy, double E, 
                    double R, double gamma, double rtol) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).ComputeT(density,momentumx,momentumy,E,R,gamma, rtol);
            else
                return ComputeT<I+1>(density,momentumx,momentumy,E,R,gamma, rtol);
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
    double ComputeT(double pressure, double density, 
                    double R) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).ComputeT(pressure, density, R);
            else
                return ComputeT<I+1>(pressure, density, R);
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

    template<int I = 0>
    AMREX_FORCE_INLINE
    AMREX_GPU_HOST_DEVICE
    double ComputeE(double density, double momentumx, double momentumy, double T, 
                    double R, double gamma) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj).ComputeE(density, momentumx, momentumy, T, R, gamma);
            else
                return ComputeE<I+1>(density, momentumx, momentumy, T, R, gamma);
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
