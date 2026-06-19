#ifndef MODEL_GAS_TRANSPORT_H_
#define MODEL_GAS_TRANSPORT_H_

#include <vector>
#include <memory>
#include <cmath>
#include "Set/Set.H"
#include "IO/ParmParse.H"

namespace Model {
namespace Gas {
namespace Transport {

template<typename... OBJ>
class Transport {
private:
    friend IO::ParmParse;
    std::tuple<OBJ...> obj;
    static constexpr const char* names[] = { OBJ::name... };
    int selected = -1;
public:

    Transport() {} ;

    template<int I = 0>
    AMREX_FORCE_INLINE
    AMREX_GPU_HOST_DEVICE
    double dynamic_viscosity(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::get<I>(obj).dynamic_viscosity(T,X,i,j,k);
            else return dynamic_viscosity<I+1>(T,X,i,j,k);
        }
        else {
            Util::Abort(INFO, selected);
            return NAN;
        }
    }


    template<int I = 0>
    AMREX_FORCE_INLINE
    AMREX_GPU_HOST_DEVICE
    double thermal_conductivity(double T, Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::get<I>(obj).thermal_conductivity(T,X,i,j,k);
            else return thermal_conductivity<I+1>(T,X,i,j,k);
        }
        else {
            Util::Abort(INFO, selected);
            return NAN;
        }
    }

    template<int I = 0>
    AMREX_FORCE_INLINE
    AMREX_GPU_HOST_DEVICE
    void diffusion_coeffs(  Set::Patch<Set::Scalar>& DKM, double T, double P,
                            Set::Patch<const Set::Scalar>& X, int i, int j, int k) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::get<I>(obj).diffusion_coeffs(DKM,T,P,X,i,j,k);
            else return diffusion_coeffs<I+1>(DKM,T,P,X,i,j,k);
        }
        else {
            Util::Abort(INFO, selected);
        }
    }

    ~Transport() = default;
}; // class Transport

} // namespace Transport
} // namespace Gas
} // namespace Model

#endif
