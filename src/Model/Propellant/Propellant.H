#ifndef MODEL_PROPELLANT_H
#define MODEL_PROPELLANT_H

#include "AMReX_Extension.H"
#include "IO/ParmParse.H"

namespace Model
{
namespace Propellant
{
template<typename... OBJ>
class Propellant
{
private:
    friend IO::ParmParse;
    std::tuple<OBJ...> obj;
    static constexpr const char* names[] = { OBJ::name... };
    int selected = -1;

public:

    template<int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    void set_pressure (Set::Scalar P)
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::get<I>(obj).set_pressure(P);
            else return set_pressure<I+1>(P);
        }
        else Util::Abort(INFO,selected);
    }

    template<int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar get_K (const Set::Scalar phi)
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::get<I>(obj).get_K(phi);
            else return get_K<I+1>(phi);
        }
        else Util::Abort(INFO,selected);
        return NAN;
    }

    template<int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar get_rho (const Set::Scalar phi)
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::get<I>(obj).get_rho(phi);
            else return get_rho<I+1>(phi);
        }
        else Util::Abort(INFO,selected);
        return NAN;
    }

    template<int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar get_cp (const Set::Scalar phi)
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::get<I>(obj).get_cp(phi);
            else return get_cp<I+1>(phi);
        }
        else Util::Abort(INFO,selected);
        return NAN;
    }

    template<int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar get_qdot (const Set::Scalar mdot, const Set::Scalar phi) 
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::get<I>(obj).get_qdot(mdot,phi);
            else return get_qdot<I+1>(mdot,phi);
        }
        else Util::Abort(INFO,selected);
        return NAN;
    }

    template<int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar get_L (const Set::Scalar phi, const Set::Scalar T) 
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::get<I>(obj).get_L(phi,T);
            else return get_L<I+1>(phi, T);
        }
        else Util::Abort(INFO,selected);
        return NAN;
    }
};
}
}


#endif


