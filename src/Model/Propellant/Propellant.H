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

    std::string get_name()
    {
        return names[selected];
    }

    template<int I = 0>
    AMREX_GPU_HOST
    Set::Scalar get (std::string property)
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::get<I>(obj).get(property);
            else return get<I+1>(property);
        }
        else Util::Abort(INFO,selected);
        return NAN;
    }


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


