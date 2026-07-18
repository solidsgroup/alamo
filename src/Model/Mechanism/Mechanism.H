#ifndef MODEL_MECHANISM_MECHANISM_H
#define MODEL_MECHANISM_MECHANISM_H

#include <string>
#include <tuple>
#include <utility>

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
    template<typename... Args>
    void Parse(IO::ParmParse& pp, const std::string& id, Args&&... args)
    {
        std::string type;
        pp.query_required(id + ".type", type);
        pp.pushPrefix(id);
        try
        {
            ParseSelected(type, pp, std::forward<Args>(args)...);
        }
        catch (...)
        {
            pp.popPrefix();
            throw;
        }
        pp.popPrefix();
    }

    template<int I = 0>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Set::Scalar VolumeSource(const State& state,
                             const int i, const int j, const int k,
                             const Set::Scalar DX[AMREX_SPACEDIM]) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(objects).VolumeSource(
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
                std::get<I>(objects).Apply(component_density_rhs, state, i, j, k, DX);
                return;
            }
            Apply<I + 1>(component_density_rhs, state, i, j, k, DX);
            return;
        }
        Util::Abort(INFO, selected);
    }

    template<int I = 0>
    Set::Scalar StabilityRate(const Set::Scalar dx) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::get<I>(objects).StabilityRate(dx);
            return StabilityRate<I + 1>(dx);
        }
        Util::Abort(INFO, selected);
        return NAN;
    }

private:
    template<int I = 0, typename... Args>
    void ParseSelected(const std::string& type, IO::ParmParse& pp, Args&&... args)
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (type == names[I])
            {
                selected = I;
                auto& object = std::get<I>(objects);
                object.Parse(object, pp, std::forward<Args>(args)...);
                return;
            }
            ParseSelected<I + 1>(type, pp, std::forward<Args>(args)...);
        }
        else
            Util::Exception(INFO, type, " is not a valid mechanism type");
    }

    std::tuple<OBJ...> objects;
    static constexpr const char* names[] = {OBJ::name...};
    int selected = -1;
};
}
}

#endif
