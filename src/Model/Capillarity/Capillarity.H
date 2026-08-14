#ifndef MODEL_CAPILLARITY_CAPILLARITY_H
#define MODEL_CAPILLARITY_CAPILLARITY_H

#include <tuple>
#include <type_traits>

#include "IO/ParmParse.H"
#include "Util/Util.H"

namespace Model
{
namespace Capillarity
{
// Input-selected static dispatch for capillarity and interface kinetics.
// Concrete models contain their parameters and local constitutive equations;
// LowMach only orchestrates the required composite AMR operators.
template<typename... Models>
class Capillarity
{
    friend IO::ParmParse;

public:
    const char* ModelName() const
    {
        return selected >= 0 ? names[selected] : "";
    }

    template<typename Model>
    bool Is() const
    {
        return selected == TupleIndex<Model>();
    }

    template<typename Model>
    const Model& Get() const
    {
        Util::Assert(INFO, TEST(Is<Model>()));
        return std::get<Model>(obj);
    }

private:
    template<typename Model, std::size_t I = 0>
    static consteval int TupleIndex()
    {
        if constexpr (I == sizeof...(Models))
            return -1;
        else if constexpr (std::is_same_v<
            Model, std::tuple_element_t<I,std::tuple<Models...>>>)
            return static_cast<int>(I);
        else
            return TupleIndex<Model,I + 1>();
    }

    std::tuple<Models...> obj;
    static constexpr const char* names[] = {Models::name...};
    int selected = -1;
};
}
}

#endif
