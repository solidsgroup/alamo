#ifndef MODEL_SURFACECOMBUSTION_SURFACECOMBUSTION_H
#define MODEL_SURFACECOMBUSTION_SURFACECOMBUSTION_H

#include "IO/ParmParse.H"

namespace Model
{
namespace SurfaceCombustion
{
template<typename... OBJ>
class SurfaceCombustion
{
private:
    friend IO::ParmParse;
    std::tuple<OBJ...> obj;
    static constexpr const char* names[] = { OBJ::name... };
    int selected = -1;

public:

    template<int I = 0>
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
    Set::Scalar get_qdot (const Set::Scalar mdot, const Set::Scalar phi)
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I) return std::get<I>(obj).get_qdot(mdot,phi);
            else return get_qdot<I+1>(mdot,phi);
        }
        else Util::Abort(INFO,selected);
    }


};
}
}


#endif


