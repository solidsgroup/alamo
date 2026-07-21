#ifndef NUMERIC_ADVECT_LIMITER_LIMITER_H
#define NUMERIC_ADVECT_LIMITER_LIMITER_H


#include <tuple>
#include "Set/Set.H"
#include "IO/ParmParse.H"

namespace Numeric
{
namespace Advect
{
namespace Limiter
{
template<typename... OBJ>
class Limiter
{
private:
    friend IO::ParmParse;
    std::tuple<OBJ...> obj;
    static constexpr const char* names[] = { OBJ::name... };
    int selected = -1;

public:
    template<int I = 0>
    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar operator() (Set::Scalar lo, Set::Scalar mid, Set::Scalar hi) const
    {
        if constexpr (I < sizeof...(OBJ))
        {
            if (selected == I)
                return std::get<I>(obj)(lo, mid, hi);
            else
                return this->operator ()<I + 1>(lo, mid, hi);
        }
        else
        {
            Util::Abort(INFO, selected);
            return NAN;
        }
    }

};
}
}
}



#endif
