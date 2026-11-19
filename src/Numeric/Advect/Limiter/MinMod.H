#ifndef NUMERIC_ADVECT_LIMITER_MINMOD_H
#define NUMERIC_ADVECT_LIMITER_MINMOD_H

#include "IO/ParmParse.H"
#include "Util/Util.H"

namespace Numeric
{
namespace Advect
{
namespace Limiter
{
class MinMod
{
public:
    static constexpr const char* name = "minmod";

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar operator() (Set::Scalar lo, Set::Scalar mid, Set::Scalar hi) const
    {
        Set::Scalar dl = mid - lo;
        Set::Scalar dr = hi - mid;
        return Util::MinMod(dl, dr);
    }

    void Parse(MinMod & /*value*/, IO::ParmParse &/*pp*/) {}
};
}
}
}

#endif
