#ifndef NUMERIC_ADVECT_LIMITER_SUPERBEE_H
#define NUMERIC_ADVECT_LIMITER_SUPERBEE_H

#include "IO/ParmParse.H"
#include "Util/Util.H"

namespace Numeric
{
namespace Advect
{
namespace Limiter
{
class Superbee
{
public:
    static constexpr const char* name = "superbee";

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar operator() (Set::Scalar lo, Set::Scalar mid, Set::Scalar hi) const
    {
        Set::Scalar dl = mid - lo;
        Set::Scalar dr = hi - mid;
        Set::Scalar slope_a = Util::MinMod(Set::Scalar(2.0) * dl, dr);
        Set::Scalar slope_b = Util::MinMod(dl, Set::Scalar(2.0) * dr);
        return Util::Abs(slope_a) > Util::Abs(slope_b) ? slope_a : slope_b;
    }

    void Parse(Superbee & /*value*/, IO::ParmParse &/*pp*/) {}
};
}
}
}

#endif
