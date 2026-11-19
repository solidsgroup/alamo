#ifndef NUMERIC_ADVECT_LIMITER_VANLEER_H
#define NUMERIC_ADVECT_LIMITER_VANLEER_H

#include "IO/ParmParse.H"
#include "Util/Util.H"

namespace Numeric
{
namespace Advect
{
namespace Limiter
{
class VanLeer
{
public:
    static constexpr const char* name = "vanleer";

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar operator() (Set::Scalar lo, Set::Scalar mid, Set::Scalar hi) const
    {
        Set::Scalar dl = mid - lo;
        Set::Scalar dr = hi - mid;
        if (dl * dr <= Set::Scalar(0.0)) return Set::Scalar(0.0);
        Set::Scalar r = dl / dr;
        Set::Scalar phi = (r + Util::Abs(r)) / (Set::Scalar(1.0) + Util::Abs(r));
        return phi * dr;
    }

    void Parse(VanLeer & /*value*/, IO::ParmParse &/*pp*/) {}
};
}
}
}

#endif
