#ifndef NUMERIC_ADVECT_LIMITER_KOREN_H
#define NUMERIC_ADVECT_LIMITER_KOREN_H

#include "IO/ParmParse.H"
#include "Util/Util.H"

namespace Numeric
{
namespace Advect
{
namespace Limiter
{
class Koren
{
public:
    static constexpr const char* name = "koren";

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar operator() (Set::Scalar lo, Set::Scalar mid, Set::Scalar hi) const
    {
        Set::Scalar dl = mid - lo;
        Set::Scalar dr = hi - mid;
        if (dl * dr <= Set::Scalar(0.0)) return Set::Scalar(0.0);
        Set::Scalar r = dl / dr;
        Set::Scalar phi = Util::Min(
            Util::Min(Set::Scalar(2.0) * r, (Set::Scalar(1.0) + Set::Scalar(2.0) * r) / Set::Scalar(3.0)),
            Set::Scalar(2.0));
        return Util::Max(Set::Scalar(0.0), phi) * dr;
    }

    void Parse(Koren & /*value*/, IO::ParmParse &/*pp*/) {}
};
}
}
}

#endif
