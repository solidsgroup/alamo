#ifndef NUMERIC_ADVECT_LIMITER_VANALBADA_H
#define NUMERIC_ADVECT_LIMITER_VANALBADA_H

#include "IO/ParmParse.H"
#include "Util/Util.H"

namespace Numeric
{
namespace Advect
{
namespace Limiter
{
class VanAlbada
{
public:
    static constexpr const char* name = "vanalbada";

    AMREX_FORCE_INLINE AMREX_GPU_HOST_DEVICE
    Set::Scalar operator() (Set::Scalar lo, Set::Scalar mid, Set::Scalar hi) const
    {
        Set::Scalar dl = mid - lo;
        Set::Scalar dr = hi - mid;
        if (dl * dr <= Set::Scalar(0.0)) return Set::Scalar(0.0);
        Set::Scalar r = dl / dr;
        Set::Scalar phi = (r * r + r) / (r * r + Set::Scalar(1.0));
        return phi * dr;
    }

    void Parse(VanAlbada & /*value*/, IO::ParmParse &/*pp*/) {}
};
}
}
}

#endif
