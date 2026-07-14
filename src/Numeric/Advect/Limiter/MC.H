#ifndef NUMERIC_ADVECT_LIMITER_MC_H
#define NUMERIC_ADVECT_LIMITER_MC_H


#include "Util/Util.H"
#include "IO/ParmParse.H"

namespace Numeric
{
namespace Advect
{
namespace Limiter
{
class MC
{
public:
    static constexpr const char* name = "mc";
    Set::Scalar operator() (Set::Scalar lo, Set::Scalar mid, Set::Scalar hi) const
    {
        double dl = mid - lo;
        double dr = hi - mid;
        return Util::MinMod(double(0.5) * (dl + dr), Util::MinMod(double(2.0) * dl, double(2.0) * dr));
    }

    void Parse(MC & /*value*/, IO::ParmParse &/*pp*/) {} 
};
}
}
}



#endif
