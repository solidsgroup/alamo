//
// This namespace has evolved to contain some general Interpolator routines that
// are actually unrelated.
// The :ref:`Numeric::Interpolator::Linear` is an I/O utility.
// The :ref:`Numeric::Interpolator::NodeBilinear` is a low-level AMReX override that
// the average user (or even coder) will probably not need to touch.
//

#ifndef NUMERIC_INTERPOLATOR_H_
#define NUMERIC_INTERPOLATOR_H_

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include "Set/Set.H"

namespace Numeric
{
namespace Interpolator
{
template<class T>
class Interpolator
{
public:
    Interpolator() {};
};
}
}
#endif
