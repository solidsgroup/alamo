#ifndef SET_SET_
#define SET_SET_

#include <iomanip>

#include "Util/Util.H"

#include "Set/Base.H"
/// \brief A collection of data types and symmetry-reduced data structures


namespace Util
{
Set::Scalar Random();
}

namespace Set
{
namespace Constant
{
static const Set::Scalar Pi = 3.14159265359;
}
}

#include "Set/Matrix3.H"
#include "Set/Matrix4.H"


#endif
