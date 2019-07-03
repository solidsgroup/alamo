#include "Set.H"
namespace Util
{
Set::Scalar Random()
{
	return ((Set::Scalar) rand()) / ((Set::Scalar) RAND_MAX);
}
}
