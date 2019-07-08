#include "Set.H"
namespace Set
{
std::ostream&
operator<< (std::ostream& os, const Matrix4<3,Sym::Full>& b)
{
   	b.Print(os);
	return os;
}
}

namespace Util
{
Set::Scalar Random()
{
	return ((Set::Scalar) rand()) / ((Set::Scalar) RAND_MAX);
}
}
