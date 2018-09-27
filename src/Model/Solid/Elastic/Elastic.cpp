#include "Elastic.H"
#include "Set/Set.H"
#include "Util/Color.H"

namespace Model
{
namespace Solid
{
namespace Elastic
{
int Elastic::operator () (int i, int j, int k, int l)
{
	if (i==j && j==k && j==l) return 0;
	else if (i==j && k==l) return 1;
	else if (i==k && j==l) return 2;
	else if (i==l && j==k) return 2;
	else return -1;
}
Set::Scalar
Elastic::C(int i, int j, int k, int l)
{
	Set::Scalar ret = 0.0;
	if (i==j && k==l) ret += mu;
	if (i==l && j==k) ret += mu;
	if (i==j && k==l) ret += lambda;
	return ret;
}
	
int
Elastic::NComp()
{
	return 3;
}

}
}
}
