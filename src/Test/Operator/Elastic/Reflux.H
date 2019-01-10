#ifndef OPERATOR_TEST_H
#define OPERATOR_TEST_H

namespace Test
{
namespace Operator
{
namespace Elastic
{
template<class Op>
class Reflux
{
public:
	int RefluxTest(int)
	{
		Util::Abort(INFO,"No general template for this test");
		return std::numeric_limits<int>::max();
	};
	int SpatiallyVaryingCTest(int)
	{
		Util::Abort(INFO,"Spatially Varying C test is not implemented yet");
		return std::numeric_limits<int>::max();	
	};
	
private:

};
}
}
}





#endif


