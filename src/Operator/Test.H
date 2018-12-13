#ifndef OPERATOR_TEST_H
#define OPERATOR_TEST_H

namespace Operator
{
template<class Op>
class Test
{
public:
	int RefluxTest(int)
	{
		Util::Abort(INFO,"No general template for this test");
		return std::numeric_limits<int>::max();
	};
	
private:

};
}





#endif


