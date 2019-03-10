/* example.i */
%module alamo
%{
	extern int test ();
#include "Test/Operator/Elastic.H"
#include "Util/Util.H"
	%}

extern int test ();
 
namespace Util
{
	void Initialize ();
}


%rename("Test_Operator_Elastic") Test::Operator::Elastic;

namespace Test
{
	namespace Operator
	{
		class List {
		public:
			List();
			~List();
			int search(char *item);
			void insert(char *item);
			void remove(char *item);
			char *get(int n);
			int length;
		};
		class Elastic{
		public:
			enum Grid;
			Elastic();
			~Elastic();
			int hello(int i);
			void Define(int _ncells, int _nlevels, int _dim = AMREX_SPACEDIM, Grid config = Grid::XYZ);
			int TrigTest(int verbose, int component = 0, int n = 1, std::string plotfile = "");


		};
	}
}
