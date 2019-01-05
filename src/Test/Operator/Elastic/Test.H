#ifndef TEST_OPERATOR_ELASTIC_TEST_H
#define TEST_OPERATOR_ELASTIC_TEST_H

namespace Test
{
namespace Operator
{
namespace Elastic
{
class Test
{
	Test() {};
	void Define(int _ncells,
		    int _nlevels);
	int TrigTest(bool verbose,
		     int component, int n);
private:


	int nlevels = 1;
	int ncells = 16;
	const int ref_ratio = 2;
 	amrex::Vector<amrex::Geometry> 			geom;
 	amrex::Vector<amrex::BoxArray> 			cgrids, ngrids;
 	amrex::Vector<amrex::DistributionMapping>       dmap;
 	amrex::Vector<amrex::MultiFab>                  u;
 	amrex::Vector<amrex::MultiFab>                  res;
 	amrex::Vector<amrex::MultiFab>                  rhs;
 	amrex::Vector<amrex::MultiFab>                  exact;

};
}
}
}




#endif


