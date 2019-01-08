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
public:
	Test() {};
	void Define(int _ncells,
		    int _nlevels);
	int TrigTest(bool verbose,
		     int component = 0,
		     int n = 1,
		     std::string plotfile = "");
private:


	int nlevels = 1;
	int ncells = 16;
	const int ref_ratio = 2;
 	amrex::Vector<amrex::Geometry> 		   geom;
 	amrex::Vector<amrex::BoxArray> 		   cgrids, ngrids;
 	amrex::Vector<amrex::DistributionMapping>  dmap;
	
 	amrex::Vector<amrex::MultiFab>             solution_exact;
 	amrex::Vector<amrex::MultiFab>             solution_numeric;
 	amrex::Vector<amrex::MultiFab>             solution_error;

 	amrex::Vector<amrex::MultiFab>             rhs_prescribed;
 	amrex::Vector<amrex::MultiFab>             rhs_exact;
 	amrex::Vector<amrex::MultiFab>             rhs_numeric;

 	amrex::Vector<amrex::MultiFab>             res_exact;
 	amrex::Vector<amrex::MultiFab>             res_numeric;

};
}
}
}




#endif


