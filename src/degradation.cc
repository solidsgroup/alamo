#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include "Util/Util.H"
//#if AMREX_SPACEDIM == 1

#include "Integrator/PolymerDegradation.H"

int main (int argc, char* argv[])
{
	if (argc < 2)
	{
		std::cout << "Missing input file" << std::endl;
		exit(-1);
	}
	Util::Initialize(argc,argv);

	srand(1.0*amrex::ParallelDescriptor::MyProc());
	{
		Integrator::PolymerDegradation model;
		model.InitData();
		model.Evolve();
	}
  
	Util::Finalize();
} 

//#else
/*int main()
{
	std::cout << "This program works in 1D only, but AMREX_SPACEDIM=" << AMREX_SPACEDIM << std::endl;
}*/

//#endif
