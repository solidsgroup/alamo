#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include <AMReX_AmrCore.H>

#include "Util/Util.H"
#include "Integrator/HeatConduction.H"


int main (int argc, char* argv[])
{
	Util::Initialize(argc,argv);

	Integrator::HeatConduction heatconduction;
	heatconduction.InitData();
	heatconduction.Evolve();

	Util::Finalize();
	return 0;
} 
