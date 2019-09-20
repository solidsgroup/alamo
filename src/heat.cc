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

	Integrator::Integrator *heatconduction = new Integrator::HeatConduction();
	heatconduction->InitData();
	heatconduction->Evolve();
	delete heatconduction;

	Util::Finalize();
	return 0;
} 
