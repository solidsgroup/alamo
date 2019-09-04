#include <iostream>
#include <fstream>
#include <iomanip>

#include "AMReX.H"
#include <AMReX_AmrCore.H>

#include "Util/Util.H"
#include "Integrator/Mobility.H"


int main (int argc, char* argv[])
{
	Util::Initialize(argc,argv);

	Integrator::Mobility mobility;
	mobility.InitData();
	mobility.Evolve();

	Util::Finalize();
} 
