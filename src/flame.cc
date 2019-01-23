#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "Integrator/Flame.H"

int main (int argc, char* argv[])
{

	Util::Initialize(argc,argv);

	Integrator::Integrator *flame =
		new Integrator::Flame();
	flame->InitData();
	flame->Evolve();
	delete flame;

	Util::Finalize();
} 
