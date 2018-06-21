#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "Integrator/Flame/Flame.H"

int main (int argc, char* argv[])
{

	Util::Initialize(argc,argv);

	Integrator::Flame flame;
	flame.InitData();
	flame.Evolve();
  
	Util::Finalize();
} 
