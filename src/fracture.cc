#include <iostream>
#include <fstream>
#include <iomanip>
#include <omp.h>

#include "Util/Util.H"
#include "Integrator/CahnHilliard.H"
#include "Integrator/Fracture.H"
#include "IO/FileNameParse.H"
#include "IO/WriteMetaData.H"

int main (int argc, char* argv[])
{
	//omp_set_num_threads(1);
	Util::Initialize(argc,argv);

	srand(1);
	Integrator::Integrator *fracturemodel =
		//new Integrator::PhaseFieldMicrostructure();
		new Integrator::Fracture();
	fracturemodel->InitData();
	fracturemodel->Evolve();
	delete fracturemodel;

	Util::Finalize();
} 
