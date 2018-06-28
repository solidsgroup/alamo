#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "Integrator/PhaseFieldMicrostructure/PhaseFieldMicrostructure.H"
#include "IO/FileNameParse.H"
#include "IO/WriteMetaData.H"

#if AMREX_SPACEDIM == 2

int main (int argc, char* argv[])
{
	Util::Initialize(argc,argv);

	srand(1);
	Integrator::Integrator *model =
		new Integrator::PhaseFieldMicrostructure();
	model->InitData();
	model->Evolve();
	delete model;

	Util::Finalize();
} 
#else
int main()
{
	std::cout << "This program works in 2D only, but AMREX_SPACEDIM=" << AMREX_SPACEDIM << std::endl;
}
#endif
