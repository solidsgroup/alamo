#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "Integrator/CahnHilliard.H"
#include "Integrator/PhaseFieldMicrostructure.H"
#include "IO/FileNameParse.H"
#include "IO/WriteMetaData.H"
#include "AMReX_ParmParse.H"

#include "Integrator/Mobility.H"

int main (int argc, char* argv[])
{
	Util::Initialize(argc,argv);

	std::string program = "microstructure";
	amrex::ParmParse pp("alamo");
	pp.query("program",program);

	if (program == "microstructure")
	{
		srand(1);
		Integrator::PhaseFieldMicrostructure pfm;
		pfm.InitData();
		pfm.Evolve();
	}
	else if (program == "mobility")
	{	
		Integrator::Mobility mobility;
		mobility.InitData();
		mobility.Evolve();
	}

	Util::Finalize();
} 
