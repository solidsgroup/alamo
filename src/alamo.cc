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
#include "Integrator/Eshelby.H"

int main (int argc, char* argv[])
{
	Util::Initialize(argc,argv);

	std::string program = "microstructure";
	amrex::ParmParse pp("alamo");
	pp.query("program",program);

	if (program == "microstructure")
	{
		srand(1);
		Integrator::Integrator *pfm = new Integrator::PhaseFieldMicrostructure();
		//Integrator::PhaseFieldMicrostructure pfm;
		pfm->InitData();
		pfm->Evolve();
		delete pfm;
	}
	else if (program == "mobility")
	{	
		Integrator::Mobility mobility;
		mobility.InitData();
		mobility.Evolve();
	}
	else if (program == "eshelby")
	{
		Integrator::Eshelby eshelby;
		eshelby.InitData();
		eshelby.Evolve();		
	}
	else
	{
		Util::Abort(INFO,"Error: ",program," is not a valid program.");
	}

	Util::Finalize();
} 
