#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "IO/FileNameParse.H"
#include "IO/WriteMetaData.H"
#include "AMReX_ParmParse.H"

#include "Integrator/CahnHilliard.H"
#include "Integrator/PhaseFieldMicrostructure.H"
#include "Integrator/Mobility.H"
#include "Integrator/Eshelby.H"
#include "Integrator/EshelbyPlastic.H"
#include "Integrator/FiniteKinematics.H"
#include "Integrator/Flame.H"
#include "Integrator/PolymerDegradation.H"
#include "Integrator/HeatConduction.H"
#include "Model/Solid/Elastic/Elastic.H"

int main (int argc, char* argv[])
{
	Util::Initialize(argc,argv);

	std::string program = "microstructure";
	amrex::ParmParse pp("alamo");
	pp.query("program",program);

	if (program == "microstructure")
	{
		srand(2);
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
		Integrator::Integrator *eshelby = new Integrator::Eshelby();
		eshelby->InitData();
		eshelby->Evolve();		
		delete eshelby;
	}
	else if (program == "crystalPlastic")
	{
		//Integrator::Integrator *cp = new Integrator::EshelbyPlastic();
		//cp->InitData();
		//cp->Evolve();		
		//delete cp;
	}
	else if (program == "finitekinematics")
	{
		Integrator::Integrator *fk = new Integrator::FiniteKinematics();
		fk->InitData();
		fk->Evolve();		
		delete fk;
	}
	else if (program == "flame")
	{
		Integrator::Integrator *flame = new Integrator::Flame();
		flame->InitData();
		flame->Evolve();
		delete flame;
	}
	else if (program == "heat")
	{
		Integrator::Integrator *heatconduction = new Integrator::HeatConduction();
		heatconduction->InitData();
		heatconduction->Evolve();
		delete heatconduction;
	}
	else if (program == "degradation")
	{
		srand(1.0*amrex::ParallelDescriptor::MyProc());
		Integrator::PolymerDegradation model;
		model.InitData();
		model.Evolve();
	}
	else
	{
		Util::Abort(INFO,"Error: \"",program,"\" is not a valid program.");
	}

	Util::Finalize();
} 
