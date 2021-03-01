#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "IO/FileNameParse.H"
#include "IO/WriteMetaData.H"
#include "AMReX_ParmParse.H"

#include "Model/Solid/Affine/Isotropic.H"
#include "Model/Solid/Elastic/NeoHookean.H"

#include "Integrator/CahnHilliard.H"
#include "Integrator/PhaseFieldMicrostructure.H"
#include "Integrator/TensionTest.H"
#include "Integrator/FiniteKinematics.H"
#include "Integrator/Flame.H"
#include "Integrator/PolymerDegradation.H"
#include "Integrator/HeatConduction.H"
#include "Integrator/Fracture.H"

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
	else if (program == "eshelby")
	{
		Integrator::Integrator *eshelby = new Integrator::TensionTest<Model::Solid::Affine::Isotropic>();
		eshelby->InitData();
		eshelby->Evolve();		
		delete eshelby;
	}
	else if (program == "finitekinematics")
	{
		//Integrator::Integrator *fk = new Integrator::FiniteKinematics();
		Integrator::Integrator *fk = new Integrator::TensionTest<Model::Solid::Elastic::NeoHookean>();
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
		//delete model;
	}
	else if (program == "fracture")
	{
		srand(1.0*amrex::ParallelDescriptor::MyProc());
		Integrator::Fracture model;
		model.InitData();
		model.Evolve();
		//delete model;
	}
	else if (program == "trigtest")
	{
		Test::Operator::Elastic test;
		test.Define(32,1);
		test.TrigTest(0,0,1,Util::GetFileName());
	}
	else
	{
		Util::Abort(INFO,"Error: \"",program,"\" is not a valid program.");
	}

	Util::Finalize();
} 
