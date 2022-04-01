#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "IO/ParmParse.H"
#include "IO/FileNameParse.H"
#include "IO/WriteMetaData.H"
#include "AMReX_ParmParse.H"

#include "Model/Solid/Affine/Isotropic.H"
#include "Model/Solid/Elastic/NeoHookean.H"

#include "Integrator/CahnHilliard.H"
#include "Integrator/PhaseFieldMicrostructure.H"
#include "Integrator/Mechanics.H"
#include "Integrator/FiniteKinematics.H"
#include "Integrator/Flame.H"
#include "Integrator/PolymerDegradation.H"
#include "Integrator/HeatConduction.H"
#include "Integrator/Fracture.H"
#include "Integrator/ThermoElastic.H"

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    std::string program = "microstructure";
    IO::ParmParse pp("alamo");
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
    else if (program == "mechanics")
    {
        Integrator::Integrator *integrator;
        std::string model = "linear.elastic";
        pp.query("program.mechanics.model",model);
        if (model == "linear.elastic") 
        {
            integrator = new Integrator::Mechanics<Model::Solid::Linear::Isotropic>();
            pp.queryclass(dynamic_cast<Integrator::Mechanics<Model::Solid::Linear::Isotropic>*>(integrator));
        }
    }
    else if (program == "eshelby")
    {
        IO::ParmParse pp;
        Integrator::Mechanics<Model::Solid::Affine::Isotropic> eshelby;
        pp.queryclass(eshelby);
        eshelby.InitData();
        eshelby.Evolve();        
    }
    else if (program == "finitekinematics")
    {
        //Integrator::Integrator *fk = new Integrator::FiniteKinematics();
        Integrator::Integrator *fk = new Integrator::Mechanics<Model::Solid::Elastic::NeoHookean>();
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
        IO::ParmParse pp;
        Integrator::HeatConduction heatconduction;
        pp.queryclass(heatconduction);
        heatconduction.InitData();
        heatconduction.Evolve();
    }
    else if (program == "degradation")
    {
        srand(1.0*amrex::ParallelDescriptor::MyProc());
        Integrator::PolymerDegradation model;
        model.InitData();
        model.Evolve();
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
    else if (program == "thermoelastic")
    {
        IO::ParmParse pp;
        Integrator::ThermoElastic te;
        pp.queryclass(te);
        te.InitData();
        te.Evolve();
    }
    else
    {
        Util::Abort(INFO,"Error: \"",program,"\" is not a valid program.");
    }

    Util::Finalize();
} 
