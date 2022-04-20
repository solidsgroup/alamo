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
#include "Model/Solid/Linear/Laplacian.H"

#include "Integrator/CahnHilliard.H"
#include "Integrator/PhaseFieldMicrostructure.H"
#include "Integrator/Mechanics.H"
#include "Integrator/Flame.H"
#include "Integrator/PolymerDegradation.H"
#include "Integrator/HeatConduction.H"
#include "Integrator/Fracture.H"
#include "Integrator/ThermoElastic.H"
#include "Integrator/Dynamics.H"

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    std::string program = "microstructure";
    IO::ParmParse pp;
    pp.query("alamo.program",program);

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
        std::string model = "linear.isotropic";
        pp.query("alamo.program.mechanics.model",model);
        if (model == "linear.isotropic") 
        {
            integrator = new Integrator::Mechanics<Model::Solid::Linear::Isotropic>();
            pp.queryclass(dynamic_cast<Integrator::Mechanics<Model::Solid::Linear::Isotropic>*>(integrator));
        }
        else if (model == "affine.isotropic") 
        {
            integrator = new Integrator::Mechanics<Model::Solid::Affine::Isotropic>();
            pp.queryclass(dynamic_cast<Integrator::Mechanics<Model::Solid::Affine::Isotropic>*>(integrator));
        }
        else if (model == "linear.laplacian") 
        {
            integrator = new Integrator::Mechanics<Model::Solid::Linear::Laplacian>();
            pp.queryclass(dynamic_cast<Integrator::Mechanics<Model::Solid::Linear::Laplacian>*>(integrator));
        }
        else
        {
            Util::Abort(INFO,model," is not a valid model");
        }

        integrator->InitData();
        integrator->Evolve();
        delete integrator;
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
    else if (program == "thermoelastic")
    {
        IO::ParmParse pp;
        Integrator::ThermoElastic te;
        pp.queryclass(te);
        te.InitData();
        te.Evolve();
    }
    else if (program == "dynamics")
    {
        IO::ParmParse pp;
        Integrator::Dynamics integrator;
        pp.queryclass(integrator);
        integrator.InitData();
        integrator.Evolve();
    }
    else
    {
        Util::Abort(INFO,"Error: \"",program,"\" is not a valid program.");
    }

    Util::Finalize();
} 
