#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "IO/ParmParse.H"
#include "IO/FileNameParse.H"
#include "IO/WriteMetaData.H"
#include "AMReX_ParmParse.H"

#include "Model/Solid/Affine/Cubic.H"
#include "Model/Solid/Affine/Hexagonal.H"

#include "Integrator/AllenCahn.H"
#include "Integrator/CahnHilliard.H"
#include "Integrator/PhaseFieldMicrostructure.H"
#include "Integrator/Mechanics.H"
#include "Integrator/Flame.H"
#include "Integrator/HeatConduction.H"
#include "Integrator/Fracture.H"
#include "Integrator/ThermoElastic.H"
#include "Integrator/Dendrite.H"
#include "Integrator/ChemReduction.H"

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    std::string program = "microstructure";
    IO::ParmParse pp;
    pp_query("alamo.program",program);
    srand(2);

    Integrator::Integrator *integrator = nullptr;
    if (program == "microstructure")
    {
        std::string model = "affine.cubic";
        pp_query("alamo.program.microstructure.model",model);
        if      (model == "affine.cubic")       integrator = new Integrator::PhaseFieldMicrostructure<Model::Solid::Affine::Cubic>(pp);
        else if (model == "affine.hexagonal")   integrator = new Integrator::PhaseFieldMicrostructure<Model::Solid::Affine::Hexagonal>(pp);
        else Util::Abort(INFO,model," is not a valid model");
    }
    else if (program == "flame")                integrator = new Integrator::Flame(pp);
    else if (program == "heat")                 integrator = new Integrator::HeatConduction(pp);
    else if (program == "thermoelastic")        integrator = new Integrator::ThermoElastic(pp);
    else if (program == "fracture")             integrator = new Integrator::Fracture();
    else if (program == "dendrite")             integrator = new Integrator::Dendrite(pp);
    else if (program == "allencahn")            integrator = new Integrator::AllenCahn(pp);
    else if (program == "chemreduction")        integrator = new Integrator::ChemReduction(pp);
    else Util::Abort(INFO,"Error: \"",program,"\" is not a valid program.");

    integrator->InitData();
    integrator->Evolve();
    delete integrator;

    Util::Finalize();
} 
