#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "IO/ParmParse.H"
#include "IO/FileNameParse.H"
#include "IO/WriteMetaData.H"
#include "AMReX_ParmParse.H"

#include "Model/Solid/Affine/Cubic.H"
#include "Integrator/ThermoMechanics.H"
#include "Integrator/PhaseFieldMicrostructure.H"


int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    std::string program = "thermomechanics";
    IO::ParmParse pp;
    pp_query("alamo.program",program);
    srand(2);

    Integrator::Integrator *integrator = nullptr;
    pp.select_only<Integrator::ThermoMechanics>(integrator);
    integrator->InitData();
    integrator->Evolve();
    
    delete integrator;
    
    Util::Finalize();
} 

