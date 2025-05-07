#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "IO/ParmParse.H"
#include "IO/FileNameParse.H"
#include "IO/WriteMetaData.H"
#include "AMReX_ParmParse.H"

#include "Model/Solid/Linear/Isotropic.H"
#include "Integrator/TopOp.H"

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    std::string program = "microstructure";
    IO::ParmParse pp;
    srand(2);

    Integrator::Integrator *integrator = nullptr;
    pp.select_only<Integrator::TopOp<Model::Solid::Linear::Isotropic>>(integrator);
    integrator->InitData();
    integrator->Evolve();
    delete integrator;

    Util::Finalize();
} 
