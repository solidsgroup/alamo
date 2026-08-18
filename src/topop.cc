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
#include "Integrator/TopOpLowMach.H"

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    std::string program;
    IO::ParmParse pp;
    // which integrator to use
    pp.query_validate("alamo.program", program, {"topop", "topop_lowmach"});
    srand(2);

    Integrator::Integrator *integrator = nullptr;
    pp.query_switch("alamo.program", {
        {"topop", [&]() {
            pp.select_only<Integrator::TopOp<Model::Solid::Linear::Isotropic>>(integrator);
        }},
        {"topop_lowmach", [&]() {
            pp.select_only<Integrator::TopOpLowMach>(integrator);
        }}
    });
    integrator->InitData();
    integrator->Evolve();
    delete integrator;

    Util::Finalize();
} 
