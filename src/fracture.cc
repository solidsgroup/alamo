// This is the main entry point for alamo and is a general-purpose launcher for
// many of the main integrators.
// Check the possible values for :code:`alamo.program` below to see the possible
// integrators that can be launched.
//

#include <iostream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "IO/ParmParse.H"
#include "Integrator/Fracture.H"
#include "Integrator/Fracture_PFCZM.H"

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    std::string program;
    IO::ParmParse pp;
    // This input determines which integrator is used.
    pp.query_validate(  "alamo.program", program,
                        {"fracture","fracture_pfczm"});
    srand(2);

    Integrator::Integrator *integrator = nullptr;
    if       (program == "fracture")             pp.select_only<Integrator::Fracture>(integrator);
    else if (program == "fracture_pfczm")       pp.select_only<Integrator::Fracture_PFCZM>(integrator);
    else Util::Abort(INFO,"Error: \"",program,"\" is not a valid program.");

    integrator->InitData();
    integrator->Evolve();
    delete integrator;

    Util::Finalize();
} 
