//
// This initializes and runs the Alamo fracture solver implemented in the
// :ref:`Integrator::Fracture` integrator.
//

#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "IO/ParmParse.H"
#include "Integrator/Fracture.H"
#include "Integrator/DuctileFracture.H"

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    std::string program;
    IO::ParmParse pp;
    // This input determines which integrator is used.
    pp.query_validate(  "alamo.program", program,
                        {"fracture", "ductilefracture"});
    srand(2);

    Integrator::Integrator *integrator = nullptr;
    if       (program == "fracture")             pp.select_only<Integrator::Fracture>(integrator);
    if       (program == "ductilefracture")      pp.select_only<Integrator::DuctileFracture>(integrator);
    else Util::Abort(INFO,"Error: \"",program,"\" is not a valid program.");

    integrator->InitData();
    integrator->Evolve();
    delete integrator;

    Util::Finalize();
} 
