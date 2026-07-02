//
// This initializes and runs the Alamo low-Mach compressible flow solver
// implemented in the :ref:`Integrator::LowMach` integrator.
//

#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "IO/ParmParse.H"
#include "IO/FileNameParse.H"
#include "IO/WriteMetaData.H"
#include "AMReX_ParmParse.H"

#include "Integrator/LowMach.H"

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    IO::ParmParse pp;
    srand(2);

    Integrator::Integrator *integrator = nullptr;
    pp.select_only<Integrator::LowMach>(integrator);

    integrator->InitData();
    integrator->Evolve();
    delete integrator;

    Util::Finalize();
}
