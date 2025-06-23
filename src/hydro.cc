//
// This initializes and runs the Alamo hydrodynamic solver implemented in the
// :ref:`Integrator::Hydro` integrator.
//

#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "IO/ParmParse.H"
#include "IO/FileNameParse.H"
#include "IO/WriteMetaData.H"
#include "AMReX_ParmParse.H"

#if AMREX_SPACEDIM==2
#include "Integrator/Hydro.H"
#endif

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    #if AMREX_SPACEDIM==2
    IO::ParmParse pp;
    srand(2);

    Integrator::Integrator *integrator = nullptr;
    pp.select_only<Integrator::Hydro>(integrator);

    integrator->InitData();
    integrator->Evolve();
    delete integrator;
    #else

    Util::Abort(INFO,"hydro currently works only in 2d");

    #endif

    Util::Finalize();
} 
