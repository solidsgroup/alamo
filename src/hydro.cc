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
    std::string program = "hydro";
    IO::ParmParse pp;
    pp.query("alamo.program",program);
    srand(2);

    Integrator::Integrator *integrator = nullptr;
    if (program == "hydro")                integrator = new Integrator::Hydro(pp);
    else Util::Abort(INFO,"Error: \"",program,"\" is not a valid program.");

    integrator->InitData();
    integrator->Evolve();
    delete integrator;
    #else

    Util::Abort(INFO,"hydro currently works only in 2d");

    #endif

    Util::Finalize();
} 
