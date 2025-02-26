#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "IO/ParmParse.H"
#include "IO/FileNameParse.H"
#include "IO/WriteMetaData.H"
#include "Integrator/Dendrite.H"
#include "Integrator/AllenCahn.H"
#include "Integrator/Hydro.H"

#if AMREX_SPACEDIM==2
#include "Integrator/SFI.H"
#endif

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    #if AMREX_SPACEDIM==2
    IO::ParmParse pp;
    std::string program;
    pp.query_default("alamo.program",program,"allencahn");
    srand(2);

    Integrator::Integrator *integrator = nullptr;
    if (program == "allencahn")     integrator = new Integrator::SFI<Integrator::AllenCahn>(pp);
    else if (program == "dendrite") integrator = new Integrator::SFI<Integrator::Dendrite>(pp);
    else Util::Abort(INFO,"Invalid program ",program);
    
    integrator->InitData();
    integrator->Evolve();

    delete integrator;
    #else

    Util::Abort(INFO,"This integrator works in 2D only");

    #endif

    
    Util::Finalize();
} 
