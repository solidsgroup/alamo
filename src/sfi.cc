#include <iostream>
#include <fstream>
#include <iomanip>

#include "Util/Util.H"
#include "IO/ParmParse.H"
#include "IO/FileNameParse.H"
#include "IO/WriteMetaData.H"
#include "AMReX_ParmParse.H"

#if AMREX_SPACEDIM==2
#include "Integrator/SFI.H"
#endif

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    #if AMREX_SPACEDIM==2
    IO::ParmParse pp;
    std::string program;
    pp.query("alamo.program",program);
    srand(2);

    Integrator::SFI integrator(pp);
    integrator.InitData();
    integrator.Evolve();
    #else

    Util::Abort(INFO,"This integrator works in 2D only");

    #endif

    
    Util::Finalize();
} 
