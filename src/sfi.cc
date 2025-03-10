#include "Util/Util.H"
#include "IO/ParmParse.H"
#include "Integrator/Dendrite.H"
#include "Integrator/AllenCahn.H"

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

    if      (program == "allencahn") pp.select_only<Integrator::SFI<Integrator::AllenCahn>>(integrator);
    else if (program == "dendrite")  pp.select_only<Integrator::SFI<Integrator::Dendrite>>(integrator);
    
    integrator->InitData();
    integrator->Evolve();

    delete integrator;
    #else

    Util::Abort(INFO,"This integrator works in 2D only");

    #endif

    
    Util::Finalize();
} 
