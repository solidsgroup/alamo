#include "IO/ParmParse.H"
#include "Integrator/AllenCahn.H"
#include "Integrator/Dendrite.H"
#include "Integrator/Flame.H"
#include "Util/Util.H"

#if AMREX_SPACEDIM==2
#include "Integrator/SFI.H"
#endif

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    #if AMREX_SPACEDIM==2
    IO::ParmParse pp;
    std::string program;

    srand(2);

    Integrator::Integrator *integrator = nullptr;

    // Validate/make sure the correct Alamo program/Inetrgator is used
    pp.query_switch("alamo.program",{
            {"allencahn", [&](){
                pp.select_only<Integrator::SFI<Integrator::AllenCahn>>(integrator);
            }},
            {"dendrite", [&](){
                pp.select_only<Integrator::SFI<Integrator::Dendrite>>(integrator);
            }},
            {"flame", [&]() {
                pp.select_only<Integrator::SFI<Integrator::Flame>>(integrator);
            }}
        });

    integrator->InitData();
    integrator->Evolve();

    delete integrator;
    #else
    if (!IO::ParmParse::InTraversalMode())
        Util::Abort(INFO,"This integrator works in 2D only");
    #endif

    
    Util::Finalize();
} 
