#include "IO/ParmParse.H"
#include "Integrator/AllenCahn.H"
#include "Integrator/Dendrite.H"
#include "Integrator/Flame.H"
#include "Util/Util.H"

#if AMREX_SPACEDIM == 2
#include "Integrator/SFI.H"
#endif

int
main(int argc, char *argv[])
{
    Util::Initialize(argc, argv);

#if AMREX_SPACEDIM == 2
    IO::ParmParse pp;
    std::string program;
    pp.query_validate("alamo.program", program, { "allencahn", "dendrite", "flame" });
    srand(2);

    Integrator::Integrator *integrator = nullptr;

    if (program == "allencahn")
        pp.select_only<Integrator::SFI<Integrator::AllenCahn>>(integrator);
    else if (program == "dendrite")
        pp.select_only<Integrator::SFI<Integrator::Dendrite>>(integrator);
    else if (program == "flame")
        pp.select_only<Integrator::SFI<Integrator::Flame>>(integrator);
    else
    {
        Util::Abort(INFO, "Invalid program option: " + program);
        return 1; // This line won't execute, but it tells the compiler the following lines won't execute.
    }


    integrator->InitData();
    integrator->Evolve();

    delete integrator;
#else

    Util::Abort(INFO, "This integrator works in 2D only");

#endif

    Util::Finalize();
}
