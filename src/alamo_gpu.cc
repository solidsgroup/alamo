// Minimal entry point for the chamber-gpu (CUDA) build. alamo.cc is a
// general-purpose launcher that compile-time includes every integrator
// (AllenCahn, CahnHilliard, Dendrite, PhaseFieldMicrostructure, ...); none of
// those are needed for chamber/Flame runs, and several have never been made
// nvcc-clean. This main only pulls in Integrator::Flame so the CUDA build's
// object closure stays limited to what Flame actually depends on.

#include <string>

#include "Util/Util.H"
#include "IO/ParmParse.H"

#include "Integrator/Flame.H"

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

    IO::ParmParse pp;
    srand(2);

    std::string program;
    pp_query_validate("alamo.program", program, {"flame"}); // Program/integrator selector for the CUDA launcher

    Integrator::Integrator *integrator = nullptr;
    pp.select_only<Integrator::Flame>(integrator);

    integrator->InitData();
    integrator->Evolve();
    delete integrator;

    Util::Finalize();
}
