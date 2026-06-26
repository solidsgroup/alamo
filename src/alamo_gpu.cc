// Minimal entry point for the chamber-gpu (CUDA) build. alamo.cc is a
// general-purpose launcher that compile-time includes every integrator
// (AllenCahn, CahnHilliard, Dendrite, PhaseFieldMicrostructure, ...); none of
// those are needed for chamber/Flame runs, and several have never been made
// nvcc-clean. This main only pulls in Integrator::Flame so the CUDA build's
// object closure stays limited to what Flame actually depends on.

#include <string>
#include <cstdlib>

#ifdef AMREX_USE_GPU
#include <cuda_runtime.h>
#include <AMReX_Print.H>
#endif

#include "Util/Util.H"
#include "IO/ParmParse.H"

#include "Integrator/Flame.H"

int main (int argc, char* argv[])
{
    Util::Initialize(argc,argv);

#ifdef AMREX_USE_GPU
    // DIAGNOSTIC PROBE: the 3D elastic prepareForSolve kernel calls a
    // non-inlined device function (NeoHookean DDW + Matrix4<3>) whose call
    // tree can exceed the default 1 KB device stack. Bump it to test whether
    // the "unspecified launch failure" at Newton.H:196 is a stack overflow.
    {
        std::size_t stk = 0;
        if (const char* e = std::getenv("ALAMO_GPU_STACK")) stk = std::strtoul(e, nullptr, 10);
        if (stk > 0)
        {
            cudaError_t rc = cudaDeviceSetLimit(cudaLimitStackSize, stk);
            amrex::Print() << "[probe] cudaDeviceSetLimit(stack=" << stk << ") rc=" << (int)rc << "\n";
        }
    }
#endif

    IO::ParmParse pp;
    srand(2);

    std::string program;
    pp.query_validate("alamo.program", program, {"flame"});

    Integrator::Integrator *integrator = nullptr;
    pp.select_only<Integrator::Flame>(integrator);

    integrator->InitData();
    integrator->Evolve();
    delete integrator;

    Util::Finalize();
}
