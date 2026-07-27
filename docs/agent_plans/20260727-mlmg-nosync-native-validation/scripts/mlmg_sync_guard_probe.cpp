#include "Solver/Nonlocal/MLMGSyncStateGuard.H"

#include <AMReX.H>

#include <iostream>
#include <stdexcept>

namespace
{
bool check_state(bool single_stream, bool no_sync)
{
    return amrex::Gpu::inSingleStreamRegion() == single_stream
        && amrex::Gpu::inNoSyncRegion() == no_sync;
}

void force_exception(bool inner_single_stream, bool inner_no_sync)
{
    Solver::Nonlocal::MLMGSyncStateGuard guard(true);
    (void)amrex::Gpu::setSingleStreamRegion(inner_single_stream);
    (void)amrex::Gpu::setNoSyncRegion(inner_no_sync);
    throw std::runtime_error("intentional synchronization-state probe");
}
}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    bool ok = true;
    const bool original_single_stream = amrex::Gpu::inSingleStreamRegion();
    const bool original_no_sync = amrex::Gpu::inNoSyncRegion();

    (void)amrex::Gpu::setSingleStreamRegion(false);
    (void)amrex::Gpu::setNoSyncRegion(false);
    try
    {
        force_exception(true, true);
    }
    catch (const std::runtime_error&)
    {
        ok = ok && check_state(false, false);
    }

    (void)amrex::Gpu::setSingleStreamRegion(true);
    (void)amrex::Gpu::setNoSyncRegion(true);
    try
    {
        force_exception(false, false);
    }
    catch (const std::runtime_error&)
    {
        ok = ok && check_state(true, true);
    }

    (void)amrex::Gpu::setSingleStreamRegion(original_single_stream);
    (void)amrex::Gpu::setNoSyncRegion(original_no_sync);
    amrex::Finalize();

    std::cout << (ok ? "PASS" : "FAIL")
              << ": MLMGSyncStateGuard restored both outer states after exceptions\n";
    return ok ? 0 : 1;
}
