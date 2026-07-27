#ifndef SOLVER_NONLOCAL_MLMGSYNCSTATEGUARD
#define SOLVER_NONLOCAL_MLMGSYNCSTATEGUARD

#include <AMReX_GpuControl.H>

namespace Solver
{
namespace Nonlocal
{

/// \brief Restore AMReX's GPU synchronization-region state if MLMG throws.
///
/// `MLMG::setNoGpuSync(true)` enters AMReX single-stream and no-sync regions
/// inside `MLMG::solve()`. AMReX restores them on normal return, but its manual
/// restoration is bypassed by an exception. Alamo supports recoverable MLMG
/// exceptions, so preserve the caller's state with an outer RAII guard.
class MLMGSyncStateGuard
{
public:
    explicit MLMGSyncStateGuard(bool enabled)
        : m_enabled(enabled)
    {
        if (m_enabled)
        {
            m_prev_single_stream = amrex::Gpu::inSingleStreamRegion();
            m_prev_no_sync = amrex::Gpu::inNoSyncRegion();
        }
    }
    ~MLMGSyncStateGuard()
    {
        if (m_enabled)
        {
            (void)amrex::Gpu::setSingleStreamRegion(m_prev_single_stream);
            (void)amrex::Gpu::setNoSyncRegion(m_prev_no_sync);
        }
    }
    MLMGSyncStateGuard(const MLMGSyncStateGuard&) = delete;
    MLMGSyncStateGuard& operator=(const MLMGSyncStateGuard&) = delete;

private:
    bool m_enabled = false;
    bool m_prev_single_stream = false;
    bool m_prev_no_sync = false;
};

}
}
#endif
