# GPU-012: Report kernel errors through device flags
Status: draft
Class: correctness
Recognizer: regex: `AMREX_GPU_DEVICE[\s\S]{0,800}(?:Util::Abort|Util::Message)\s*\(`
Applies: A device kernel detects non-finite data but attempts host abort/logging directly.
Transform:
  Before:
    Device body calls host-only `Util::Abort` or emits diagnostics on NaN/Inf.
  After:
    Set an `Util::DeviceErrorFlag`/device integer in the kernel, synchronize, then call `Util::AbortIfDeviceError` on host.
Constraints: Mandatory when recognizer matches. Keep flag writes race-safe and check every relevant kernel; do not suppress or clamp the bad value.
Verify: `make -j4`; compare CPU/GPU outputs at identical timestep and dimensions; inject NaN and expect deterministic host diagnostic after the launch.
Failure modes: Any mismatch or compile diagnostic is a failed conversion. nvcc host/device error, silent NaN propagation, or abort before asynchronous writes complete. Check and clear the flag at a defined host boundary so one stale error does not poison later steps. Validation must include the smallest representative input and a CPU reference; do not waive a failure as numerical noise.
Evidence: commits 54a941433b7582578cb5d56db794b1d648fb03cc
