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
Verify: `DIM=2 CUDA_FP=strict bash benchmark/build_alamo_local_gpu.sh`; inject NaN in a non-last box and run a multi-box case; expect one deterministic host diagnostic after all streams synchronize.
Failure modes: Any mismatch or diagnostic is a failed conversion. Host/device compile errors, silent NaN propagation, a false-negative from checking only the current MFIter stream, or a stale flag poisoning later steps all fail.
Evidence: commit 54a941433b7582578cb5d56db794b1d648fb03cc; `docs/llm/changelog/2026-07-02-gpu-audit.md`; `docs/llm/BUG_PATTERNS.md`
