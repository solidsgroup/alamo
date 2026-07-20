# GPU-008: Stage storage before its lifetime ends and fence the stream
Status: draft
Class: correctness
Recognizer: manual: Host-owned data feeds a device launch without device staging, or asynchronous work can outlive local GPU-backed storage.
Applies: A host buffer or DeviceVector is staged/copied for asynchronous GPU work, or a local operator/MultiFab is destroyed before that work completes.
Transform:
  Before:
    A kernel captures host storage directly, or a local DeviceVector/operator/MultiFab reaches destruction while asynchronous work still references it.
  After:
    Stage host buffers into retained `DeviceVector`/POD `GpuArray` storage and synchronize the copy before launch; separately fence outstanding work with `amrex::Gpu::streamSynchronizeAll()` at a local-storage lifetime boundary.
Constraints: Mandatory when recognizer matches. Fence at the lifetime boundary, not arbitrarily in inner loops; do not hide a missing ownership fix with repeated global syncs.
Verify: `TIERS=2 benchmark/local_a100_gate.sh`; then repeat the strict multi-box elastic case five times. Expect sanitizer `ERROR SUMMARY: 0 errors`, no CUDA-700/UAF, and stable golden output.
Failure modes: Any mismatch or diagnostic is a failed conversion. Intermittent illegal address, host segfault, timing-sensitive output, or a device flag missed on a non-current MFIter stream indicates an ownership/fence defect. Single-box success is insufficient.
Evidence: commits f296eae5d31d39d2bf0bbfefdfddc5ea8142b0da, 54a941433b7582578cb5d56db794b1d648fb03cc, dd4054056aa5ef24bde948ac1b04013b5b106eb6; `docs/llm/BUG_PATTERNS.md`; `docs/llm/changelog/2026-07-02-gpu-audit.md`
