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
Verify: `make -j4`; compare CPU/GPU outputs at identical timestep and dimensions; run repeated elastic solves; expected no CUDA-700/UAF and stable output.
Failure modes: Any mismatch or compile diagnostic is a failed conversion. intermittent illegal address, host segfault, or results changing with stream timing. Validation must include the smallest representative input and a CPU reference; do not waive a failure as numerical noise.
Evidence: commits f296eae5d31d39d2bf0bbfefdfddc5ea8142b0da, 54a941433b7582578cb5d56db794b1d648fb03cc, dd4054056aa5ef24bde948ac1b04013b5b106eb6
