# GPU-006: Hoist kernel captures into device-safe locals
Status: draft
Class: correctness
Recognizer: regex: `AMREX_GPU_DEVICE[^\n]*\n[^\n]*this->`
Applies: A kernel body reads integrator members through `this` or host-owned pointers.
Transform:
  Before:
    `ParallelFor(..., [=] AMREX_GPU_DEVICE { use this->thermal.hc; })`
  After:
    Copy scalar/aggregate values before launch (`const auto thermal_hc = thermal.hc;`) and capture those locals by value.
Constraints: Mandatory when recognizer matches. Capture only trivially copyable/device-safe state; never capture allocators, polymorphic objects, or host references.
Verify: `make -j4`; compare CPU/GPU outputs at identical timestep and dimensions; expected GPU compile succeeds and a one-step run produces finite fields without illegal-address errors.
Failure modes: Any mismatch or compile diagnostic is a failed conversion. nvcc rejects `this`, device faults on host pointer, or stale member values are observed. Validation must include the smallest representative input and a CPU reference; do not waive a failure as numerical noise. Include launch-region behavior, ownership, and dimensional assumptions in review; the transform is complete only when those remain explicit.
Evidence: commits d522e1ac08729b306a215ad896d8a304983d55de
