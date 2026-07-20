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
Verify: `DIM=2 CUDA_FP=strict bash benchmark/build_alamo_local_gpu.sh && benchmark/golden_compare_flame.sh input 1`; expect CUDA compilation and finite output matching the CPU golden case.
Failure modes: Any mismatch or diagnostic is a failed conversion. nvcc can reject `this`/inaccessible captured types; a real GPU can fault on host pointers even when HMM masks the defect locally. Validate the smallest representative input and a CPU reference; review capture ownership and dimensional assumptions explicitly.
Evidence: commit d522e1ac08729b306a215ad896d8a304983d55de; `docs/gpu_device_capture_conventions.md`; `docs/llm/BUG_PATTERNS.md`
