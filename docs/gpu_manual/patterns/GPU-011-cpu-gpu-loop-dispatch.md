# GPU-011: Dispatch loops explicitly between CPU and GPU
Status: draft
Class: correctness
Recognizer: manual: A data-parallel host loop or backend-dependent loop lacks an explicit CPU/GPU dispatch boundary.
Applies: A loop has separate CPU/GPU launch behavior or an unconditional host loop.
Transform:
  Before:
    Host-only `for` loop or launch that cannot select the execution backend.
  After:
    Use AMReX launch-region/ParallelFor dispatch with a device-callable body and preserve a CPU-safe path.
Constraints: Mandatory when the recognizer matches and the loop is data-parallel. Do not move reductions, I/O, or ordering-dependent logic without proving semantics; preserve a tested CPU path.
Verify: `make -j4`; compare CPU/GPU outputs at identical timestep and dimensions; expected CPU and GPU runs complete with matching fields; GPU profiler shows kernel dispatch where applicable.
Failure modes: Any mismatch or compile diagnostic is a failed conversion. race, divergent reduction, or CPU/GPU numerical mismatch. An apparently faster launch that changes ordering or omits the CPU branch is incorrect. Validation must include the smallest representative input and a CPU reference; do not waive a failure as numerical noise.
Evidence: commits d522e1ac08729b306a215ad896d8a304983d55de
