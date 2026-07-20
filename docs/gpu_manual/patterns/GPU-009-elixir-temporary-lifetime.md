# GPU-009: Keep GPU temporaries alive through the launch
Status: draft
Class: correctness
Recognizer: manual: A local FArrayBox temporary is consumed asynchronously without a retained Elixir.
Applies: A GPU launch uses temporary Fab data without retaining an AMReX Elixir for its lifetime.
Transform:
  Before:
    Device path comments out or omits `amrex::Elixir` ownership.
  After:
    Retain `amrex::Gpu::Elixir tmpfab_eli = tmpfab.elixir()` when running on device in a launch region.
Constraints: Mandatory when recognizer matches. Retain the Elixir until all dependent asynchronous work completes; do not create one for CPU-only execution. The Elixir belongs to the temporary Fab lifetime, not merely the lexical launch statement.
Verify: `make -j4`; compare CPU/GPU outputs at identical timestep and dimensions; expected interpolator run has no use-after-free and identical CPU/GPU values.
Failure modes: Any mismatch or compile diagnostic is a failed conversion. intermittent invalid memory access or corrupted interpolation when temporary storage is reclaimed early. Synchronization may mask the bug but does not replace ownership. Validation must include the smallest representative input and a CPU reference; do not waive a failure as numerical noise.
Evidence: commits c00f69086c6d1dc67bdf71e2bd174dbdcc953c85
