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
Verify: `benchmark/ci_golden_compare.sh && TIERS=2 benchmark/local_a100_gate.sh`; force a multi-box MLMG case and repeat it. Expect identical interpolation values and no sanitizer/UAF failure.
Failure modes: Any mismatch or diagnostic is a failed conversion. MFIter's stream pool can make this a multi-box-only silent interpolation corruption or intermittent invalid access; synchronization may mask the bug but does not replace ownership.
Evidence: commit c00f69086c6d1dc67bdf71e2bd174dbdcc953c85; `benchmark/archive/elixir_race_audit.md`; `docs/llm/BUG_PATTERNS.md`
