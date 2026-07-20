# GPU-013: Use ReduceOps for aggregate kernel results
Status: draft
Class: correctness
Recognizer: regex: `AMREX_GPU_DEVICE[\s\S]{0,800}\b(?:trac_hi|disp_hi|massflux|mdot|volume)\b[^;\n]*(?:\+=|=)`
Applies: A kernel accumulates into host/member state or performs several independent reductions.
Transform:
  Before:
    Device lambda writes `trac_hi[face] += ...`, `disp_hi[face] = ...`, or another host/member aggregate.
  After:
    Use `ReduceOps`/`ReduceData` for tuple reductions, then assign results on host; fused tuple reductions are a performance variant.
Constraints: Mandatory when the recognizer exposes unsafe device accumulation. Preserve reduction identities and ordering tolerance. Fusing otherwise-correct independent reductions is an optional performance variant requiring profiling justification and unchanged golden results; never mix incompatible operators.
Verify: `benchmark/ci_golden_compare.sh`; repeat a strict real-GPU multi-box run and compare aggregates to CPU. Expect deterministic values and no CUDA-700/race under decomposition changes.
Failure modes: Any mismatch or diagnostic is a failed conversion. Direct host writes can appear valid under HMM yet fault with CUDA error 700 on an A100; wrong identities/types yield race-dependent totals. Confirm dimensional and reduction identities explicitly.
Evidence: commits 76ae550e0388c96fd315724628942c2c5394426b, 02772a36997539b183c70ba195e469478b894007; `docs/llm/BUG_PATTERNS.md`; `docs/llm/changelog/2026-07-02-gpu-audit.md`
