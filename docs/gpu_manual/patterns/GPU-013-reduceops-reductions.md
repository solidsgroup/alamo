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
Verify: `make -j4`; compare CPU/GPU outputs at identical timestep and dimensions; expected deterministic aggregate values and no races under repeated GPU runs.
Failure modes: Any mismatch or compile diagnostic is a failed conversion. race-dependent totals, stale host values, or reduction type mismatch. Validation must include the smallest representative input and a CPU reference; do not waive a failure as numerical noise. Include launch-region behavior, ownership, and dimensional assumptions in review; the transform is complete only when those remain explicit. Confirm dimensional and reduction identities explicitly.
Evidence: commits 76ae550e0388c96fd315724628942c2c5394426b, 02772a36997539b183c70ba195e469478b894007
