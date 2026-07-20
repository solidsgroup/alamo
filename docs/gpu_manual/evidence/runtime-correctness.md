# Runtime correctness evidence

- `docs/llm/BUG_PATTERNS.md`: a local `FArrayBox` consumed across the MFIter stream pool can be reclaimed early; multi-box interpolation exposed silent corruption. This supports GPU-009 and requires multi-box verification.
- `docs/llm/changelog/2026-07-02-gpu-audit.md`: reading a device error flag after synchronizing only the current stream produced silent false negatives. `streamSynchronizeAll()` at the host check/lifetime boundary supports GPU-008 and GPU-012.
- `docs/llm/BUG_PATTERNS.md`: direct host-member accumulation failed with CUDA error 700 on an A100 while local HMM masked it. Device reduction followed by host assignment supports GPU-013.
- `docs/llm/BUG_PATTERNS.md`: chained Eigen inverse/transpose evaluation failed with CUDA error 719 in 3-D. Materializing the intermediate supports GPU-017.
- `docs/gpu_safe_ic_bc_matrix.md`: host-loop IC/BC paths remain guarded outside the supported device closure. This is evidence for temporary GPU-010/GPU-025 scaffolding, not permission to expand it.

Verification evidence: `benchmark/ci_golden_compare.sh`, `benchmark/golden_compare_flame.sh`, repeated strict multi-box runs, and compute-sanitizer `ERROR SUMMARY: 0 errors` are recorded in `benchmark/GPU_BRANCH_GUIDE.md` and the 2026-07-02 audit changelog.
