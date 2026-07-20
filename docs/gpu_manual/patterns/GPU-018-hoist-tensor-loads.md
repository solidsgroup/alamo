# GPU-018: Hoist reused field and tensor loads
Status: draft
Class: performance
Recognizer: regex: `AMREX_GPU_DEVICE\s*\([^)]*\)\s*\{(?:(?!\}\s*\);)[\s\S]){0,6000}?\bDDW\s*\(\s*i\s*,\s*j\s*,\s*k\s*\)(?:(?!\}\s*\);)[\s\S]){0,6000}?\bDDW\s*\(\s*i\s*,\s*j\s*,\s*k\s*\)`
Applies: Device kernels repeatedly loading the same model/tensor or member state for one cell.
Transform:
  Before:
    Repeated `DDW(i,j,k)`, `m_psi_set`, or member reads in several branches.
  After:
    Load once into a local (`auto ddw = DDW(i,j,k); auto psi = m_psi_set;`) and reuse it.
Constraints: Apply only with profiling justification; preserve evaluation order, lifetime, and correctness/golden results. Hoist immutable values only; do not cache data whose value changes during the kernel.
Verify: `make -j4`; run elastic correctness/golden tests and a profiler comparison; expect identical checksums/residuals and reduced loads or kernel time. Repeat in the configured 2-D and 3-D modes where applicable, and retain the command/output in the build log.
Failure modes: Hoisting a mutable or out-of-bounds value produces stale physics or illegal access; unproven changes may increase registers and reduce performance. Review the failing signature before changing the pattern; do not broaden its scope to silence an unrelated failure.
Evidence: chamber-gpu commit 9470889b14f10a902dab6dd9573deefc09972e8d; `src/Operator/Elastic.cpp:173-264#5`, `317-347#2`; `docs/agent_plans/20260709-fapply-kernel-surgery/results/RESULT.md`.
