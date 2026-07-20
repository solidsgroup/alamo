# GPU-020: Specialize only the consumed tensor contraction
Status: draft
Class: performance
Recognizer: regex: `Matrix4.*operator\*|for\s*\(.*AMREX_SPACEDIM.*\).*operator\(.*\)`
Applies: Generic Matrix4×Matrix3 contractions whose full operator surface is used in a hot device kernel.
Transform:
  Before:
    Generic nested loops/operator dispatch for a contraction.
  After:
    A dedicated contraction for the consumed symmetry/layout, with explicit terms or a narrow helper.
Constraints: Keep scope narrow. Apply only with profiling justification; preserve exact accumulation order and correctness/golden results. Do not specialize unused paths or silently change symmetry assumptions.
Verify: `make -j4`; run elastic golden tests and compare register/instruction metrics; expect bitwise or accepted-tolerance agreement and a measured improvement. Repeat in the configured 2-D and 3-D modes where applicable, and retain the command/output in the build log.
Failure modes: Wrong index mapping or accumulation order corrupts stress/flux; specialization can increase code size or registers without benefit. Review the failing signature before changing the pattern; do not broaden its scope to silence an unrelated failure.
Evidence: chamber-gpu commit 9470889b14f10a902dab6dd9573deefc09972e8d; `src/Set/Matrix4_Major.H:531-550`; `docs/agent_plans/20260709-fapply-kernel-surgery/results/RESULT.md`.
