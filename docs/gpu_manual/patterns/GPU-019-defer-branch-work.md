# GPU-019: Defer branch-only expensive work
Status: draft
Class: performance
Recognizer: regex: `if\s*\(.*boundary.*\).*DDW|if\s*\(.*conservative.*\)`
Applies: Kernels computing expensive tensor/stencil quantities before a branch that may not consume them.
Transform:
  Before:
    Compute `ddw`, stress, or boundary tensors unconditionally, then test whether the row/face needs them.
  After:
    Test the cheap branch predicate first; construct expensive values only in the consuming branch.
Constraints: Apply only with profiling justification; preserve branch semantics, numerical ordering, and correctness/golden results. Do not defer values required by multiple branches without retaining equivalent initialization.
Verify: `make -j4`; run elastic golden/regression checks plus profiler comparison; expect unchanged residuals/outputs and lower instruction or kernel time. Repeat in the configured 2-D and 3-D modes where applicable, and retain the command/output in the build log.
Failure modes: Missing a branch assignment leaves undefined values; moving computation can change synchronization or constitutive ordering. Review the failing signature before changing the pattern; do not broaden its scope to silence an unrelated failure.
Evidence: chamber-gpu commit 332ecffdd93a83def543d0348d592f815950cc5f; `src/Operator/Elastic.cpp:173-264#5`.
