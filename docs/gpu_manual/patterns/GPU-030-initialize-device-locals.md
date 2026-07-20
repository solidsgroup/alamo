# GPU-030: Initialize device-local aggregates before use
Status: draft
Class: correctness
Recognizer: regex: `Set::Vector\s+ximg\s*;|[A-Za-z]+\s+[A-Za-z_]+\s*;\s*//.*device`
Applies: Local vectors/aggregates populated conditionally inside a device-callable sampler or kernel.
Transform:
  Before:
    `Set::Vector ximg;` followed by partial component assignments.
  After:
    `Set::Vector ximg = Set::Vector::Zero();` before conditional writes.
Constraints: Keep scope narrow. Mandatory when recognizer matches. Initialize only to the mathematically neutral value and preserve all subsequent clamps/assignments. Initialize before every control-flow path; do not rely on a branch to assign every active dimension.
Verify: `make -j4 && scripts/runtests.py --dim=2 --serial tests/PNG`; expect no uninitialized-value diagnostic and unchanged PNG field comparison.
Failure modes: Uninitialized components yield nondeterministic device results, NaNs, or sanitizer failures; wrong zeroing changes boundary sampling. Review the failing signature before changing the pattern; do not broaden its scope to silence an unrelated failure.
Evidence: chamber-gpu commit f296eae5d31d39d2bf0bbfefdfddc5ea8142b0da; `src/Util/PNG.H:168-174`.
