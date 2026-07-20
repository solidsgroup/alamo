# GPU-021: Use dimension-aware compile-time constants
Status: draft
Class: correctness
Recognizer: manual: Shared 2-D/3-D code uses an unguarded fixed-dimensional factor, component, or index that is inactive in one configured dimension.
Applies: Device stencil/constants code that indexes dimensions not present in the configured build.
Transform:
  Before:
    Unconditionally define or access 3-D constants/terms in 2-D and 3-D code.
  After:
    Guard terms with `AMREX_SPACEDIM`/`AMREX_D_DECL` and provide dimension-valid constants.
Constraints: Keep scope narrow. Mandatory when recognizer matches. Keep 2-D and 3-D formulas equivalent in their active dimensions; do not hide required components behind runtime branches.
Verify: `./configure --dim=2 --debug && make -j4`, then repeat for 3-D; expect both builds and stencil tests to pass. Repeat in the configured 2-D and 3-D modes where applicable, and retain the command/output in the build log.
Failure modes: Unprotected z terms cause compile errors or invalid indexing in 2-D; wrong constants change gradients and golden results. Review the failing signature before changing the pattern; do not broaden its scope to silence an unrelated failure.
Evidence: chamber-gpu commit 5efce7cff26a13c2741b011d08735d2a679ef036; `src/Numeric/Stencil.H:1402-1455`.
