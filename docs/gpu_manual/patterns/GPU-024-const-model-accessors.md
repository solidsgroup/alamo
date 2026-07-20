# GPU-024: Const-qualify captured model accessors
Status: draft
Class: correctness
Recognizer: manual: A named model accessor such as `get_K`, `get_rho`, `get_cp`, `get_qdot`, or `get_L` is invoked through a const object in device code but lacks a trailing `const`.
Applies: Model methods called through const objects captured by device lambdas.
Transform:
  Before:
    A device-reachable accessor lacks trailing `const` or mutates implicit host state.
  After:
    Declare the accessor `AMREX_GPU_HOST_DEVICE ... const` and keep its implementation read-only.
Constraints: Keep scope narrow. Mandatory when recognizer matches. Do not add const by casting away mutation or change model semantics; stage mutable state separately.
Verify: `make -j4 && scripts/runtests.py --dim=3 --serial tests/GPU/F1_smoke_flame_only`; expect nvcc const-capture compilation and unchanged Flame source-term checks.
Failure modes: Missing const causes device compile errors or forces unsafe `this` capture; fake constness can race shared model data. Review the failing signature before changing the pattern; do not broaden its scope to silence an unrelated failure.
Evidence: chamber-gpu commit 587853d79432c0394928edb457683961ef5f08da; `src/Model/Propellant/Propellant.H:49-108`, `src/Util/PNG.H:155-161#2`.
