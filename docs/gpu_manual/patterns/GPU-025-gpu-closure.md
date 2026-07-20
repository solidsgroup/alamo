# GPU-025: Declare and enforce the supported GPU closure
Status: draft
Class: scaffolding
Recognizer: manual: inspect GPU build policy and integrator source/object closure
Applies: GPU builds that accidentally compile unsupported integrators or omit required device-safe sources.
Transform:
  Before:
    General launcher/object closure includes every integrator and relies on accidental nvcc compatibility.
  After:
    Define a narrow supported closure in `IntegratorPolicy.mk`, select supported sources, and fail unsupported paths explicitly.
Constraints: Temporary incremental-porting scaffolding whose objective is to shrink the closure footprint, never grow it. Add a source only after a correctness pass; do not use closure exclusions to hide a required production path.
Verify: `make -n ALAMO_GPU_INTEGRATOR=flame` then `make -j4`; expect only the declared Flame closure in the object list and a green GPU build; `make -n ALAMO_GPU_INTEGRATOR=bogus` must fail with the policy diagnostic.
Failure modes: Closure growth reintroduces non-nvcc-clean code; over-pruning causes link errors or silently removes required physics. Review the failing signature before changing the pattern; do not broaden its scope to silence an unrelated failure.
Evidence: chamber-gpu commit 54a941433b7582578cb5d56db794b1d648fb03cc; `src/GPU/IntegratorPolicy.mk`, `src/alamo_gpu.cc`; `docs/gpu_safe_ic_bc_matrix.md`; `benchmark/GPU_BRANCH_GUIDE.md`.
