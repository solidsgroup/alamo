# GPU-005: Make captured parameter types visible to extended lambdas
Status: draft
Class: correctness
Recognizer: manual: An extended device lambda is defined in protected/private scope or captures a protected/private unnamed parameter type.
Applies: An extended device lambda captures members through `this` or inaccessible protected/private parameter types.
Transform:
  Before:
    `[=] AMREX_GPU_DEVICE (...) { use this->pf; }`
  After:
    Expose/copy the parameter struct, capture it by value, and reference local fields inside the lambda.
Constraints: Mandatory when recognizer matches. This is a temporary visibility/type fix, not permission to enlarge public API unnecessarily; retain encapsulation outside the launch site.
Verify: `make -j4`; compare CPU/GPU outputs at identical timestep and dimensions; expected nvcc accepts the extended lambda and generated code has no inaccessible-type diagnostic.
Failure modes: Any mismatch or compile diagnostic is a failed conversion. nvcc reports private/protected captured type, or lambda still dereferences host `this`. Validation must include the smallest representative input and a CPU reference; do not waive a failure as numerical noise. Include launch-region behavior, ownership, and dimensional assumptions in review; the transform is complete only when those remain explicit.
Evidence: commits d522e1ac08729b306a215ad896d8a304983d55de, 4afea70999397683597dc2fd8dbf4ad004e86bb1
