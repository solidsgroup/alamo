# GPU-001: Make the complete device call chain callable
Status: draft
Class: correctness
Recognizer: build-error: `calling a __host__ function from a __device__ function`
Applies: A GPU kernel reaches a helper or model method that lacks a host/device annotation.
Transform:
  Before:
    `AMREX_FORCE_INLINE Set::Scalar operator()(...)` is reached from an `AMREX_GPU_DEVICE` lambda.
  After:
    Add `AMREX_GPU_HOST_DEVICE` to that helper and every device-safe transitive callee; keep host-only operations outside the chain.
Constraints: Mandatory when the recognizer matches. Annotate the full transitive value path, including inline helpers and operators. Do not annotate routines that perform I/O, allocation, virtual dispatch, or other host-only work; split those paths instead and keep their host boundary explicit.
Verify: `make -j4`; compare CPU/GPU outputs at identical timestep and dimensions; expected nvcc compiles the target and no host/device call diagnostic remains.
Failure modes: Any mismatch or compile diagnostic is a failed conversion. nvcc reports the quoted diagnostic with the inserted function name, or runtime fails after a host-only call executes on device. A partial annotation can compile one translation unit while failing in a deeper template instantiation. Validation must include the smallest representative input and a CPU reference; do not waive a failure as numerical noise.
Evidence: commits f296eae5d31d39d2bf0bbfefdfddc5ea8142b0da, 3d382a5655c715c7483a15d53c3f8e5c908330af, 4afea70999397683597dc2fd8dbf4ad004e86bb1
