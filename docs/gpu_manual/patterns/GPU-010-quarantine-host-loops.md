# GPU-010: Quarantine host-only loops during incremental porting
Status: draft
Class: scaffolding
Recognizer: manual: Host loop or virtual/plugin path reached from a GPU integrator entry point.
Applies: A conversion cannot yet make a loop/device dependency GPU-safe.
Transform:
  Before:
    GPU entry point executes an unconverted host-only loop or polymorphic callback.
  After:
    Temporarily isolate the loop behind a host-only branch/dispatch boundary and keep the GPU path explicit.
Constraints: Temporary incremental-porting scaffolding. Objective is to shrink the quarantine, never grow it. Dodging a hard conversion requires user approval; do not silently expand the host fallback.
Verify: `make -j4`; compare CPU/GPU outputs at identical timestep and dimensions; expected GPU target builds and the quarantined path is observable in logs/tests.
Failure modes: Any mismatch or compile diagnostic is a failed conversion. hidden CPU fallback, missing updates, or unsupported host call in device code. Validation must include the smallest representative input and a CPU reference; do not waive a failure as numerical noise.
Evidence: commits d522e1ac08729b306a215ad896d8a304983d55de
