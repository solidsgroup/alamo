# GPU Pattern Index

Load this file first. Use `COVERAGE.csv` to select candidate files, then load only the matched Tier 1 files. Manual and build-error recognizers still require inspection.

## Invariants

- Apply every correctness pattern whose recognizer and Applies condition match; a clean regex scan does not waive manual or compiler-diagnostic gates.
- Apply performance patterns only with profiling justification and preserve CPU/golden correctness.
- Scaffolding is temporary incremental-porting machinery: shrink its footprint, never grow it. GPU-010 requires user approval when used to avoid a hard conversion.
- Device closures contain only device-callable code and device-safe value state: no virtual dispatch, host pointer, host logging, I/O, allocation, or hidden `this` capture.
- Preserve execution dimension, component ranges, reduction identities, ownership, synchronization, and the tested CPU path.
- Treat `FEATURES.md` as an explicit do-not-import list. A FEATURE requires task-level user opt-in.
- Surface every `[NUM]` ONEOFF or FEATURE to the user; never apply numerical behavior changes implicitly.

## Triage

### Correctness

- nvcc host/device call diagnostic -> GPU-001; runtime model pointer or virtual BC call -> GPU-002, GPU-004; static parser drops context -> GPU-003.
- Extended-lambda access failure -> GPU-005; kernel reaches member/host state -> GPU-006; raw `CellSize()` pointer -> GPU-007.
- Async storage/lifetime fault -> GPU-008, GPU-009; backend loop lacks a safe split -> GPU-011.
- Kernel logs/aborts or writes host aggregates -> GPU-012, GPU-013; direct host FFT path -> GPU-016.
- Chained fixed-matrix temporary -> GPU-017; fixed-dimension assumption -> GPU-021; host sentinel -> GPU-022.
- Implicit Fab execution mode -> GPU-023; non-const model getter -> GPU-024; uninitialized device local -> GPU-030.

### Performance

- Component loop launches one kernel per component -> GPU-015.
- Repeated per-cell tensor/member load -> GPU-018; expensive value built before a selective branch -> GPU-019.
- Generic hot tensor contraction computes more than the consumer needs -> GPU-020.

### Scaffolding

- Unconverted host/plugin loop blocks incremental progress -> GPU-010, with user approval for hard-conversion avoidance.
- GPU target builds an unclosed integrator object graph -> GPU-025.

## Patterns

- GPU-001: annotate the complete device-safe transitive call chain.
- GPU-002: replace Gas-family runtime base pointers with concrete value dispatch.
- GPU-003: forward constructor/context arguments through static selection recursion.
- GPU-004: replace kernel virtual BC calls with static POD-state dispatch.
- GPU-005: expose and name parameter types captured by extended lambdas.
- GPU-006: hoist device-safe member values into by-value launch locals.
- GPU-007: copy geometry cell sizes into a device-safe value.
- GPU-008: stage host storage and fence asynchronous lifetime boundaries.
- GPU-009: retain an Elixir while launches consume temporary Fab storage.
- GPU-010: temporarily quarantine an approved host-only loop.
- GPU-011: make CPU/GPU loop dispatch explicit.
- GPU-012: report kernel failures through a device flag checked on host.
- GPU-013: use ReduceOps for aggregate results; fuse only as a profiled variant.
- GPU-015: fuse equivalent component launches into a component-aware launch.
- GPU-016: route transforms through the GPU-capable FFT wrapper.
- GPU-017: materialize chained fixed-size matrix expressions before reuse.
- GPU-018: hoist immutable repeated tensor and member loads.
- GPU-019: defer expensive work until its consuming branch is selected.
- GPU-020: specialize only the tensor contraction result actually consumed.
- GPU-021: guard constants and indices by configured dimension.
- GPU-022: replace host-oriented sentinel storage with device-callable values.
- GPU-023: run Fab initialization and algebra explicitly on device.
- GPU-024: make read-only model accessors device-callable and `const`.
- GPU-025: close the GPU build over the selected integrator source graph.
- GPU-030: initialize device-local aggregates before conditional writes.

## Numerical changes

No NUM Tier 1 transform survives device-port classification. Numerical material is opt-in evidence or a `[NUM]` entry in `FEATURES.md`/`ONE_OFFS.md`.

## One-offs

See `ONE_OFFS.md`. Consult `FEATURES.md` separately as the upstream roadmap and do-not-import ledger.
