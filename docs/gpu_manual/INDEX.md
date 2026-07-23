# GPU Port Manual Index

Load this file first. `[I]` is invariant CUDA/AMReX semantics, `[P]` is a
port-supplied contract or decision, and `[C]` is historical corpus evidence.
Regex is advisory: compiler diagnostics and recorded inspection outrank scan
results. Start a port with `ONBOARDING.md`; load only matching Tier 1 patterns
and the contract/policy they reference.

## Status

- `[C]` Transform evidence: 2/26 file-verified (GPU-007, GPU-016), 24 draft,
  0 transfer-verified, 0 cross-family. The repaired authoring targets do not
  demonstrate transfer.
- `[P]` Expanded onboarding and validation gate: `pending-pilot` until the next
  authorized non-Flame port completes it.
- `[P]` Port status is reported independently for scope, closure, inspection,
  validation, GPU safety, GPU-native shape, baseline efficiency, and harvest;
  see `STATUS.md`.

## Invariants

- `[I]` A device-reachable call chain contains only device-callable code and
  device-safe state; no hidden host pointer, virtual call, I/O, allocation, or
  logging path. Grounding: `evidence/primary-sources.md`.
- `[I]` Preserve ownership and asynchronous lifetime until every consuming
  stream is complete; synchronize at the required boundary, not reflexively in
  an inner loop. Grounding:
  `evidence/primary-sources.md#asynchrony-ownership-and-synchronization`.
- `[I]` Preserve dimensions, launch bounds, component ranges, reduction
  operators, and identities. Grounding:
  `evidence/primary-sources.md#amrex-launches-dimensions-reductions-and-explicit-execution`.
- `[P]` Preserve the CPU path declared by the port's scope and validation
  contract.
- `[P]` Scope, closure, physics oracles, tolerances, and performance evidence
  belong to the port. Flame commands and chamber-gpu diffs are not procedure.
- `[P]` FEATURE/[NUM] choices require the owner workflow in
  `ARCHITECTURE_POLICIES.md`; never import numerical behavior implicitly.

## Triage

### Correctness

- Host/device call diagnostic -> GPU-001; GPU-reachable runtime model or BC
  dispatch -> GPU-002, GPU-003, GPU-004.
- Extended-lambda/capture fault -> GPU-005, GPU-006; geometry pointer -> GPU-007.
- Storage lifetime or temporary Fab fault -> GPU-008, GPU-009; unsafe execution
  boundary -> GPU-011, GPU-023.
- Device logging/error or aggregate write -> GPU-012, GPU-013; direct FFT path
  -> GPU-016.
- Host-oriented iterative, branch-heavy, or container-backed numerical call
  chain -> GPU-031; annotation alone does not close it.
- Expression temporary, dimension, sentinel, accessor, or local initialization
  -> GPU-017, GPU-021, GPU-022, GPU-024, GPU-030.

### Scaffolding and optimization

- Temporary host quarantine or incomplete source closure -> GPU-010, GPU-025;
  both require lifecycle/retirement records.
- After `BASELINE_EFFICIENCY=pass`, profiled launch/load/branch/contraction
  hypotheses -> GPU-015, GPU-018, GPU-019, GPU-020.

## Patterns

- GPU-001 `[I]`: make the complete device-safe call chain callable.
- GPU-002 `[P]`: replace GPU-reachable runtime model pointers with value dispatch.
- GPU-003 `[P]`: forward port context through static selection.
- GPU-004 `[P]`: move virtual boundary-condition dispatch outside kernels.
- GPU-005 `[I]`: expose device-lambda parameter types safely.
- GPU-006 `[I]`: hoist device-safe launch state; do not capture host `this`.
- GPU-007 `[I]`: copy geometry values instead of capturing host pointers.
- GPU-008 `[I]`: stage storage and fence its asynchronous lifetime boundary.
- GPU-009 `[I]`: retain temporary Fab storage through consuming launches.
- GPU-010 `[P]`: quarantine only named, owned, retiring host scaffolding.
- GPU-011 `[I]`: make CPU/GPU loop execution explicit.
- GPU-012 `[P]`: propagate device errors to one owned host observation boundary.
- GPU-013 `[I]`: produce aggregate results with device reductions.
- GPU-015 `[P]`: optionally fuse equivalent component launches after profiling.
- GPU-016 `[P]`: route transforms through the selected GPU-capable wrapper.
- GPU-017 `[I]`: materialize unsafe fixed-size expression intermediates.
- GPU-018 `[P]`: optionally hoist repeated immutable loads after profiling.
- GPU-019 `[P]`: optionally defer work to its consuming branch after profiling.
- GPU-020 `[P]`: optionally specialize a consumed contraction after profiling.
- GPU-021 `[I]`: guard constants and indices by configured dimension.
- GPU-022 `[I]`: use device-callable sentinel values.
- GPU-023 `[I]`: make Fab execution location explicit where AMReX requires it.
- GPU-024 `[I]`: make read-only device-reached accessors callable and const.
- GPU-025 `[P]`: declare and retire the selected integrator closure scaffold.
- GPU-030 `[I]`: initialize device-local aggregates before conditional writes.
- GPU-031 `[P]`: port a host-only numerical algorithm as an explicit sub-port.

## Contracts and evidence

- `[P]` Correctness: `VALIDATION.md`; device workload shape:
  `GPU_NATIVE_SHAPE.md`; baseline efficiency: `PERFORMANCE.md`; recognizers:
  `RECOGNIZERS.md`; status/harvest: `STATUS.md`.
- `[P]` Templates: `templates/`; known limits: `BLIND_SPOTS.md`.
- `[C]` `FEATURES.md`, `ONE_OFFS.md`, `BUILD_LOG.md`, and `evidence/` retain the
  chamber-gpu corpus as evidence, not a port plan.
