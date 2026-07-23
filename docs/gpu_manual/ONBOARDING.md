# Port onboarding protocol

This protocol is the ordered entry point for a new integrator. It requires no
`chamber-gpu` access. Store filled artifacts under the port's own evidence
directory and give every artifact the same `port_id`.

## 1. Declare scope

Copy `templates/SCOPE.md`. Record supported dimensions and backends, in-scope
entry points and physics, explicitly excluded paths, CPU-path obligations, and
every FEATURE/[NUM] decision. Scope is frozen before conversion; changes require
a dated owner decision and re-running affected gates.

Gate: no `undecided` scope item and no implicit FEATURE/[NUM] import.

## 2. Derive the source closure compiler-first

Copy `templates/CLOSURE.csv`. Start from the selected executable and integrator
entry point, build with strict device diagnostics, and record each required
translation unit, header, template instantiation, model, operator, and runtime
selection as the compiler exposes the transitive graph. Inspection may add a
missed runtime edge; regex may only suggest one. Never broaden the closure just
to silence a diagnostic, and never exclude a required path to claim success.

Gate: every strict-compiler diagnostic is a closure or worklist row with an
owner and disposition; the exact declared closure links in every supported
dimension.

## 3. Build the per-file inspection ledger

Copy `templates/INSPECTION_LEDGER.csv`. For every file in the closure, create
and disposition each applicable taxonomy row:

- `host-loop`: hot or GPU-reachable serial loops and host fallbacks;
- `launch`: bounds, components, dimension, backend, and synchronization;
- `receiver-type`: iterator, Fab/MultiFab, pointer, and dispatch receiver type;
- `diagnostic`: abort, logging, exception, and non-finite paths;
- `capture`: `this`, members, references, ownership, and value safety;
- `lifetime`: staging, temporary storage, Elixir/arena, and async destruction;
- `reduction`: identity, operator, destination, determinism, and host result;
- `dispatch`: virtual/plugin/runtime selection in a device-reachable chain.
- `numerical-kernel`: host-oriented solver/model calls, iteration/convergence,
  exceptions, branching/divergence, and transitive device callability;
- `field-layout`: component ordering, access locality, ghost use, and residency;
- `kernel-graph`: MFIter/box/tile granularity, launch multiplicity,
  dependencies, transfers, synchronization, reductions/atomics, and allocation.

Generate coverage for this port and revision; the manual does not ship an
operational snapshot:

```text
python3 docs/gpu_manual/recognizers/scan.py \
  --root <repository-or-isolated-closure> \
  --table docs/gpu_manual/recognizers/table.csv \
  --port-id <port-id> --source-revision <revision> \
  --out <port-evidence>/COVERAGE.csv \
  [--previous <prior-port-coverage.csv>]
```

Scanner sites join the ledger as advisory findings. Every hit is linked to a
row or dispositioned `false-positive`/`not-applicable`; a clean scan never
closes an inspection row. Never reuse another port's coverage file.

Gate: zero `open` ledger rows and zero unexplained scanner `candidate` rows.

## 4. Define validation and GPU-native shape

Instantiate `VALIDATION.md` before changing code. Inventory fields and the
proposed timestep/kernel graph with `GPU_NATIVE_SHAPE.md`; this must expose
layout, component order, MFIter granularity, transfers, numerical call-chain
complexity, resource evidence to collect, shared-memory/atomic decisions, and
arena/allocation lifetime. A host-only numerical algorithm is a GPU-031 sub-port,
not a batch of annotations.

Gate: the validation deck is runnable; every pre-conversion shape row has an
owner and proposed disposition. Missing hardware evidence remains `open` until
post-conversion measurement.

## 5. Convert under oracle gating

Apply Tier 1 patterns only when their invariant and port contract match a ledger
row. Re-run the smallest relevant oracle after each bounded conversion; run the
full contract after the closure is clean. Complete the compiler resource report,
device timeline, field map, and kernel graph. Numerical or capability questions
enter the decision workflow in `ARCHITECTURE_POLICIES.md` rather than
hitchhiking on the port.

Gate: validation passes, then `GPU_NATIVE_SHAPE=pass`, then
`BASELINE_EFFICIENCY=pass`. Optional optimization starts only afterward.

## 6. Close and harvest

Retire or explicitly time-bound every closure/quarantine scaffold. Copy
`templates/HARVEST.md` and `templates/PORT_STATUS.csv`, update
`BLIND_SPOTS.md`, add newly reproduced evidence, and update at least one manual
artifact: recognizer disposition/rule, pattern evidence/status, policy
clarification, or blind-spot entry. A port is not closed until this harvest is
reviewed.
