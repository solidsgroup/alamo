# Agent Brief: GPU Port and Manual Maintainer v3

Mission: use the three-tier GPU manual to start, validate, and harvest any ALAMO
integrator port without access to `chamber-gpu`. The manual separates universal
device semantics, port-supplied decisions and oracles, and historical corpus
evidence. Every completed port must improve the manual.

## 1. Authority and labels

Evidence authority, strongest first:

1. Executable results reproduced for the selected port and declared closure:
   strict compiler, analytic/conservation/golden/restart/multi-box oracles,
   sanitizer, and recorded timeline evidence.
2. Primary CUDA and AMReX documentation for universal execution semantics.
3. Executable results from other completed ports.
4. `chamber-gpu` diffs and markdown, explicitly labeled corpus evidence.
5. Inference, which may propose a finding but may not close one.

Tier 0 uses `[I]` for invariant semantics, `[P]` for a port contract, and `[C]`
for corpus evidence. Tier 1 keeps the same separation in `Invariant`, `Port
contract`, and `Corpus example` fields. Scanner results are always advisory;
compiler diagnostics and recorded inspection outrank regex hits.

## 2. Layout and token budget

```text
docs/gpu_manual/
  INDEX.md                         Tier 0, <= 1500 approximate tokens
  patterns/GPU-NNN-*.md           Tier 1, 200-400 approximate tokens each
  evidence/*.md                   Tier 2, never loaded by default
  VALIDATION.md                   reusable correctness contract
  ONBOARDING.md                   ordered port protocol and gates
  ARCHITECTURE_POLICIES.md        single homes for recurring decisions
  GPU_NATIVE_SHAPE.md             mandatory layout/kernel/resource contract
  PERFORMANCE.md                  baseline-efficiency contract
  RECOGNIZERS.md                  advisory recognizer lifecycle
  STATUS.md                       pattern, port, closed-book, harvest status
  BLIND_SPOTS.md                  known limitations
  templates/*                     schemas copied by a port
  recognizers/table.csv           versioned candidate/converted rules
  recognizers/scan.py             stateful site scanner
  templates/COVERAGE.csv          per-port advisory status schema
  FEATURES.md / ONE_OFFS.md       chamber-gpu corpus ledgers
  BUILD_LOG.md                    append-only manual decisions
  build/                           manual construction evidence
```

The three-tier layout and token budgets are unchanged. Contract and policy
documents are loaded only when the corresponding Tier 0 or Tier 1 entry points
to them.

## 3. Frozen v3 schemas

### Tier 1 pattern

```text
# GPU-NNN: <imperative name>
Transform status: draft | file-verified | transfer-verified | cross-family
Class: correctness | optimization | scaffolding
Detection: <advisory regex, manual inspection, or compiler diagnostic>
Invariant: <universal CUDA/AMReX semantic claim>
Port contract: <decision, state, or proof this port must supply>
Transform:
  Before: <bad semantic shape>
  After: <device-safe semantic shape>
Corpus example: <clearly labeled chamber-gpu example; never procedure>
Constraints: <when not to apply and preserved semantics>
Verify: <VALIDATION.md or PERFORMANCE.md categories and observation>
Failure modes: <compile/runtime/review signature>
Evidence: Primary: <primary-sources anchor>; Corpus: <hash/path/result>
```

`file-verified` means the exercise used a target that authored or repaired the
pattern. `transfer-verified` requires a frozen, unseen target; a repair makes
that target part of the authoring corpus and a new target is required.
`cross-family` additionally requires completed-port harvests from two distinct
physics families. Recognizer state never upgrades transform status.

### Recognizer table and coverage

`recognizers/table.csv` schema v3:

```text
schema_version,pattern_id,type,candidate_expression,converted_expression,exclude_expression,confirmation,notes
```

Each port writes its own coverage report using schema v3:

```text
schema_version,port_id,source_revision,site_id,file,line,pattern_id,state,evidence
```

Allowed states are `candidate`, `converted`, `not-applicable`, and
`false-positive`. `port_id` and `source_revision` are required; dispositions
carry only within the same port and unchanged revision. A rerun recomputes candidate/converted sites
and preserves reviewed dispositions. It never hides a candidate because the
same file also contains a converted site. There is no operational root coverage
snapshot. Schema changes increment the version and receive a `BUILD_LOG.md`
entry.

### Port artifacts

- `SCOPE.md`: supported execution, physics/capability scope, FEATURE/[NUM]
  decisions, owner approval.
- `CLOSURE.csv`: compiler-first transitive source graph and scaffolding state.
- `INSPECTION_LEDGER.csv`: per-file taxonomy questions and dispositions.
- `VALIDATION.csv`: oracle category, exact command/reference/result, tolerance
  rationale, evidence, owner.
- `EFFICIENCY.md`: fixed baseline checklist plus trace reproduction.
- `FIELD_LAYOUT.csv`, `KERNEL_GRAPH.csv`, and `SHAPE_PROFILE.md`: mandatory
  workload-shape, transfer, call-chain, and compiler/profiler evidence.
- `COVERAGE.csv`: revision-bound advisory findings for this port.
- `FEATURE_DECISIONS.csv`: surfaced capability/numerical decisions.
- `HARVEST.md`: misses, recognizer feedback, status evidence, mandatory manual
  update, future work.
- `PORT_STATUS.csv`: independent gate axes and closed-book outcomes used to
  compare evidence across completed ports.

Templates under `templates/` are the schema definitions; filled instances live
with the consuming port's evidence.

## 4. Port and maintenance phases

### Phase 0: Intake

Read INDEX, STATUS, ONBOARDING, and BLIND_SPOTS. Select a stable `port_id`, owner,
entry point, and toolchain. Do not use `chamber-gpu` as procedure.

Exit: scope artifact exists and contains no undecided FEATURE/[NUM] item.

### Phase 1: Compiler-first closure

Derive the transitive closure from the entry point with strict device compiler
diagnostics. Record every translation unit, header, template/runtime edge,
exclusion, and scaffold. Inspection can add an edge; regex cannot remove one.

Exit: the exact closure links in all supported dimensions and every diagnostic
has a closure/worklist disposition.

### Phase 2: Inspection worklist

For every closure file, instantiate the fixed taxonomy in ONBOARDING: host
loops, launches, receiver types, diagnostics, captures, lifetimes, reductions,
dispatch, host-only numerical kernels, field layout, and kernel graph. Generate
revision-bound scanner sites and merge them without treating them as proof.

Exit: bounded worklist exists; no required taxonomy or scan row is unexplained.

### Phase 3: Oracle-gated conversion

Fill the validation contract before conversion. Apply only matching Tier 1
patterns, preserve scope and existing numerical behavior, and gate each bounded
change with the smallest relevant oracle. Surface architecture or numerical
choices through the policy workflow.

Exit: zero open correctness rows and every applicable validation category
passes with a written tolerance rationale.

### Phase 4: GPU-native shape and baseline

Fill `GPU_NATIVE_SHAPE.md` and run the physics-agnostic procedure in PERFORMANCE.
A port may report a reasoned failure, but it is not GPU-native until layout,
kernel graph, residency/transfers, numerical call-chain complexity, compiler
resources, and all applicable baseline items pass. Optimization patterns remain
optional, measured future work.

Exit: `GPU_NATIVE_SHAPE=pass` and `BASELINE_EFFICIENCY=pass` with named-GPU
compiler/timeline evidence, or honest `fail`/`blocked` statuses.

### Phase 5: Closed-book gate and harvest

Run the frozen unseen-target gate in STATUS on this port's material: a fresh session first
produces scope/closure/inspection artifacts, then converts under its validation
contract. Update only statuses actually exercised. Fill HARVEST, update blind
spots, add executable evidence, and make at least one manual change.

Exit: reviewed harvest exists, scanner candidates are zero for the completed
port, and the manual contains the port's new evidence or correction.

## 5. Recognizer feedback loop

After compiler and inspection closure, promote any missed recurring shape to a
candidate recognizer rule. Link every scan hit to an inspection row; disposition
unmatched hits as a new finding, `false-positive`, or `not-applicable` with
evidence. Rerun until a completed port has zero `candidate` rows. Do not maintain
formal precision/recall or a held-out corpus for this single-developer manual.

## 6. Stop conditions

- A compiler diagnostic or required inspection row is absent from the plan:
  add it before continuing.
- A scan is clean while manual/compiler rows remain open: the port is incomplete.
- A physics/capability/numerical choice appears mid-port: preserve behavior and
  run the FEATURE/[NUM] decision workflow.
- Required hardware/toolchain evidence is unavailable: mark the affected gate
  `blocked`, never pass.
- Shape evidence is missing for layout, launch decomposition, transfer bytes,
  numerical call chains, or compiler resources: fail GPU-native shape.
- Baseline trace shows a hot host fallback, per-tile transfer/synchronization,
  runtime kernel dispatch, or avoidable component launch multiplication: fail
  baseline efficiency.
- A performance change lacks a named-GPU before/after result or changes a
  correctness oracle: keep it future work.
- Completing a port would leave the manual unchanged: the harvest is incomplete.

## 7. Definition of done

- A reader can distinguish invariant, port contract, and corpus evidence at the
  point of use.
- A new port can start from scope, compiler-first closure, inspection,
  validation, shape, efficiency, decision, coverage, and harvest templates
  without chamber-gpu.
- Correctness Verify fields reference the validation contract; optimization
  fields are gated behind a passing baseline.
- Scanner output is versioned, stateful, advisory, and fully dispositioned for
  every completed port.
- Architecture policies have one home, an owner decision, and a worked example.
- Status is per pattern and per port axis; no partial gate becomes a global pass.
- Every completed port produces a harvest and changes the manual.

The next authorized non-Flame port is the pilot for the contracts and expanded
closed-book gate. Until it finishes, pilot-dependent acceptance remains
`pending-pilot` rather than inferred from Flame or file-only exercises.
