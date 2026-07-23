# GPU-native workload shape contract

This mandatory contract prevents a CPU algorithm from passing merely because it
wears device annotations. Fill it twice: first during inspection to choose a
credible device shape, then after conversion with compiler and profiler
evidence. It is physics-agnostic and does not prescribe an optimization.

Copy `templates/FIELD_LAYOUT.csv`, `templates/KERNEL_GRAPH.csv`, and
`templates/SHAPE_PROFILE.md` into the port evidence directory.

## 1. Field layout and component ordering

Inventory every steady-state field: storage type, number and meaning of
components, component order, ghost region, kernel consumers, and access pattern.
Record whether adjacent threads access adjacent data, whether a kernel traverses
unused components, and whether separate fields or components are repeatedly
gathered. Give a coalescing/locality rationale. AoS, SoA, and component ordering
are decisions to measure, not universal answers.

Gate: no unexplained layout-hostile hot access, whole-field traversal for a
small subset, or component-serialized launch structure.

## 2. Kernel graph and iteration granularity

Draw the steady-state timestep as kernels, library calls, reductions, transfers,
and synchronization edges. Record each `MFIter`/box/tile iteration space,
launches per step, components per launch, dependencies, and fusion/splitting
rationale. Compare per-box and whole-level AMReX launch forms where semantics
allow; never assume CPU tiling is the right GPU granularity.

Gate: no unexplained launch explosion, repeated traversal, per-tile global
synchronization, or host stage between device producers and consumers.

## 3. Residency and transfer accounting

For each graph edge, record memory owner, arena, host/device residence, transfer
direction, frequency, and bytes per timestep. Separate required bounded scalar
results and I/O from bulk field traffic. Managed-memory success does not waive
the accounting.

Gate: hot fields remain device-resident and all recurring transfers have an
owner-approved necessity and timeline evidence.

## 4. Numerical call-chain complexity

For every hot kernel, record device-call depth, runtime dispatch, host/container
state, exceptions/diagnostics, dynamic allocation, fixed or data-dependent
iteration bounds, early exits, and branch/divergence hotspots. A Riemann solver,
constitutive update, or failure-surface evaluation is a sub-port: annotation
alone is not a disposition. GPU-031 defines that correctness pattern.

Gate: no unresolved host-only edge; convergence/termination and error
propagation are explicit; divergence and loop bounds have an inspection result.

## 5. Resource snapshot

On the named GPU and exact binary, record block size, registers per thread,
local-memory spills, static/dynamic shared memory, achieved occupancy/active
warps, and bandwidth where meaningful. Relate resource pressure to inlined
matrix math, branches, and kernel boundaries. No universal occupancy or register
threshold applies, but missing or unexplained limiting resources block the gate.

## 6. Shared memory, atomics, and allocation

Record each shared-memory use or explicit decision not to use it, including
reuse and bank-conflict rationale. Record every atomic/reduction target,
operator/identity, contention expectation, ordering requirement, and bounded
host result. Record arena choice, allocation frequency, lifetime, and whether
allocation occurs inside a timestep, tile, or cell path.

Gate: no per-cell device allocation, unexplained hot allocation, unsafe atomic,
or allocation/lifetime mismatch.

## Result

`GPU_NATIVE_SHAPE=pass` requires every row to be `pass` or owner-approved
`not-applicable`, with the field map, kernel graph, compiler resource report,
timeline, and validation links present. `blocked` and `fail` are valid honest
outcomes, not passes. This gate precedes `BASELINE_EFFICIENCY` and optional
measured optimization.

Grounding: `evidence/primary-sources.md#workload-shape-and-resource-evidence`.

