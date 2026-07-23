# GPU-native baseline efficiency contract

Baseline efficiency is mandatory and physics-agnostic. It follows a passing
correctness contract and `GPU_NATIVE_SHAPE=pass`, then precedes optional
optimization. Copy `templates/EFFICIENCY.md`; record `pass`, `fail`, `blocked`,
or `not-run` for every item. A port may fail for a stated reason, but missing
evidence is not a pass.

## Profiling procedure

Record the immutable binary/configuration, backend and named GPU, dimensions,
mesh/box layout, MPI ranks, input, warm-up interval, steady-state interval, and
profiling command. Capture a device timeline plus an enclosing AMReX/TinyProfiler
region. Map the intended steady-state kernels to the scope and retain the trace
and summary beside the filled checklist. Link the field map, kernel graph, and
compiler resource report from `GPU_NATIVE_SHAPE.md`. Repeat after any change
claimed to affect this gate.

## Pass/fail checklist

- No hot host loop or quarantine executes in the steady-state path.
- Field state remains device-resident; there is no per-cell or per-tile host
  transfer.
- There is no per-tile device-wide/global synchronization; any synchronization
  has a documented dependency and sits at the narrowest valid boundary.
- Device-reachable calls use static/value dispatch; no runtime virtual/plugin
  dispatch occurs inside kernels.
- Equivalent work is not launched once per component when a component-aware
  launch preserves dependencies and semantics.
- Intended compute/library kernels dominate the steady-state device timeline;
  avoidable transfers, synchronization, and launch overhead do not dominate it.
- No hot per-cell/per-tile allocation occurs; arena choice and temporary
  lifetime match the kernel graph.
- Field layout, component ordering, kernel/MFIter granularity, transfer bytes,
  divergence/iteration, block size, registers/spills, occupancy, shared memory,
  and atomics all have dispositions in the mandatory shape profile.

All applicable items must pass for `BASELINE_EFFICIENCY=pass`. No speedup,
occupancy, or bandwidth threshold is implied.

## Optional optimization

GPU-015, GPU-018, GPU-019, and GPU-020 are optimization hypotheses, not baseline
port procedure. Attempt one only after the baseline passes and record exact
before/after commands, binary, GPU, dimensions, box layout, kernel/device time,
registers, occupancy, meaningful bandwidth, and preserved validation results.
Unattempted work goes in the per-port future-work section; it does not block a
GPU-native baseline.

The same rule applies to broader architecture hypotheses: AoS/SoA or component
reordering, kernel fusion/splitting or whole-level launch changes, transfer
elimination, block-size changes, register-pressure surgery, divergence
restructuring, shared-memory caching, atomic contention changes, and arena or
allocation reuse. The manual requires these decisions to be measured but does
not assume any one strategy wins. Record unattempted ideas in the per-port
future-work section rather than hiding them behind device annotations.
