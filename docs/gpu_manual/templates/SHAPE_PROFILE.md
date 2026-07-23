# GPU-native shape profile

Schema version: 1
Port ID: `<stable-id>`
Owner: `<name>`

## Reproduction

- Source revision and binary/configuration:
- Backend, compiler, and named GPU:
- Dimensions, mesh/box layout, and MPI ranks:
- Input, warm-up, and steady-state interval:
- Compiler resource-report command/evidence:
- Device timeline command/evidence:
- Linked field map and kernel graph:
- Linked correctness validation result:

## Per-kernel resources

| Kernel ID | Block size | Registers/thread | Local spill | Shared memory | Occupancy/active warps | Bandwidth | Limiter/rationale | Result |
|-----------|------------|------------------|-------------|---------------|------------------------|-----------|-------------------|--------|
| | | | | | | | | `open` |

## Architecture decisions

| Topic | Decision and evidence | Result (`pass`, `fail`, `blocked`, `not-applicable`) | Owner |
|-------|-----------------------|---------------------------------------------------|-------|
| Data layout/component order | | | |
| Kernel/MFIter granularity | | | |
| Residency/transfers per step | | | |
| Iteration/divergence | | | |
| Shared memory | | | |
| Atomics/reductions | | | |
| Arena/allocation/lifetime | | | |

`GPU_NATIVE_SHAPE=`

