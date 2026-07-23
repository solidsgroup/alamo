# Baseline efficiency record

Schema version: 1
Port ID: `<stable-id>`
Owner: `<name>`

## Reproduction

- Binary/configuration:
- Backend and named GPU:
- Dimensions and box layout:
- MPI ranks:
- Input, warm-up, and measured interval:
- Timeline command and evidence path:
- TinyProfiler command and evidence path:
- GPU-native shape profile, field map, kernel graph, and resource evidence:

## Checklist

| Item | Result (`not-run`, `blocked`, `fail`, `pass`) | Evidence/rationale |
|------|------------------------------------------------|--------------------|
| No hot host loop/quarantine | | |
| No per-cell/per-tile host transfer | | |
| No per-tile global synchronization | | |
| No runtime dispatch inside kernels | | |
| No avoidable per-component launch multiplication | | |
| Intended kernels dominate device timeline | | |
| No hot per-cell/per-tile allocation | | |
| `GPU_NATIVE_SHAPE=pass` | | |

`BASELINE_EFFICIENCY=`

## Deferred optimization

List hypotheses, owners, and future-work references. These do not change the
baseline result.
