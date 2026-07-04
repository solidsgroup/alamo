# Phase 2.1 Box/Grid Sweep

Elastic solve is parsed but skipped with `elastic.tstart=1000000000.0` so this isolates the phase-field/thermal path.

| Case | n_cell | max_level | blocking | max_grid | grid_eff | CPU np8 s | GPU s | GPU/CPU | Correctness |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| base_512_bf8_mgs32 | 64 64 64 | 3 | 8 | 32 | 0.7 | 3.833 | 3.748 | 0.978 | ok |
| wide_512_bf16_mgs64 | 64 64 64 | 3 | 16 | 64 | 0.9 | 2.610 | 2.795 | 1.071 | ok |
| wide_512_bf32_mgs128 | 64 64 64 | 3 | 32 | 128 | 0.9 | 1.418 | 1.835 | 1.294 | ok |
| wide_1024_bf32_mgs128 | 128 128 128 | 3 | 32 | 128 | 0.9 | 2.082 | 2.591 | 1.244 | ok |
| shallow_1024_bf32_mgs128 | 256 256 256 | 2 | 32 | 128 | 0.9 | 1.503 | 1.818 | 1.209 | ok |

CSV: `benchmark/phase2_box_sweep/summary.csv`
