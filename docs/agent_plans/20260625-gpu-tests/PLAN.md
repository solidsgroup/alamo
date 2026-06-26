# Plan: GPU Test Suite (10 tests)
Date: 2026-06-25

## User goal
Implement 10 GPU tests for the ALAMO flame simulator covering functionality (F1/F2),
performance (P1/P2/P3), and correctness regression (C1–C4).

## Current architecture
- `tests/` — test cases with `input` (AMReX ParmParse format, `#@` metadata) and `test.py`
- `scripts/runtests.py` — discovers `#@` blocks, builds binary name from `exe`+`dim`+compiler
- `scripts/testlib.py` — shared Python helpers (validate, readContours, etc.)
- `benchmark/baseline_suite.py` — independent GPU/CPU runner with thermo comparison
- Binaries (under `bin/`):
  - `alamo_gpu-2d-cuda86-g++`         — GPU fast (--use_fast_math), 2D
  - `alamo_gpu-2d-nofast-cuda86-g++`  — GPU strict (no fast-math), 2D
  - `alamo_gpu-3d-cuda86-g++`         — GPU fast, 3D
  - `alamo-2d-g++`                    — CPU reference, 2D
- Thermo output: tab-separated `thermo.dat` at `<plot_file>/thermo.dat`
- Elastic fix (2026-06-25): `tmpfab.elixir()` in `Operator.cpp` fixes GPU cross-stream UAF

## Desired architecture
New `tests/GPU/` directory:
```
tests/GPU/
├── run_gpu_tests.py              # CLI driver (auto-discovers test dirs)
├── testlib_gpu.py                # shared utilities
├── F1_smoke_flame_only/          # input + test.py
├── F2_smoke_elastic/
├── P1_perf_2d_hiRes_noAMR/
├── P2_perf_2d_hiRes_AMR3/
├── P3_perf_3d_256/
├── C1_correctness_elastic/
├── C2_restart_roundtrip/
├── C3_multibox_elastic_stress/
└── C4_amr_correctness/
```

Each `test.py` receives `outdir` as argv[1] and exits 0 on pass, nonzero on fail.
Driver `run_gpu_tests.py` handles binary selection, timing, and summary reporting.

## Invariants and constraints
1. GPU-safe ICs only: `IC::Constant`, `IC::Expression` (scalar), `IC::BMP` are device-safe.
   `IC::PSRead`, `IC::PNG`, `IC::Trig` are guarded/abort on GPU — do NOT use.
2. GPU-safe BCs: `BC::Constant` is device-safe. `BC::Operator::Elastic::Expression` is not.
3. For elastic tests use `elastic.type = static` (not disable) and `BC::Constant` elastic BCs.
4. With `elastic.type = disable` the elastic model/BC blocks must be OMITTED — strict parser aborts.
5. `allow_unused=1` on the command line suppresses unused-parameter errors (use for overrides).
6. In 2D builds `amr.n_cell` takes 3 values (z is ignored); `geometry.is_periodic` similarly 3.
7. Checkpoint write: `amr.check_int = N`, `amr.check_file = chk`. Restart: `restart=chkNNNNNN`.
8. `plot_file` must be a relative path inside the test's output dir (passed by driver).
9. Workers must NOT modify any file under `src/`, `benchmark/`, `scripts/`, or `tests/<existing>`.
10. 3D binary suffix: `alamo_gpu-3d-cuda86-g++`. 2D binary: `alamo_gpu-2d-cuda86-g++`.

## Files involved
### Created (new)
- `tests/GPU/run_gpu_tests.py`
- `tests/GPU/testlib_gpu.py`
- `tests/GPU/F1_smoke_flame_only/input`
- `tests/GPU/F1_smoke_flame_only/test.py`
- `tests/GPU/F2_smoke_elastic/input`
- `tests/GPU/F2_smoke_elastic/test.py`
- `tests/GPU/P1_perf_2d_hiRes_noAMR/input`
- `tests/GPU/P1_perf_2d_hiRes_noAMR/test.py`
- `tests/GPU/P2_perf_2d_hiRes_AMR3/input`
- `tests/GPU/P2_perf_2d_hiRes_AMR3/test.py`
- `tests/GPU/P3_perf_3d_256/input`
- `tests/GPU/P3_perf_3d_256/test.py`
- `tests/GPU/C1_correctness_elastic/input`
- `tests/GPU/C1_correctness_elastic/test.py`
- `tests/GPU/C2_restart_roundtrip/input`
- `tests/GPU/C2_restart_roundtrip/test.py`
- `tests/GPU/C3_multibox_elastic_stress/input`
- `tests/GPU/C3_multibox_elastic_stress/test.py`
- `tests/GPU/C4_amr_correctness/input`
- `tests/GPU/C4_amr_correctness/test.py`
- `.claude/agents/refactor-worker.md`

### Read (existing, never modified)
- `tests/SCPSandwich/input`
- `tests/SCPChamber/input`
- `tests/SCPSpheresElastic/input`
- `tests/SCPThermalSandwich/input`
- `tests/SCPThermalVoid/input`
- `input_3d_flame`
- `scripts/testlib.py`
- `benchmark/baseline_suite.py`
- `scripts/runtests.py`

## Build and test commands
```bash
# Run all GPU tests (requires CUDA GPU):
cd /home/jackplum/Projects/alamo
python3 tests/GPU/run_gpu_tests.py

# Run a single test:
python3 tests/GPU/run_gpu_tests.py --test F1_smoke_flame_only

# Dry-run (shows what would run, no GPU needed):
python3 tests/GPU/run_gpu_tests.py --dry-run
```

## Task graph
All tasks are PARALLEL-SAFE — each owns different files with no overlap.

| Task | Slug                      | Files owned                                          | Status          |
|------|---------------------------|------------------------------------------------------|-----------------|
| 001  | gpu-test-framework        | testlib_gpu.py, run_gpu_tests.py                     | parallel-safe   |
| 002  | smoke-tests               | F1_smoke_flame_only/*, F2_smoke_elastic/*             | parallel-safe   |
| 003  | perf-2d-tests             | P1_perf_2d_hiRes_noAMR/*, P2_perf_2d_hiRes_AMR3/*   | parallel-safe   |
| 004  | perf-3d-test              | P3_perf_3d_256/*                                     | parallel-safe   |
| 005  | correctness-elastic       | C1_correctness_elastic/*                             | parallel-safe   |
| 006  | restart-roundtrip         | C2_restart_roundtrip/*                               | parallel-safe   |
| 007  | multibox-elastic-stress   | C3_multibox_elastic_stress/*                         | parallel-safe   |
| 008  | amr-correctness           | C4_amr_correctness/*                                 | parallel-safe   |

## Parallelization strategy
All 8 tasks can be dispatched simultaneously — they write to disjoint subdirectories
under `tests/GPU/`. The driver (Task 001) auto-discovers test dirs at runtime, so
it does not depend on Tasks 002-008 completing first.

## Integration strategy
After all workers complete, read RESULT.md files and verify:
1. All task directories exist under tests/GPU/
2. Each input file has the required keys (alamo.program, timestep, etc.)
3. run_gpu_tests.py --dry-run exits 0
4. If a GPU is available, run the two cheapest smoke tests as a final check

## Risks
- Input parameter errors (wrong units, missing required field) — workers must validate
  by reading the reference inputs carefully and following the same format
- GPU not available at test-write time — write tests so they report SKIP gracefully
- Restart checkpoint path assumptions — use a relative path inside outdir

## Rollback plan
All changes are new files under `tests/GPU/` and `docs/agent_plans/`. Nothing existing
is modified. Delete `tests/GPU/` to roll back completely.

---

## Post-run fixes (2026-06-26) — 5P/4F → 9P/0F

First execution of this suite against the current parser was **5 passed /
4 failed** (snapshot in `docs/llm/changelog/2026-06-26-3d-elastic-gpu-fix.md`).
The 4 failures were analyzed and fixed — full writeup in
**`benchmark/GPU_TEST_SUITE_FIXES.md`**, perf baseline in
**`benchmark/GPU_TEST_PERF_TRACKING.md`**, changelog
`docs/llm/changelog/2026-06-26-gpu-test-suite-fixes.md`.

Key correction to several **assumptions/risks above**:
- "Restart checkpoint path assumptions" (Risk) was the big one: there is **no**
  `amr.check_int`/`amr.check_file` in ALAMO. Checkpoints are embedded in
  plotfile dirs (`plot/<step>{cell,node}/Checkpoint`), gated by `amr.plot_int`;
  restart needs `restart_cell=`+`restart_node=` (Flame has a nodal fab). The C2
  deck/test were rewritten accordingly, and a **real `Integrator::Restart`
  node-fab out-of-bounds segfault** + a **headerless-restart-`thermo.dat`** bug
  were fixed in `src/Integrator/Integrator.cpp`.
- Invariant 9 ("Workers must NOT modify any file under `src/`…") held for the
  test-authoring tasks, but closing out the suite **required** two `src/` fixes
  (the restart defects above) — those are genuine framework bugs the tests
  surfaced, not test-harness changes.
- The decks (F2/C1/C3) used the **old** Homogenize parameter schema and had
  never actually run; they were re-authored to the validated single-`*_prop`
  physics, with dt lowered to 1e-5_s for thermal CFL on the 0.001 m grid and
  the elastic `model_casing` softened to 500 MPa to keep MLMG conditioned.
