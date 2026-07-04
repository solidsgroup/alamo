# Task 005: Correctness Test C1 (CPU vs GPU Parity with Elastic)

## Goal
Create `tests/GPU/C1_correctness_elastic/input` and `tests/GPU/C1_correctness_elastic/test.py`.

Run the SAME elastic-enabled input on:
1. CPU binary (`alamo-2d-g++`) — reference
2. GPU strict binary (`alamo_gpu-2d-nofast-cuda86-g++`) — candidate

Compare `thermo.dat` outputs within tolerance. This guards against regressions of the
elastic cross-stream UAF fix (2026-06-25).

**Pass criteria:** both runs complete, thermo outputs match within 5% relative tolerance
on all columns (except `time`, which should be bit-for-bit identical).

## Context
Project root: `/home/jackplum/Projects/alamo`
Plan: `docs/agent_plans/20260625-gpu-tests/PLAN.md`

The elastic fix root cause was a GPU cross-stream use-after-free in `interpolation()`.
The fix (`tmpfab.elixir()`) defers freeing until the stream completes. This test
runs enough elastic solves (at least 5) to detect if the race recurs.

Key constraints for GPU-safe elastic input:
- ICs: `Constant` or `Expression` (scalar) or `BMP` — GPU-safe
- Elastic BCs: `disp` type (maps to `BC::Operator::Elastic::Constant`) — GPU-safe
- `elastic.type = static`
- Small grid for speed, but enough resolution to stress the MLMG solver

## Files to read first
- `tests/SCPChamber/input`            — reference elastic input (has elastic.on=1)
- `benchmark/elastic_sensitivity_20260621/fix_notes.md`  (if readable) — context
- `docs/agent_plans/20260625-gpu-tests/PLAN.md`

## Files allowed to modify
- `tests/GPU/C1_correctness_elastic/input`   (CREATE)
- `tests/GPU/C1_correctness_elastic/test.py` (CREATE)

## Files NOT allowed to modify
Everything else.

## Implementation steps

### C1 input (`tests/GPU/C1_correctness_elastic/input`)

Base from `tests/SCPChamber/input`. Modifications:
- Remove all `#@` metadata lines
- `plot_file = output`
- `amr.n_cell = 8 8 8`  (small, fast, exercises MLMG with multiple boxes)
- `amr.max_level = 2`
- `amr.blocking_factor = 4`
- `amr.max_grid_size = 64`   (forces multi-box — this is the critical condition)
- `amr.plot_int = -1`
- `amr.thermo.int = 1`
- `elastic.type = static`
- `elastic.interval = 5`    (elastic solve every 5 steps → 6 solves in 30-step run)
- `elastic.solver.fixed_iter = 100`
- `elastic.solver.nriters = 200`
- `elastic.solver.nrtolerance = 1e-5`
- `elastic.solver.verbose = 0`
- `elastic.print_model = 0`
- `elastic.print_residual = 0` (if this param exists)
- `variable_pressure = 0`
- `chamber.pressure = 4.0_MPa`
- Keep all elastic.bc.*, model_ap.*, model_htpb.*, model_void.* from SCPChamber
- Keep pf.*, thermal.*, temp.ic.*, laser.ic.* from SCPChamber
- `amr.thermo.plot_int = 1`  (write thermo every step for comparison)

### C1 test.py (`tests/GPU/C1_correctness_elastic/test.py`)

```python
#!/usr/bin/env python3
import os, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests/GPU"))
import testlib_gpu

outdir = Path(sys.argv[1])
input_file = Path(__file__).parent / "input"

# Find binaries
cpu_bin = testlib_gpu.find_binary(ROOT, "bin/alamo-2d-g++")
if cpu_bin is None:
    env = os.environ.get("ALAMO_CPU_BIN")
    cpu_bin = Path(env) if env else None

gpu_bin = testlib_gpu.find_binary(ROOT, "bin/alamo_gpu-2d-nofast-cuda86-g++")
if gpu_bin is None:
    env = os.environ.get("ALAMO_GPU_STRICT_BIN")
    gpu_bin = Path(env) if env else None

if cpu_bin is None:
    print("SKIP: CPU binary not found (bin/alamo-2d-g++)")
    sys.exit(77)
if gpu_bin is None:
    print("SKIP: GPU strict binary not found (bin/alamo_gpu-2d-nofast-cuda86-g++)")
    sys.exit(77)

NUM_STEPS = 30   # 30 steps → 6 elastic solves at interval=5

cpu_out = outdir / "cpu"
gpu_out = outdir / "gpu"

# Run CPU reference
print("Running CPU reference...")
rc_cpu, log_cpu, elapsed_cpu = testlib_gpu.run_alamo(
    cpu_bin, input_file, cpu_out,
    extra_args=["max_step={}".format(NUM_STEPS)],
    timeout=300)
if rc_cpu != 0:
    print("FAIL: CPU run exited", rc_cpu)
    print(log_cpu[-2000:])
    sys.exit(1)

# Run GPU strict
print("Running GPU strict...")
rc_gpu, log_gpu, elapsed_gpu = testlib_gpu.run_alamo(
    gpu_bin, input_file, gpu_out,
    extra_args=["max_step={}".format(NUM_STEPS)],
    timeout=300)
if rc_gpu != 0:
    print("FAIL: GPU run exited", rc_gpu)
    print(log_gpu[-2000:])
    sys.exit(1)

# Compare thermo
cpu_thermo = testlib_gpu.parse_thermo(cpu_out / "plot" / "thermo.dat")
gpu_thermo = testlib_gpu.parse_thermo(gpu_out / "plot" / "thermo.dat")

ok_cpu, cpu_issues = testlib_gpu.check_sanity(cpu_thermo)
ok_gpu, gpu_issues = testlib_gpu.check_sanity(gpu_thermo)

if not ok_cpu:
    print("FAIL: CPU thermo insane:", cpu_issues)
    sys.exit(1)
if not ok_gpu:
    print("FAIL: GPU thermo insane:", gpu_issues)
    sys.exit(1)

ok, diff_issues = testlib_gpu.compare_thermo(cpu_thermo, gpu_thermo,
                                              rel_tol=0.05, abs_tol=1e-10)
if not ok:
    print("FAIL: CPU/GPU thermo mismatch:")
    for i in diff_issues:
        print(" ", i)
    sys.exit(1)

print(f"PASS: C1 correctness elastic parity")
print(f"  CPU elapsed: {elapsed_cpu:.1f}s, GPU elapsed: {elapsed_gpu:.1f}s")
print(f"  Rows compared: {len(cpu_thermo.get('time',[]))}")
sys.exit(0)
```

## Invariants
- `max_grid_size = 64` with `n_cell = 8 8 8` and `max_level = 2` creates
  effective grid 32×32 at level 2. With blocking_factor=4 and max_grid_size=64,
  this should produce multiple boxes — which is the critical condition.
- Both runs use the SAME input file (not copies), ensuring parameter parity.
- GPU run uses `nofast` binary for exact comparison.
- `elastic.interval = 5` with `NUM_STEPS = 30` → at least 6 elastic solves.

## Build and test commands
```bash
cd /home/jackplum/Projects/alamo
python3 -m py_compile tests/GPU/C1_correctness_elastic/test.py
python3 -c "
from pathlib import Path
t = Path('tests/GPU/C1_correctness_elastic/input').read_text()
assert 'elastic.type = static' in t
assert 'elastic.interval = 5' in t
assert 'max_grid_size = 64' in t
print('C1 input OK')
"
```

## Expected result
- 2 files created.
- Input has elastic.type=static, max_grid_size=64, elastic.interval=5.
- test.py compiles cleanly.

## Non-goals
- Do not compare full field data (thermo comparison is sufficient for this test).
- Do not run the test.

## Stop conditions
Stop if SCPChamber/input cannot be read.

## Final report
Write `results/005-RESULT.md` with: files created, key parameter values confirmed,
syntax check results.
