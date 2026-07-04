# Task 003: Performance Tests P1 and P2

## Goal
Create input files and test scripts for:
- **P1** `tests/GPU/P1_perf_2d_hiRes_noAMR/` — 2D, 256×256 base, max_level=0, no AMR
- **P2** `tests/GPU/P2_perf_2d_hiRes_AMR3/`  — 2D, 64×64 base, max_level=3 (eff. 512 cells)

Both run for 300 steps and report wall-clock ms/step. There is no pass/fail threshold —
these are measurement tests. The test.py exits 0 as long as the run completes without crash.

## Context
Project root: `/home/jackplum/Projects/alamo`
Plan: `docs/agent_plans/20260625-gpu-tests/PLAN.md`

Performance tests use the GPU-fast binary (`alamo_gpu-2d-cuda86-g++`) for meaningful
throughput numbers. `--use_fast_math` is fine here since we care about speed not exact match.

The 300-step run should take ~1-5 minutes on an A1000 GPU. Thermo-only output
(no plot dumps) keeps I/O from distorting the timing.

**Key metric:** `ms/step = elapsed_seconds * 1000 / steps_completed`

For AMR (P2), the effective finest-level cell count grows during the run. The
`ms/step` still measures wall time per timestep, which is what we care about.

## Files to read first
- `tests/SCPSandwich/input`   — reference for a typical 2D flame input structure
- `input_3d_flame`            — reference for large-grid GPU-optimized settings
  (blocking_factor=32, max_grid_size=128, grid_eff=0.7, amr.node.all=1 etc.)
- `benchmark/phase2_box_sweep.py` — see COMMON_OVERRIDES for useful suppression flags
- `docs/agent_plans/20260625-gpu-tests/PLAN.md`

## Files allowed to modify
- `tests/GPU/P1_perf_2d_hiRes_noAMR/input`   (CREATE)
- `tests/GPU/P1_perf_2d_hiRes_noAMR/test.py` (CREATE)
- `tests/GPU/P2_perf_2d_hiRes_AMR3/input`    (CREATE)
- `tests/GPU/P2_perf_2d_hiRes_AMR3/test.py`  (CREATE)

## Files NOT allowed to modify
Everything else.

## Implementation steps

### P1 input (`tests/GPU/P1_perf_2d_hiRes_noAMR/input`)

Goal: large 2D grid, no AMR, GPU-optimized box layout.

```
alamo.program = flame
plot_file = output

system.length = m
system.time = s

# AMR — no refinement
amr.plot_int = -1
amr.thermo.plot_int = -1
amr.thermo.int = 1
amr.max_level = 0
amr.n_cell = 256 256 4
amr.blocking_factor = 32
amr.max_grid_size = 128
amr.grid_eff = 0.7
amr.node.all = 1

# Geometry: 4cm × 4cm × 0.5mm (thin in z for 2D build)
geometry.prob_lo = 0.0_m 0.0_m 0.0_m
geometry.prob_hi = 0.04_m 0.04_m 0.0005_m
geometry.is_periodic = 0 0 0

timestep = 1.0e-5_s
stop_time = 1e9_s    # controlled by max_step override

# Phase field
pf.eps = 5.0e-5_m
pf.lambda = 0.001_J/m^2
pf.kappa = 1.0_J/m^2
pf.relax_steps = 0
pf.w1 = 1.0_1
pf.w12 = 2.0_1
pf.w0 = 0.0_1
small = 1E-4

# IC: central bore (eta=0 inside circle, eta=1 solid)
pf.eta.ic.type = expression
pf.eta.ic.expression.constant.cx = 0.02
pf.eta.ic.expression.constant.cy = 0.02
pf.eta.ic.expression.constant.rb = 0.008
pf.eta.ic.expression.constant.w  = 0.001
pf.eta.ic.expression.region0 = "0.5 + 0.5*tanh((sqrt((x-cx)^2 + (y-cy)^2) - rb)/w)"

# phi IC: propellant disk
phi.ic.type = expression
phi.ic.expression.constant.cx = 0.02
phi.ic.expression.constant.cy = 0.02
phi.ic.expression.constant.rd = 0.016
phi.ic.expression.constant.w  = 0.001
phi.ic.expression.region0 = "0.5 + 0.5*tanh((rd - sqrt((x-cx)^2 + (y-cy)^2))/w)"

# eta BCs (all Dirichlet = solid wall)
pf.eta.bc.type = constant
pf.eta.bc.constant.type.xlo = dirichlet
pf.eta.bc.constant.type.xhi = dirichlet
pf.eta.bc.constant.type.ylo = dirichlet
pf.eta.bc.constant.type.yhi = dirichlet
pf.eta.bc.constant.type.zlo = dirichlet
pf.eta.bc.constant.type.zhi = dirichlet
pf.eta.bc.constant.val.xlo = 1.0
pf.eta.bc.constant.val.xhi = 1.0
pf.eta.bc.constant.val.ylo = 1.0
pf.eta.bc.constant.val.yhi = 1.0
pf.eta.bc.constant.val.zlo = 1.0
pf.eta.bc.constant.val.zhi = 1.0

# Propellant model
propellant.type = homogenize
propellant.homogenize.dispersion1 = 0.025_W/m/K
propellant.homogenize.dispersion2 = 1.2_kg/m^3
propellant.homogenize.dispersion3 = 1000_J/kg/K
propellant.homogenize.h1 = 1.5e6_W/m^2
propellant.homogenize.h2 = 8.0e6_W/m^2
propellant.homogenize.P_reference = 1_MPa
propellant.homogenize.rho_prop = 1745.0_kg/m^3
propellant.homogenize.k_prop = 1.5_W/m/K
propellant.homogenize.cp_prop = 1476.0_J/kg/K
propellant.homogenize.m_prop = 8000_cm/s
propellant.homogenize.mlocal_prop = 0.1_kg/m^2/s
propellant.homogenize.E_prop = 3500.0_K
propellant.homogenize.mob_prop = true
propellant.homogenize.pressure_exponent = 0.372531
propellant.homogenize.bound = 0_K
propellant.homogenize.bound_width = 50_K

# Thermal
thermal.on = 1
thermal.Tref = 300
thermal.Tfluid = 300
thermal.temp.bc.type = constant
thermal.temp.bc.constant.type.xlo = dirichlet
thermal.temp.bc.constant.type.xhi = dirichlet
thermal.temp.bc.constant.type.ylo = dirichlet
thermal.temp.bc.constant.type.yhi = dirichlet
thermal.temp.bc.constant.type.zlo = dirichlet
thermal.temp.bc.constant.type.zhi = dirichlet
thermal.temp.bc.constant.val.xlo = 300_K
thermal.temp.bc.constant.val.xhi = 300_K
thermal.temp.bc.constant.val.ylo = 300_K
thermal.temp.bc.constant.val.yhi = 300_K
thermal.temp.bc.constant.val.zlo = 300_K
thermal.temp.bc.constant.val.zhi = 300_K

temp.ic.type = constant
temp.ic.constant.value = 300_K

laser.ic.type = expression
laser.ic.expression.region0 = "(t < 0.025) * (150000000)"
thermal.hc = 1.0

# Pressure
variable_pressure = 1
chamber.pressure = 0.101325_MPa
chamber.ballistic.At = 8.5e-4_m^2
chamber.ballistic.T0 = 3200_K
chamber.ballistic.R = 287_J/kg/K
chamber.ballistic.gamma = 1.25
chamber.ballistic.pressure = 4.0e6_Pa

# Elastic DISABLED
elastic.type = disable
```

### P2 input (`tests/GPU/P2_perf_2d_hiRes_AMR3/input`)

Copy P1 input and change only:
- `amr.max_level = 3`
- `amr.n_cell = 64 64 4`  (start coarse; effective resolution with L3 = 512 cells)
- `amr.base_regrid_int = 10`
- `amr.regrid_int = 10`
- `amr.refinement_criterion = 0.1`
- `amr.refinement_criterion_temp = 10.0_K`
- `amr.phi_refinement_criterion = 0.5`
- Everything else identical to P1.

### P1 test.py and P2 test.py

Same structure for both (adjust test name in PASS message):

```python
#!/usr/bin/env python3
import os, sys, time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests/GPU"))
import testlib_gpu

outdir = Path(sys.argv[1])
input_file = Path(__file__).parent / "input"
test_name = Path(__file__).parent.name

binary = testlib_gpu.find_binary(ROOT, "bin/alamo_gpu-2d-cuda86-g++")
if binary is None:
    env_bin = os.environ.get("ALAMO_GPU_BIN")
    binary = Path(env_bin) if env_bin else None
if binary is None:
    print("SKIP: no GPU binary found")
    sys.exit(77)

NUM_STEPS = 300

rc, log, elapsed = testlib_gpu.run_alamo(
    binary, input_file, outdir,
    extra_args=["max_step={}".format(NUM_STEPS), "allow_unused=1"],
    timeout=600)

print(log[-3000:])

if rc != 0:
    print(f"FAIL: alamo exited {rc}")
    sys.exit(1)

thermo = testlib_gpu.parse_thermo(outdir / "plot" / "thermo.dat")
ok, issues = testlib_gpu.check_sanity(thermo)
if not ok:
    for i in issues:
        print("FAIL:", i)
    sys.exit(1)

steps_done = len(thermo.get("time", [])) - 1
ms_per_step = testlib_gpu.wall_per_step(elapsed, thermo)

print(f"PASS: {test_name}")
print(f"  Steps completed : {steps_done}")
print(f"  Elapsed         : {elapsed:.1f}s")
print(f"  ms/step         : {ms_per_step:.2f}" if ms_per_step else "  ms/step: N/A")
sys.exit(0)
```

## Invariants
- No elastic blocks in either input (elastic.type = disable, all model/BC blocks omitted).
- `amr.plot_int = -1` and `amr.thermo.plot_int = -1` — suppress field dumps, keep thermo.
- Expression IC must be scalar (single `region0`, no vector component syntax) — GPU-safe.
- Use `allow_unused=1` in extra_args since perf runs may skip unused parameters.

## Build and test commands
```bash
cd /home/jackplum/Projects/alamo
python3 -m py_compile tests/GPU/P1_perf_2d_hiRes_noAMR/test.py
python3 -m py_compile tests/GPU/P2_perf_2d_hiRes_AMR3/test.py
# Check input files parse:
python3 -c "
from pathlib import Path
for p in ['tests/GPU/P1_perf_2d_hiRes_noAMR/input', 'tests/GPU/P2_perf_2d_hiRes_AMR3/input']:
    t = Path(p).read_text()
    assert 'elastic.type = disable' in t, f'{p}: missing elastic.type=disable'
    assert 'alamo.program = flame' in t, f'{p}: missing program'
    print(p, 'OK')
"
```

## Expected result
- 4 files created, all syntactically valid.
- P1 has max_level=0, P2 has max_level=3.
- Both have elastic.type=disable with no model/BC elastic blocks.

## Non-goals
- Do not add a perf pass/fail threshold (timing varies by hardware).
- Do not run the tests.

## Stop conditions
Stop if you cannot read `tests/SCPSandwich/input` or `input_3d_flame`.

## Final report
Write `results/003-RESULT.md` with: files created, key parameter values confirmed,
syntax check results.
