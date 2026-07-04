# Task 006: Correctness Test C2 (Restart/Checkpoint Roundtrip)

## Goal
Create `tests/GPU/C2_restart_roundtrip/input` and `tests/GPU/C2_restart_roundtrip/test.py`.

Procedure:
1. **Run A**: run 20 steps with checkpointing every 10 steps → produces `chk000010`
2. **Run B**: restart from `chk000010`, run 10 more steps → produces final state at step 20
3. Compare final `thermo.dat` from Run A (step 20) vs Run B (step 10 from restart = step 20)

A clean restart should produce bit-for-bit identical (or near-identical) thermo output.
If the GPU serialize/deserialize path has bugs, this test catches them.

**Pass criteria:** final thermo row of Run A and Run B match within 1E-6 relative tolerance.

## Context
Project root: `/home/jackplum/Projects/alamo`
Plan: `docs/agent_plans/20260625-gpu-tests/PLAN.md`

Checkpoint parameters (AMReX standard, used in ALAMO):
- `amr.check_int = N`   — write checkpoint every N steps
- `amr.check_file = chk` — prefix for checkpoint directories (creates `chk000010/`)

Restart: pass `restart=<path>` on the command line (not in the input file).

The recon found (from FlowRestart): `cmdargs += " restart=" + restartfile`
So ALAMO accepts `restart=path/to/chk000010` as a command-line argument.

Use the GPU-strict binary for this test (correctness, not performance).

## Files to read first
- `tests/FlowRestart/input`     — reference for checkpoint/restart parameters
- `tests/GPU/testlib_gpu.py`    — run_alamo interface (extra_args)
- `docs/agent_plans/20260625-gpu-tests/PLAN.md`

## Files allowed to modify
- `tests/GPU/C2_restart_roundtrip/input`   (CREATE)
- `tests/GPU/C2_restart_roundtrip/test.py` (CREATE)

## Files NOT allowed to modify
Everything else.

## Implementation steps

### C2 input (`tests/GPU/C2_restart_roundtrip/input`)

Simple flame input WITHOUT elastic (restart with elastic adds complexity; keep focused).
Base it on the P1 perf input pattern but smaller:

```
alamo.program = flame
plot_file = output

system.length = m
system.time = s

amr.plot_int = -1
amr.thermo.plot_int = -1
amr.thermo.int = 1
amr.max_level = 0
amr.n_cell = 32 32 4
amr.blocking_factor = 8
amr.max_grid_size = 32
amr.grid_eff = 0.7
amr.node.all = 1

# Checkpoint configuration
amr.check_int = 10
amr.check_file = chk

geometry.prob_lo = 0.0_m 0.0_m 0.0_m
geometry.prob_hi = 0.04_m 0.04_m 0.0005_m
geometry.is_periodic = 0 0 0

timestep = 1.0e-5_s
stop_time = 1e9_s

pf.eps = 5.0e-5_m
pf.lambda = 0.001_J/m^2
pf.kappa = 1.0_J/m^2
pf.relax_steps = 0
pf.w1 = 1.0_1
pf.w12 = 2.0_1
pf.w0 = 0.0_1
small = 1E-4

pf.eta.ic.type = expression
pf.eta.ic.expression.constant.cx = 0.02
pf.eta.ic.expression.constant.cy = 0.02
pf.eta.ic.expression.constant.rb = 0.008
pf.eta.ic.expression.constant.w  = 0.001
pf.eta.ic.expression.region0 = "0.5 + 0.5*tanh((sqrt((x-cx)^2 + (y-cy)^2) - rb)/w)"

phi.ic.type = expression
phi.ic.expression.constant.cx = 0.02
phi.ic.expression.constant.cy = 0.02
phi.ic.expression.constant.rd = 0.016
phi.ic.expression.constant.w  = 0.001
phi.ic.expression.region0 = "0.5 + 0.5*tanh((rd - sqrt((x-cx)^2 + (y-cy)^2))/w)"

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

variable_pressure = 1
chamber.pressure = 0.101325_MPa
chamber.ballistic.At = 8.5e-4_m^2
chamber.ballistic.T0 = 3200_K
chamber.ballistic.R = 287_J/kg/K
chamber.ballistic.gamma = 1.25
chamber.ballistic.pressure = 4.0e6_Pa

elastic.type = disable
```

### C2 test.py (`tests/GPU/C2_restart_roundtrip/test.py`)

```python
#!/usr/bin/env python3
import os, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests/GPU"))
import testlib_gpu

outdir = Path(sys.argv[1])
input_file = Path(__file__).parent / "input"

gpu_bin = testlib_gpu.find_binary(ROOT, "bin/alamo_gpu-2d-nofast-cuda86-g++")
if gpu_bin is None:
    env = os.environ.get("ALAMO_GPU_STRICT_BIN")
    gpu_bin = Path(env) if env else None
if gpu_bin is None:
    print("SKIP: no GPU strict binary found")
    sys.exit(77)

run_a_out = outdir / "run_a"
run_b_out = outdir / "run_b"

# Run A: 20 steps, checkpoint at step 10
print("Run A: 20 steps with checkpoint...")
rc_a, log_a, _ = testlib_gpu.run_alamo(
    gpu_bin, input_file, run_a_out,
    extra_args=["max_step=20"],
    timeout=120)
if rc_a != 0:
    print("FAIL: Run A exited", rc_a)
    print(log_a[-2000:])
    sys.exit(1)

# Locate the checkpoint written at step 10
chk_dir = run_a_out / "chk000010"
if not chk_dir.is_dir():
    # AMReX zero-pads to 6 digits
    candidates = sorted(run_a_out.glob("chk*"))
    if not candidates:
        print("FAIL: no checkpoint directory found in", run_a_out)
        sys.exit(1)
    chk_dir = candidates[0]
    print(f"  Using checkpoint: {chk_dir.name}")

# Run B: restart from step 10, run 10 more steps
print(f"Run B: restart from {chk_dir}, 10 more steps...")
rc_b, log_b, _ = testlib_gpu.run_alamo(
    gpu_bin, input_file, run_b_out,
    extra_args=["max_step=20", f"restart={chk_dir}"],
    timeout=120)
if rc_b != 0:
    print("FAIL: Run B exited", rc_b)
    print(log_b[-2000:])
    sys.exit(1)

# Compare final thermo rows
thermo_a = testlib_gpu.parse_thermo(run_a_out / "plot" / "thermo.dat")
thermo_b = testlib_gpu.parse_thermo(run_b_out / "plot" / "thermo.dat")

if not thermo_a or not thermo_b:
    print("FAIL: missing thermo.dat from one or both runs")
    sys.exit(1)

# Trim thermo_b to only rows after the restart point (time > checkpoint time)
# Both should end at the same simulation time; compare last rows only
def last_row(thermo):
    return {k: [v[-1]] for k, v in thermo.items() if v}

final_a = last_row(thermo_a)
final_b = last_row(thermo_b)

ok, issues = testlib_gpu.compare_thermo(final_a, final_b, rel_tol=1e-6, abs_tol=1e-12)
if not ok:
    print("FAIL: restart final-row mismatch:")
    for i in issues:
        print(" ", i)
    sys.exit(1)

print("PASS: C2 restart roundtrip")
print(f"  Run A steps: {len(thermo_a.get('time',[]))-1}")
print(f"  Run B steps: {len(thermo_b.get('time',[]))-1}")
sys.exit(0)
```

## Invariants
- `amr.check_int = 10` and `amr.check_file = chk` must be in the input.
- The restart path passed to run_alamo must be an absolute or relative-from-cwd path.
  `testlib_gpu.run_alamo` runs from the project root, so `run_a_out / "chk000010"` as
  an absolute path should work.
- The restart run still uses the same input file (IC parameters are ignored when
  restarting; AMReX restores state from the checkpoint).
- `max_step=20` applies to both runs: Run A stops at 20; Run B restarts from step 10
  and runs until step 20 (10 more steps).

## Build and test commands
```bash
cd /home/jackplum/Projects/alamo
python3 -m py_compile tests/GPU/C2_restart_roundtrip/test.py
python3 -c "
from pathlib import Path
t = Path('tests/GPU/C2_restart_roundtrip/input').read_text()
assert 'amr.check_int = 10' in t
assert 'amr.check_file = chk' in t
assert 'elastic.type = disable' in t
print('C2 input OK')
"
```

## Expected result
- 2 files created.
- Input has `amr.check_int = 10`, `amr.check_file = chk`.
- test.py compiles cleanly.

## Non-goals
- Do not test elastic restart (separate concern from checkpoint correctness).
- Do not run the test.

## Stop conditions
Stop if FlowRestart/input cannot be read (use it for reference on checkpoint format).

## Final report
Write `results/006-RESULT.md` with: files created, checkpoint parameters confirmed,
syntax check results.
