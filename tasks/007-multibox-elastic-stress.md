# Task 007: Correctness Test C3 (Multi-box Elastic Stress Test)

## Goal
Create `tests/GPU/C3_multibox_elastic_stress/input` and
`tests/GPU/C3_multibox_elastic_stress/test.py`.

This test directly guards the `tmpfab.elixir()` fix in `Operator.cpp`. It forces
the exact conditions that previously caused the GPU MLMG divergence:
- **Multi-box decomposition** (`max_grid_size = 32` or `64`, so many small boxes)
- **Elastic MLMG solve** running for 200+ steps with many solves
- **GPU strict binary** (no fast-math, exactness matters)

**Pass criteria:**
- All 200 steps complete without crash or abort.
- `thermo.dat` values remain finite throughout.
- `chamber_pressure` stays within [0.05 MPa, 50 MPa] (physically sane range).
- No "MLMG failed to converge" or "diverged" messages in the log.

## Context
Project root: `/home/jackplum/Projects/alamo`
Plan: `docs/agent_plans/20260625-gpu-tests/PLAN.md`

Root-cause recap (from GPU_BRANCH_GUIDE.md): The elastic MLMG diverged when
`max_grid_size` forced multiple small boxes because `tmpfab` (a per-box temp in
`Operator::interpolation()`) was freed at iteration end and its device memory
reused by a different CUDA stream before the kernels finished. Fix: `amrex::Gpu::Elixir
tmpfab_eli = tmpfab.elixir()`. The test regresses that fix by:
- Using `max_grid_size = 32` (forces many boxes)
- Running 200 steps with `elastic.interval = 5` → 40 MLMG solves

## Files to read first
- `tests/SCPChamber/input`          — reference for elastic parameters
- `benchmark/GPU_BRANCH_GUIDE.md`   — elastic fix description (D1 section)
- `docs/agent_plans/20260625-gpu-tests/PLAN.md`

## Files allowed to modify
- `tests/GPU/C3_multibox_elastic_stress/input`   (CREATE)
- `tests/GPU/C3_multibox_elastic_stress/test.py` (CREATE)

## Files NOT allowed to modify
Everything else.

## Implementation steps

### C3 input (`tests/GPU/C3_multibox_elastic_stress/input`)

Base from `tests/SCPChamber/input`. Key parameters:

```
alamo.program = flame
plot_file = output

system.length = m
system.time = s

amr.plot_int = -1
amr.thermo.plot_int = -1
amr.thermo.int = 1
amr.max_level = 1
amr.n_cell = 8 8 8
amr.blocking_factor = 4
amr.max_grid_size = 32   # <-- CRITICAL: forces multi-box to stress the fix
amr.grid_eff = 0.7
amr.node.all = 1

geometry.prob_lo = 0.0 0.0 0.0
geometry.prob_hi = 0.001 0.001 0.001
geometry.is_periodic = 0 0 0

timestep = 1.0e-4
stop_time = 1e9

# Phase field (from SCPChamber)
pf.eps = 0.00008
pf.lambda = 0.001
pf.kappa = 1.0
pf.w1 = 1.0
pf.w12 = 2.0
pf.w0 = 0.0
small = 1E-4

pf.eta.ic.type = bmp
pf.eta.ic.bmp.filename = simple_circle.bmp
pf.eta.ic.bmp.fit = stretch
pf.eta.ic.bmp.channel = g

phi.ic.type = bmp
phi.ic.bmp.filename = base_circle0.bmp
phi.ic.bmp.fit = stretch
phi.ic.bmp.channel = g

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
```

For propellant, thermal, laser, chamber: copy from SCPChamber/input verbatim.

For elastic:
```
elastic.type = static
elastic.interval = 5
elastic.on = 1
elastic.print_model = 0
elastic.solver.fixed_iter = 200
elastic.solver.nriters = 200
elastic.solver.nrtolerance = 1e-5
elastic.solver.verbose = 0
elastic.traction = 4.0
```

Elastic BCs — copy from SCPChamber/input verbatim:
```
elastic.bc.type = constant
elastic.bc.constant.type.xhi = disp disp
elastic.bc.constant.type.xlo = disp disp
elastic.bc.constant.type.xloyhi = disp disp
elastic.bc.constant.type.xloylo = disp disp
elastic.bc.constant.type.ylo = disp disp
```

Model parameters — copy from SCPChamber/input verbatim:
```
model_ap.eps0 = 2.217e-5 0.0 0.0 0.0 2.217e-5 0.0 0.0 0.0 2.217e-5
model_ap.kappa = 150
model_ap.mu = 140
model_htpb.eps0 = 51e-7 0.0 0.0 0.0 51e-7 0.0 0.0 0.0 51e-7
model_htpb.kappa = 210
model_htpb.mu = 8
model_void.eps0 = 51e-7 0.0 0.0 0.0 51e-7 0.0 0.0 0.0 51e-7
model_void.kappa = 2000
model_void.mu = 8
```

NOTE: SCPChamber does NOT have system.length/time — the input uses raw SI values.
Follow that convention (no `_m`, `_K` unit suffixes in SCPChamber style; or add
`system.length = m` and use `_m` suffixes — be consistent within one input).
Look at what SCPChamber/input actually does and match it.

### C3 test.py (`tests/GPU/C3_multibox_elastic_stress/test.py`)

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

NUM_STEPS = 200  # 200 steps → ~40 elastic solves at interval=5

print(f"Running C3: {NUM_STEPS} steps, elastic every 5...")
rc, log, elapsed = testlib_gpu.run_alamo(
    gpu_bin, input_file, outdir,
    extra_args=["max_step={}".format(NUM_STEPS)],
    timeout=600)

print(log[-4000:])

if rc != 0:
    print(f"FAIL: alamo exited {rc}")
    sys.exit(1)

# Check for divergence signals in log
lower_log = log.lower()
diverge_signals = ["failed to converge", "diverged", "nan detected", "inf detected",
                   "mlmg failed", "elastic solver failed"]
for sig in diverge_signals:
    if sig in lower_log:
        print(f"FAIL: divergence signal found in log: '{sig}'")
        sys.exit(1)

thermo = testlib_gpu.parse_thermo(outdir / "plot" / "thermo.dat")
ok, issues = testlib_gpu.check_sanity(thermo)
if not ok:
    for i in issues:
        print("FAIL:", i)
    sys.exit(1)

# Check step count
steps_done = len(thermo.get("time", [])) - 1
if steps_done < NUM_STEPS * 0.95:   # allow 5% for any early-stop weirdness
    print(f"FAIL: only {steps_done}/{NUM_STEPS} steps completed")
    sys.exit(1)

# Check pressure stayed in physical range
pressure = thermo.get("chamber_pressure", [])
if pressure:
    min_p = min(pressure)
    max_p = max(pressure)
    if min_p < 5e4 or max_p > 5e7:   # 0.05 MPa to 50 MPa
        print(f"FAIL: chamber_pressure out of range [{min_p:.3e}, {max_p:.3e}] Pa")
        sys.exit(1)

print(f"PASS: C3 multi-box elastic stress ({elapsed:.1f}s, {steps_done} steps)")
print(f"  Elastic solves expected: {steps_done // 5}")
sys.exit(0)
```

## Invariants
- `max_grid_size = 32` is critical — do not increase it (the fix must be tested under multi-box)
- `elastic.type = static` (not disable)
- `elastic.interval = 5` with 200 steps → many solves
- Use GPU strict binary (no fast-math)
- BMP ICs are GPU-safe for 2D builds

## Build and test commands
```bash
cd /home/jackplum/Projects/alamo
python3 -m py_compile tests/GPU/C3_multibox_elastic_stress/test.py
python3 -c "
from pathlib import Path
t = Path('tests/GPU/C3_multibox_elastic_stress/input').read_text()
assert 'elastic.type = static' in t
assert 'max_grid_size = 32' in t
assert 'elastic.interval = 5' in t
print('C3 input OK')
"
```

## Expected result
- 2 files created.
- Input has max_grid_size=32, elastic.type=static, elastic.interval=5.
- test.py compiles cleanly.

## Non-goals
- Do not reduce NUM_STEPS — 200 steps with 40 elastic solves is intentionally stressful.

## Stop conditions
Stop if SCPChamber/input cannot be read (need to copy BMP IC paths from there).

## Final report
Write `results/007-RESULT.md` with: files created, critical parameters confirmed,
note about BMP IC paths (are they relative to project root?), syntax check results.
