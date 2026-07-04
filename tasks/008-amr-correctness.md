# Task 008: Correctness Test C4 (GPU vs CPU AMR Correctness)

## Goal
Create `tests/GPU/C4_amr_correctness/input` and `tests/GPU/C4_amr_correctness/test.py`.

Run the same AMR-enabled input on CPU and GPU strict. Compare the eta and temp
fields at a cross-section using the existing `scripts/testlib.py` infrastructure.

AMR grid transfers (`interpolation()`, `restriction()`) were the exact site of the
elastic UAF fix. Even for flame-only (no elastic), AMR-path correctness on GPU
should be explicitly verified.

**Pass criteria:**
- Both runs complete without crash.
- `eta` field along a center cross-section: GPU vs CPU relative error < 5%.
- `temp` field along the same cross-section: GPU vs CPU relative error < 5%.

## Context
Project root: `/home/jackplum/Projects/alamo`
Plan: `docs/agent_plans/20260625-gpu-tests/PLAN.md`

`scripts/testlib.py` provides `validate()` and `readContours()` which load AMReX
plot directories using `yt` and compare against a reference. We adapt this for
a GPU vs CPU comparison by running CPU first (as reference) and GPU second.

The test uses a medium AMR case (max_level=2) with a non-trivial IC so AMR
actually triggers refinement.

## Files to read first
- `tests/SCPThermalSandwich/input`   — reference for a medium AMR 2D flame input
- `scripts/testlib.py`               — readContours, validate API
- `tests/SCPThermalVoid/test`        — example of how test.py uses testlib.validate
- `docs/agent_plans/20260625-gpu-tests/PLAN.md`

## Files allowed to modify
- `tests/GPU/C4_amr_correctness/input`   (CREATE)
- `tests/GPU/C4_amr_correctness/test.py` (CREATE)

## Files NOT allowed to modify
Everything else.

## Implementation steps

### C4 input (`tests/GPU/C4_amr_correctness/input`)

Use a simplified version of SCPThermalSandwich with:
- `amr.max_level = 2`  (enough to exercise AMR grid transfer paths)
- `amr.n_cell = 64 8 8`  (elongated for a sandwich-like geometry)
- `amr.blocking_factor = 8`
- `amr.max_grid_size = 32`
- `elastic.type = disable`  (no elastic — focus on AMR path correctness)
- Small domain, short run: 50 steps total

Model the input on SCPThermalSandwich closely (same propellant, thermal, pf params)
but reduce the grid and drop elastic. Key parameters to copy from SCPThermalSandwich:
- `system.length = mm`
- geometry: `prob_lo = 0.0_mm -0.25_mm -0.1_mm`, `prob_hi = 4.0_mm 0.25_mm 0.1_mm`
- `geometry.is_periodic = 0 1 1`
- `amr.n_cell = 64 8 8`  (reduce x from 512 to 64 for speed)
- All propellant, thermal, pf params from SCPThermalSandwich
- `amr.plot_int = 50`  (write one plot at the end for field comparison)
- `amr.thermo.int = 1`

Remove ALL elastic.* lines (with type=disable, model/BC blocks must be absent).
Add `elastic.type = disable` as the only elastic line.

The `amr.plot_int = 50` ensures a plot directory is written at step 50, which
testlib.py can load.

### C4 test.py (`tests/GPU/C4_amr_correctness/test.py`)

```python
#!/usr/bin/env python3
import os, sys, glob
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests/GPU"))
sys.path.insert(0, str(ROOT / "scripts"))
import testlib_gpu
import testlib

outdir = Path(sys.argv[1])
input_file = Path(__file__).parent / "input"

cpu_bin = testlib_gpu.find_binary(ROOT, "bin/alamo-2d-g++")
if cpu_bin is None:
    env = os.environ.get("ALAMO_CPU_BIN")
    cpu_bin = Path(env) if env else None

gpu_bin = testlib_gpu.find_binary(ROOT, "bin/alamo_gpu-2d-nofast-cuda86-g++")
if gpu_bin is None:
    env = os.environ.get("ALAMO_GPU_STRICT_BIN")
    gpu_bin = Path(env) if env else None

if cpu_bin is None:
    print("SKIP: CPU binary not found")
    sys.exit(77)
if gpu_bin is None:
    print("SKIP: GPU strict binary not found")
    sys.exit(77)

NUM_STEPS = 50

cpu_out = outdir / "cpu"
gpu_out = outdir / "gpu"

print("Running CPU...")
rc_cpu, log_cpu, _ = testlib_gpu.run_alamo(
    cpu_bin, input_file, cpu_out,
    extra_args=["max_step={}".format(NUM_STEPS)],
    timeout=300)
if rc_cpu != 0:
    print("FAIL: CPU run exited", rc_cpu)
    print(log_cpu[-2000:])
    sys.exit(1)

print("Running GPU strict...")
rc_gpu, log_gpu, _ = testlib_gpu.run_alamo(
    gpu_bin, input_file, gpu_out,
    extra_args=["max_step={}".format(NUM_STEPS)],
    timeout=300)
if rc_gpu != 0:
    print("FAIL: GPU run exited", rc_gpu)
    print(log_gpu[-2000:])
    sys.exit(1)

# Find the final plot directories
def find_latest_plot(base: Path):
    candidates = sorted(base.glob("plot/*/"), key=lambda p: p.name)
    if not candidates:
        # flat structure: plot/ IS the plot dir
        if (base / "plot" / "Header").exists():
            return base / "plot"
    return candidates[-1] if candidates else None

cpu_plot = find_latest_plot(cpu_out)
gpu_plot = find_latest_plot(gpu_out)

if cpu_plot is None or gpu_plot is None:
    print("FAIL: could not find plot directories")
    print("  CPU plot:", cpu_plot)
    print("  GPU plot:", gpu_plot)
    sys.exit(1)

# Use testlib to read cross-sections and compare
# Cross-section: x axis through the middle of the domain
try:
    cpu_df = testlib.readContours(str(cpu_plot),
                                   start=[0.0, 0, 0], end=[4.0, 0, 0],
                                   vars=["eta", "temp"])
    gpu_df = testlib.readContours(str(gpu_plot),
                                   start=[0.0, 0, 0], end=[4.0, 0, 0],
                                   vars=["eta", "temp"])
except Exception as e:
    print("FAIL: could not read plot data:", e)
    sys.exit(1)

import numpy as np
tolerance = 0.05

all_ok = True
for var in ["eta", "temp"]:
    try:
        cpu_vals = np.array(cpu_df[var])
        gpu_vals = np.array(gpu_df[var])
    except KeyError:
        print(f"WARN: variable '{var}' not in output, skipping")
        continue

    # Interpolate GPU to CPU x-points for fair comparison
    cpu_x = np.array(cpu_df["x"])
    gpu_x = np.array(gpu_df["x"])
    cpu_sorted = np.argsort(cpu_x)
    gpu_sorted = np.argsort(gpu_x)
    gpu_interp = np.interp(cpu_x[cpu_sorted], gpu_x[gpu_sorted], gpu_vals[gpu_sorted])
    cpu_s = cpu_vals[cpu_sorted]

    denom = np.maximum(np.abs(cpu_s) + np.abs(gpu_interp), 1e-30)
    rel_err = np.max(np.abs(cpu_s - gpu_interp) / denom)
    status = "OK" if rel_err <= tolerance else "FAIL"
    print(f"  {var}: max rel error = {rel_err:.4f} ({status})")
    if rel_err > tolerance:
        all_ok = False

if not all_ok:
    print("FAIL: C4 AMR correctness — GPU/CPU field mismatch")
    sys.exit(1)

print("PASS: C4 AMR correctness")
sys.exit(0)
```

NOTE: `testlib.readContours` expects string paths and unit-aware coordinates.
Adjust `start`/`end` to match the geometry in the input. With `system.length = mm`
and domain `[0, 4]_mm × [-0.25, 0.25]_mm`, the cross-section is
`start=[0.0, 0, 0], end=[4.0, 0, 0]` (in mm).

If `testlib` import fails (yt not available), catch the ImportError and `sys.exit(77)`.

## Invariants
- AMR must actually trigger refinement for this test to be meaningful:
  use `amr.refinement_criterion = 0.1` (eta gradient-based refinement).
- `amr.plot_int = 50` so a plot directory exists for field comparison.
- `elastic.type = disable` with no model/BC elastic blocks.
- The cross-section coordinates must match the actual domain in the input.
- `system.length = mm` must be in the input (from SCPThermalSandwich convention).

## Build and test commands
```bash
cd /home/jackplum/Projects/alamo
python3 -m py_compile tests/GPU/C4_amr_correctness/test.py
python3 -c "
from pathlib import Path
t = Path('tests/GPU/C4_amr_correctness/input').read_text()
assert 'amr.max_level = 2' in t
assert 'amr.plot_int = 50' in t
assert 'elastic.type = disable' in t
print('C4 input OK')
"
```

## Expected result
- 2 files created.
- Input has max_level=2, amr.plot_int=50, elastic.type=disable.
- test.py compiles cleanly.
- testlib import is guarded with try/except → exits 77 if yt is unavailable.

## Non-goals
- Do not compare full 2D field arrays (1D cross-section is sufficient and avoids
  grid interpolation complexity).
- Do not run the test.

## Stop conditions
Stop if SCPThermalSandwich/input cannot be read.

## Final report
Write `results/008-RESULT.md` with: files created, key parameter values confirmed,
note about yt availability guard, syntax check results.
