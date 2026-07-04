# Task 004: Performance Test P3 (3D 256×256×128)

## Goal
Create `tests/GPU/P3_perf_3d_256/input` and `tests/GPU/P3_perf_3d_256/test.py`.

This is a 3D performance test: 256×256×128 base grid, max_level=1, wide-shallow box
strategy matching the existing `input_3d_flame` production input. Run 150 steps,
report ms/step.

## Context
Project root: `/home/jackplum/Projects/alamo`
Plan: `docs/agent_plans/20260625-gpu-tests/PLAN.md`

**Differences from P1/P2:**
- Uses the 3D GPU binary: `alamo_gpu-3d-cuda86-g++`
- `input_3d_flame` is the canonical reference — copy it almost verbatim, just add
  test-specific overrides (plot_file, thermo.int, max_step).
- 3D runs are slower: 150 steps is sufficient.
- `elastic.type = disable` is already in `input_3d_flame` (the elastic fix is not yet
  ported to 3D; the 3D perf test stays flame-only).
- Expression ICs and BCs in `input_3d_flame` are already GPU-safe.

## Files to read first
- `input_3d_flame`     — PRIMARY reference; copy this almost exactly
- `input_3d_flame_256` — verify it matches input_3d_flame (they should be identical)
- `docs/agent_plans/20260625-gpu-tests/PLAN.md`

## Files allowed to modify
- `tests/GPU/P3_perf_3d_256/input`   (CREATE)
- `tests/GPU/P3_perf_3d_256/test.py` (CREATE)

## Files NOT allowed to modify
Everything else.

## Implementation steps

### P3 input (`tests/GPU/P3_perf_3d_256/input`)

Read `input_3d_flame` and copy it nearly verbatim. Apply these changes:
1. Change `plot_file = output_3d_flame` → `plot_file = output`
2. Change `amr.plot_int = -1` (already -1, confirm)
3. Add `amr.thermo.plot_int = -1` if not already -1
4. Add `amr.thermo.int = 1`
5. Remove the large comment block at the top (##...## header) — keep parameters only
6. All other parameters: copy exactly from `input_3d_flame`

The result should be a clean copy of `input_3d_flame` with just the plot_file renamed
and thermo.int set. Everything else (n_cell, max_level, ICs, BCs, propellant, thermal,
elastic.type=disable) should be identical.

### P3 test.py (`tests/GPU/P3_perf_3d_256/test.py`)

```python
#!/usr/bin/env python3
import os, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests/GPU"))
import testlib_gpu

outdir = Path(sys.argv[1])
input_file = Path(__file__).parent / "input"
test_name = Path(__file__).parent.name

# 3D binary
binary = testlib_gpu.find_binary(ROOT, "bin/alamo_gpu-3d-cuda86-g++")
if binary is None:
    env_bin = os.environ.get("ALAMO_GPU_3D_BIN")
    binary = Path(env_bin) if env_bin else None
if binary is None:
    print("SKIP: no 3D GPU binary found (bin/alamo_gpu-3d-cuda86-g++)")
    sys.exit(77)

NUM_STEPS = 150

rc, log, elapsed = testlib_gpu.run_alamo(
    binary, input_file, outdir,
    extra_args=["max_step={}".format(NUM_STEPS), "allow_unused=1"],
    timeout=900)   # 15 min timeout for 3D

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
- Input must be a clean copy of `input_3d_flame` (no new physics parameters added).
- `elastic.type = disable` must be present (and no elastic model/BC blocks).
- `amr.n_cell = 256 256 128` (matching the production input).
- `amr.max_level = 1` (matching the production input).
- Binary lookup uses pattern `bin/alamo_gpu-3d-cuda86-g++` (3D, not 2D).

## Build and test commands
```bash
cd /home/jackplum/Projects/alamo
python3 -m py_compile tests/GPU/P3_perf_3d_256/test.py
python3 -c "
from pathlib import Path
t = Path('tests/GPU/P3_perf_3d_256/input').read_text()
assert '256 256 128' in t, 'n_cell mismatch'
assert 'elastic.type = disable' in t, 'elastic not disabled'
assert 'alamo.program = flame' in t
print('P3 input OK')
"
```

## Expected result
- 2 files created.
- Input matches `input_3d_flame` except `plot_file = output` and `amr.thermo.int = 1`.

## Non-goals
- Do not port elastic to 3D — 3D elastic is out of scope for this test suite.
- Do not reduce grid size — 256×256×128 is the intended production test grid.

## Stop conditions
Stop if `input_3d_flame` cannot be read.

## Final report
Write `results/004-RESULT.md` with: files created, diff summary vs input_3d_flame,
syntax check results.
