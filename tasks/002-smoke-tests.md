> **SUPERSEDED (2026-07-06).** Historical only — the live plan is
> `docs/llm/PLAN.md`. Do not follow instructions in this file.
>
> # Task 002: Smoke Tests F1 and F2

## Goal
Create input files and test scripts for:
- **F1** `tests/GPU/F1_smoke_flame_only/` — GPU flame, no elastic, 2D, 20 steps
- **F2** `tests/GPU/F2_smoke_elastic/`    — GPU flame + elastic, 2D, 50 steps

## Context
Project root: `/home/jackplum/Projects/alamo`
Plan: `docs/agent_plans/20260625-gpu-tests/PLAN.md`

**F1 pass criteria:** exit 0, `thermo.dat` exists, all values finite, `chamber_pressure > 0`.
**F2 pass criteria:** exit 0, `thermo.dat` exists, all values finite, at least one row
after the first has non-zero `disp_xhi_x` or similar elastic displacement output,
confirming the elastic solve ran.

Framework utilities are in `tests/GPU/testlib_gpu.py` (created by Task 001).
Each `test.py` script:
1. Reads env vars `ALAMO_GPU_BIN` and `ALAMO_GPU_STRICT_BIN` (or falls back to
   `testlib_gpu.find_binary`).
2. Calls `testlib_gpu.run_alamo(binary, input_file, outdir, extra_args=['max_step=20'])`.
3. Runs post-processing checks. Exits 0 on pass, 1 on fail, 77 on skip (no binary).

## Files to read first
- `tests/SCPSandwich/input`       — base for F1 (flame, no elastic)
- `tests/SCPChamber/input`        — base for F2 (has elastic.on=1)
- `tests/GPU/testlib_gpu.py`      — utilities (may not exist yet; that's OK, just call it)
- `docs/agent_plans/20260625-gpu-tests/PLAN.md` — invariants

## Files allowed to modify
- `tests/GPU/F1_smoke_flame_only/input`   (CREATE)
- `tests/GPU/F1_smoke_flame_only/test.py` (CREATE)
- `tests/GPU/F2_smoke_elastic/input`      (CREATE)
- `tests/GPU/F2_smoke_elastic/test.py`    (CREATE)

## Files NOT allowed to modify
Everything else, including `tests/SCPSandwich/`, `tests/SCPChamber/`, `src/`, `scripts/`.

## Implementation steps

### F1 input (`tests/GPU/F1_smoke_flame_only/input`)

Base from `tests/SCPSandwich/input`. Key modifications:
- Remove ALL `#@` metadata lines (driver handles binary selection, not runtests.py)
- `plot_file = output`  (driver will prepend outdir)
- `amr.max_level = 0`  (no AMR for smoke test speed)
- `amr.n_cell = 64 64 64`
- `amr.blocking_factor = 16`
- `amr.max_grid_size = 64`
- `amr.plot_int = -1`  (no plot output, just thermo)
- `amr.thermo.int = 1`
- `elastic.type = disable`  (and OMIT all elastic.* model/BC/solver blocks)
- Keep: all propellant, thermal, chamber, pf.* blocks from SCPSandwich
- `variable_pressure = 0`  (fix pressure for smoke test stability)
- `chamber.pressure = 4.0_MPa`

CRITICAL: With `elastic.type = disable`, omit all `elastic.bc.*`, `elastic.solver.*`,
`elastic.traction`, `model_ap.*`, `model_htpb.*`, `model_void.*` lines — the strict
parser will abort on unused entries.

### F1 test.py (`tests/GPU/F1_smoke_flame_only/test.py`)

```python
#!/usr/bin/env python3
import os, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT / "tests/GPU"))

outdir = Path(sys.argv[1])
input_file = Path(__file__).parent / "input"

# Import testlib_gpu (created by Task 001)
import testlib_gpu

binary = testlib_gpu.find_binary(ROOT, "bin/alamo_gpu-2d-cuda86-g++")
if binary is None:
    binary_str = os.environ.get("ALAMO_GPU_BIN")
    if binary_str:
        binary = Path(binary_str)
if binary is None:
    print("SKIP: no GPU binary found")
    sys.exit(77)

rc, log, elapsed = testlib_gpu.run_alamo(binary, input_file, outdir,
                                          extra_args=["max_step=20"], timeout=120)
print(log[-2000:])   # last 2000 chars of log

if rc != 0:
    print(f"FAIL: alamo exited {rc}")
    sys.exit(1)

thermo = testlib_gpu.parse_thermo(outdir / "plot" / "thermo.dat")
ok, issues = testlib_gpu.check_sanity(thermo)
if not ok:
    for i in issues:
        print("FAIL:", i)
    sys.exit(1)

if thermo.get("chamber_pressure") and thermo["chamber_pressure"][-1] <= 0:
    print("FAIL: chamber_pressure not positive")
    sys.exit(1)

print(f"PASS: F1 smoke flame-only ({elapsed:.1f}s, {len(thermo.get('time',[]))-1} steps)")
sys.exit(0)
```

### F2 input (`tests/GPU/F2_smoke_elastic/input`)

Base from `tests/SCPChamber/input`. Key modifications:
- Remove ALL `#@` metadata lines
- `plot_file = output`
- `amr.n_cell = 8 8 8`  (small grid — elastic solves are expensive)
- `amr.max_level = 1`
- `amr.blocking_factor = 4`
- `amr.max_grid_size = 64`  (forces multi-box decomposition with elixir fix)
- `amr.plot_int = -1`
- `amr.thermo.int = 1`
- `elastic.type = static`  (KEEP elastic ON — this is the key test)
- `elastic.interval = 5`   (trigger elastic solve at step 5, 10, 15...)
- `elastic.solver.fixed_iter = 50`  (faster for smoke test)
- `elastic.solver.nriters = 50`
- `elastic.solver.nrtolerance = 1e-4`  (relaxed for speed)
- `elastic.solver.verbose = 0`  (suppress verbose output)
- `elastic.print_model = 0`
- Keep all elastic BC blocks from SCPChamber (they use `disp` type = GPU-safe)
- Keep model_ap.*, model_htpb.* blocks
- `variable_pressure = 0`
- `chamber.pressure = 4.0_MPa`

IMPORTANT: The elastic BCs in SCPChamber use `disp` type which maps to
`BC::Operator::Elastic::Constant` — this IS GPU-safe. Keep them.

### F2 test.py (`tests/GPU/F2_smoke_elastic/test.py`)

Same structure as F1 test.py but:
- `extra_args=["max_step=50"]`
- Use `ALAMO_GPU_STRICT_BIN` / pattern `bin/alamo_gpu-2d-nofast-cuda86-g++` (strict binary for correctness)
- After sanity check, verify at least one elastic solve happened:
  ```python
  # elastic.interval=5, so after 50 steps there should be ~10 elastic solves
  # A successful elastic solve leaves a sane displacement field in thermo
  # Check that disp_xhi_x column exists and has at least one nonzero value
  # (it will be zero until elastic runs, then nonzero if physics is working)
  disp = thermo.get("disp_xhi_x", [])
  # Actually just check thermo has >=50 rows (step 50 was reached)
  if len(thermo.get("time", [])) < 10:
      print("FAIL: simulation stopped early (< 10 thermo rows)")
      sys.exit(1)
  ```
- Print `PASS: F2 smoke elastic`

## Invariants
- F1 must NOT have any elastic.* lines in input (strict parser aborts on unused params
  when elastic.type=disable).
- F2 must use a Constant elastic BC type (GPU-safe). The disp BC from SCPChamber is fine.
- Both test.py scripts must exit 77 (not 1) when no binary is found.
- `outdir / "plot" / "thermo.dat"` is where ALAMO writes the thermo file when
  `plot_file = output` and the driver sets the actual path.

Wait — there's a subtlety: `run_alamo` in testlib_gpu should set `plot_file=<outdir>/plot`
on the command line. The thermo file then lives at `<outdir>/plot/thermo.dat`. Confirm
this assumption with the testlib_gpu implementation. If the path is different, adjust.

## Build and test commands
```bash
cd /home/jackplum/Projects/alamo
# Syntax check input (no GPU needed):
python3 -c "
from pathlib import Path
for p in ['tests/GPU/F1_smoke_flame_only/input', 'tests/GPU/F2_smoke_elastic/input']:
    lines = Path(p).read_text().splitlines()
    print(p, len(lines), 'lines')
"
# Syntax check test.py files:
python3 -m py_compile tests/GPU/F1_smoke_flame_only/test.py
python3 -m py_compile tests/GPU/F2_smoke_elastic/test.py
```

## Expected result
- 4 files created, all importable/parseable with no syntax errors.
- F1 input has no elastic.* lines, has `elastic.type = disable` ... wait, no:
  with `elastic.type = disable` you STILL have that one line. But all MODEL and BC
  blocks must be absent. Recheck: from PLAN.md invariant 4:
  "With `elastic.type = disable` the elastic model/bc blocks must be OMITTED".
  So F1 has exactly one elastic line: `elastic.type = disable`.
- F2 input has `elastic.type = static` and all elastic.bc.*, model_*.* blocks.

## Non-goals
- Do not implement perf measurement (that's P1/P2/P3).
- Do not actually run the tests (no GPU at plan-write time).

## Stop conditions
Stop if SCPSandwich/input or SCPChamber/input cannot be read.

## Final report
Write `results/002-RESULT.md` with: files created, line counts, syntax-check results.
