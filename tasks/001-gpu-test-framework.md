> **SUPERSEDED (2026-07-06).** Historical only — the live plan is
> `docs/llm/PLAN.md`. Do not follow instructions in this file.
>
> # Task 001: GPU Test Framework

## Goal
Create `tests/GPU/testlib_gpu.py` (shared utilities) and `tests/GPU/run_gpu_tests.py`
(CLI driver). These are the framework files — individual tests live in sibling tasks.

## Context
Project root: `/home/jackplum/Projects/alamo`
Plan: `docs/agent_plans/20260625-gpu-tests/PLAN.md`

The GPU test suite lives under `tests/GPU/`. Each sub-directory has an `input` file
and a `test.py` that accepts `outdir` as argv[1] and exits 0 on pass. The driver
auto-discovers all such subdirectories and runs them.

Available binaries (under `bin/`):
- `alamo_gpu-2d-cuda86-g++`         — GPU fast, 2D
- `alamo_gpu-2d-nofast-cuda86-g++`  — GPU strict (no fast-math), 2D
- `alamo_gpu-3d-cuda86-g++`         — GPU fast, 3D
- `alamo-2d-g++`                    — CPU reference, 2D

Thermo output: `<outdir>/thermo.dat`, tab-separated, first row is header,
subsequent rows are one-per-step data.

## Files to read first
- `benchmark/baseline_suite.py`   — see how it runs alamo and reads thermo
- `scripts/runtests.py`           — see how it constructs binary paths and runs tests
- `scripts/testlib.py`            — see existing Python test helpers

## Files allowed to modify
- `tests/GPU/testlib_gpu.py`      (CREATE)
- `tests/GPU/run_gpu_tests.py`    (CREATE)

## Files NOT allowed to modify
Everything else. Do not touch `scripts/`, `src/`, `benchmark/`, or existing `tests/`.

## Implementation steps

### Step 1: Create `tests/GPU/testlib_gpu.py`

Provide these functions:

```python
def find_binary(root: Path, pattern: str) -> Path | None:
    """Glob for a binary matching pattern under root/bin/, return newest match or None."""

def run_alamo(binary: Path, input_file: Path, output_dir: Path,
              extra_args: list[str] = (), timeout: int = 600) -> tuple[int, str, float]:
    """Run alamo binary. Returns (returncode, combined_stdout+stderr, elapsed_seconds).
    Sets plot_file=<output_dir>/plot in extra_args if not already set.
    Captures stdout+stderr together."""

def parse_thermo(thermo_path: Path) -> dict[str, list[float]]:
    """Parse tab-separated thermo.dat. Returns dict of column_name -> list of float values.
    Returns empty dict if file does not exist."""

def check_sanity(thermo: dict) -> tuple[bool, list[str]]:
    """Check all values are finite and at least one row exists.
    Returns (ok: bool, issues: list of strings describing failures)."""

def compare_thermo(ref: dict, cand: dict,
                   rel_tol: float = 0.05, abs_tol: float = 1e-12
                   ) -> tuple[bool, list[str]]:
    """Compare two thermo dicts column-by-column. Passes if for every (ref, cand) pair:
    abs(cand-ref) <= abs_tol  OR  abs(cand-ref)/max(abs(ref),1e-30) <= rel_tol.
    Returns (ok, list of failure descriptions)."""

def wall_per_step(elapsed_sec: float, thermo: dict) -> float | None:
    """Return ms per step = elapsed_sec*1000 / (len(thermo['time'])-1).
    Returns None if thermo has fewer than 2 rows."""
```

Use only stdlib + pathlib. No numpy/yt dependency here.

### Step 2: Create `tests/GPU/run_gpu_tests.py`

CLI: `python3 tests/GPU/run_gpu_tests.py [--test NAME] [--dry-run] [--timeout SEC]`

Logic:
1. `ROOT = Path(__file__).resolve().parents[2]`  (points to alamo project root)
2. Discover test dirs: all subdirs of `tests/GPU/` that contain both `input` and `test.py`.
3. Select GPU-fast 2D binary: `bin/alamo_gpu-2d-cuda86-g++`; GPU-strict: `bin/alamo_gpu-2d-nofast-cuda86-g++`; 3D: `bin/alamo_gpu-3d-cuda86-g++`; CPU: `bin/alamo-2d-g++`.
4. If `--dry-run`, print discovered tests and binary paths, exit 0.
5. For each discovered test (filtered by `--test` if given):
   a. Create `outdir = ROOT / "tests/GPU" / test_name / f"output_{timestamp}"`
   b. Run `test.py` with `outdir` as argv[1] — but FIRST the driver must run the alamo
      binary; `test.py` only does post-processing checks. Actually: let test.py handle
      everything (it imports testlib_gpu). The driver just calls:
      `subprocess.run([sys.executable, test_py, str(outdir)], ...)`
   c. Capture result (exit code, stdout/stderr, elapsed)
   d. Print PASS / FAIL / SKIP with elapsed time
6. Print final summary table: test name | result | time
7. Exit nonzero if any test failed.

For binary not found: print WARNING and SKIP the test (exit 0 from driver perspective
for that test, but note it in summary as SKIP).

Environment variable overrides:
- `ALAMO_GPU_BIN`        — override the 2D GPU-fast binary path
- `ALAMO_GPU_STRICT_BIN` — override the 2D GPU-strict binary path
- `ALAMO_GPU_3D_BIN`     — override the 3D GPU binary path
- `ALAMO_CPU_BIN`        — override the CPU binary path

Each `test.py` gets the output dir and also needs to know binary paths. Pass them
as environment variables before calling subprocess: set the four `ALAMO_*_BIN` env
vars so test.py scripts can import testlib_gpu and call `find_binary()` or read env.

Actually simpler: `test.py` scripts read `ALAMO_GPU_BIN` etc. from environment, or
call `testlib_gpu.find_binary(ROOT, pattern)` as fallback.

## Invariants
- `tests/GPU/` directory must exist before writing files (create it with `mkdir -p`).
- testlib_gpu.py must not import yt/numpy (keep it stdlib-only for portability).
- run_gpu_tests.py must exit 0 on `--dry-run` even if no GPU is present.
- All output directories created under `tests/GPU/<test_name>/output_<timestamp>/`.

## Build and test commands
```bash
cd /home/jackplum/Projects/alamo
python3 -c "import tests.GPU.testlib_gpu"   # syntax check
python3 tests/GPU/run_gpu_tests.py --dry-run  # must exit 0 and list tests
```

## Expected result
- `tests/GPU/testlib_gpu.py` exists and imports cleanly.
- `tests/GPU/run_gpu_tests.py` exists, `--dry-run` exits 0.
- Running `--dry-run` shows all 9 test directories (F1, F2, P1, P2, P3, C1, C2, C3, C4).

## Non-goals
- Do not implement the individual test input files or test.py scripts (that's Tasks 002-008).
- Do not wire into `scripts/runtests.py`.

## Stop conditions
Stop if you cannot find the alamo project root or if required reference files are missing.

## Final report
Write `results/001-RESULT.md` with: summary, files changed, import test result,
dry-run output, any issues.
