"""Shared utilities for the ALAMO GPU test suite.

Stdlib-only (no numpy/yt) so it can be imported on any machine regardless of
whether the heavy scientific Python stack is installed.
"""

from __future__ import annotations

import math
import os
import subprocess
import time
from pathlib import Path


# ---------------------------------------------------------------------------
# Binary discovery
# ---------------------------------------------------------------------------

def find_binary(root: Path, pattern: str) -> Path | None:
    """Glob for a binary matching *pattern* under ``root/bin/``.

    Returns the lexicographically newest match (most specific arch suffix
    sorts last alphabetically), or ``None`` if no match is found.
    """
    bin_dir = root / "bin"
    candidates = [
        p for p in bin_dir.glob(pattern)
        if p.is_file() and os.access(p, os.X_OK)
    ]
    if not candidates:
        return None
    return sorted(candidates)[-1]


# ---------------------------------------------------------------------------
# Running alamo
# ---------------------------------------------------------------------------

def run_alamo(
    binary: Path,
    input_file: Path,
    output_dir: Path,
    extra_args: list[str] = (),
    timeout: int = 600,
) -> tuple[int, str, float]:
    """Run an alamo binary.

    Parameters
    ----------
    binary:
        Absolute path to the alamo executable.
    input_file:
        Path to the AMReX ParmParse input file.
    output_dir:
        Directory where alamo should write its output.  ``plot_file`` will be
        set to ``<output_dir>/plot`` unless already present in *extra_args*.
    extra_args:
        Additional ``key=value`` overrides forwarded on the command line.
    timeout:
        Wall-clock timeout in seconds (default: 600).

    Returns
    -------
    (returncode, combined_stdout_stderr, elapsed_seconds)
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Inject plot_file if not already overridden by the caller.
    extra = list(extra_args)
    if not any(a.startswith("plot_file=") for a in extra):
        extra.append(f"plot_file={output_dir}/plot")

    cmd = [str(binary), str(input_file)] + extra

    started = time.monotonic()
    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            timeout=timeout,
            text=True,
        )
    except subprocess.TimeoutExpired as exc:
        elapsed = time.monotonic() - started
        partial = (exc.output or b"").decode("utf-8", errors="replace") if isinstance(exc.output, bytes) else (exc.output or "")
        return -1, f"TIMEOUT after {timeout}s\n{partial}", elapsed

    elapsed = time.monotonic() - started
    return result.returncode, result.stdout, elapsed


# ---------------------------------------------------------------------------
# Thermo parsing
# ---------------------------------------------------------------------------

def parse_thermo(thermo_path: Path) -> dict[str, list[float]]:
    """Parse a tab-separated ``thermo.dat`` file.

    The first row is treated as column headers; subsequent rows are data.

    Returns a dict mapping ``column_name -> list[float]``, or an empty dict
    if the file does not exist or is malformed.
    """
    if not thermo_path.exists():
        return {}
    try:
        with thermo_path.open("r", encoding="utf-8") as fh:
            header_line = fh.readline()
            if not header_line.strip():
                return {}
            headers = header_line.split()
            result: dict[str, list[float]] = {h: [] for h in headers}
            for line in fh:
                if not line.strip():
                    continue
                parts = line.split()
                if len(parts) != len(headers):
                    continue  # skip malformed rows
                for h, v in zip(headers, parts):
                    try:
                        result[h].append(float(v))
                    except ValueError:
                        result[h].append(float("nan"))
        return result
    except OSError:
        return {}


# ---------------------------------------------------------------------------
# Sanity checking
# ---------------------------------------------------------------------------

def check_sanity(thermo: dict) -> tuple[bool, list[str]]:
    """Check that all values are finite and at least one data row exists.

    Returns ``(ok, issues)`` where *issues* is a list of human-readable
    descriptions of any failures found.
    """
    issues: list[str] = []

    if not thermo:
        issues.append("thermo dict is empty (file missing or unreadable)")
        return False, issues

    # Check we have at least one row
    first_col = next(iter(thermo.values()))
    if len(first_col) == 0:
        issues.append("thermo file has a header but no data rows")
        return False, issues

    # Check all values are finite
    for col_name, values in thermo.items():
        for i, v in enumerate(values):
            if not math.isfinite(v):
                issues.append(
                    f"column '{col_name}' row {i}: non-finite value {v!r}"
                )

    return len(issues) == 0, issues


# ---------------------------------------------------------------------------
# Thermo comparison
# ---------------------------------------------------------------------------

def compare_thermo(
    ref: dict,
    cand: dict,
    rel_tol: float = 0.05,
    abs_tol: float = 1e-12,
) -> tuple[bool, list[str]]:
    """Compare two thermo dicts column-by-column.

    A (ref, cand) value pair passes if::

        abs(cand - ref) <= abs_tol
        OR
        abs(cand - ref) / max(abs(ref), 1e-30) <= rel_tol

    Returns ``(ok, failures)`` where *failures* is a list of descriptions.
    """
    if not ref:
        return False, ["reference thermo is empty"]
    if not cand:
        return False, ["candidate thermo is empty"]

    failures: list[str] = []

    # Only compare columns present in both dicts
    common_cols = [c for c in ref if c in cand]
    if not common_cols:
        return False, ["no common columns between reference and candidate"]

    for col in common_cols:
        ref_vals = ref[col]
        cand_vals = cand[col]
        if len(ref_vals) != len(cand_vals):
            failures.append(
                f"column '{col}': row count mismatch "
                f"ref={len(ref_vals)} cand={len(cand_vals)}"
            )
            continue
        for i, (r, c) in enumerate(zip(ref_vals, cand_vals)):
            diff = abs(c - r)
            rel = diff / max(abs(r), 1e-30)
            if diff > abs_tol and rel > rel_tol:
                failures.append(
                    f"column '{col}' row {i}: ref={r} cand={c} "
                    f"abs_err={diff:.3e} rel_err={rel:.3e} "
                    f"(abs_tol={abs_tol} rel_tol={rel_tol})"
                )

    return len(failures) == 0, failures


# ---------------------------------------------------------------------------
# Performance helper
# ---------------------------------------------------------------------------

def wall_per_step(elapsed_sec: float, thermo: dict) -> float | None:
    """Return milliseconds per simulation step.

    Computed as ``elapsed_sec * 1000 / (rows - 1)`` so that the first
    (initialisation) step is excluded from the average.

    Returns ``None`` if *thermo* has fewer than 2 rows.
    """
    if not thermo:
        return None
    first_col = next(iter(thermo.values()))
    n_rows = len(first_col)
    if n_rows < 2:
        return None
    return elapsed_sec * 1000.0 / (n_rows - 1)
