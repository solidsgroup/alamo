#!/usr/bin/env python3
"""Shared parsing/classification helpers for alamo's MLMG/Newton solver
diagnostics.

Part of the diagnostics-first pass on MLMG convergence at high modulus
contrast (see ~/.claude/plans/improve-MLMG-solver-ideas.md, item F6).
Before this module, iteration counts and convergence status were recovered
ad hoc by regexing stdout in individual test scripts (e.g.
tests/ElasticSoftVoidAMR/test). This module gives every consumer -- the
benchmark harness in solver_benchmark.py, and eventually that test script --
one shared, tested parser.

Two layers, preferred in this order:

1. The "SOLVER_STATS" line emitted by Solver::Nonlocal::Linear
   (src/Solver/Nonlocal/Linear.H, reportSolverStats) whenever
   elastic.solver.verbose >= 1. This gives exact iteration counts,
   residuals, and a status already classified as converged/stalled/
   diverged by the code that actually ran the solve.
2. A fallback regex over amrex::MLMG's own stdout text (the
   "MLMG: Final Iter." / "MLMG: Failed to converge" / "MLMG: Failing to
   converge" family), for logs captured before SOLVER_STATS existed, or
   for executables/configurations where verbose was left at 0.

Layer 2 is intentionally conservative: it cannot recover convfactor (no
per-iteration residual history without SOLVER_STATS or verbose>=2), so
records built from it leave that field as NaN.
"""
import dataclasses
import math
import re
from typing import List, Optional


# ---------------------------------------------------------------------------
# Regexes
# ---------------------------------------------------------------------------

# src/Solver/Nonlocal/Linear.H reportSolverStats():
#   "SOLVER_STATS iters=18 resid0=2449.19 resid=2.23597e-05 rhs0=2449.19
#    convfactor=0.383003 status=converged"
SOLVER_STATS_RE = re.compile(
    r"SOLVER_STATS\s+iters=(?P<iters>\S+)\s+resid0=(?P<resid0>\S+)\s+"
    r"resid=(?P<resid>\S+)\s+rhs0=(?P<rhs0>\S+)\s+"
    r"convfactor=(?P<convfactor>\S+)\s+status=(?P<status>\S+)"
)

# ext/AMReX-Codes/amrex/Src/LinearSolvers/MLMG/AMReX_MLMG.H
MLMG_FINAL_RE = re.compile(
    r"MLMG: Final Iter\.\s+(?P<iters>\d+)\s+resid, resid/\S+\s*=\s*"
    r"(?P<resid>\S+),\s*(?P<ratio>\S+)"
)
MLMG_STALL_RE = re.compile(r"MLMG: Failed to converge after (?P<iters>\d+) iterations")
MLMG_DIVERGE_RE = re.compile(r"MLMG: Failing to converge after (?P<iters>\d+) iterations")

# src/Solver/Nonlocal/Newton.H
NEWTON_ITER_RE = re.compile(
    r"NR iteration (?P<nriter>\d+), alpha = (?P<alpha>\S+), "
    r"full max norm\(ddisp\) = (?P<full>\S+), "
    r"accepted max norm\(ddisp\) = (?P<accepted>\S+)"
)
NEWTON_ABORT_RE = re.compile(
    r"Newton line search failed after (?P<backtracks>\d+) backtracks: "
    r"initial residual=(?P<r0>\S+), final trial residual=(?P<r1>\S+)"
)
NONFINITE_RE = re.compile(r"non-finite|\bnan\b|\binfinity\b", re.IGNORECASE)


@dataclasses.dataclass
class SolveRecord:
    iters: int
    resid0: float
    resid: float
    rhs0: float
    convfactor: float
    status: str  # "converged" | "stalled" | "diverged" | "unknown"


@dataclasses.dataclass
class NewtonRecord:
    nriter: int
    alpha: float
    full_update: float
    accepted_update: float


def parse_solver_stats(stdout: str) -> List[SolveRecord]:
    """One SolveRecord per SOLVER_STATS line, in order."""
    records = []
    for m in SOLVER_STATS_RE.finditer(stdout):
        records.append(SolveRecord(
            iters=int(m.group("iters")),
            resid0=float(m.group("resid0")),
            resid=float(m.group("resid")),
            rhs0=float(m.group("rhs0")),
            convfactor=float(m.group("convfactor")),
            status=m.group("status"),
        ))
    return records


def parse_solver_stats_fallback(stdout: str) -> List[SolveRecord]:
    """Recover what we can from AMReX's own prose when SOLVER_STATS is
    absent (older logs, or verbose=0). Order is not reliable across mixed
    converged/stalled/diverged solves in one run -- use parse_solver_stats
    when available."""
    records = []
    for m in MLMG_FINAL_RE.finditer(stdout):
        records.append(SolveRecord(
            iters=int(m.group("iters")), resid0=math.nan,
            resid=float(m.group("resid")), rhs0=math.nan,
            convfactor=math.nan, status="converged",
        ))
    for m in MLMG_STALL_RE.finditer(stdout):
        records.append(SolveRecord(
            iters=int(m.group("iters")), resid0=math.nan, resid=math.nan,
            rhs0=math.nan, convfactor=math.nan, status="stalled",
        ))
    for m in MLMG_DIVERGE_RE.finditer(stdout):
        records.append(SolveRecord(
            iters=int(m.group("iters")), resid0=math.nan, resid=math.nan,
            rhs0=math.nan, convfactor=math.nan, status="diverged",
        ))
    return records


def parse_newton_iters(stdout: str) -> List[NewtonRecord]:
    return [
        NewtonRecord(
            nriter=int(m.group("nriter")), alpha=float(m.group("alpha")),
            full_update=float(m.group("full")), accepted_update=float(m.group("accepted")),
        )
        for m in NEWTON_ITER_RE.finditer(stdout)
    ]


def newton_line_search_failed(stdout: str) -> bool:
    return NEWTON_ABORT_RE.search(stdout) is not None


# ---------------------------------------------------------------------------
# Stable "what happened" classification, matching amrex::MLMG's own
# divergence criterion (composite_norminf > 1e20*max_norm, AMReX_MLMG.H)
# rather than inventing a new one.
# ---------------------------------------------------------------------------

def classify_run(returncode: int, stdout: str, timed_out: bool = False) -> str:
    """Returns one of: OK, STALL, DIVERGE, NR_ABORT, TIMEOUT, CRASH.

    Priority matters: a run can converge in MLMG and still fail overall
    (tests/ElasticSoftVoidAMR does exactly this -- MLMG converges at 1055
    iterations, then Newton's line search aborts), so NR_ABORT is checked
    before treating a converged SolveRecord as OK.
    """
    if timed_out:
        return "TIMEOUT"
    if newton_line_search_failed(stdout):
        return "NR_ABORT"

    records = parse_solver_stats(stdout)
    if records:
        last = records[-1]
        if last.status == "diverged":
            return "DIVERGE"
        if last.status == "stalled":
            return "STALL"
        if last.status == "converged":
            return "OK" if returncode == 0 else "CRASH"
        # status == "unknown": fall through to the raw-text checks below.
    else:
        if MLMG_DIVERGE_RE.search(stdout):
            return "DIVERGE"
        if MLMG_STALL_RE.search(stdout):
            return "STALL"
        if MLMG_FINAL_RE.search(stdout):
            return "OK" if returncode == 0 else "CRASH"

    if NONFINITE_RE.search(stdout):
        return "CRASH"
    return "OK" if returncode == 0 else "CRASH"


def convergence_factor(records: List[SolveRecord]) -> float:
    """Geometric mean of the reported per-solve convfactors, ignoring NaNs.
    Prefer SolveRecord.convfactor from a single solve's residual history
    (Linear::getConvergenceFactor) when available; this is only useful as a
    summary across multiple solves (e.g. several Newton iterations)."""
    vals = [r.convfactor for r in records if r.convfactor == r.convfactor]  # drop NaN
    if not vals:
        return math.nan
    logsum = sum(math.log(v) for v in vals if v > 0)
    n = sum(1 for v in vals if v > 0)
    return math.exp(logsum / n) if n else math.nan


# ---------------------------------------------------------------------------
# #@ test-header reader (the idiom scripts/runtests.py uses to discover
# exe/dim/args from a tests/<Name>/input file), so the benchmark harness
# reads its base configuration from the same source of truth as the test
# rather than duplicating it.
# ---------------------------------------------------------------------------

def read_test_config(input_path: str, section: Optional[str] = None) -> dict:
    """Parses the '#@ [section]' / '#@ key = value' header block at the top
    of a tests/<Name>/input file (the same syntax scripts/runtests.py
    consumes) and returns the merged key/value dict for the requested
    section (first section in the file if None)."""
    import configparser
    lines = []
    with open(input_path) as f:
        for line in f:
            if line.startswith("#@"):
                lines.append(line[2:].strip())
            elif lines:
                break  # header block ends at the first non-#@ line
    cp = configparser.ConfigParser()
    cp.read_string("\n".join(lines))
    sect = section or cp.sections()[0]
    return dict(cp[sect])


# ---------------------------------------------------------------------------
# Error columns for tests/LowMachElasticPressure, factored out of
# tests/LowMachElasticPressure/test so solver_benchmark.py can reuse the
# exact same ray sample and closed-form comparison instead of duplicating it.
# Imports testlib lazily so this module stays import-cheap (no yt/pandas
# dependency) for callers that only need the regex/classification helpers.
# ---------------------------------------------------------------------------

def lowmach_pressure_errors(plotfile: str, pressure: float,
                             mu: float = 140.0, kappa: float = 150.0,
                             x_mid: float = 1.6e-3,
                             y_lo: float = 0.5e-3, y_hi: float = 2.4e-3,
                             y_top: float = 6.4e-3):
    """Returns (disp_y_relL2_err, sigma_yy_relerr, sigma_ratio_relerr).

    Mirrors tests/LowMachElasticPressure/test's analytic comparison:
        sigma_yy(y) = -P
        disp_y(y)   = -P*y/M,  M = kappa + 4*mu/3
        sigma_xx/sigma_yy = lambda/M,  lambda = kappa - 2*mu/3
    sampled on the same vertical ray (x=1.6e-3, y in [0.5e-3, 2.4e-3]) used
    by that test, so results are directly comparable to it.
    """
    import os
    import sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import testlib  # noqa: local import, needs yt/pandas

    lam = kappa - 2.0 * mu / 3.0
    M = lam + 2.0 * mu

    df = testlib.readContours(
        path=plotfile,
        start=[x_mid, 0.0, 0.0],
        end=[x_mid, y_top, 0.0],
        vars=["disp_x", "disp_y", "stress_xx", "stress_yy"],
    )
    y = df["y"].to_numpy()
    mask = (y >= y_lo) & (y <= y_hi)
    if mask.sum() < 4:
        raise RuntimeError(f"too few sample points in [{y_lo},{y_hi}]: {mask.sum()}")

    y_s = y[mask]
    disp_y = df["disp_y"].to_numpy()[mask]
    stress_xx = df["stress_xx"].to_numpy()[mask]
    stress_yy = df["stress_yy"].to_numpy()[mask]

    disp_y_exact = -pressure * y_s / M
    err = testlib.numpy.sqrt(testlib.integrate(y_s, (disp_y - disp_y_exact) ** 2))
    mag = testlib.numpy.sqrt(testlib.integrate(y_s, disp_y_exact ** 2))
    disp_y_relerr = float(err / mag) if mag > 1.0e-300 else float(err)

    sigma_yy_mean = float(stress_yy.mean())
    sigma_yy_expected = -pressure
    sigma_yy_relerr = abs(sigma_yy_mean - sigma_yy_expected) / abs(sigma_yy_expected)

    ratio = float((stress_xx / stress_yy).mean())
    ratio_expected = lam / M
    ratio_relerr = abs(ratio - ratio_expected) / abs(ratio_expected)

    return disp_y_relerr, sigma_yy_relerr, ratio_relerr


if __name__ == "__main__":
    # Minimal self-check against a synthetic log, so this module can be
    # validated without running the executable.
    sample = """
MLMG: Initial rhs               = 2449.186624
MLMG: Initial residual (resid0) = 2449.186624
MLMG: Final Iter. 18 resid, resid/bnorm = 2.235970385e-05, 9.129440619e-09
MESSAGE: ./src/Solver/Nonlocal/Linear.H:472 (reportSolverStats) SOLVER_STATS iters=18 resid0=2449.19 resid=2.23597e-05 rhs0=2449.19 convfactor=0.383003 status=converged
NR iteration 1, alpha = 1, full max norm(ddisp) = 1.19164e-05, accepted max norm(ddisp) = 1.19164e-05
"""
    recs = parse_solver_stats(sample)
    assert len(recs) == 1 and recs[0].iters == 18 and recs[0].status == "converged", recs
    assert classify_run(0, sample) == "OK"
    nr = parse_newton_iters(sample)
    assert len(nr) == 1 and nr[0].nriter == 1
    print("solverlib self-check OK")
