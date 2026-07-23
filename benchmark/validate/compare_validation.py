#!/usr/bin/env python3
"""Physics-budget comparator for the GPU Roadmap v3 validation suite.

Roadmap task 1.D (`benchmark/GPU_ROADMAP_V3.md` Phase 1). Compares two output
bundles (schema: `benchmark/validate/README.md`, task 1.C) against the
observable registry in `physics_budget.yaml` (task 1.A) and reports a
per-observable PASS/FAIL table plus a roll-up verdict.

This tool reads each bundle's per-case `metrics.json` only -- it does not
re-parse `run.log`/`thermo.dat`/plotfiles itself. That extraction (regex over
`run.log`, thermo.dat column selection, plotfile field norms) happens once at
record time per `physics_budget.yaml`'s `pattern`/`column`/`field` entries and
is baked into `metrics.json` (the bundle's "canonical comparison surface",
per the 1.C schema doc) -- the responsibility of the validation runner
(task 1.B), not this comparator. Keeping that one-directional means there is
exactly one place that knows how to turn raw run output into an observable
value.

`metrics.json` shape (one entry per physics_budget observable):
    {
      "<observable_name>": {"class": "correctness", "value": <scalar>},
      "<observable_name>": {"class": "correctness", "value": {"x": <scalar>, "y": <scalar>}},   # per_component
      "<observable_name>": {"class": "engineering_trajectory", "value": [[t0, v0], [t1, v1], ...]},
      "<observable_name>": {"class": "engineering_trajectory", "value": <scalar>, "endpoint_only": true},  # derived_field source
      "<observable_name>": {"class": "solver_health", "value": <scalar-or-list>}
    }

The PASS/FAIL logic (value passes if within abs_tol OR within rel_tol) mirrors
the existing convention in `benchmark/baseline_suite.py`'s `compare_rows` /
`benchmark/compare_thermo.py` -- this tool extends that convention to
field-norm and trajectory observables, it does not reinvent it.
"""

from __future__ import annotations

import argparse
import json
import math
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any

import yaml

VALIDATE_DIR = Path(__file__).resolve().parent
DEFAULT_BUDGET = VALIDATE_DIR / "physics_budget.yaml"

EPS = 1.0e-300  # avoids div-by-zero on exact-zero references without masking real misses


# ---------------------------------------------------------------------------
# Budget loading
# ---------------------------------------------------------------------------

def load_budget(path: Path) -> dict[str, list[dict[str, Any]]]:
    data = yaml.safe_load(path.read_text(encoding="utf-8"))
    for key in ("correctness", "engineering_trajectory", "solver_health"):
        data.setdefault(key, [])
    return data


def load_metrics(bundle_dir: Path, case: str) -> dict[str, Any]:
    path = bundle_dir / case / "metrics.json"
    if not path.exists():
        raise SystemExit(f"missing metrics.json: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def discover_cases(bundle_dir: Path) -> list[str]:
    if not bundle_dir.is_dir():
        raise SystemExit(f"not a bundle directory: {bundle_dir}")
    return sorted(
        p.name for p in bundle_dir.iterdir()
        if p.is_dir() and (p / "metrics.json").exists()
    )


def check_manifest_compatibility(reference: Path, candidate: Path) -> None:
    """Reject A/B bundles whose execution/oracle identity is not comparable."""
    rp, cp = reference / "manifest.json", candidate / "manifest.json"
    if not rp.is_file() or not cp.is_file():
        raise SystemExit(f"manifest.json required in both bundles: {rp}, {cp}")
    a, b = json.loads(rp.read_text()), json.loads(cp.read_text())
    required = ("host", "device", "profile", "gpu_name", "gpu_uuid", "driver_version", "cuda_compute_cap", "build_flags", "oracle_scripts")
    for key in required:
        if a.get(key) in (None, "") or b.get(key) in (None, ""):
            raise SystemExit(f"manifest identity field missing/null: {key}")
    mismatches = []
    for key in required:
        if a.get(key) != b.get(key):
            mismatches.append(key)
    # Per-case inputs and overrides are part of the command identity.
    def case_identity(m):
        return [(c.get("id"), c.get("input"), c.get("input_sha256"), c.get("overrides"), c.get("overrides_sha256"),
                 c.get("command_shape"), c.get("max_step"), c.get("np")) for c in m.get("cases", [])]
    if case_identity(a) != case_identity(b):
        mismatches.append("case_inputs_overrides")
    if a.get("command_shape") != b.get("command_shape"):
        mismatches.append("command_shape")
    if mismatches:
        raise SystemExit("incompatible validation manifests: " + ", ".join(mismatches))


# ---------------------------------------------------------------------------
# Per-class comparison
# ---------------------------------------------------------------------------

def _passes(delta_abs: float, ref_mag: float, abs_tol: float | None, rel_tol: float | None) -> bool:
    rel = delta_abs / max(ref_mag, EPS)
    ok_abs = abs_tol is not None and delta_abs <= abs_tol
    ok_rel = rel_tol is not None and rel <= rel_tol
    return ok_abs or ok_rel


def _scalar_row(name: str, cls: str, a: float, b: float, abs_tol: float | None, rel_tol: float | None,
                 extra: dict[str, Any] | None = None) -> dict[str, Any]:
    delta = abs(b - a)
    ref_mag = max(abs(a), abs(b))
    rel = delta / max(ref_mag, EPS)
    ok = _passes(delta, ref_mag, abs_tol, rel_tol)
    row = {
        "name": name, "class": cls, "A": a, "B": b,
        "delta_abs": delta, "delta_rel": rel,
        "abs_tol": abs_tol, "rel_tol": rel_tol,
        "status": "PASS" if ok else "FAIL",
    }
    if extra:
        row.update(extra)
    return row


def compare_correctness(spec: dict[str, Any], a_entry: dict[str, Any] | None,
                         b_entry: dict[str, Any] | None) -> list[dict[str, Any]]:
    name = spec["name"]
    abs_tol = spec.get("abs_tol")
    rel_tol = spec.get("rel_tol")
    if a_entry is None and b_entry is None:
        # Consistently absent from both sides (e.g. a solver mode that never emits
        # this observable, like MLMG under fixed_iter) is not a divergence -- N/A,
        # non-gating. Missing from only one side means the two runs took genuinely
        # different paths, which stays a FAIL below.
        return [{"name": name, "class": "correctness", "status": "N/A",
                  "note": "observable absent from both bundles (consistent, non-gating)"}]
    if a_entry is None or b_entry is None:
        return [{"name": name, "class": "correctness", "status": "FAIL",
                  "note": "observable missing from one bundle only (divergent solver path)"}]
    a_val, b_val = a_entry["value"], b_entry["value"]
    if spec.get("per_component"):
        if not isinstance(a_val, dict) or not isinstance(b_val, dict):
            return [{"name": name, "class": "correctness", "status": "FAIL",
                      "note": "expected per-component dict value"}]
        components = sorted(set(a_val) | set(b_val))
        rows = []
        for comp in components:
            if comp not in a_val or comp not in b_val:
                rows.append({"name": f"{name}.{comp}", "class": "correctness", "status": "FAIL",
                              "note": "component missing from one bundle"})
                continue
            rows.append(_scalar_row(f"{name}.{comp}", "correctness", a_val[comp], b_val[comp], abs_tol, rel_tol))
        return rows
    return [_scalar_row(name, "correctness", a_val, b_val, abs_tol, rel_tol)]


def _interp_series(series: list[list[float]], t: float) -> float:
    """Linear-interpolate a [[t, v], ...] series (sorted by t) at time t, clamped at ends."""
    if t <= series[0][0]:
        return series[0][1]
    if t >= series[-1][0]:
        return series[-1][1]
    for (t0, v0), (t1, v1) in zip(series, series[1:]):
        if t0 <= t <= t1:
            if t1 == t0:
                return v0
            frac = (t - t0) / (t1 - t0)
            return v0 + frac * (v1 - v0)
    return series[-1][1]


def _max_rel_dev(a_series: list[list[float]], b_series: list[list[float]]) -> float:
    """Max relative deviation of B vs A, B resampled onto A's time grid."""
    worst = 0.0
    for t, a_v in a_series:
        b_v = _interp_series(b_series, t)
        worst = max(worst, abs(b_v - a_v) / max(abs(a_v), abs(b_v), EPS))
    return worst


def _phase_lag_samples(a_series: list[list[float]], b_series: list[list[float]], max_shift: int = 10) -> int:
    """Integer-sample lag (in A's index space) that maximizes cross-correlation of B against A.

    Returns the shift k (can be negative) such that B(t_{i-k}) best matches A(t_i).
    Samples are resampled onto A's own time grid first so unequal-length series compare cleanly.
    """
    a_vals = [v for _, v in a_series]
    n = len(a_vals)
    if n < 3:
        return 0
    b_vals = [_interp_series(b_series, t) for t, _ in a_series]
    a_mean = sum(a_vals) / n
    b_mean = sum(b_vals) / n
    a_c = [v - a_mean for v in a_vals]
    b_c = [v - b_mean for v in b_vals]
    best_k, best_score = 0, -math.inf
    span = min(max_shift, n - 1)
    for k in range(-span, span + 1):
        score = 0.0
        count = 0
        for i in range(n):
            j = i - k
            if 0 <= j < n:
                score += a_c[i] * b_c[j]
                count += 1
        if count == 0:
            continue
        score /= count
        if score > best_score:
            best_score, best_k = score, k
    return best_k


def compare_engineering(spec: dict[str, Any], a_entry: dict[str, Any] | None,
                         b_entry: dict[str, Any] | None) -> dict[str, Any]:
    name = spec["name"]
    cls = "engineering_trajectory"
    if a_entry is None and b_entry is None:
        return {"name": name, "class": cls, "status": "N/A",
                "note": "observable absent from both bundles (consistent, non-gating)"}
    if a_entry is None or b_entry is None:
        return {"name": name, "class": cls, "status": "FAIL",
                "note": "observable missing from one bundle only (divergent solver path)"}
    a_val, b_val = a_entry["value"], b_entry["value"]
    endpoint_only = bool(a_entry.get("endpoint_only") or spec.get("source") == "derived_field")
    max_rel_dev_tol = spec.get("max_rel_dev")
    max_abs_dev_tol = spec.get("max_abs_dev")

    if spec.get("columns"):
        # multi-column observable (e.g. corner_displacement): value is {col: [[t,v],...]}
        if not isinstance(a_val, dict) or not isinstance(b_val, dict):
            return {"name": name, "class": cls, "status": "FAIL",
                      "note": "expected per-column dict value"}
        worst_rel = 0.0
        for col, a_series in a_val.items():
            b_series = b_val.get(col, [])
            if not a_series or not b_series:
                return {"name": name, "class": cls, "status": "FAIL",
                          "note": f"column {col} missing or empty series"}
            worst_rel = max(worst_rel, _max_rel_dev(a_series, b_series))
        ok = max_rel_dev_tol is not None and worst_rel <= max_rel_dev_tol
        return {
            "name": name, "class": cls, "max_rel_dev": worst_rel, "rel_tol": max_rel_dev_tol,
            "endpoint_only": False, "status": "PASS" if ok else "FAIL",
        }

    if endpoint_only or not isinstance(a_val, list):
        # scalar endpoint comparison only (no trajectory data available for this observable)
        delta = abs(b_val - a_val)
        ref_mag = max(abs(a_val), abs(b_val))
        rel = delta / max(ref_mag, EPS)
        if max_abs_dev_tol is not None:
            ok = delta <= max_abs_dev_tol
        else:
            ok = max_rel_dev_tol is not None and rel <= max_rel_dev_tol
        return {
            "name": name, "class": cls, "A": a_val, "B": b_val,
            "max_rel_dev": rel, "max_abs_dev": delta,
            "rel_tol": max_rel_dev_tol, "abs_tol": max_abs_dev_tol,
            "endpoint_only": True,
            "status": "PASS" if ok else "FAIL",
        }

    if not a_val or not b_val:
        return {"name": name, "class": cls, "status": "FAIL", "note": "empty trajectory series"}

    rel = _max_rel_dev(a_val, b_val)
    abs_dev = max(abs(_interp_series(b_val, t) - v) for t, v in a_val)
    lag = spec.get("phase_lag_max_intervals")
    measured_lag = _phase_lag_samples(a_val, b_val) if lag is not None else None

    if max_abs_dev_tol is not None:
        ok = abs_dev <= max_abs_dev_tol
    else:
        ok = max_rel_dev_tol is not None and rel <= max_rel_dev_tol
    if lag is not None and measured_lag is not None and abs(measured_lag) > lag:
        ok = False

    row = {
        "name": name, "class": cls, "max_rel_dev": rel, "max_abs_dev": abs_dev,
        "rel_tol": max_rel_dev_tol, "abs_tol": max_abs_dev_tol,
        "endpoint_only": False, "status": "PASS" if ok else "FAIL",
    }
    if lag is not None:
        row["phase_lag_samples"] = measured_lag
        row["phase_lag_tol_intervals"] = lag
    return row


def compare_solver_health(spec: dict[str, Any], a_entry: dict[str, Any] | None,
                           b_entry: dict[str, Any] | None) -> dict[str, Any]:
    name = spec["name"]
    if a_entry is None or b_entry is None:
        return {"name": name, "class": "solver_health", "status": "N/A",
                "note": "observable missing from one or both bundles"}
    a_val, b_val = a_entry["value"], b_entry["value"]

    def summarize(v: Any) -> Any:
        if isinstance(v, list):
            nums = [x for x in v if isinstance(x, (int, float))]
            if not nums:
                return v
            return {"mean": sum(nums) / len(nums), "max": max(nums), "n": len(nums)}
        return v

    a_summary, b_summary = summarize(a_val), summarize(b_val)
    row = {"name": name, "class": "solver_health", "status": "INFO",
           "A": a_summary, "B": b_summary}
    if isinstance(a_summary, dict) and isinstance(b_summary, dict):
        row["delta_mean"] = b_summary["mean"] - a_summary["mean"]
    elif isinstance(a_val, (int, float)) and isinstance(b_val, (int, float)):
        row["delta"] = b_val - a_val
    return row


# ---------------------------------------------------------------------------
# Case / bundle level orchestration
# ---------------------------------------------------------------------------

def compare_case(budget: dict[str, list[dict[str, Any]]], a_metrics: dict[str, Any],
                  b_metrics: dict[str, Any]) -> dict[str, Any]:
    correctness_rows: list[dict[str, Any]] = []
    for spec in budget["correctness"]:
        correctness_rows.extend(compare_correctness(spec, a_metrics.get(spec["name"]), b_metrics.get(spec["name"])))

    engineering_rows = [
        compare_engineering(spec, a_metrics.get(spec["name"]), b_metrics.get(spec["name"]))
        for spec in budget["engineering_trajectory"]
    ]

    solver_rows = [
        compare_solver_health(spec, a_metrics.get(spec["name"]), b_metrics.get(spec["name"]))
        for spec in budget["solver_health"]
    ]

    correctness_fail = any(r["status"] == "FAIL" for r in correctness_rows)
    engineering_fail = any(r["status"] == "FAIL" for r in engineering_rows)

    return {
        "correctness": correctness_rows,
        "engineering_trajectory": engineering_rows,
        "solver_health": solver_rows,
        "correctness_fail": correctness_fail,
        "engineering_fail": engineering_fail,
        "verdict": "FAIL" if correctness_fail else "PASS",
        "gate_verdict": "FAIL" if (correctness_fail or engineering_fail) else "PASS",
    }


def render_markdown(report: dict[str, Any], reference: Path, candidate: Path) -> str:
    lines = [
        "# Physics validation compare report",
        "",
        f"- reference: `{reference}`",
        f"- candidate: `{candidate}`",
        f"- overall verdict: **{report['overall_verdict']}**"
        f" (gate verdict: **{report['overall_gate_verdict']}**)",
        "",
    ]
    for case, case_report in report["cases"].items():
        lines.append(f"## {case} -- {case_report['verdict']} (gate: {case_report['gate_verdict']})")
        lines.append("")
        lines.append("### CORRECTNESS")
        lines.append("")
        lines.append("| observable | A | B | abs Δ | rel Δ | tol | status |")
        lines.append("| --- | --- | --- | --- | --- | --- | --- |")
        for r in case_report["correctness"]:
            if "A" not in r:
                lines.append(f"| {r['name']} | -- | -- | -- | -- | -- | {r['status']} ({r.get('note', '')}) |")
                continue
            tol = r.get("abs_tol") if r.get("abs_tol") is not None else r.get("rel_tol")
            lines.append(
                f"| {r['name']} | {r['A']:.6g} | {r['B']:.6g} | {r['delta_abs']:.3e} | "
                f"{r['delta_rel']:.3e} | {tol:.3e} | {r['status']} |"
            )
        lines.append("")
        lines.append("### ENGINEERING-TRAJECTORY")
        lines.append("")
        lines.append("| observable | max rel Δ | max abs Δ | phase-lag | tol | endpoint-only | status |")
        lines.append("| --- | --- | --- | --- | --- | --- | --- |")
        for r in case_report["engineering_trajectory"]:
            if "max_rel_dev" not in r:
                lines.append(f"| {r['name']} | -- | -- | -- | -- | -- | {r['status']} ({r.get('note', '')}) |")
                continue
            lag = r.get("phase_lag_samples")
            lag_str = f"{lag} (tol {r.get('phase_lag_tol_intervals')})" if lag is not None else "--"
            tol = r.get("rel_tol") if r.get("rel_tol") is not None else r.get("abs_tol")
            lines.append(
                f"| {r['name']} | {r['max_rel_dev']:.3e} | {r.get('max_abs_dev', float('nan')):.3e} | "
                f"{lag_str} | {tol} | {r.get('endpoint_only')} | {r['status']} |"
            )
        lines.append("")
        lines.append("### SOLVER-HEALTH (non-gating)")
        lines.append("")
        lines.append("| observable | A | B | delta |")
        lines.append("| --- | --- | --- | --- |")
        for r in case_report["solver_health"]:
            delta = r.get("delta_mean", r.get("delta", ""))
            lines.append(f"| {r['name']} | {r.get('A')} | {r.get('B')} | {delta} |")
        lines.append("")
    return "\n".join(lines)


def run_compare(budget_path: Path, reference: Path, candidate: Path,
                 cases: list[str] | None, require_compatible_manifest: bool = False) -> dict[str, Any]:
    if require_compatible_manifest:
        check_manifest_compatibility(reference, candidate)
    budget = load_budget(budget_path)
    ref_cases = set(discover_cases(reference))
    cand_cases = set(discover_cases(candidate))
    selected = sorted(cases) if cases else sorted(ref_cases & cand_cases)
    missing = (set(cases) - (ref_cases & cand_cases)) if cases else (ref_cases ^ cand_cases)
    if not selected:
        raise SystemExit(
            f"no common cases between {reference} ({sorted(ref_cases)}) and "
            f"{candidate} ({sorted(cand_cases)})"
        )

    case_reports = {}
    for case in selected:
        a_metrics = load_metrics(reference, case)
        b_metrics = load_metrics(candidate, case)
        case_reports[case] = compare_case(budget, a_metrics, b_metrics)

    overall_verdict = "FAIL" if any(c["verdict"] == "FAIL" for c in case_reports.values()) else "PASS"
    overall_gate_verdict = (
        "FAIL" if any(c["gate_verdict"] == "FAIL" for c in case_reports.values()) else "PASS"
    )

    return {
        "reference": str(reference),
        "candidate": str(candidate),
        "cases": case_reports,
        "cases_skipped_no_match": sorted(missing) if not cases else [],
        "overall_verdict": overall_verdict,
        "overall_gate_verdict": overall_gate_verdict,
    }


# ---------------------------------------------------------------------------
# Self-test (task 1.D Done-when: known-good -> PASS, deliberately-detuned -> FAIL)
# ---------------------------------------------------------------------------

def _write_metrics(path: Path, metrics: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(metrics, indent=2), encoding="utf-8")


def _selftest_budget() -> dict[str, Any]:
    return {
        "correctness": [
            {"name": "eta_field_l2", "rel_tol": 1.0e-6},
            # Regression coverage for the "missing observable" bug: absent from
            # both sides (e.g. fixed_iter never printing MLMG's Final Iter. line)
            # must be N/A/non-gating, not FAIL; absent from only one side (a real
            # divergent solver path) must still FAIL.
            {"name": "absent_both", "rel_tol": 1.0e-6},
            {"name": "absent_one_side", "rel_tol": 1.0e-6},
        ],
        "engineering_trajectory": [
            {"name": "chamber_pressure", "max_rel_dev": 0.02, "phase_lag_max_intervals": 1, "source": "thermo_dat"},
        ],
        "solver_health": [
            {"name": "newton_iters_per_solve"},
        ],
    }


def selftest() -> int:
    tmp = Path(tempfile.mkdtemp(prefix="compare_validation_selftest_"))
    try:
        budget_path = tmp / "budget.yaml"
        budget_path.write_text(yaml.safe_dump(_selftest_budget()), encoding="utf-8")

        reference = tmp / "ref"
        good = tmp / "good"
        bad = tmp / "bad"

        series_ref = [[float(i), 1.0 + 0.001 * i] for i in range(10)]
        series_good = [[float(i), 1.0 + 0.0011 * i] for i in range(10)]  # within 2% everywhere
        series_bad = [[float(i), 1.0 + 0.05 * i] for i in range(10)]  # blows past 2%

        # absent_both is omitted from all three bundles (N/A, non-gating).
        # absent_one_side is present in reference+good (both sides agree ->
        # PASS-eligible) but dropped from bad only (a genuine divergence -> FAIL).
        _write_metrics(reference / "case1" / "metrics.json", {
            "eta_field_l2": {"class": "correctness", "value": 1.234567},
            "absent_one_side": {"class": "correctness", "value": 5.0},
            "chamber_pressure": {"class": "engineering_trajectory", "value": series_ref},
            "newton_iters_per_solve": {"class": "solver_health", "value": [2, 2, 3]},
        })
        _write_metrics(good / "case1" / "metrics.json", {
            "eta_field_l2": {"class": "correctness", "value": 1.234567 + 1.0e-9},
            "absent_one_side": {"class": "correctness", "value": 5.0},
            "chamber_pressure": {"class": "engineering_trajectory", "value": series_good},
            "newton_iters_per_solve": {"class": "solver_health", "value": [2, 3, 3]},
        })
        _write_metrics(bad / "case1" / "metrics.json", {
            "eta_field_l2": {"class": "correctness", "value": 1.34},  # ~8.5% off, well past 1e-6
            "chamber_pressure": {"class": "engineering_trajectory", "value": series_bad},
            "newton_iters_per_solve": {"class": "solver_health", "value": [9, 9, 9]},
        })

        good_report = run_compare(budget_path, reference, good, None)
        bad_report = run_compare(budget_path, reference, bad, None)

        # Manifest preflight: one valid pair passes; each identity mismatch is rejected.
        base_manifest = {"host": "host", "device": "a1000_sm86_strict", "profile": "gpu_strict",
                         "gpu_name": "A1000", "gpu_uuid": "GPU-1", "driver_version": "1",
                         "cuda_compute_cap": "8.6", "build_flags": "flags",
                         "oracle_scripts": {"extract": "abc"}, "command_shape": {"launcher": ["mpiexec"], "profile": "gpu_strict"},
                         "cases": [{"id": "case1", "input": "input", "input_sha256": "i", "overrides_sha256": "o", "max_step": 1, "np": 1}]}
        (reference / "manifest.json").write_text(json.dumps(base_manifest))
        (good / "manifest.json").write_text(json.dumps(base_manifest))
        check_manifest_compatibility(reference, good)
        ok = True
        for key, value in (("device", "other"), ("profile", "gpu_fast"),
                           ("oracle_scripts", {"extract": "different"}),
                           ("cases", [{"id": "case1", "input": "other", "input_sha256": "i", "overrides_sha256": "o", "max_step": 1, "np": 1}])):
            altered = dict(base_manifest)
            altered[key] = value
            (bad / "manifest.json").write_text(json.dumps(altered))
            try:
                check_manifest_compatibility(reference, bad)
                print(f"FAIL: manifest mismatch {key} was accepted", file=sys.stderr)
                ok = False
            except SystemExit:
                pass

        def _row(report: dict[str, Any], name: str) -> dict[str, Any]:
            return next(r for r in report["cases"]["case1"]["correctness"] if r["name"] == name)

        if good_report["overall_verdict"] != "PASS":
            print(f"FAIL: expected good bundle to PASS, got {good_report['overall_verdict']}", file=sys.stderr)
            ok = False
        if bad_report["overall_verdict"] != "FAIL":
            print(f"FAIL: expected detuned bundle to FAIL, got {bad_report['overall_verdict']}", file=sys.stderr)
            ok = False
        if bad_report["cases"]["case1"]["correctness"][0]["status"] != "FAIL":
            print("FAIL: expected eta_field_l2 to fail on detuned bundle", file=sys.stderr)
            ok = False
        if _row(good_report, "absent_both")["status"] != "N/A":
            print("FAIL: expected absent_both to be N/A (non-gating) when missing from both sides", file=sys.stderr)
            ok = False
        if _row(bad_report, "absent_both")["status"] != "N/A":
            print("FAIL: expected absent_both to stay N/A on the detuned bundle too", file=sys.stderr)
            ok = False
        if _row(good_report, "absent_one_side")["status"] != "PASS":
            print("FAIL: expected absent_one_side to PASS when present+matching on both sides", file=sys.stderr)
            ok = False
        if _row(bad_report, "absent_one_side")["status"] != "FAIL":
            print("FAIL: expected absent_one_side to FAIL when missing from only the candidate", file=sys.stderr)
            ok = False
        if bad_report["overall_gate_verdict"] != "FAIL":
            print("FAIL: expected gate verdict FAIL on detuned bundle", file=sys.stderr)
            ok = False

        if ok:
            print("compare_validation.py selftest passed: known-good -> PASS, detuned -> FAIL")
            return 0
        return 1
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reference", type=Path, nargs="?", help="reference bundle directory")
    parser.add_argument("candidate", type=Path, nargs="?", help="candidate bundle directory")
    parser.add_argument("--budget", type=Path, default=DEFAULT_BUDGET)
    parser.add_argument("--case", action="append", default=[], help="restrict to specific case id(s)")
    parser.add_argument("--gate", action="store_true",
                         help="exit non-zero on any ENGINEERING regression too, not just CORRECTNESS")
    parser.add_argument("--require-compatible-manifest", action="store_true",
                        help="reject bundles with different device/profile/input/oracle identity")
    parser.add_argument("--out-md", type=Path, default=None)
    parser.add_argument("--out-json", type=Path, default=None)
    parser.add_argument("--quiet", action="store_true")
    parser.add_argument("--selftest", action="store_true",
                         help="run the known-good/detuned self-test and exit (no bundles needed)")
    args = parser.parse_args()

    if args.selftest:
        return selftest()

    if args.reference is None or args.candidate is None:
        parser.error("reference and candidate bundle directories are required (unless --selftest)")

    report = run_compare(args.budget, args.reference, args.candidate, args.case or None,
                         args.require_compatible_manifest)

    md = render_markdown(report, args.reference, args.candidate)
    out_md = args.out_md or (args.candidate / "compare_report.md")
    out_json = args.out_json or (args.candidate / "compare_report.json")
    out_md.parent.mkdir(parents=True, exist_ok=True)
    out_md.write_text(md, encoding="utf-8")
    out_json.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")

    if not args.quiet:
        print(md)
        print(f"\nwrote {out_md}")
        print(f"wrote {out_json}")
        if report["cases_skipped_no_match"]:
            print(f"NOTE: cases present in only one bundle, skipped: {report['cases_skipped_no_match']}", file=sys.stderr)

    verdict = report["overall_gate_verdict"] if args.gate else report["overall_verdict"]
    return 0 if verdict == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())
