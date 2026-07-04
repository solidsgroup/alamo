#!/usr/bin/env python3
"""Turn one ALAMO run's raw output into metrics.json + field_norms.json.

Roadmap task 1.B support module (`benchmark/GPU_ROADMAP_V3.md` Phase 1). This
is the one place that knows how to turn `run.log` / `thermo.dat` / the final
AMReX plotfile into the observable values `physics_budget.yaml` (task 1.A)
defines and `compare_validation.py` (task 1.D) compares. Called by
`run_validation_local.sh` (and, on NOVA, `run_validation_nova.slurm`) once per
case, after the run completes, with the case's bundle directory as input.

Expects the case directory to already contain (this is what
`run_validation_local.sh` lays down):
    run.log                  -- full stdout of the alamo binary
    thermo.dat                -- copied straight from the run's plot_file dir
    <N>node/, <N>cell/        -- one or more AMReX plotfile pairs (ALAMO's own
                                 node/cell split), `<N>` zero-padded step number.
                                 The highest `<N>` is treated as "the final
                                 plotfile" for CORRECTNESS / derived-field norms.

Writes metrics.json (the canonical comparison surface, schema in
benchmark/validate/README.md) and field_norms.json (raw per-field norms,
kept separately for human/regression inspection) into the same directory.

Known scoping limits (documented, not silently swept under the rug):
  - SOLVER-HEALTH aggregates are computed over the *whole* run.log, not
    segmented per individual elastic solve -- acceptable because that class is
    explicitly non-gating regression telemetry (physics_budget.md), and a
    validation case is short (1-2 solves) so the distinction rarely matters.
  - field_norms.json stores norms of the *final* plotfile only, so
    ENGINEERING-TRAJECTORY observables sourced from `derived_field` (von Mises
    / principal stress / elastic energy) are endpoint-only comparisons, not
    true trajectories -- this is the same limitation physics_budget.md flags
    for those three observables.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path
from typing import Any

import numpy as np
import yaml

VALIDATE_DIR = Path(__file__).resolve().parent
DEFAULT_BUDGET = VALIDATE_DIR / "physics_budget.yaml"

# plotfile field name -> (registered name, component suffix) per src/Integrator/*.H
# RegisterGeneralFab/RegisterNodalFab field naming: scalar "name", vector
# "name_x"/"name_y"/"name_z", matrix "name_xx"/"name_xy"/.../"name_zz".
SCALAR_FIELDS = ["eta", "phi", "temp"]
VECTOR_FIELDS = ["disp"]
MATRIX_FIELDS = ["stress", "strain"]
AXES = ["x", "y", "z"]


# ---------------------------------------------------------------------------
# thermo.dat
# ---------------------------------------------------------------------------

def read_thermo(path: Path) -> dict[str, list[list[float]]]:
    with path.open(encoding="utf-8") as fh:
        header = fh.readline().split()
        columns: dict[str, list[list[float]]] = {name: [] for name in header}
        time_idx = header.index("time")
        for line in fh:
            if not line.strip():
                continue
            row = [float(v) for v in line.split()]
            t = row[time_idx]
            for name, value in zip(header, row):
                columns[name].append([t, value])
    return columns


# ---------------------------------------------------------------------------
# run.log
# ---------------------------------------------------------------------------

def _aggregate(kind: str, matches: list[str]) -> Any:
    if not matches:
        return None
    if kind == "count_per_solve":
        return len(matches)
    if kind == "last_per_solve":
        return float(matches[-1])
    if kind == "max_per_solve":
        return max(int(m) for m in matches)
    if kind == "list_per_vcycle":
        return [int(m) for m in matches]
    if kind == "monotonic_nonincreasing_per_solve":
        floats = [float(m) for m in matches]
        return all(a >= b - 1.0e-300 for a, b in zip(floats, floats[1:]))
    raise ValueError(f"unknown aggregate kind: {kind}")


def parse_run_log(path: Path, budget: dict[str, list[dict[str, Any]]]) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8", errors="replace")
    results: dict[str, Any] = {}
    for spec in budget["solver_health"]:
        pattern = re.compile(spec["pattern"])
        matches = pattern.findall(text)
        results[spec["name"]] = _aggregate(spec["aggregate"], matches)
    for spec in budget["correctness"]:
        if spec.get("source") != "run_log":
            continue
        pattern = re.compile(spec["pattern"])
        matches = pattern.findall(text)
        results[spec["name"]] = float(matches[-1]) if matches else None
    return results


# ---------------------------------------------------------------------------
# plotfile field norms (final step only)
# ---------------------------------------------------------------------------

def find_final_plotfile(case_dir: Path) -> Path | None:
    node_dirs = sorted(
        (p for p in case_dir.glob("[0-9]*node") if p.is_dir() and (p / "Header").exists()),
        key=lambda p: p.name,
    )
    return node_dirs[-1] if node_dirs else None


def _l2_linf(arr: np.ndarray, vol: np.ndarray) -> dict[str, float]:
    l2 = float(np.sqrt(np.sum((arr.astype(np.float64) ** 2) * vol)))
    linf = float(np.max(np.abs(arr))) if arr.size else 0.0
    return {"l2": l2, "linf": linf}


def compute_field_norms(node_plotfile: Path) -> dict[str, Any]:
    import yt  # heavy import, kept local so --selftest/help don't need it

    yt.set_log_level(50)  # suppress yt's INFO chatter; this script's own output should be the signal
    ds = yt.load(str(node_plotfile))
    ad = ds.all_data()
    vol = ad[("index", "cell_volume")].to_ndarray()
    available = {name for _, name in ds.field_list}

    def get(field: str) -> np.ndarray | None:
        if field not in available:
            return None
        return ad[("boxlib", field)].to_ndarray()

    norms: dict[str, Any] = {"dimensionality": int(ds.dimensionality)}

    for name in SCALAR_FIELDS:
        arr = get(name)
        if arr is not None:
            norms[name] = _l2_linf(arr, vol)

    vector_arrays: dict[str, dict[str, np.ndarray]] = {}
    for name in VECTOR_FIELDS:
        comps = {}
        for axis in AXES:
            arr = get(f"{name}_{axis}")
            if arr is not None:
                comps[axis] = arr
        if comps:
            norms[name] = {axis: _l2_linf(arr, vol) for axis, arr in comps.items()}
            vector_arrays[name] = comps

    matrix_arrays: dict[str, dict[str, np.ndarray]] = {}
    for name in MATRIX_FIELDS:
        comps = {}
        for a in AXES:
            for b in AXES:
                arr = get(f"{name}_{a}{b}")
                if arr is not None:
                    comps[f"{a}{b}"] = arr
        if comps:
            norms[name] = {comp: _l2_linf(arr, vol) for comp, arr in comps.items()}
            matrix_arrays[name] = comps

    norms["derived"] = compute_derived(matrix_arrays, vol)
    return norms


def _sym(comps: dict[str, np.ndarray], a: str, b: str) -> np.ndarray:
    """Average comps[ab] and comps[ba] if both present (ALAMO stores the full matrix,
    not just the symmetric part, so xy/yx can differ at FP noise level)."""
    ab, ba = comps.get(a + b), comps.get(b + a)
    if ab is not None and ba is not None:
        return 0.5 * (ab + ba)
    if ab is not None:
        return ab
    if ba is not None:
        return ba
    return None


def compute_derived(matrix_arrays: dict[str, dict[str, np.ndarray]], vol: np.ndarray) -> dict[str, Any]:
    stress = matrix_arrays.get("stress")
    if not stress:
        return {}
    n = next(iter(stress.values())).shape[0]
    zero = np.zeros(n)
    sxx = stress.get("xx", zero)
    syy = stress.get("yy", zero)
    szz = stress.get("zz", zero)
    sxy = _sym(stress, "x", "y")
    sxy = sxy if sxy is not None else zero
    syz = _sym(stress, "y", "z")
    syz = syz if syz is not None else zero
    szx = _sym(stress, "z", "x")
    szx = szx if szx is not None else zero

    von_mises = np.sqrt(0.5 * ((sxx - syy) ** 2 + (syy - szz) ** 2 + (szz - sxx) ** 2
                                + 6.0 * (sxy ** 2 + syz ** 2 + szx ** 2)))
    total_vol = float(np.sum(vol))
    derived: dict[str, Any] = {
        "von_mises": {
            "max": float(np.max(von_mises)) if von_mises.size else 0.0,
            "mean": float(np.sum(von_mises * vol) / total_vol) if total_vol > 0 else 0.0,
        }
    }

    # principal stress: largest eigenvalue of the per-cell 3x3 symmetric stress tensor
    # (2D cases naturally have zz/xz/yz == 0, so this degrades correctly to the 2D case)
    mats = np.zeros((n, 3, 3))
    mats[:, 0, 0], mats[:, 1, 1], mats[:, 2, 2] = sxx, syy, szz
    mats[:, 0, 1] = mats[:, 1, 0] = sxy
    mats[:, 1, 2] = mats[:, 2, 1] = syz
    mats[:, 2, 0] = mats[:, 0, 2] = szx
    eigvals = np.linalg.eigvalsh(mats)  # ascending order per cell
    max_principal = eigvals[:, -1]
    derived["max_principal"] = {"max": float(np.max(max_principal)) if max_principal.size else 0.0}

    strain = matrix_arrays.get("strain")
    if strain:
        exx = strain.get("xx", zero)
        eyy = strain.get("yy", zero)
        ezz = strain.get("zz", zero)
        exy = _sym(strain, "x", "y")
        exy = exy if exy is not None else zero
        eyz = _sym(strain, "y", "z")
        eyz = eyz if eyz is not None else zero
        ezx = _sym(strain, "z", "x")
        ezx = ezx if ezx is not None else zero
        # stress:strain double contraction, off-diagonals counted twice (symmetric tensor)
        contraction = (sxx * exx + syy * eyy + szz * ezz
                        + 2.0 * (sxy * exy + syz * eyz + szx * ezx))
        derived["elastic_energy"] = float(0.5 * np.sum(contraction * vol))

    return derived


# ---------------------------------------------------------------------------
# metrics.json assembly (matches the schema compare_validation.py expects)
# ---------------------------------------------------------------------------

def build_metrics(budget: dict[str, list[dict[str, Any]]], thermo: dict[str, list[list[float]]],
                   run_log: dict[str, Any], field_norms: dict[str, Any]) -> dict[str, Any]:
    metrics: dict[str, Any] = {}

    for spec in budget["correctness"]:
        name = spec["name"]
        if spec.get("source") == "run_log":
            value = run_log.get(name)
            if value is not None:
                metrics[name] = {"class": "correctness", "value": value}
            continue
        field = spec["field"]
        norm = spec["norm"]
        entry = field_norms.get(field)
        if entry is None:
            continue
        if spec.get("per_component"):
            value = {comp: data[norm] for comp, data in entry.items() if isinstance(data, dict)}
        else:
            value = entry[norm]
        metrics[name] = {"class": "correctness", "value": value}

    for spec in budget["engineering_trajectory"]:
        name = spec["name"]
        if spec.get("source") == "thermo_dat":
            column = spec.get("column")
            if column:
                series = thermo.get(column)
                if series is not None:
                    metrics[name] = {"class": "engineering_trajectory", "value": series}
            else:
                cols = spec.get("columns", [])
                value = {c: thermo[c] for c in cols if c in thermo}
                if value:
                    metrics[name] = {"class": "engineering_trajectory", "value": value}
        elif spec.get("source") == "derived_field":
            derivation = spec["derivation"]
            reduce_ = spec["reduce"]
            derived = field_norms.get("derived", {})
            if derivation == "von_mises":
                value = derived.get("von_mises", {}).get(reduce_)
            elif derivation == "max_principal":
                value = derived.get("max_principal", {}).get(reduce_)
            elif derivation == "strain_energy":
                value = derived.get("elastic_energy")
            else:
                value = None
            if value is not None:
                metrics[name] = {"class": "engineering_trajectory", "value": value, "endpoint_only": True}

    for spec in budget["solver_health"]:
        name = spec["name"]
        value = run_log.get(name)
        if value is not None:
            metrics[name] = {"class": "solver_health", "value": value}

    return metrics


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("case_dir", type=Path, help="bundle case directory (contains run.log, thermo.dat, <N>node/)")
    parser.add_argument("--budget", type=Path, default=DEFAULT_BUDGET)
    parser.add_argument("--skip-field-norms", action="store_true",
                         help="skip the yt plotfile pass (faster iteration on thermo/log-only changes)")
    args = parser.parse_args()

    budget = yaml.safe_load(args.budget.read_text(encoding="utf-8"))
    for key in ("correctness", "engineering_trajectory", "solver_health"):
        budget.setdefault(key, [])

    thermo_path = args.case_dir / "thermo.dat"
    thermo = read_thermo(thermo_path) if thermo_path.exists() else {}

    run_log_path = args.case_dir / "run.log"
    run_log = parse_run_log(run_log_path, budget) if run_log_path.exists() else {}

    field_norms: dict[str, Any] = {}
    if not args.skip_field_norms:
        final_plotfile = find_final_plotfile(args.case_dir)
        if final_plotfile is not None:
            field_norms = compute_field_norms(final_plotfile)
            (args.case_dir / "field_norms.json").write_text(
                json.dumps(field_norms, indent=2, sort_keys=True), encoding="utf-8"
            )
        else:
            print(f"WARNING: no <N>node plotfile found under {args.case_dir}, skipping field norms",
                  file=sys.stderr)

    metrics = build_metrics(budget, thermo, run_log, field_norms)
    (args.case_dir / "metrics.json").write_text(
        json.dumps(metrics, indent=2, sort_keys=True), encoding="utf-8"
    )
    print(f"wrote {args.case_dir / 'metrics.json'} ({len(metrics)} observables)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
