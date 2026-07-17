#!/usr/bin/env python3
"""Run a non-mutating ElasticSoftVoid diagnostic matrix.

All run products are placed outside the repository by default.  This is a
diagnostic harness: the acceptance oracle in tests/ElasticSoftVoid/test is
not changed or bypassed.
"""

from __future__ import annotations

import argparse
import json
import re
import subprocess
import time
from pathlib import Path

import numpy as np
import yt

yt.set_log_level(50)

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "tests/ElasticSoftVoid/input"


def cases():
    result = []
    for dim in (2, 3):
        for modulus, value in (("soft", "0.2_MPa"), ("hard", "10_MPa")):
            for thermal in ("on", "off"):
                for amr in ("on", "off"):
                    result.append({
                        "dim": dim, "modulus": modulus, "modulus_value": value,
                        "thermal": thermal, "amr": amr,
                    })
    return result


def metric_rows(stdout):
    pat = re.compile(
        r"NR convergence metrics:.*?max_update=([^,\s]+).*?"
        r"nonlinear_resid=([^,\s]+),\s*nonlinear_resid_rel=([^,\s]+),\s*"
        r"converged=([01])")
    return [
        {"update": float(a), "nonlinear_residual": float(b),
         "nonlinear_ratio": float(c), "converged": bool(int(d))}
        for a, b, c, d in pat.findall(stdout)
    ]


def mlmg_rows(stdout):
    pat = re.compile(
        r"MLMG: Final Iter\.\s+(\d+)\s+resid,\s+"
        r"resid/(?:bnorm|resid0)\s*=\s*([^,\s]+),\s*([^,\s]+)")
    return [{"iterations": int(n), "absolute": float(a), "relative": float(r)}
            for n, a, r in pat.findall(stdout)]


def classification(coord, level, level_lo, level_hi, domain_lo, domain_hi, dim):
    eps = 2.0e-9 * max(1.0, float(np.max(domain_hi - domain_lo)))
    physical_axes = []
    outer_axes = []
    # z is a periodic extrusion in the 3-D sections, so only x/y are
    # physical boundaries for this chamber.
    for axis in range(min(dim, 2)):
        if abs(coord[axis] - domain_lo[axis]) <= eps or abs(coord[axis] - domain_hi[axis]) <= eps:
            physical_axes.append(axis)
        if level > 0 and (abs(coord[axis] - level_lo[axis]) <= eps or
                          abs(coord[axis] - level_hi[axis]) <= eps):
            outer_axes.append(axis)
    if physical_axes:
        return "physical-boundary"
    if outer_axes:
        return "outer coarse/fine-union row"
    return "interior"


def plot_diagnostics(plot):
    ds = yt.load(str(plot))
    dim = ds.dimensionality
    domain_lo = np.asarray(ds.domain_left_edge, dtype=float)
    domain_hi = np.asarray(ds.domain_right_edge, dtype=float)
    residual_fields = [f for _, f in ds.field_list if f.startswith("res_")]
    rhs_fields = [f for _, f in ds.field_list if f.startswith("rhs_")]
    best = None
    rhs_max = 0.0
    for grid in ds.index.grids:
        level = int(grid.Level)
        lo = np.asarray(grid.LeftEdge, dtype=float)
        dd = np.asarray(grid.dds, dtype=float)
        for field in rhs_fields:
            rhs_max = max(rhs_max, float(np.max(np.abs(np.asarray(grid[("boxlib", field)])))))
        for field in residual_fields:
            values = np.asarray(grid[("boxlib", field)], dtype=float)
            index = np.unravel_index(int(np.nanargmax(np.abs(values))), values.shape)
            value = float(values[index])
            candidate = (abs(value), field, value, grid, index, lo, dd)
            if best is None or candidate[0] > best[0]:
                best = candidate

    by_level = {}
    for grid in ds.index.grids:
        by_level.setdefault(int(grid.Level), []).append(grid)
    level = best[3].Level
    level_grids = by_level[int(level)]
    level_lo = np.min([np.asarray(g.LeftEdge, dtype=float) for g in level_grids], axis=0)
    level_hi = np.max([np.asarray(g.RightEdge, dtype=float) for g in level_grids], axis=0)
    _, field, value, grid, index, lo, dd = best
    coord = lo + np.asarray(index, dtype=float) * dd
    # Plotfile nodal fields use the left edge plus integer grid indices.
    coord = coord[:dim]
    mu = float("nan")
    kap = float("nan")
    for name in ("model_mu", "_mu"):
        if ("boxlib", name) in grid.ds.field_list:
            mu = float(np.min(np.asarray(grid[("boxlib", name)])))
            break
    for name in ("model_kappa", "_kappa"):
        if ("boxlib", name) in grid.ds.field_list:
            kap = float(np.min(np.asarray(grid[("boxlib", name)])))
            break
    return {
        "plot": str(plot), "max_residual": best[0], "residual_field": field,
        "signed_residual": value, "level": int(level), "index": list(map(int, index)),
        "coordinates": coord.tolist(),
        "classification": classification(coord, int(level), level_lo[:dim], level_hi[:dim],
                                           domain_lo[:dim], domain_hi[:dim], dim),
        "rhs_max": rhs_max, "ratio": best[0] / rhs_max if rhs_max else float("inf"),
        "level_union_lo": level_lo[:dim].tolist(), "level_union_hi": level_hi[:dim].tolist(),
        "mu_min_on_max_grid": mu, "kappa_min_on_max_grid": kap,
        "amr_max_level": int(ds.index.max_level),
    }


def run_case(case, root, timeout, max_step):
    label = "{dim}d-{modulus}-{thermal}-thermal-{amr}-amr".format(**case)
    out = root / label
    out.mkdir(parents=True, exist_ok=True)
    dim = case["dim"]
    args = [
        str(ROOT / f"bin/alamo-{dim}d-g++"), str(INPUT),
        f"plot_file={out / 'plot'}", f"max_step={max_step}", "allow_unused=true",
        # Keep thermal transport enabled in both branches.  Flame's phase
        # kinetics require a finite temperature even when transport is off;
        # thermal.on=0 is therefore not a controlled elastic-loading toggle.
        "thermal.on=1",
        "pf.eta.ic.expression.constant.w=0.002",
        f"model_void.kappa={case['modulus_value']}",
        f"model_void.mu={case['modulus_value']}",
    ]
    if case["thermal"] == "on":
        # Flame interprets F0-I as thermal expansion per kelvin before it
        # multiplies by T-Telastic.  alpha=1e-5/K and a 300 K offset impose a
        # controlled 0.3% eigenstrain without changing the thermal evolution.
        args += ["Telastic=0"]
        for material in ("model_prop", "model_void", "model_casing"):
            if dim == 2:
                args += [f"{material}.F0=1.00001 0 0 1.00001"]
            else:
                args += [f"{material}.F0=1.00001 0 0 0 1.00001 0 0 0 1.00001"]
    else:
        args += ["Telastic=300"]
    if dim == 2:
        args += ["amr.n_cell=64 64 8", "amr.max_level=2"]
        if case["amr"] == "on":
            args += ["explicitmesh.lo1=32 32 0", "explicitmesh.hi1=95 95 0",
                     "explicitmesh.lo2=72 72 0", "explicitmesh.hi2=183 183 0"]
    else:
        args += ["amr.n_cell=64 64 8", "amr.max_level=1"]
        # The 3-D diagnostic is an exact periodic extrusion of the 2-D
        # chamber.  Project linear corrections onto that declared invariant
        # subspace so roundoff in the redundant z planes is not reported as a
        # physical symmetry error.
        args += ["elastic.solver.invariant_periodic=0 0 1"]
        if case["amr"] == "on":
            args += ["explicitmesh.lo1=32 32 0", "explicitmesh.hi1=95 95 15"]
    if case["amr"] == "off":
        args += ["amr.max_level=0", "explicitmesh.on=0"]
    started = time.monotonic()
    try:
        proc = subprocess.run(args, cwd=ROOT, text=True, capture_output=True,
                              timeout=timeout)
        timed_out = False
    except subprocess.TimeoutExpired as exc:
        proc = exc
        timed_out = True
    stdout = getattr(proc, "stdout", "") or ""
    stderr = getattr(proc, "stderr", "") or ""
    (out / "stdout").write_text(stdout)
    (out / "stderr").write_text(stderr)
    result = dict(case, label=label, returncode=None if timed_out else proc.returncode,
                  timeout=timed_out, seconds=time.monotonic() - started,
                  newton=metric_rows(stdout), mlmg=mlmg_rows(stdout),
                  status="timeout" if timed_out else ("complete" if proc.returncode == 0 else "abort"))
    plots = sorted(out.glob("plot/*node"))
    if plots:
        try:
            result["plot_diagnostics"] = plot_diagnostics(plots[-1])
        except Exception as exc:
            result["plot_diagnostics_error"] = repr(exc)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path("/tmp/alamo-elastic-void-matrix"))
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--max-step", type=int, default=2)
    parser.add_argument("--case", action="append", help="run only labels containing this text")
    args = parser.parse_args()
    args.root.mkdir(parents=True, exist_ok=True)
    selected = [c for c in cases() if not args.case or any(x in
               f"{c['dim']}d-{c['modulus']}-{c['thermal']}-thermal-{c['amr']}-amr"
               for x in args.case)]
    results = [run_case(c, args.root, args.timeout, args.max_step) for c in selected]
    report = args.root / "matrix.json"
    report.write_text(json.dumps(results, indent=2, allow_nan=True) + "\n")
    print("| case | status | max residual / rhs | Newton ratio | MLMG iters | location |")
    print("|---|---|---:|---:|---:|---|")
    for result in results:
        diag = result.get("plot_diagnostics", {})
        ratio = diag.get("ratio", float("nan"))
        nr = result.get("newton", [])
        nr_ratio = nr[-1]["nonlinear_ratio"] if nr else float("nan")
        iters = ",".join(str(x["iterations"]) for x in result.get("mlmg", [])) or "-"
        loc = (f"L{diag.get('level','-')} {diag.get('classification','-')} "
               f"{diag.get('coordinates','-')}")
        print(f"| {result['label']} | {result['status']} | {ratio:.6g} | "
              f"{nr_ratio:.6g} | {iters} | {loc} |")
    print(f"\nJSON report: {report}")


if __name__ == "__main__":
    main()
