#!/usr/bin/env python3
"""Grid sweep over HTPB fullfeedback (pre_exponential, activation_temperature)
to find a good initial guess for optimize_htpb_fullfeedback_combined.py,
meant to run as a Slurm job array on an HPC cluster (see
scripts/htpb_sweep.slurm) rather than on the local machine.

Evaluates the SAME combined objective (both geometries -- Chorpening 200um
and the internal-group 100um lamina -- at their respective fit pressures,
relative-error residuals concatenated) as optimize_htpb_fullfeedback_combined.py,
just over a fixed grid instead of via least_squares, so the grid point with
the lowest residual norm is a principled --pre-exponential0/
--activation-temperature0 for that script, rather than its current untuned
defaults (0.001, 4000 -- picked before either template was smoke-tested).

Modes
-----
--mode list       print the grid size (use for the Slurm array range,
                   e.g. --array=0-$(( $(python3 scripts/htpb_grid_sweep.py --mode list) - 1 ))
--mode eval       evaluate ONE grid point (selected by --task-id, meant to be
                   $SLURM_ARRAY_TASK_ID) and write its result as JSON under
                   <workdir>/grid_results/task_XXXX.json
--mode aggregate  after all array tasks finish, collect every
                   <workdir>/grid_results/task_*.json, rank by residual norm,
                   and write <workdir>/best_grid_point.json

Grid
----
Default grid is 9 log-spaced pre_exponential values x 7 activation_temperature
values = 63 points, each a full evaluation against both geometries' fit
pressures (11 sims/point, 693 sims total) -- sized for many small, easily
scheduled Slurm array tasks rather than a few large ones (see htpb_sweep.slurm
for the --cpus-per-task=8 rationale). Override with --pre-exponential-grid /
--activation-temperature-grid for a coarser/finer or differently-ranged sweep.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import optimize_htpb_fullfeedback_combined as opt  # noqa: E402

DEFAULT_PRE_EXPONENTIAL_GRID = [1e-5, 3e-5, 1e-4, 3e-4, 1e-3, 3e-3, 1e-2, 3e-2, 1e-1]
DEFAULT_ACTIVATION_TEMPERATURE_GRID = [1000.0, 2000.0, 3000.0, 4000.0, 5000.0, 6000.0, 7000.0]


def build_grid(pre_exp_grid: list[float],
                act_temp_grid: list[float]) -> list[tuple[float, float]]:
    return [(pe, at) for at in act_temp_grid for pe in pre_exp_grid]


def evaluate_point(pre_exponential: float, activation_temperature: float,
                    geometries: list[opt.Geometry],
                    targets_by_geom: dict[str, list[float]],
                    task_dir: Path, lowmach_bin: str | None,
                    cores: int = opt.DEFAULT_CORES) -> dict:
    residuals = []
    rates_by_geom: dict[str, dict[float, float | None]] = {}
    for geom in geometries:
        geom_dir = task_dir / geom.name
        rates = opt.run_sweep(pre_exponential, activation_temperature,
                               geom.fit_pressures, geom_dir, lowmach_bin, geom.template,
                               cores=cores)
        rates_by_geom[geom.name] = rates
        for p, target in zip(geom.fit_pressures, targets_by_geom[geom.name]):
            sim = rates.get(p)
            if sim is None:
                residuals.append(opt.PENALTY_RESIDUAL)
            else:
                residuals.append((sim - target) / target)
    residuals = np.array(residuals)
    return {
        "pre_exponential": pre_exponential,
        "activation_temperature": activation_temperature,
        "rates_mm_s": rates_by_geom,
        "residual_norm": float(np.linalg.norm(residuals)),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--mode", choices=["list", "eval", "aggregate"], required=True)
    parser.add_argument("--workdir", type=Path,
                         help="directory for grid_results/ and best_grid_point.json "
                              "(required for eval/aggregate)")
    parser.add_argument("--task-id", type=int,
                         help="grid point index to evaluate (required for eval, "
                              "e.g. pass $SLURM_ARRAY_TASK_ID)")
    parser.add_argument("--lowmach-bin", help="override LOWMACH_BIN")
    parser.add_argument("--cores", type=int, default=opt.DEFAULT_CORES,
                         help="total MPI ranks for one grid-point job, split across "
                              "whatever pressures run concurrently for it "
                              f"(default: {opt.DEFAULT_CORES})")
    parser.add_argument("--all-pressures", action="store_true",
                         help="evaluate every fit pressure per geometry instead of just "
                              "the cheap min/max endpoints (default: endpoints only)")
    parser.add_argument("--pre-exponential-grid", type=float, nargs="+",
                         default=DEFAULT_PRE_EXPONENTIAL_GRID)
    parser.add_argument("--activation-temperature-grid", type=float, nargs="+",
                         default=DEFAULT_ACTIVATION_TEMPERATURE_GRID)
    parser.add_argument("--points-file", type=Path,
                         help="JSON list of {pre_exponential, activation_temperature} "
                              "objects (e.g. from generate_next_batch.py) to evaluate "
                              "instead of the cartesian pre/activation grid")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.points_file:
        with open(args.points_file) as fh:
            points = json.load(fh)
        grid = [(float(p["pre_exponential"]), float(p["activation_temperature"])) for p in points]
    else:
        grid = build_grid(args.pre_exponential_grid, args.activation_temperature_grid)

    if args.mode == "list":
        print(len(grid))
        return

    if args.workdir is None:
        raise SystemExit("--workdir is required for --mode eval/aggregate")

    geometries = opt.DEFAULT_GEOMETRIES if args.all_pressures \
        else opt.endpoint_geometries(opt.DEFAULT_GEOMETRIES)
    data_by_geom = {g.name: opt.load_experimental_data(g.data) for g in geometries}
    targets_by_geom = {g.name: opt.nearest_experimental(data_by_geom[g.name], g.fit_pressures)
                       for g in geometries}

    if args.mode == "eval":
        if args.task_id is None:
            raise SystemExit("--task-id is required for --mode eval")
        if not (0 <= args.task_id < len(grid)):
            raise SystemExit(f"--task-id {args.task_id} out of range [0, {len(grid)})")
        pre_exponential, activation_temperature = grid[args.task_id]

        results_dir = args.workdir / "grid_results"
        results_dir.mkdir(parents=True, exist_ok=True)
        task_dir = args.workdir / "grid_points" / f"task_{args.task_id:04d}"
        task_dir.mkdir(parents=True, exist_ok=True)

        print(f"[task {args.task_id}] pre_exponential={pre_exponential:.6g} "
              f"activation_temperature={activation_temperature:.6g}", flush=True)
        record = evaluate_point(pre_exponential, activation_temperature,
                                geometries, targets_by_geom, task_dir, args.lowmach_bin,
                                cores=args.cores)
        record["task_id"] = args.task_id
        out_path = results_dir / f"task_{args.task_id:04d}.json"
        with open(out_path, "w") as fh:
            json.dump(record, fh, indent=2)
        print(f"[task {args.task_id}] residual_norm={record['residual_norm']:.6g} "
              f"-> wrote {out_path}", flush=True)
        return

    if args.mode == "aggregate":
        results_dir = args.workdir / "grid_results"
        records = []
        for path in sorted(results_dir.glob("task_*.json")):
            with open(path) as fh:
                records.append(json.load(fh))
        if not records:
            raise SystemExit(f"no task_*.json files found under {results_dir}")
        records.sort(key=lambda r: r["residual_norm"])

        print(f"{'task':>5} {'pre_exponential':>16} {'activation_temperature':>22} "
              f"{'residual_norm':>14}")
        for r in records[:15]:
            print(f"{r['task_id']:>5} {r['pre_exponential']:>16.6g} "
                  f"{r['activation_temperature']:>22.6g} {r['residual_norm']:>14.6g}")

        best = records[0]
        best_out = {
            "pre_exponential0": best["pre_exponential"],
            "activation_temperature0": best["activation_temperature"],
            "residual_norm": best["residual_norm"],
            "task_id": best["task_id"],
        }
        out_path = args.workdir / "best_grid_point.json"
        with open(out_path, "w") as fh:
            json.dump(best_out, fh, indent=2)
        print(f"\nBest grid point: pre_exponential={best['pre_exponential']:.6g} "
              f"activation_temperature={best['activation_temperature']:.6g} "
              f"residual_norm={best['residual_norm']:.6g}")
        print(f"Wrote {out_path}")
        print(f"\nResume the real optimization with this as the initial guess:\n"
              f"  python scripts/optimize_htpb_fullfeedback_combined.py "
              f"--workdir <workdir> --pre-exponential0 {best['pre_exponential']:.6g} "
              f"--activation-temperature0 {best['activation_temperature']:.6g}")
        return


if __name__ == "__main__":
    main()
