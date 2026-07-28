#!/usr/bin/env python3
"""Fit ONE HTPB fullfeedback (pre_exponential, activation_temperature) pair
against BOTH AP/HTPB sandwich datasets/geometries at once.

Mirrors the single-geometry driver (optimize_htpb_fullfeedback.py), but each
objective evaluation runs pressure sweeps against *two* geometries -- the
Chorpening et al. 2000 ~200 um-lamina case (input.lm.ap_htpb_fullfeedback.template,
HTPB_sandwich_reg_rate_chorpening2000_summary.csv) and the internal-group
100 um-lamina case (input.lm.ap_htpb_fullfeedback_100um.template,
HTPB_sandwich_reg_rate_group_100um.csv) -- and concatenates both datasets'
relative residuals into one vector for least_squares. This is the same
"one parameter pair across the whole curve" pattern used for
AP_decomposition against AP_reg_rate.csv (FULLFEEDBACK_CALIBRATION.md),
just extended across two physical geometries instead of one pressure range.

AP_decomposition's fullfeedback parameters are baked into both templates
already (held fixed at the values calibrated against pure AP monopropellant
data) and are not touched here.

Two-stage fit: since both datasets are smooth in P, stage 1 fits only the
lowest and highest pressure per geometry (4 sims/iteration instead of 11) --
a fast, cheap optimization that should land close to the real optimum given
the smoothness. Stage 2 then runs the real fit over every pressure, starting
from stage 1's result instead of the untuned defaults. Pass --skip-stage1
(with --pre-exponential0/--activation-temperature0) to go straight to
stage 2 from an already-known-good starting point.

Optimizer: Gaussian-process-based Bayesian optimization (scikit-optimize's
gp_minimize) instead of a gradient/finite-difference method. Each objective
call is one expensive multi-sim evaluation with no wasted finite-difference
probes -- a GP surrogate models residual_norm(log10(pre_exponential),
activation_temperature) from every point evaluated so far, and Expected
Improvement picks the next point to try. Much more sample-efficient than
least_squares/trf here, since every evaluation counts and the surrogate
can route around parameter regions that turned out to be numerically
expensive (e.g. very stiff/slow high-pressure sims) once they're observed.

This is a *long-running, expensive* driver -- each stage-2 objective
evaluation runs 11 sims total (4 Chorpening + 7 group, sequentially by
geometry to avoid oversubscribing the machine, each internally parallel
across its own pressures) and each sim is a full 2D sandwich run. Expect
multi-day wall-clock for a real optimization; run under nohup/tmux/sbatch.
Progress streams to stdout and to <workdir>/iterations_stage1.jsonl /
<workdir>/iterations.jsonl.

Example
-------
    python scripts/optimize_htpb_fullfeedback_combined.py --workdir /tmp/htpb_combined_calib
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from skopt import gp_minimize
from skopt.space import Real

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent

DEFAULT_PRE_EXPONENTIAL0 = 0.001
DEFAULT_ACTIVATION_TEMPERATURE0 = 4000.0

PENALTY_RESIDUAL = 5.0  # relative-error stand-in for a failed/non-igniting sim

# Total MPI ranks to spread across whatever sims are running concurrently in
# a given stage (leaves a couple cores free for the OS/driver process).
DEFAULT_CORES = max((os.cpu_count() or 4) - 2, 2)


@dataclass
class Geometry:
    name: str
    template: Path
    data: Path
    fit_pressures: list[float]


# Only fitting the 100um internal-group dataset for now (per latest
# instructions) -- chorpening_200um is left defined but unused so it's easy
# to re-add later.
CHORPENING_200UM = Geometry(
    name="chorpening_200um",
    template=REPO_ROOT / "input.lm.ap_htpb_fullfeedback.template",
    data=REPO_ROOT / "HTPB_sandwich_reg_rate_chorpening2000_summary.csv",
    fit_pressures=[0.213, 0.467, 1.498, 3.098],
)

GROUP_100UM = Geometry(
    name="group_100um",
    template=REPO_ROOT / "input.lm.ap_htpb_fullfeedback_100um.template",
    data=REPO_ROOT / "HTPB_sandwich_reg_rate_group_100um.csv",
    fit_pressures=[0.8, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0],
)

DEFAULT_GEOMETRIES = [GROUP_100UM]


def endpoint_geometries(geometries: list[Geometry]) -> list[Geometry]:
    """Same template/data, but fit_pressures trimmed to each geometry's
    min and max pressure only -- for the cheap stage-1 fit."""
    return [
        Geometry(name=g.name, template=g.template, data=g.data,
                 fit_pressures=sorted({min(g.fit_pressures), max(g.fit_pressures)}))
        for g in geometries
    ]


def load_experimental_data(csv_path: Path) -> dict[float, float]:
    """Average duplicate/near-duplicate pressures."""
    rows = []
    with open(csv_path) as fh:
        reader = csv.reader(fh)
        next(reader)  # header
        for pressure, rate in reader:
            rows.append((float(pressure), float(rate)))
    data: dict[float, list[float]] = {}
    for pressure, rate in rows:
        data.setdefault(pressure, []).append(rate)
    return {p: float(np.mean(v)) for p, v in data.items()}


def nearest_experimental(data: dict[float, float], pressures: list[float]) -> list[float]:
    known_p = np.array(sorted(data))
    known_r = np.array([data[p] for p in known_p])
    targets = []
    for p in pressures:
        if p in data:
            targets.append(data[p])
        else:
            targets.append(float(np.interp(p, known_p, known_r)))
    return targets


def run_sweep(pre_exponential: float, activation_temperature: float,
              pressures: list[float], workdir: Path,
              lowmach_bin: str | None, template: Path,
              cores: int = DEFAULT_CORES) -> dict[float, float | None]:
    results_csv = workdir / "results.csv"
    env = os.environ.copy()
    if lowmach_bin:
        env["LOWMACH_BIN"] = lowmach_bin
    env["TEMPLATE"] = str(template)
    # Spread all available cores across however many pressures run
    # concurrently in this sweep, so a 2-pressure (endpoints) stage uses
    # much bigger MPI jobs per sim than a 7-pressure (full) stage.
    env["MPI_NP"] = str(max(1, cores // max(1, len(pressures))))
    cmd = [
        str(SCRIPT_DIR / "run_htpb_sandwich_pressure_sweep.sh"),
        repr(pre_exponential), repr(activation_temperature),
        str(workdir), str(results_csv),
        *[repr(p) for p in pressures],
    ]
    subprocess.run(cmd, check=True, env=env)
    rates: dict[float, float | None] = {}
    with open(results_csv) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            p = float(row["pressure_mpa"])
            r = row["reg_rate_mm_s"]
            rates[p] = float(r) if r else None
    return rates


def make_objective(geometries: list[Geometry], targets_by_geom: dict[str, list[float]],
                    workdir: Path, lowmach_bin: str | None, log_path: Path,
                    stage_label: str = "iter", cores: int = DEFAULT_CORES):
    iteration = [0]

    def objective(x: np.ndarray) -> float:
        iteration[0] += 1
        pre_exponential = float(10.0 ** x[0])
        activation_temperature = float(x[1])
        eval_dir = workdir / f"{stage_label}_{iteration[0]:03d}"
        eval_dir.mkdir(parents=True, exist_ok=True)

        residuals = []
        rates_by_geom: dict[str, dict[float, float | None]] = {}
        for geom in geometries:
            geom_dir = eval_dir / geom.name
            rates = run_sweep(pre_exponential, activation_temperature,
                               geom.fit_pressures, geom_dir, lowmach_bin, geom.template,
                               cores=cores)
            rates_by_geom[geom.name] = rates
            for p, target in zip(geom.fit_pressures, targets_by_geom[geom.name]):
                sim = rates.get(p)
                if sim is None:
                    residuals.append(PENALTY_RESIDUAL)
                else:
                    residuals.append((sim - target) / target)
        residuals = np.array(residuals)
        residual_norm = float(np.linalg.norm(residuals))

        record = {
            "iteration": iteration[0],
            "pre_exponential": pre_exponential,
            "activation_temperature": activation_temperature,
            "rates_mm_s": rates_by_geom,
            "residual_norm": residual_norm,
        }
        with open(log_path, "a") as fh:
            fh.write(json.dumps(record) + "\n")
        print(f"[{stage_label} {iteration[0]}] pre_exponential={pre_exponential:.6g} "
              f"activation_temperature={activation_temperature:.6g} "
              f"residual_norm={residual_norm:.6g}", flush=True)
        return residual_norm

    return objective


def load_prior_points(prior_dir: Path) -> tuple[list[list[float]], list[float]]:
    """Load externally-evaluated (pre_exponential, activation_temperature) ->
    residual_norm points (e.g. from an HPC grid sweep's grid_results/*.json)
    to seed gp_minimize's x0/y0, so it starts with real data instead of
    burning local evaluations rediscovering the same terrain. Only valid as
    a prior for a stage whose objective matches the same fit pressures the
    prior points were evaluated against (residual_norm isn't comparable
    across different pressure sets)."""
    x0, y0 = [], []
    for path in sorted(prior_dir.glob("task_*.json")):
        with open(path) as fh:
            record = json.load(fh)
        x0.append([float(np.log10(record["pre_exponential"])),
                   float(record["activation_temperature"])])
        y0.append(float(record["residual_norm"]))
    return x0, y0


def run_stage(stage_label: str, geometries: list[Geometry], workdir: Path,
              lowmach_bin: str | None, pre_exponential0: float,
              activation_temperature0: float, n_calls: int,
              n_initial_points: int, cores: int,
              random_state: int = 0,
              prior_x0: list[list[float]] | None = None,
              prior_y0: list[float] | None = None) -> tuple[float, float]:
    data_by_geom = {g.name: load_experimental_data(g.data) for g in geometries}
    targets_by_geom = {g.name: nearest_experimental(data_by_geom[g.name], g.fit_pressures)
                       for g in geometries}
    for g in geometries:
        print(f"[{stage_label}] Fitting {g.name} to {len(g.fit_pressures)} pressures: "
              f"{list(zip(g.fit_pressures, targets_by_geom[g.name]))}")

    n_pressures = max(len(g.fit_pressures) for g in geometries)
    print(f"[{stage_label}] {n_pressures} sim(s) run concurrently -> "
          f"MPI_NP={max(1, cores // n_pressures)} ranks/sim ({cores} cores total)")

    log_path = workdir / f"iterations_{stage_label}.jsonl"
    objective = make_objective(geometries, targets_by_geom, workdir, lowmach_bin,
                               log_path, stage_label=stage_label, cores=cores)

    dimensions = [Real(-6.0, 1.0, name="log_pre_exponential"),
                  Real(0.0, 15000.0, name="activation_temperature")]
    seed = [float(np.log10(pre_exponential0)), float(activation_temperature0)]
    x0 = [seed]
    y0 = None
    if prior_x0:
        # Skip re-evaluating the seed if a prior point already covers it
        # (within float rendering tolerance) -- e.g. the seed is often
        # exactly a grid point an HPC sweep already ran.
        duplicate = next((i for i, p in enumerate(prior_x0)
                          if abs(p[0] - seed[0]) < 1e-6 and abs(p[1] - seed[1]) < 1e-6), None)
        print(f"[{stage_label}] seeding GP with {len(prior_x0)} prior evaluation(s) "
              f"(no local sim cost) before {n_calls} new call(s)")
        if duplicate is not None:
            x0, y0 = prior_x0, list(prior_y0)
            print(f"[{stage_label}] seed point already covered by prior "
                  f"(residual_norm={y0[duplicate]:.6g}), not re-evaluating")
        else:
            x0 = prior_x0 + x0
            y0 = prior_y0 + [objective(seed)]
        # n_calls is "new evaluations of func" per skopt's own docs --
        # points supplied via x0/y0 don't count against it, so it's left
        # as-is (the caller's requested number of genuinely new sims).
        n_initial_points = 0

    result = gp_minimize(objective, dimensions, x0=x0, y0=y0,
                         n_calls=n_calls,
                         n_initial_points=max(n_initial_points - len(x0), 0),
                         acq_func="EI", random_state=random_state)

    pre_exponential = float(10.0 ** result.x[0])
    activation_temperature = float(result.x[1])
    print(f"\n=== [{stage_label}] Best point found (GP Bayesian optimization) ===")
    print(f"htpb_pre_exponential      = {pre_exponential:.6g} 1/Pa/s")
    print(f"htpb_activation_temperature = {activation_temperature:.6g} K")
    print(f"best residual norm    = {result.fun:.6g}")

    best = {
        "htpb_pre_exponential": pre_exponential,
        "htpb_activation_temperature": activation_temperature,
    }
    best_path = workdir / f"best_fit_{stage_label}.json"
    with open(best_path, "w") as fh:
        json.dump(best, fh, indent=2)
    print(f"Wrote {best_path}\n")

    return pre_exponential, activation_temperature


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--workdir", type=Path, required=True,
                         help="directory for per-iteration sim outputs and logs")
    parser.add_argument("--pre-exponential0", type=float,
                         default=DEFAULT_PRE_EXPONENTIAL0,
                         help="initial guess for HTPB pre_exponential [1/Pa/s]")
    parser.add_argument("--activation-temperature0", type=float,
                         default=DEFAULT_ACTIVATION_TEMPERATURE0,
                         help="initial guess for HTPB activation_temperature [K]")
    parser.add_argument("--lowmach-bin",
                         help="override LOWMACH_BIN for run_htpb_sandwich_pressure_sweep.sh")
    parser.add_argument("--cores", type=int, default=DEFAULT_CORES,
                         help=f"total MPI ranks to spread across whatever sims run "
                              f"concurrently in a stage (default: {DEFAULT_CORES}, "
                              f"i.e. cpu_count()-2)")
    parser.add_argument("--n-calls", type=int, default=15,
                         help="total gp_minimize evaluations for the real (stage 2) fit "
                              "(default: 15)")
    parser.add_argument("--n-initial-points", type=int, default=5,
                         help="random/seed evaluations before the GP surrogate starts "
                              "picking points, for stage 2 (default: 5)")
    parser.add_argument("--stage1-n-calls", type=int, default=10,
                         help="total gp_minimize evaluations for the cheap endpoints-only "
                              "stage 1 fit (default: 10)")
    parser.add_argument("--stage1-n-initial-points", type=int, default=4,
                         help="random/seed evaluations before the GP surrogate starts "
                              "picking points, for stage 1 (default: 4)")
    parser.add_argument("--stage1-prior-dir", type=Path,
                         help="directory of grid_results/task_*.json (e.g. from an HPC "
                              "htpb_grid_sweep.py --mode eval run) to seed stage 1's GP "
                              "with real, already-evaluated points at no extra sim cost "
                              "-- only valid if those points were evaluated against the "
                              "same endpoint pressures as stage 1")
    parser.add_argument("--stage2-prior-dir", type=Path,
                         help="directory of grid_results/task_*.json evaluated against "
                              "every fit pressure (htpb_grid_sweep.py --all-pressures) to "
                              "seed stage 2's GP with real, already-evaluated points at no "
                              "extra sim cost")
    parser.add_argument("--skip-stage1", action="store_true",
                         help="go straight to the full (stage 2) fit using "
                              "--pre-exponential0/--activation-temperature0 as its x0, "
                              "skipping the endpoints-only warm-start stage")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.workdir.mkdir(parents=True, exist_ok=True)

    geometries = DEFAULT_GEOMETRIES

    if args.skip_stage1:
        pre_exponential0, activation_temperature0 = (
            args.pre_exponential0, args.activation_temperature0)
    else:
        print("=== Stage 1: cheap fit against lowest/highest pressure per geometry ===\n")
        prior_x0 = prior_y0 = None
        if args.stage1_prior_dir:
            prior_x0, prior_y0 = load_prior_points(args.stage1_prior_dir)
            print(f"Loaded {len(prior_x0)} prior point(s) from {args.stage1_prior_dir}")
        pre_exponential0, activation_temperature0 = run_stage(
            "stage1", endpoint_geometries(geometries), args.workdir, args.lowmach_bin,
            args.pre_exponential0, args.activation_temperature0,
            args.stage1_n_calls, args.stage1_n_initial_points, args.cores,
            prior_x0=prior_x0, prior_y0=prior_y0)
        print(f"Stage 1 result becomes stage 2's initial guess: "
              f"pre_exponential0={pre_exponential0:.6g}, "
              f"activation_temperature0={activation_temperature0:.6g}\n")

    print("=== Stage 2: real fit against every pressure, all fitted geometries ===\n")
    prior_x0 = prior_y0 = None
    if args.stage2_prior_dir:
        prior_x0, prior_y0 = load_prior_points(args.stage2_prior_dir)
        print(f"Loaded {len(prior_x0)} prior point(s) from {args.stage2_prior_dir}")
    pre_exponential, activation_temperature = run_stage(
        "stage2", geometries, args.workdir, args.lowmach_bin,
        pre_exponential0, activation_temperature0,
        args.n_calls, args.n_initial_points, args.cores,
        prior_x0=prior_x0, prior_y0=prior_y0)

    best = {
        "htpb_pre_exponential": pre_exponential,
        "htpb_activation_temperature": activation_temperature,
    }
    with open(args.workdir / "best_fit.json", "w") as fh:
        json.dump(best, fh, indent=2)
    print(f"Wrote {args.workdir / 'best_fit.json'}")

    data_by_geom = {g.name: load_experimental_data(g.data) for g in geometries}

    print("\nRunning full validation sweep over all experimental pressures, both geometries...")
    val_dir = args.workdir / "validation"
    val_dir.mkdir(parents=True, exist_ok=True)
    val_rates_by_geom = {}
    all_pressures_by_geom = {}
    for geom in geometries:
        all_pressures = sorted(data_by_geom[geom.name])
        all_pressures_by_geom[geom.name] = all_pressures
        geom_val_dir = val_dir / geom.name
        geom_val_dir.mkdir(parents=True, exist_ok=True)
        val_rates_by_geom[geom.name] = run_sweep(
            pre_exponential, activation_temperature, all_pressures,
            geom_val_dir, args.lowmach_bin, geom.template)

        calibrated_input = args.workdir / f"input.lm.ap_htpb_fullfeedback_{geom.name}"
        render_cmd = [
            sys.executable, str(SCRIPT_DIR / "render_htpb_sandwich_input.py"),
            "--template", str(geom.template),
            "--pressure-mpa", "3.0",
            "--htpb-pre-exponential", repr(pre_exponential),
            "--htpb-activation-temperature", repr(activation_temperature),
            "--out", str(calibrated_input),
        ]
        subprocess.run(render_cmd, check=True)

    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        styles = {
            "chorpening_200um": dict(exp_fmt="ko", sim_fmt="r^-",
                                      label="Chorpening et al. 2000 (~200 um)"),
            "group_100um": dict(exp_fmt="bs", sim_fmt="gD--",
                                 label="Internal group (100 um)"),
        }

        plt.figure(figsize=(8, 5))
        for geom in geometries:
            all_pressures = all_pressures_by_geom[geom.name]
            data = data_by_geom[geom.name]
            val_rates = val_rates_by_geom[geom.name]
            style = styles[geom.name]
            exp_p = np.array(all_pressures)
            exp_r = np.array([data[p] for p in all_pressures])
            sim_p = np.array([p for p in all_pressures if val_rates.get(p) is not None])
            sim_r = np.array([val_rates[p] for p in sim_p])
            plt.plot(exp_p, exp_r, style["exp_fmt"], label=f"{style['label']} - Experiment")
            plt.plot(sim_p, sim_r, style["sim_fmt"], label=f"{style['label']} - Sim")

        plt.xlabel("Pressure (MPa)")
        plt.ylabel("Sandwich Regression Rate (mm/s)")
        plt.xscale("log")
        plt.yscale("log")
        plt.title(f"HTPB fullfeedback combined fit: pre_exponential={pre_exponential:.4g}, "
                  f"activation_temperature={activation_temperature:.4g} K")
        plt.grid(True, alpha=0.4, which="both")
        plt.legend(fontsize=8)
        plt.tight_layout()
        plot_path = args.workdir / "htpb_fullfeedback_combined_fit.png"
        plt.savefig(plot_path, dpi=300)
        print(f"Wrote {plot_path}")
    except ImportError:
        print("matplotlib not available; skipping htpb_fullfeedback_combined_fit.png")

    for geom in geometries:
        print(f"\nValidation rates for {geom.name} (mm/s):")
        for p in all_pressures_by_geom[geom.name]:
            sim = val_rates_by_geom[geom.name].get(p)
            exp = data_by_geom[geom.name][p]
            sim_str = f"{sim:.4f}" if sim is not None else "FAILED"
            print(f"  P={p:>6.3f} MPa: sim={sim_str}  exp={exp:.4f}")


if __name__ == "__main__":
    main()
