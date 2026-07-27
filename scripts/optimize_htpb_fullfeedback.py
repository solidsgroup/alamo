#!/usr/bin/env python3
"""Fit HTPB fullfeedback (pre_exponential, activation_temperature) to sandwich data.

Drives scripts/run_htpb_sandwich_pressure_sweep.sh with
scipy.optimize.least_squares to find a single (pre_exponential,
activation_temperature) pair for HTPB_pyrolysis whose simulated AP/HTPB
sandwich regression rate best matches the experimental r(P) curve in
HTPB_sandwich_reg_rate_chorpening2000_summary.csv, across the whole pressure
range at once (not a per-pressure fit).

AP_decomposition's fullfeedback parameters are NOT touched by this script --
they are baked into input.lm.ap_htpb_fullfeedback.template at the values
already calibrated against pure AP monopropellant data (see
FULLFEEDBACK_CALIBRATION.md). Only HTPB's own fullfeedback parameters are
fit here. pressure_dependence stays linear (P/reference_pressure) for HTPB,
mirroring AP; no C++ change is involved.

Each objective evaluation runs one pressure sweep (all pressures in
parallel, via run_htpb_sandwich_pressure_sweep.sh) and compares the
resulting steady-state regression rates (mm/s) to the experimental values at
the same pressures, using relative residuals.

By default, optimization runs in two stages: first fitting against only the
lowest and highest --fit-pressures (cheapest possible sweep, since it still
brackets the full range) to get (pre_exponential, activation_temperature)
into the right region, then fitting against the full --fit-pressures set
starting from that result. Pass --no-bracket-stage to skip straight to the
full-pressure fit.

CAVEAT: the experimental data only spans 0.2-3.2 MPa (Chorpening, Knott &
Brewster 2000), well below the 2.76-13.79 MPa range AP_decomposition's own
parameters were calibrated against. Evaluating AP's fixed fullfeedback law
at these lower pressures is an extrapolation below its own calibration
floor; treat the resulting HTPB fit as approximate for this reason, not
just because of the usual sandwich-vs-monopropellant modeling gap.

This is a *long-running, expensive* driver -- each sweep is a 2D run and
takes longer than the monopropellant sweeps. Run it as a background job; it
logs every iteration's parameters and residual norm as it goes, and writes
the best-fit input + a validation plot when it converges (or is interrupted).

Examples
--------
Optimize using all 4 summary pressures, then write the calibrated input and
a validation plot::

    python scripts/optimize_htpb_fullfeedback.py --workdir /tmp/htpb_calib

Resume/extend an interrupted run using its saved iteration log (rerun the
same command; results.csv/iterations.jsonl accumulate under --workdir)::

    python scripts/optimize_htpb_fullfeedback.py --workdir /tmp/htpb_calib
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent

# All 4 points in the summary CSV are used during optimization by default --
# unlike the AP fit, there are only 4 experimental points here (see the
# digitization notes in HTPB_sandwich_reg_rate_chorpening2000.csv), so there
# is no cheaper representative subset to fall back on.
DEFAULT_FIT_PRESSURES = [0.213, 0.467, 1.498, 3.098]

# Initial guess: HTPB_pyrolysis's old allencahn mobility was a constant
# 0.01 1/Pa/s with a separate mechanism-level activation_temperature of
# 7500 K; used as a starting point for the fullfeedback pre_exponential search
# (log10-space) and activation_temperature.
DEFAULT_PRE_EXPONENTIAL0 = 0.001
DEFAULT_ACTIVATION_TEMPERATURE0 = 4000.0

PENALTY_RESIDUAL = 5.0  # relative-error stand-in for a failed/non-igniting sim


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
    """Look up (or interpolate) experimental rates at the requested pressures."""
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
              lowmach_bin: str | None, template: Path | None,
              min_time: float = 0.0) -> dict[float, float | None]:
    results_csv = workdir / "results.csv"
    import os
    env = os.environ.copy()
    if lowmach_bin:
        env["LOWMACH_BIN"] = lowmach_bin
    if template:
        env["TEMPLATE"] = str(template)
    if min_time:
        env["MIN_TIME"] = repr(min_time)
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


def make_objective(pressures: list[float], targets: list[float], workdir: Path,
                    lowmach_bin: str | None, template: Path | None, log_path: Path,
                    min_time: float = 0.0, stage: str = "full"):
    iteration = [0]

    def objective(x: np.ndarray) -> np.ndarray:
        iteration[0] += 1
        pre_exponential = float(10.0 ** x[0])
        activation_temperature = float(x[1])
        eval_dir = workdir / f"{stage}_iter_{iteration[0]:03d}"
        eval_dir.mkdir(parents=True, exist_ok=True)
        rates = run_sweep(pre_exponential, activation_temperature, pressures,
                           eval_dir, lowmach_bin, template, min_time)

        residuals = []
        for p, target in zip(pressures, targets):
            sim = rates.get(p)
            if sim is None:
                residuals.append(PENALTY_RESIDUAL)
            else:
                residuals.append((sim - target) / target)
        residuals = np.array(residuals)

        record = {
            "stage": stage,
            "iteration": iteration[0],
            "pre_exponential": pre_exponential,
            "activation_temperature": activation_temperature,
            "rates_mm_s": rates,
            "residual_norm": float(np.linalg.norm(residuals)),
        }
        with open(log_path, "a") as fh:
            fh.write(json.dumps(record) + "\n")
        print(f"[{stage} iter {iteration[0]}] pre_exponential={pre_exponential:.6g} "
              f"activation_temperature={activation_temperature:.6g} "
              f"residual_norm={record['residual_norm']:.6g}", flush=True)
        return residuals

    return objective


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", type=Path,
                         default=REPO_ROOT / "HTPB_sandwich_reg_rate_chorpening2000_summary.csv",
                         help="experimental data CSV (default: "
                              "HTPB_sandwich_reg_rate_chorpening2000_summary.csv)")
    parser.add_argument("--workdir", type=Path, required=True,
                         help="directory for per-iteration sim outputs and logs")
    parser.add_argument("--fit-pressures", type=float, nargs="+",
                         default=DEFAULT_FIT_PRESSURES,
                         help=f"pressures (MPa) used during optimization "
                              f"(default: {DEFAULT_FIT_PRESSURES})")
    parser.add_argument("--pre-exponential0", type=float,
                         default=DEFAULT_PRE_EXPONENTIAL0,
                         help="initial guess for HTPB pre_exponential [1/Pa/s]")
    parser.add_argument("--activation-temperature0", type=float,
                         default=DEFAULT_ACTIVATION_TEMPERATURE0,
                         help="initial guess for HTPB activation_temperature [K]")
    parser.add_argument("--lowmach-bin",
                         help="override LOWMACH_BIN for run_htpb_sandwich_pressure_sweep.sh")
    parser.add_argument("--template", type=Path,
                         help="override TEMPLATE for run_htpb_sandwich_pressure_sweep.sh")
    parser.add_argument("--min-time", type=float, default=0.0,
                         help="seconds of simulated time to exclude from the start "
                              "of each run before measuring the regression rate, to "
                              "skip the startup transient (default: 0.0, no "
                              "exclusion; forwarded to regression_rate.py --min-time "
                              "via MIN_TIME)")
    parser.add_argument("--xtol", type=float, default=1.0e-3,
                         help="least_squares xtol (default: 1e-3)")
    parser.add_argument("--max-nfev", type=int, default=30,
                         help="cap on objective evaluations for the full-pressure-set "
                              "stage (default: 30)")
    parser.add_argument("--no-bracket-stage", dest="bracket_stage",
                         action="store_false",
                         help="skip the initial lowest+highest-pressure-only stage "
                              "and optimize all --fit-pressures directly")
    parser.add_argument("--bracket-max-nfev", type=int, default=15,
                         help="cap on objective evaluations for the initial "
                              "lowest+highest-pressure bracket stage (default: 15)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.workdir.mkdir(parents=True, exist_ok=True)
    log_path = args.workdir / "iterations.jsonl"

    data = load_experimental_data(args.data)
    targets = nearest_experimental(data, args.fit_pressures)
    print(f"Fitting to {len(args.fit_pressures)} pressures: "
          f"{list(zip(args.fit_pressures, targets))}")

    x0 = np.array([np.log10(args.pre_exponential0), args.activation_temperature0])
    bounds = ([-6.0, 0.0], [1.0, 15000.0])

    if args.bracket_stage and len(args.fit_pressures) > 2:
        # Optimize against just the lowest and highest fit pressures first --
        # two sims per iteration instead of the full set -- to get close to
        # the right (pre_exponential, activation_temperature) region cheaply
        # before paying for every pressure.
        bracket_pressures = [min(args.fit_pressures), max(args.fit_pressures)]
        bracket_targets = nearest_experimental(data, bracket_pressures)
        print(f"\n=== Stage 1: bracket fit to {bracket_pressures} ===")
        bracket_objective = make_objective(
            bracket_pressures, bracket_targets, args.workdir,
            args.lowmach_bin, args.template, log_path, args.min_time,
            stage="bracket")
        bracket_result = least_squares(
            bracket_objective, x0, bounds=bounds, xtol=args.xtol,
            max_nfev=args.bracket_max_nfev, diff_step=0.05)
        print(f"Stage 1 result: pre_exponential={10.0 ** bracket_result.x[0]:.6g} "
              f"activation_temperature={bracket_result.x[1]:.6g} "
              f"residual_norm={np.linalg.norm(bracket_result.fun):.6g}")
        x0 = bracket_result.x

    print(f"\n=== Stage 2: full fit to {args.fit_pressures} ===")
    objective = make_objective(args.fit_pressures, targets, args.workdir,
                               args.lowmach_bin, args.template, log_path,
                               args.min_time, stage="full")

    result = least_squares(objective, x0, bounds=bounds, xtol=args.xtol,
                           max_nfev=args.max_nfev, diff_step=0.05)

    pre_exponential = float(10.0 ** result.x[0])
    activation_temperature = float(result.x[1])
    print("\n=== Converged (or hit max_nfev) ===")
    print(f"htpb_pre_exponential      = {pre_exponential:.6g} 1/Pa/s")
    print(f"htpb_activation_temperature = {activation_temperature:.6g} K")
    print(f"final residual norm   = {np.linalg.norm(result.fun):.6g}")
    print(f"success={result.success} status={result.status}: {result.message}")

    best = {
        "htpb_pre_exponential": pre_exponential,
        "htpb_activation_temperature": activation_temperature,
    }
    with open(args.workdir / "best_fit.json", "w") as fh:
        json.dump(best, fh, indent=2)
    print(f"Wrote {args.workdir / 'best_fit.json'}")

    # Final validation across every experimental pressure (same 4 points --
    # there is no larger held-out set for this dataset).
    print("\nRunning full validation sweep over all experimental pressures...")
    all_pressures = sorted(data)
    val_dir = args.workdir / "validation"
    val_dir.mkdir(parents=True, exist_ok=True)
    val_rates = run_sweep(pre_exponential, activation_temperature, all_pressures,
                          val_dir, args.lowmach_bin, args.template, args.min_time)

    calibrated_input = args.workdir / "input.lm.ap_htpb_fullfeedback"
    render_cmd = [
        sys.executable, str(SCRIPT_DIR / "render_htpb_sandwich_input.py"),
        "--template", str(args.template or REPO_ROOT / "input.lm.ap_htpb_fullfeedback.template"),
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

        exp_p = np.array(all_pressures)
        exp_r = np.array([data[p] for p in all_pressures])
        sim_p = np.array([p for p in all_pressures if val_rates.get(p) is not None])
        sim_r = np.array([val_rates[p] for p in sim_p])

        plt.figure(figsize=(8, 5))
        plt.plot(exp_p, exp_r, "ko", label="Experiment (Chorpening et al. 2000)")
        plt.plot(sim_p, sim_r, "r^-", label="HTPB fullfeedback fit")
        plt.xlabel("Pressure (MPa)")
        plt.ylabel("Sandwich Regression Rate (mm/s)")
        plt.xscale("log")
        plt.yscale("log")
        plt.title(f"HTPB fullfeedback fit: pre_exponential={pre_exponential:.4g}, "
                  f"activation_temperature={activation_temperature:.4g} K")
        plt.grid(True, alpha=0.4, which="both")
        plt.legend()
        plt.tight_layout()
        plot_path = args.workdir / "htpb_fullfeedback_fit.png"
        plt.savefig(plot_path, dpi=300)
        print(f"Wrote {plot_path}")
    except ImportError:
        print("matplotlib not available; skipping htpb_fullfeedback_fit.png")

    print("\nValidation rates (mm/s):")
    for p in all_pressures:
        sim = val_rates.get(p)
        exp = data[p]
        sim_str = f"{sim:.4f}" if sim is not None else "FAILED"
        print(f"  P={p:>6.3f} MPa: sim={sim_str}  exp={exp:.4f}")


if __name__ == "__main__":
    main()
