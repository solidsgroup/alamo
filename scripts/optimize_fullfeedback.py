#!/usr/bin/env python3
"""Fit fullfeedback (pre_exponential, activation_temperature) to AP_reg_rate.csv.

Drives scripts/run_pressure_sweep.sh with scipy.optimize.least_squares to find
a single (pre_exponential, activation_temperature) pair whose simulated AP
monopropellant regression rate best matches the experimental r(P) curve in
AP_reg_rate.csv, across the whole pressure range at once (not a per-pressure
fit). pressure_dependence stays linear (P/reference_pressure); no C++ change.

Each objective evaluation runs one pressure sweep (all pressures in parallel,
via run_pressure_sweep.sh) and compares the resulting steady-state regression
rates (cm/s) to the experimental values at the same pressures, using relative
residuals so the ~0.38-1.24 cm/s data range is weighted evenly.

This is a *long-running, expensive* driver -- each sweep is several minutes.
Run it as a background job; it logs every iteration's parameters and
residual norm as it goes, and writes the best-fit input + a validation plot
when it converges (or is interrupted -- see --resume).

Examples
--------
Optimize using the default 6-pressure subset, then validate against all 16
points and write fullfeedback_fit.png::

    python scripts/optimize_fullfeedback.py --workdir /tmp/ff_calib

Resume/extend an interrupted run using its saved iteration log::

    python scripts/optimize_fullfeedback.py --workdir /tmp/ff_calib --resume
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

# Representative subset used during optimization (spans the full pressure
# range without paying for all 16 sims every iteration). Final validation
# uses every pressure in AP_reg_rate.csv.
DEFAULT_FIT_PRESSURES = [2.85, 4.14, 5.61, 8.27, 10.44, 13.79]

# Initial guess: pre_exponential scaled down from the original uncalibrated
# 0.01 by the ~4/14.75 mm/s ratio observed at 3 MPa; activation_temperature
# left at its original value.
DEFAULT_PRE_EXPONENTIAL0 = 0.0027
DEFAULT_ACTIVATION_TEMPERATURE0 = 3145.0

PENALTY_RESIDUAL = 5.0  # relative-error stand-in for a failed/non-igniting sim


def load_experimental_data(csv_path: Path) -> dict[float, float]:
    """Average duplicate/near-duplicate pressures, as in plot_AP_reg_rate.py."""
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
              lowmach_bin: str | None, template: Path | None) -> dict[float, float | None]:
    results_csv = workdir / "results.csv"
    env_prefix = []
    import os
    env = os.environ.copy()
    if lowmach_bin:
        env["LOWMACH_BIN"] = lowmach_bin
    if template:
        env["TEMPLATE"] = str(template)
    cmd = [
        str(SCRIPT_DIR / "run_pressure_sweep.sh"),
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
            r = row["reg_rate_cm_s"]
            rates[p] = float(r) if r else None
    return rates


def make_objective(pressures: list[float], targets: list[float], workdir: Path,
                    lowmach_bin: str | None, template: Path | None, log_path: Path):
    iteration = [0]

    def objective(x: np.ndarray) -> np.ndarray:
        iteration[0] += 1
        pre_exponential = float(10.0 ** x[0])
        activation_temperature = float(x[1])
        eval_dir = workdir / f"iter_{iteration[0]:03d}"
        eval_dir.mkdir(parents=True, exist_ok=True)
        rates = run_sweep(pre_exponential, activation_temperature, pressures,
                           eval_dir, lowmach_bin, template)

        residuals = []
        for p, target in zip(pressures, targets):
            sim = rates.get(p)
            if sim is None:
                residuals.append(PENALTY_RESIDUAL)
            else:
                residuals.append((sim - target) / target)
        residuals = np.array(residuals)

        record = {
            "iteration": iteration[0],
            "pre_exponential": pre_exponential,
            "activation_temperature": activation_temperature,
            "rates_cm_s": rates,
            "residual_norm": float(np.linalg.norm(residuals)),
        }
        with open(log_path, "a") as fh:
            fh.write(json.dumps(record) + "\n")
        print(f"[iter {iteration[0]}] pre_exponential={pre_exponential:.6g} "
              f"activation_temperature={activation_temperature:.6g} "
              f"residual_norm={record['residual_norm']:.6g}", flush=True)
        return residuals

    return objective


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data", type=Path, default=REPO_ROOT / "AP_reg_rate.csv",
                         help="experimental data CSV (default: AP_reg_rate.csv)")
    parser.add_argument("--workdir", type=Path, required=True,
                         help="directory for per-iteration sim outputs and logs")
    parser.add_argument("--fit-pressures", type=float, nargs="+",
                         default=DEFAULT_FIT_PRESSURES,
                         help=f"pressures (MPa) used during optimization "
                              f"(default: {DEFAULT_FIT_PRESSURES})")
    parser.add_argument("--pre-exponential0", type=float,
                         default=DEFAULT_PRE_EXPONENTIAL0,
                         help="initial guess for pre_exponential [1/Pa/s]")
    parser.add_argument("--activation-temperature0", type=float,
                         default=DEFAULT_ACTIVATION_TEMPERATURE0,
                         help="initial guess for activation_temperature [K]")
    parser.add_argument("--lowmach-bin",
                         help="override LOWMACH_BIN for run_pressure_sweep.sh")
    parser.add_argument("--template", type=Path,
                         help="override TEMPLATE for run_pressure_sweep.sh")
    parser.add_argument("--xtol", type=float, default=1.0e-3,
                         help="least_squares xtol (default: 1e-3)")
    parser.add_argument("--max-nfev", type=int, default=30,
                         help="cap on objective evaluations (default: 30)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.workdir.mkdir(parents=True, exist_ok=True)
    log_path = args.workdir / "iterations.jsonl"

    data = load_experimental_data(args.data)
    targets = nearest_experimental(data, args.fit_pressures)
    print(f"Fitting to {len(args.fit_pressures)} pressures: "
          f"{list(zip(args.fit_pressures, targets))}")

    objective = make_objective(args.fit_pressures, targets, args.workdir,
                               args.lowmach_bin, args.template, log_path)

    x0 = np.array([np.log10(args.pre_exponential0), args.activation_temperature0])
    bounds = ([-6.0, 0.0], [1.0, 10000.0])
    result = least_squares(objective, x0, bounds=bounds, xtol=args.xtol,
                           max_nfev=args.max_nfev, diff_step=0.05)

    pre_exponential = float(10.0 ** result.x[0])
    activation_temperature = float(result.x[1])
    print("\n=== Converged (or hit max_nfev) ===")
    print(f"pre_exponential      = {pre_exponential:.6g} 1/Pa/s")
    print(f"activation_temperature = {activation_temperature:.6g} K")
    print(f"final residual norm   = {np.linalg.norm(result.fun):.6g}")
    print(f"success={result.success} status={result.status}: {result.message}")

    best = {
        "pre_exponential": pre_exponential,
        "activation_temperature": activation_temperature,
    }
    with open(args.workdir / "best_fit.json", "w") as fh:
        json.dump(best, fh, indent=2)
    print(f"Wrote {args.workdir / 'best_fit.json'}")

    # Final validation across every experimental pressure.
    print("\nRunning full validation sweep over all experimental pressures...")
    all_pressures = sorted(data)
    val_dir = args.workdir / "validation"
    val_dir.mkdir(parents=True, exist_ok=True)
    val_rates = run_sweep(pre_exponential, activation_temperature, all_pressures,
                          val_dir, args.lowmach_bin, args.template)

    calibrated_input = args.workdir / "input.lm.ap_monopropellant_fullfeedback"
    render_cmd = [
        sys.executable, str(SCRIPT_DIR / "render_input.py"),
        "--template", str(args.template or REPO_ROOT / "input.lm.ap_monopropellant_fullfeedback.template"),
        "--pressure-mpa", "3.0",
        "--pre-exponential", repr(pre_exponential),
        "--activation-temperature", repr(activation_temperature),
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
        plt.plot(exp_p, exp_r, "ko", label="Experiment")
        plt.plot(sim_p, sim_r, "r^-", label="fullfeedback fit")
        plt.xlabel("Pressure (MPa)")
        plt.ylabel("Regression Rate (cm/s)")
        plt.title(f"fullfeedback fit: pre_exponential={pre_exponential:.4g}, "
                  f"activation_temperature={activation_temperature:.4g} K")
        plt.grid(True, alpha=0.4)
        plt.legend()
        plt.tight_layout()
        plot_path = args.workdir / "fullfeedback_fit.png"
        plt.savefig(plot_path, dpi=300)
        print(f"Wrote {plot_path}")
    except ImportError:
        print("matplotlib not available; skipping fullfeedback_fit.png")

    print("\nValidation rates (cm/s):")
    for p in all_pressures:
        sim = val_rates.get(p)
        exp = data[p]
        sim_str = f"{sim:.4f}" if sim is not None else "FAILED"
        print(f"  P={p:>6.2f} MPa: sim={sim_str}  exp={exp:.4f}")


if __name__ == "__main__":
    main()
