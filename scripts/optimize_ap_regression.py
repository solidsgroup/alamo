#!/usr/bin/env python3
"""Fit (rate_multiplier, activation_temperature) to AP_reg_rate.csv.

Drives scripts/run_pressure_sweep.sh with scipy.optimize.least_squares to find
a single AP_decomposition.phase_change (rate_multiplier, activation_temperature)
pair whose simulated AP monopropellant regression rate best matches the
experimental r(P) curve in AP_reg_rate.csv, across the whole fit-pressure
range at once (not a per-pressure fit). No C++ change: pressure sensitivity
comes entirely from the Rocfire gas-phase feedback already in
input.lm.ap_monopropellant.template; allencahn.mobility stays fixed because
only rate_multiplier * mobility is identifiable (see MassSource in
src/Model/Mechanism/PhaseChange.H).

Each objective evaluation runs one pressure sweep (all pressures in parallel,
via run_pressure_sweep.sh) and compares the resulting steady-state regression
rates (mm/s) to the experimental values at the same pressures, using relative
residuals. Experimental targets that are <= 0 (the 1.0 MPa deflagration-limit
point in AP_reg_rate.csv) are dropped from the fit automatically -- a relative
residual against a zero target is undefined, and 1 MPa is checked separately
in the validation sweep as a pass/fail extinction check instead.

This is a *long-running, expensive* driver -- each sweep is several minutes.
Run it as a background job; it logs every iteration's parameters and
residual norm as it goes, and writes the best-fit input + a validation plot
when it converges.

Example
-------
Optimize using the default 5-pressure subset (2-6 MPa), then validate against
every pressure in AP_reg_rate.csv (including the 1 MPa extinction check) and
write ap_regression_fit.png::

    python scripts/optimize_ap_regression.py --workdir /tmp/ap_calib
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

# All non-zero pressures in AP_reg_rate.csv (2-6 MPa); 1 MPa is the
# deflagration-limit point (rate = 0) and is excluded from the fit.
DEFAULT_FIT_PRESSURES = [2.0, 3.0, 4.0, 5.0, 6.0]

# Initial guess: the coefficients currently in input.lm.ap_monopropellant.
DEFAULT_RATE_MULTIPLIER0 = 2450.0
DEFAULT_ACTIVATION_TEMPERATURE0 = 3145.0

RATE_UNIT = "mm/s"
RATE_UNIT_SUFFIX = "mm_s"

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


def run_sweep(rate_multiplier: float, activation_temperature: float,
              pressures: list[float], workdir: Path,
              lowmach_bin: str | None, template: Path | None,
              min_time_low: float = 0.0, min_time_high: float = 0.0
              ) -> dict[float, float | None]:
    results_csv = workdir / "results.csv"
    import os
    env = os.environ.copy()
    if lowmach_bin:
        env["LOWMACH_BIN"] = lowmach_bin
    if template:
        env["TEMPLATE"] = str(template)
    env["RATE_UNIT"] = RATE_UNIT
    if min_time_low:
        env["MIN_TIME_LOW"] = repr(min_time_low)
    if min_time_high:
        env["MIN_TIME_HIGH"] = repr(min_time_high)
    cmd = [
        str(SCRIPT_DIR / "run_pressure_sweep.sh"),
        repr(rate_multiplier), repr(activation_temperature),
        str(workdir), str(results_csv),
        *[repr(p) for p in pressures],
    ]
    subprocess.run(cmd, check=True, env=env)
    rates: dict[float, float | None] = {}
    with open(results_csv) as fh:
        reader = csv.DictReader(fh)
        for row in reader:
            p = float(row["pressure_mpa"])
            r = row[f"reg_rate_{RATE_UNIT_SUFFIX}"]
            rates[p] = float(r) if r else None
    return rates


def make_fun_jac(pressures: list[float], targets: list[float], workdir: Path,
                  lowmach_bin: str | None, template: Path | None, log_path: Path,
                  min_time_low: float = 0.0, min_time_high: float = 0.0,
                  stage: str = "full", diff_step: float = 0.05):
    """Build (fun, jac) for least_squares that run independent sweeps concurrently.

    Each gradient step needs the residual at x plus one forward-difference
    perturbation per parameter -- 3 independent sweeps here (2 params). scipy
    calls these one at a time by default; running all of them at once (each
    on its own thread, since run_sweep just blocks on a subprocess) uses the
    machine's idle cores instead of leaving them idle between sequential
    sweeps. Evaluations are cached by parameter vector so repeated x's (e.g.
    fun(x) then jac(x) at the same point) don't re-run a sweep.

    Experimental targets that are <= 0 (the 1 MPa deflagration-limit point,
    if present in `pressures`) are dropped before computing residuals.
    """
    import threading
    from concurrent.futures import ThreadPoolExecutor

    fit_pairs = [(p, t) for p, t in zip(pressures, targets) if t > 0.0]
    dropped = [p for p, t in zip(pressures, targets) if t <= 0.0]
    if dropped:
        print(f"[{stage}] dropping non-positive experimental targets at "
              f"pressures {dropped} from the fit residuals", flush=True)

    iteration = [0]
    cache: dict[tuple, np.ndarray] = {}
    lock = threading.Lock()

    def eval_at(x: np.ndarray) -> np.ndarray:
        key = tuple(round(float(v), 12) for v in x)
        with lock:
            if key in cache:
                return cache[key]
            iteration[0] += 1
            idx = iteration[0]

        rate_multiplier = float(10.0 ** x[0])
        activation_temperature = float(x[1])
        eval_dir = workdir / f"{stage}_iter_{idx:03d}"
        eval_dir.mkdir(parents=True, exist_ok=True)
        rates = run_sweep(rate_multiplier, activation_temperature, pressures,
                           eval_dir, lowmach_bin, template,
                           min_time_low, min_time_high)

        residuals = []
        for p, target in fit_pairs:
            sim = rates.get(p)
            if sim is None:
                residuals.append(PENALTY_RESIDUAL)
            else:
                residuals.append((sim - target) / target)
        residuals = np.array(residuals)

        record = {
            "stage": stage,
            "iteration": idx,
            "rate_multiplier": rate_multiplier,
            "activation_temperature": activation_temperature,
            "rates_mm_s": rates,
            "residual_norm": float(np.linalg.norm(residuals)),
        }
        with lock:
            with open(log_path, "a") as fh:
                fh.write(json.dumps(record) + "\n")
            print(f"[{stage} iter {idx}] rate_multiplier={rate_multiplier:.6g} "
                  f"activation_temperature={activation_temperature:.6g} "
                  f"residual_norm={record['residual_norm']:.6g}", flush=True)
            cache[key] = residuals
        return residuals

    def fun(x: np.ndarray) -> np.ndarray:
        return eval_at(x)

    def jac(x: np.ndarray, f0: np.ndarray | None = None) -> np.ndarray:
        steps = [diff_step * max(abs(xi), 1.0) for xi in x]
        perturbed = []
        for i, step in enumerate(steps):
            xp = x.copy()
            xp[i] = xp[i] + step
            perturbed.append(xp)

        with ThreadPoolExecutor(max_workers=len(x) + 1) as ex:
            base_future = ex.submit(eval_at, x)
            pert_futures = [ex.submit(eval_at, xp) for xp in perturbed]
            base = base_future.result()
            columns = [(fut.result() - base) / step
                       for fut, step in zip(pert_futures, steps)]
        return np.column_stack(columns)

    return fun, jac


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
    parser.add_argument("--rate-multiplier0", type=float,
                         default=DEFAULT_RATE_MULTIPLIER0,
                         help="initial guess for rate_multiplier (dimensionless)")
    parser.add_argument("--activation-temperature0", type=float,
                         default=DEFAULT_ACTIVATION_TEMPERATURE0,
                         help="initial guess for activation_temperature [K]")
    parser.add_argument("--lowmach-bin",
                         help="override LOWMACH_BIN for run_pressure_sweep.sh")
    parser.add_argument("--template", type=Path,
                         help="override TEMPLATE for run_pressure_sweep.sh")
    parser.add_argument("--min-time-low", type=float, default=0.0,
                         help="seconds of simulated time to exclude from the start "
                              "of the run at the lowest pressure in each sweep, to "
                              "skip its startup transient (default: 0.0, no "
                              "exclusion). Pressures between the lowest and highest "
                              "in a given sweep get a --min-time linearly "
                              "interpolated between --min-time-low and "
                              "--min-time-high (forwarded via MIN_TIME_LOW/"
                              "MIN_TIME_HIGH)")
    parser.add_argument("--min-time-high", type=float, default=0.0,
                         help="same as --min-time-low but for the highest pressure "
                              "in each sweep (default: 0.0)")
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

    x0 = np.array([np.log10(args.rate_multiplier0), args.activation_temperature0])
    bounds = ([-1.0, 0.0], [8.0, 10000.0])

    if args.bracket_stage and len(args.fit_pressures) > 2:
        # Optimize against just the lowest and highest fit pressures first --
        # two sims per iteration instead of the full set -- to get close to
        # the right (rate_multiplier, activation_temperature) region cheaply
        # before paying for every pressure.
        bracket_pressures = [min(args.fit_pressures), max(args.fit_pressures)]
        bracket_targets = nearest_experimental(data, bracket_pressures)
        print(f"\n=== Stage 1: bracket fit to {bracket_pressures} ===")
        bracket_fun, bracket_jac = make_fun_jac(
            bracket_pressures, bracket_targets, args.workdir,
            args.lowmach_bin, args.template, log_path,
            args.min_time_low, args.min_time_high, stage="bracket")
        bracket_result = least_squares(
            bracket_fun, x0, jac=bracket_jac, bounds=bounds, xtol=args.xtol,
            max_nfev=args.bracket_max_nfev)
        print(f"Stage 1 result: rate_multiplier={10.0 ** bracket_result.x[0]:.6g} "
              f"activation_temperature={bracket_result.x[1]:.6g} "
              f"residual_norm={np.linalg.norm(bracket_result.fun):.6g}")
        x0 = bracket_result.x

    print(f"\n=== Stage 2: full fit to {args.fit_pressures} ===")
    fun, jac = make_fun_jac(args.fit_pressures, targets, args.workdir,
                             args.lowmach_bin, args.template, log_path,
                             args.min_time_low, args.min_time_high, stage="full")

    result = least_squares(fun, x0, jac=jac, bounds=bounds, xtol=args.xtol,
                           max_nfev=args.max_nfev)

    rate_multiplier = float(10.0 ** result.x[0])
    activation_temperature = float(result.x[1])
    print("\n=== Converged (or hit max_nfev) ===")
    print(f"rate_multiplier         = {rate_multiplier:.6g}")
    print(f"activation_temperature  = {activation_temperature:.6g} K")
    print(f"final residual norm     = {np.linalg.norm(result.fun):.6g}")
    print(f"success={result.success} status={result.status}: {result.message}")

    best = {
        "rate_multiplier": rate_multiplier,
        "activation_temperature": activation_temperature,
    }
    with open(args.workdir / "best_fit.json", "w") as fh:
        json.dump(best, fh, indent=2)
    print(f"Wrote {args.workdir / 'best_fit.json'}")

    # Final validation across every experimental pressure, including the
    # 1 MPa deflagration-limit point (checked qualitatively, not fitted).
    print("\nRunning full validation sweep over all experimental pressures...")
    all_pressures = sorted(data)
    val_dir = args.workdir / "validation"
    val_dir.mkdir(parents=True, exist_ok=True)
    val_rates = run_sweep(rate_multiplier, activation_temperature, all_pressures,
                          val_dir, args.lowmach_bin, args.template,
                          args.min_time_low, args.min_time_high)

    calibrated_input = args.workdir / "input.lm.ap_monopropellant"
    render_cmd = [
        sys.executable, str(SCRIPT_DIR / "render_input.py"),
        "--template", str(args.template or REPO_ROOT / "input.lm.ap_monopropellant.template"),
        "--pressure-mpa", "3.0",
        "--rate-multiplier", repr(rate_multiplier),
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
        plt.plot(sim_p, sim_r, "r^-", label="Simulated fit")
        plt.xlabel("Pressure (MPa)")
        plt.ylabel("Regression Rate (mm/s)")
        plt.title(f"AP regression fit: rate_multiplier={rate_multiplier:.4g}, "
                  f"activation_temperature={activation_temperature:.4g} K")
        plt.grid(True, alpha=0.4)
        plt.legend()
        plt.tight_layout()
        plot_path = args.workdir / "ap_regression_fit.png"
        plt.savefig(plot_path, dpi=300)
        print(f"Wrote {plot_path}")
    except ImportError:
        print("matplotlib not available; skipping ap_regression_fit.png")

    print("\nValidation rates (mm/s):")
    for p in all_pressures:
        sim = val_rates.get(p)
        exp = data[p]
        sim_str = f"{sim:.4f}" if sim is not None else "EXTINGUISHED/FAILED"
        print(f"  P={p:>6.2f} MPa: sim={sim_str}  exp={exp:.4f}")
        if exp <= 0.0:
            note = "matches deflagration limit (extinguished)" if sim is None \
                else "DID NOT extinguish as expected"
            print(f"    -> deflagration-limit check: {note}")


if __name__ == "__main__":
    main()
