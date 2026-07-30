#!/usr/bin/env python3
"""Fit (rate_multiplier, activation_temperature) to AP_reg_rate.csv.

Drives tests/LMRFMonoAP/input directly (the same input used by the
LMRFMonoAP regression test) with command-line ParmParse overrides for
pressure, rate_multiplier and activation_temperature, and measures the
resulting regression rate with the same method as tests/LMRFMonoAP/test:
track the rigid_eta = 0.5 interface position vs. time via yt, then fit a
line to the back half of the run (back quarter for the 1 MPa case) to get
a steady-state rate in mm/s. No template files, no custom HDF5 rate
extractor -- this is exactly what the regression test itself checks,
just without the pass/fail assertions and with variable coefficients.

scipy.optimize.least_squares searches x = [log10(rate_multiplier),
activation_temperature] against the experimental r(P) curve in
AP_reg_rate.csv (matches tests/LMRFMonoAP/reference.csv), across
2-6 MPa at once. allencahn.mobility stays fixed at the test's value
(0.01_1/Pa/s) because only rate_multiplier * mobility is identifiable
(see MassSource in src/Model/Mechanism/PhaseChange.H). 1 MPa (the
deflagration-limit point, rate = 0) is excluded from the fit residuals
and checked separately in the validation sweep as a pass/fail extinction
check.

Run as a background job; it logs every iteration's parameters and
residual norm to <workdir>/iterations.jsonl as it goes, and writes
<workdir>/best_fit.json plus a validation plot when it converges.

Example
-------
    python scripts/optimize_ap_regression.py --workdir /tmp/ap_calib
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import subprocess
import sys
import threading
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
TEST_INPUT = REPO_ROOT / "tests" / "LMRFMonoAP" / "input"
DEFAULT_LOWMACH_BIN = REPO_ROOT / "bin" / "lowmach-2d-clang++"

sys.path.insert(0, str(REPO_ROOT / "scripts"))
import testlib  # noqa: E402

# All non-zero pressures in AP_reg_rate.csv (2-6 MPa); 1 MPa is the
# deflagration-limit point (rate = 0) and is excluded from the fit.
DEFAULT_FIT_PRESSURES = [2.0, 3.0, 4.0, 5.0, 6.0]

# tests/LMRFMonoAP/input's own coefficients -- a reasonable starting point
# since they already land within ~15-75% of the experimental rates at
# 4-6 MPa (measured directly against tests/LMRFMonoAP/test's own rate
# extraction), unlike earlier attempts that pushed rate_multiplier to
# ~1e7-1e8 on a modified mesh and saturated the phase-field interface
# velocity instead of moving it.
DEFAULT_RATE_MULTIPLIER0 = 8.0e5
DEFAULT_ACTIVATION_TEMPERATURE0 = 3145.0

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
    known_p = np.array(sorted(data))
    known_r = np.array([data[p] for p in known_p])
    targets = []
    for p in pressures:
        if p in data:
            targets.append(data[p])
        else:
            targets.append(float(np.interp(p, known_p, known_r)))
    return targets


def run_case(pressure_mpa: float, rate_multiplier: float,
             activation_temperature: float, outdir: Path,
             lowmach_bin: Path) -> subprocess.Popen:
    """Launch tests/LMRFMonoAP/input at the given pressure/coefficients.

    Mirrors the per-pressure #@ args blocks in tests/LMRFMonoAP/input
    (pressure substitution into the four ParmParse keys the test itself
    overrides), plus the 1 MPa case's longer stop_time/coarser plot_dt so
    its slower transient/extinction has time to show up.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    pressure_pa = pressure_mpa * 1.0e6
    args = [
        str(lowmach_bin), str(TEST_INPUT),
        f"Final.density.ic.expression.constant.P={pressure_pa!r}",
        f"component_density.bc.expression.constant.P={pressure_pa!r}",
        f"pressure.ic.constant.value={pressure_pa!r}",
        f"pressure.bc.constant.val.yhi={pressure_pa!r}",
        f"AP_decomposition.phase_change.rate_multiplier={rate_multiplier!r}",
        f"AP_decomposition.phase_change.activation_temperature={activation_temperature!r}_K",
        f"plot_file={outdir}/output",
    ]
    if pressure_mpa == 1.0:
        args += ["stop_time=3.0e-3", "amr.plot_dt=1.5e-4"]
    log_path = outdir / "run.log"
    with open(log_path, "w") as log_fh:
        return subprocess.Popen(args, stdout=log_fh, stderr=subprocess.STDOUT,
                                 cwd=str(REPO_ROOT))


def measure_rate(outdir: Path, pressure_mpa: float) -> float | None:
    """tests/LMRFMonoAP/test's own rate extraction, without the assertions."""
    plotfiles = sorted(glob.glob(f"{outdir}/output/*cell/"))
    if len(plotfiles) < 3:
        return None
    times, positions = [], []
    for path in plotfiles:
        try:
            ds = testlib.yt.load(path)
            data = ds.all_data().to_dataframe(
                [("index", "y"), ("boxlib", "rigid_eta"),
                 ("boxlib", "temperature")])
        except Exception:
            return None
        if not testlib.numpy.all(testlib.numpy.isfinite(data["temperature"].to_numpy())):
            return None
        profile = data.groupby("y", as_index=False)["rigid_eta"].mean()
        y = profile["y"].to_numpy()
        eta = profile["rigid_eta"].to_numpy()
        if not testlib.numpy.all(testlib.numpy.isfinite(eta)):
            return None
        crossing = testlib.numpy.flatnonzero((eta[:-1] >= 0.5) & (eta[1:] < 0.5))
        if not len(crossing):
            return None
        n = crossing[-1]
        position = y[n] + (0.5 - eta[n]) * (y[n + 1] - y[n]) / (eta[n + 1] - eta[n])
        times.append(float(ds.current_time))
        positions.append(position)

    times = testlib.numpy.asarray(times)
    positions = testlib.numpy.asarray(positions)
    fit = times >= (0.75 if pressure_mpa == 1.0 else 0.5) * times[-1]
    if fit.sum() < 2:
        return None
    slope = testlib.numpy.polyfit(times[fit], positions[fit], 1)[0]
    return float(-1000.0 * slope)


def run_sweep(rate_multiplier: float, activation_temperature: float,
              pressures: list[float], workdir: Path,
              lowmach_bin: Path) -> dict[float, float | None]:
    """Run all pressures concurrently, then measure each one's rate."""
    procs = {}
    for p in pressures:
        outdir = workdir / f"P{p:g}MPa"
        procs[p] = (outdir, run_case(p, rate_multiplier, activation_temperature,
                                      outdir, lowmach_bin))
    rates: dict[float, float | None] = {}
    for p, (outdir, proc) in procs.items():
        proc.wait()
        rates[p] = measure_rate(outdir, p)
    return rates


def make_fun_jac(pressures: list[float], targets: list[float], workdir: Path,
                  lowmach_bin: Path, log_path: Path, stage: str = "full",
                  diff_step: float = 0.05):
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
                           eval_dir, lowmach_bin)

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
    parser.add_argument("--lowmach-bin", type=Path, default=DEFAULT_LOWMACH_BIN,
                         help="lowmach binary to run")
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
        bracket_pressures = [min(args.fit_pressures), max(args.fit_pressures)]
        bracket_targets = nearest_experimental(data, bracket_pressures)
        print(f"\n=== Stage 1: bracket fit to {bracket_pressures} ===")
        bracket_fun, bracket_jac = make_fun_jac(
            bracket_pressures, bracket_targets, args.workdir,
            args.lowmach_bin, log_path, stage="bracket")
        bracket_result = least_squares(
            bracket_fun, x0, jac=bracket_jac, bounds=bounds, xtol=args.xtol,
            max_nfev=args.bracket_max_nfev)
        print(f"Stage 1 result: rate_multiplier={10.0 ** bracket_result.x[0]:.6g} "
              f"activation_temperature={bracket_result.x[1]:.6g} "
              f"residual_norm={np.linalg.norm(bracket_result.fun):.6g}")
        x0 = bracket_result.x

    print(f"\n=== Stage 2: full fit to {args.fit_pressures} ===")
    fun, jac = make_fun_jac(args.fit_pressures, targets, args.workdir,
                             args.lowmach_bin, log_path, stage="full")

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

    print("\nRunning full validation sweep over all experimental pressures...")
    all_pressures = sorted(data)
    val_dir = args.workdir / "validation"
    val_dir.mkdir(parents=True, exist_ok=True)
    val_rates = run_sweep(rate_multiplier, activation_temperature, all_pressures,
                          val_dir, args.lowmach_bin)

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
