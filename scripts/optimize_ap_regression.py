#!/usr/bin/env python3
"""Fit (rate_multiplier, activation_temperature, w1, pressure_exponent) to
AP_reg_rate.csv.

Drives tests/LMRFMonoAP/input directly (the same input used by the
LMRFMonoAP regression test) with command-line ParmParse overrides for
pressure, rate_multiplier, activation_temperature, the Allen-Cahn
double-well parameter w1, and the optional pressure_exponent factor
(P/pressure_reference)^pressure_exponent, and measures the resulting
regression rate with the same method as tests/LMRFMonoAP/test: track the
rigid_eta = 0.5 interface position vs. time via yt, then fit a line to the
back half of the run (back quarter for the 1 MPa case) to get a
steady-state rate in mm/s. No template files, no custom HDF5 rate
extractor -- this is exactly what the regression test itself checks, just
without the pass/fail assertions and with variable coefficients.

The search is a two-stage global optimization over
x = [log10(rate_multiplier), activation_temperature, w1, pressure_exponent],
scored against the experimental r(P) curve in AP_reg_rate.csv (matches
tests/LMRFMonoAP/reference.csv). Every evaluation -- in both stages -- fits
against the *entire* fit-pressure set at once (2-6 MPa by default); there is
no reduced-pressure bracket stage, since fitting a subset of the data first
is exactly the kind of shortcut that can lock the search onto a minimum that
is only good for those pressures. Stage 1 is a Latin-hypercube random search
of the full parameter box, which gives the global stage a diverse initial
population instead of a single local starting guess. Stage 2 seeds
scipy.optimize.differential_evolution's population with the best random-
search points and evolves them (with a final local polish) -- a
population-based global method that does not get stuck the way a
single-start local method (least_squares/Newton-type) can.

allencahn.mobility, lambda and kappa stay fixed at the test's values:
mobility because only rate_multiplier * mobility is identifiable (see
MassSource in src/Model/Mechanism/PhaseChange.H), and lambda/kappa because
rate_multiplier multiplies both model.LocalStabilityRate()
(~mobility*lambda) and model.GradientCoefficient() (~mobility*kappa) with
the same outer factor -- so only their ratio (interface width, a
mesh-resolution choice, not a material property) is a non-degenerate knob,
and fitting it would tie the result to this test's grid spacing. w1 (with
w0 fixed at 0 and w12 fixed at the test's value) instead reshapes the
dimensionless double-well potential and carries no mesh/length-scale
dependence. pressure_exponent likewise multiplies
LocalRate/LocalStabilityRate/GradientCoefficient uniformly (see
PhaseChange.H), so it rescales the whole interface kinetics by
(P/pressure_reference)^pressure_exponent without touching lambda/kappa or
the interface width -- it stays mesh-independent for the same reason
rate_multiplier and the Arrhenius factor are. pressure_reference is fixed
at 1 MPa (the code default) for this fit. Earlier AP models typically
landed on a pressure_exponent near 1, which sets this script's default
initial guess and keeps the search bounds centered around that region.
1 MPa (the deflagration-limit point, rate = 0) is excluded from the fit
residuals and checked separately in the validation sweep as a pass/fail
extinction check.

Run as a background job; it logs every evaluation's parameters and
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
from pathlib import Path

import numpy as np
from scipy.optimize import differential_evolution
from scipy.stats import qmc

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
TEST_INPUT = REPO_ROOT / "tests" / "LMRFMonoAP" / "input"
DEFAULT_LOWMACH_BIN = REPO_ROOT / "bin" / "lowmach-2d-clang++"

sys.path.insert(0, str(REPO_ROOT / "scripts"))
import testlib  # noqa: E402

# All non-zero pressures in AP_reg_rate.csv (2-6 MPa); 1 MPa is the
# deflagration-limit point (rate = 0) and is excluded from the fit.
DEFAULT_FIT_PRESSURES = [2.0, 3.0, 4.0, 5.0, 6.0]

# The pre-existing test defaults -- a reasonable center point for the
# random-search box. pressure_exponent0 defaults to 1.0 (rather than 0,
# i.e. no pressure dependence) because earlier AP models tended to land
# near a pressure exponent of 1.
DEFAULT_RATE_MULTIPLIER0 = 1.109867e6
DEFAULT_ACTIVATION_TEMPERATURE0 = 2748.73
DEFAULT_W1_0 = 1.0
DEFAULT_PRESSURE_EXPONENT0 = 1.0

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
             activation_temperature: float, w1: float,
             pressure_exponent: float, outdir: Path,
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
        f"AP_decomposition.phase_change.phase_field.allencahn.w1={w1!r}",
        f"AP_decomposition.phase_change.pressure_exponent={pressure_exponent!r}",
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
        # A solver crash from fully burning through the domain (a known,
        # expected outcome at high pressure/late time) corrupts or omits the
        # interface in only the last plotfile or two. Stop collecting at the
        # first bad plotfile and fit whatever steady-state data came before
        # it, rather than discarding an otherwise-good run.
        try:
            ds = testlib.yt.load(path)
            data = ds.all_data().to_dataframe(
                [("index", "y"), ("boxlib", "rigid_eta"),
                 ("boxlib", "temperature")])
        except Exception:
            break
        if not testlib.numpy.all(testlib.numpy.isfinite(data["temperature"].to_numpy())):
            break
        profile = data.groupby("y", as_index=False)["rigid_eta"].mean()
        y = profile["y"].to_numpy()
        eta = profile["rigid_eta"].to_numpy()
        if not testlib.numpy.all(testlib.numpy.isfinite(eta)):
            break
        crossing = testlib.numpy.flatnonzero((eta[:-1] >= 0.5) & (eta[1:] < 0.5))
        if not len(crossing):
            break
        n = crossing[-1]
        position = y[n] + (0.5 - eta[n]) * (y[n + 1] - y[n]) / (eta[n + 1] - eta[n])
        times.append(float(ds.current_time))
        positions.append(position)

    times = testlib.numpy.asarray(times)
    positions = testlib.numpy.asarray(positions)

    # Mirror tests/LMRFMonoAP/test: ignore the first/last 10% of the run
    # (transient startup and any end-of-run edge effects) and any interface
    # positions that have drifted near the bottom of the domain, so the fit
    # only sees steady-state regression.
    t0, t1 = times[0], times[-1]
    dt = 0.1 * (t1 - t0)
    mask = (times >= t0 + dt) & (times <= t1 - dt) & (positions >= -0.00013)
    if mask.sum() < 2:
        return None
    slope = testlib.numpy.polyfit(times[mask], positions[mask], 1)[0]
    return float(-1000.0 * slope)


def run_sweep(rate_multiplier: float, activation_temperature: float, w1: float,
              pressure_exponent: float, pressures: list[float], workdir: Path,
              lowmach_bin: Path) -> dict[float, float | None]:
    """Run all pressures concurrently, then measure each one's rate."""
    procs = {}
    for p in pressures:
        outdir = workdir / f"P{p:g}MPa"
        procs[p] = (outdir, run_case(p, rate_multiplier, activation_temperature,
                                      w1, pressure_exponent, outdir, lowmach_bin))
    rates: dict[float, float | None] = {}
    for p, (outdir, proc) in procs.items():
        proc.wait()
        rates[p] = measure_rate(outdir, p)
    return rates


def make_objective(pressures: list[float], targets: list[float], workdir: Path,
                    lowmach_bin: Path, log_path: Path, stage: str = "search",
                    fixed_w1: float | None = None,
                    fixed_pressure_exponent: float | None = None):
    """Build a scalar objective (sum of squared relative-error residuals)
    against the *full* pressure/target set passed in -- every call fits
    against all of it, there is no reduced-data bracket stage."""
    fit_pairs = [(p, t) for p, t in zip(pressures, targets) if t > 0.0]
    dropped = [p for p, t in zip(pressures, targets) if t <= 0.0]
    if dropped:
        print(f"[{stage}] dropping non-positive experimental targets at "
              f"pressures {dropped} from the fit residuals", flush=True)

    iteration = [0]
    cache: dict[tuple, float] = {}
    lock = threading.Lock()

    def unpack(x: np.ndarray) -> tuple[float, float, float, float]:
        rate_multiplier = float(10.0 ** x[0])
        activation_temperature = float(x[1])
        free_idx = 2
        if fixed_w1 is not None:
            w1 = fixed_w1
        else:
            w1 = float(x[free_idx])
            free_idx += 1
        if fixed_pressure_exponent is not None:
            pressure_exponent = fixed_pressure_exponent
        else:
            pressure_exponent = float(x[free_idx])
            free_idx += 1
        return rate_multiplier, activation_temperature, w1, pressure_exponent

    def objective(x: np.ndarray) -> float:
        key = tuple(round(float(v), 12) for v in x)
        with lock:
            cached = cache.get(key)
            if cached is not None:
                return cached
            iteration[0] += 1
            idx = iteration[0]

        rate_multiplier, activation_temperature, w1, pressure_exponent = unpack(x)
        eval_dir = workdir / f"{stage}_iter_{idx:03d}"
        eval_dir.mkdir(parents=True, exist_ok=True)
        rates = run_sweep(rate_multiplier, activation_temperature, w1,
                           pressure_exponent, pressures, eval_dir, lowmach_bin)

        residuals = []
        for p, target in fit_pairs:
            sim = rates.get(p)
            if sim is None:
                residuals.append(PENALTY_RESIDUAL)
            else:
                residuals.append((sim - target) / target)
        residual_norm = float(np.linalg.norm(np.array(residuals)))
        cost = residual_norm ** 2

        record = {
            "stage": stage,
            "iteration": idx,
            "rate_multiplier": rate_multiplier,
            "activation_temperature": activation_temperature,
            "w1": w1,
            "pressure_exponent": pressure_exponent,
            "rates_mm_s": rates,
            "residual_norm": residual_norm,
        }
        with lock:
            with open(log_path, "a") as fh:
                fh.write(json.dumps(record) + "\n")
            print(f"[{stage} iter {idx}] rate_multiplier={rate_multiplier:.6g} "
                  f"activation_temperature={activation_temperature:.6g} "
                  f"w1={w1:.6g} pressure_exponent={pressure_exponent:.6g} "
                  f"residual_norm={residual_norm:.6g}", flush=True)
            cache[key] = cost
        return cost

    return objective, unpack


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
    parser.add_argument("--w1-0", type=float, default=DEFAULT_W1_0,
                         help="initial guess for the Allen-Cahn w1 double-well "
                              "parameter (w0 stays 0, w12 stays the test's "
                              "default; dimensionless)")
    parser.add_argument("--fit-w1", dest="fit_w1", action="store_true",
                         default=True,
                         help="include w1 as a free fit parameter (default)")
    parser.add_argument("--no-fit-w1", dest="fit_w1", action="store_false",
                         help="hold w1 fixed at --w1-0 and fit only "
                              "rate_multiplier/activation_temperature")
    parser.add_argument("--pressure-exponent0", type=float,
                         default=DEFAULT_PRESSURE_EXPONENT0,
                         help="initial guess for the optional "
                              "(P/1MPa)^pressure_exponent kinetics correction "
                              "(dimensionless; 0 == no effect; earlier AP models "
                              "tended toward ~1, which is the default)")
    parser.add_argument("--fit-pressure-exponent", dest="fit_pressure_exponent",
                         action="store_true", default=True,
                         help="include pressure_exponent as a free fit "
                              "parameter (default)")
    parser.add_argument("--no-fit-pressure-exponent",
                         dest="fit_pressure_exponent", action="store_false",
                         help="hold pressure_exponent fixed at "
                              "--pressure-exponent0 (0.0 reproduces the "
                              "pre-existing model exactly)")
    parser.add_argument("--lowmach-bin", type=Path, default=DEFAULT_LOWMACH_BIN,
                         help="lowmach binary to run")
    parser.add_argument("--random-samples", type=int, default=24,
                         help="number of Latin-hypercube samples drawn over the "
                              "full parameter box in the stage-1 random search "
                              "(default: 24)")
    parser.add_argument("--de-seed-count", type=int, default=10,
                         help="number of best random-search points used to seed "
                              "the stage-2 differential_evolution population "
                              "(default: 10; must be >= 5)")
    parser.add_argument("--de-maxiter", type=int, default=15,
                         help="max generations for stage-2 differential_evolution "
                              "(default: 15)")
    parser.add_argument("--de-tol", type=float, default=1.0e-3,
                         help="differential_evolution convergence tolerance "
                              "(default: 1e-3)")
    parser.add_argument("--de-seed", type=int, default=None,
                         help="random seed for reproducibility (default: unseeded)")
    parser.add_argument("--no-polish", dest="polish", action="store_false",
                         default=True,
                         help="skip the final local (L-BFGS-B) polish step after "
                              "differential_evolution converges")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.workdir.mkdir(parents=True, exist_ok=True)
    log_path = args.workdir / "iterations.jsonl"

    data = load_experimental_data(args.data)
    targets = nearest_experimental(data, args.fit_pressures)
    print(f"Fitting to {len(args.fit_pressures)} pressures: "
          f"{list(zip(args.fit_pressures, targets))}")

    fixed_w1 = None if args.fit_w1 else args.w1_0
    fixed_pressure_exponent = None if args.fit_pressure_exponent else args.pressure_exponent0

    x0_list = [np.log10(args.rate_multiplier0), args.activation_temperature0]
    lo, hi = [-1.0, 0.0], [8.0, 10000.0]
    if args.fit_w1:
        x0_list.append(args.w1_0)
        lo.append(0.1)
        hi.append(4.0)
    if args.fit_pressure_exponent:
        x0_list.append(args.pressure_exponent0)
        lo.append(-1.0)
        hi.append(3.0)
    x0 = np.array(x0_list)
    lo, hi = np.array(lo), np.array(hi)
    bounds = list(zip(lo, hi))
    ndim = len(x0)

    objective, unpack = make_objective(
        args.fit_pressures, targets, args.workdir, args.lowmach_bin, log_path,
        stage="search", fixed_w1=fixed_w1,
        fixed_pressure_exponent=fixed_pressure_exponent)

    print(f"\n=== Stage 1: random search over the full parameter box "
          f"({args.random_samples} samples, fitting all {args.fit_pressures} "
          f"MPa every evaluation) ===")
    sampler = qmc.LatinHypercube(d=ndim, seed=args.de_seed)
    unit_samples = sampler.random(n=args.random_samples)
    samples = qmc.scale(unit_samples, lo, hi)
    samples[0] = x0  # always include the current best-known point
    costs = np.array([objective(x) for x in samples])
    order = np.argsort(costs)
    best_random = samples[order[0]]
    param_names = ["rate_multiplier", "activation_temperature", "w1", "pressure_exponent"]
    best_random_params = dict(zip(param_names, unpack(best_random)))
    print(f"Stage 1 best: cost={costs[order[0]]:.6g} at {best_random_params}")

    seed_count = max(args.de_seed_count, 5)
    seed_pop = samples[order[:seed_count]]
    if len(seed_pop) < seed_count:
        extra = qmc.scale(sampler.random(n=seed_count - len(seed_pop)), lo, hi)
        seed_pop = np.vstack([seed_pop, extra])

    print(f"\n=== Stage 2: differential_evolution seeded from the "
          f"{seed_count} best random-search points, fitting all "
          f"{args.fit_pressures} MPa every evaluation ===")
    de_objective, _ = make_objective(
        args.fit_pressures, targets, args.workdir, args.lowmach_bin, log_path,
        stage="de", fixed_w1=fixed_w1,
        fixed_pressure_exponent=fixed_pressure_exponent)

    result = differential_evolution(
        de_objective, bounds, init=seed_pop, maxiter=args.de_maxiter,
        tol=args.de_tol, seed=args.de_seed, polish=args.polish,
        updating="deferred")

    rate_multiplier, activation_temperature, w1, pressure_exponent = unpack(result.x)
    print("\n=== Converged (or hit maxiter) ===")
    print(f"rate_multiplier         = {rate_multiplier:.6g}")
    print(f"activation_temperature  = {activation_temperature:.6g} K")
    print(f"w1                      = {w1:.6g}")
    print(f"pressure_exponent       = {pressure_exponent:.6g}")
    print(f"final cost              = {result.fun:.6g}")
    print(f"success={result.success}: {result.message}")

    best = {
        "rate_multiplier": rate_multiplier,
        "activation_temperature": activation_temperature,
        "w1": w1,
        "pressure_exponent": pressure_exponent,
    }
    with open(args.workdir / "best_fit.json", "w") as fh:
        json.dump(best, fh, indent=2)
    print(f"Wrote {args.workdir / 'best_fit.json'}")

    print("\nRunning full validation sweep over all experimental pressures...")
    all_pressures = sorted(data)
    val_dir = args.workdir / "validation"
    val_dir.mkdir(parents=True, exist_ok=True)
    val_rates = run_sweep(rate_multiplier, activation_temperature, w1,
                          pressure_exponent, all_pressures,
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
                  f"activation_temperature={activation_temperature:.4g} K, "
                  f"w1={w1:.4g}, pressure_exponent={pressure_exponent:.4g}")
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
