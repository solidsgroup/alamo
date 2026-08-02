#!/usr/bin/env python3
"""Fit (rate_multiplier, activation_temperature, w1, pressure_exponent,
pressure_reference, latent_heat) to AP_reg_rate.csv.

Drives tests/LMRFMonoAP/input directly (the same input used by the
LMRFMonoAP regression test) with command-line ParmParse overrides for
pressure, rate_multiplier, activation_temperature, the Allen-Cahn
double-well parameter w1, the optional pressure_exponent factor
(P/pressure_reference)^pressure_exponent, and latent_heat, and measures the
resulting regression rate with the same method as tests/LMRFMonoAP/test:
track the rigid_eta = 0.5 interface position vs. time via yt, then fit a
line to the back half of the run (back quarter for the 1 MPa case) to get a
steady-state rate in mm/s. No template files, no custom HDF5 rate
extractor -- this is exactly what the regression test itself checks, just
without the pass/fail assertions and with variable coefficients.

The search is a two-stage global optimization over
x = [log10(rate_multiplier), activation_temperature, w1, pressure_exponent,
     pressure_reference, latent_heat] (any of
w1/pressure_exponent/pressure_reference/latent_heat can be held fixed
instead of fit -- see --no-fit-*), scored against the experimental
r(P) curve in AP_reg_rate.csv (matches tests/LMRFMonoAP/reference.csv).
Every evaluation -- in both stages -- fits against the *entire*
fit-pressure set at once (2-6 MPa by default); there is no reduced-pressure
bracket stage, since fitting a subset of the data first is exactly the
kind of shortcut that can lock the search onto a minimum that is only good
for those pressures. Stage 1 is a Latin-hypercube random search of the
full parameter box, which gives the global stage a diverse initial
population instead of a single local starting guess. Stage 2 seeds a
particle-swarm optimization (PSO, hand-rolled below -- no extra dependency)
with the best random-search points and evolves them (with a final local
L-BFGS-B polish) -- a population-based global method that does not get
stuck the way a single-start local method (least_squares/Newton-type) can.

allencahn.mobility, lambda and kappa stay fixed at the test's values:
mobility because only rate_multiplier * mobility is identifiable (see
MassSource in src/Model/Mechanism/PhaseChange.H -- verified exactly
degenerate: LocalRate, GradientCoefficient, and LocalStabilityRate all
scale as the simple product rate_multiplier * mobility with no term where
they appear separately), and lambda/kappa because rate_multiplier
multiplies both model.LocalStabilityRate() (~mobility*lambda) and
model.GradientCoefficient() (~mobility*kappa) with the same outer factor --
so only their ratio (interface width, a mesh-resolution choice, not a
material property) is a non-degenerate knob, and fitting it would tie the
result to this test's grid spacing. w1 (with w0 fixed at 0 and w12 fixed
at the test's default, i.e. the Allen-Cahn double-well potential is left
at its default shape) instead reshapes the dimensionless double-well
potential and carries no mesh/length-scale dependence. pressure_exponent
likewise multiplies LocalRate/LocalStabilityRate/GradientCoefficient
uniformly (see PhaseChange.H), so it rescales the whole interface kinetics
by (P/pressure_reference)^pressure_exponent without touching lambda/kappa
or the interface width -- it stays mesh-independent for the same reason
rate_multiplier and the Arrhenius factor are, and is fit by default
(--fit-pressure-exponent). pressure_reference (the P0 in
(P/P0)^pressure_exponent) is also fit by default, strictly bounded above 0
(see --pressure-reference-min) since a non-positive reference pressure is
not physical/would blow up the power law. latent_heat sets the
enthalpy cost of solid->gas conversion (HeatSource = latent_heat *
MassSource, a direct sink in the energy equation at the interface); it is
a physical material property (not mesh-dependent), but literature values
vary, so a bounded amount of tuning around the test's default (100 cal/g)
is included here -- strictly bounded above 0 (see --latent-heat-min), since
latent_heat <= 0 is not physical. 1 MPa (the deflagration-limit point,
rate = 0) is excluded from the fit residuals and checked separately in the
validation sweep as a pass/fail extinction check.

Run as a background job; it logs every evaluation's parameters and
residual norm to <workdir>/iterations.jsonl as it goes, and writes
<workdir>/best_fit.json plus a validation plot when it converges.

Example
-------
    python scripts/optimize_ap_regression.py --workdir /tmp/ap_calib
"""

from __future__ import annotations

import argparse
import concurrent.futures
import csv
import glob
import json
import subprocess
import sys
import threading
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
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
# random-search box. pressure_exponent0 defaults to 0.0 (no pressure
# effect) per this script's current configuration; latent_heat0 is the
# test's literature default (100 cal/g).
DEFAULT_RATE_MULTIPLIER0 = 1.109867e6
DEFAULT_ACTIVATION_TEMPERATURE0 = 2748.73
DEFAULT_W1_0 = 1.0
DEFAULT_PRESSURE_EXPONENT0 = 0.0
DEFAULT_PRESSURE_REFERENCE0_MPA = 1.0
DEFAULT_LATENT_HEAT0 = 100.0

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
             pressure_exponent: float, pressure_reference_mpa: float,
             latent_heat: float, outdir: Path,
             lowmach_bin: Path) -> subprocess.Popen:
    """Launch tests/LMRFMonoAP/input at the given pressure/coefficients.

    Mirrors the per-pressure #@ args blocks in tests/LMRFMonoAP/input
    (pressure substitution into the four ParmParse keys the test itself
    overrides), plus the 1 MPa case's longer stop_time/coarser plot_dt so
    its slower transient/extinction has time to show up.
    """
    outdir.mkdir(parents=True, exist_ok=True)
    pressure_pa = pressure_mpa * 1.0e6
    pressure_reference_pa = pressure_reference_mpa * 1.0e6
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
        f"AP_decomposition.phase_change.pressure_reference={pressure_reference_pa!r}_Pa",
        f"AP_decomposition.phase_change.latent_heat={latent_heat!r}_cal/g",
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
              pressure_exponent: float, pressure_reference_mpa: float,
              latent_heat: float,
              pressures: list[float], workdir: Path,
              lowmach_bin: Path) -> dict[float, float | None]:
    """Run all pressures concurrently, then measure each one's rate."""
    procs = {}
    for p in pressures:
        outdir = workdir / f"P{p:g}MPa"
        procs[p] = (outdir, run_case(p, rate_multiplier, activation_temperature,
                                      w1, pressure_exponent, pressure_reference_mpa,
                                      latent_heat, outdir, lowmach_bin))
    rates: dict[float, float | None] = {}
    for p, (outdir, proc) in procs.items():
        proc.wait()
        rates[p] = measure_rate(outdir, p)
    return rates


def make_objective(pressures: list[float], targets: list[float], workdir: Path,
                    lowmach_bin: Path, log_path: Path, stage: str = "search",
                    fixed_w1: float | None = None,
                    fixed_pressure_exponent: float | None = None,
                    fixed_pressure_reference: float | None = None,
                    fixed_latent_heat: float | None = None):
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

    def unpack(x: np.ndarray) -> tuple[float, float, float, float, float, float]:
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
        if fixed_pressure_reference is not None:
            pressure_reference = fixed_pressure_reference
        else:
            pressure_reference = float(x[free_idx])
            free_idx += 1
        if fixed_latent_heat is not None:
            latent_heat = fixed_latent_heat
        else:
            latent_heat = float(x[free_idx])
            free_idx += 1
        return (rate_multiplier, activation_temperature, w1, pressure_exponent,
                pressure_reference, latent_heat)

    def objective(x: np.ndarray) -> float:
        key = tuple(round(float(v), 12) for v in x)
        with lock:
            cached = cache.get(key)
            if cached is not None:
                return cached
            iteration[0] += 1
            idx = iteration[0]

        (rate_multiplier, activation_temperature, w1, pressure_exponent,
         pressure_reference, latent_heat) = unpack(x)
        eval_dir = workdir / f"{stage}_iter_{idx:03d}"
        eval_dir.mkdir(parents=True, exist_ok=True)
        rates = run_sweep(rate_multiplier, activation_temperature, w1,
                           pressure_exponent, pressure_reference, latent_heat,
                           pressures, eval_dir, lowmach_bin)

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
            "pressure_reference_mpa": pressure_reference,
            "latent_heat": latent_heat,
            "rates_mm_s": rates,
            "residual_norm": residual_norm,
        }
        with lock:
            with open(log_path, "a") as fh:
                fh.write(json.dumps(record) + "\n")
            print(f"[{stage} iter {idx}] rate_multiplier={rate_multiplier:.6g} "
                  f"activation_temperature={activation_temperature:.6g} "
                  f"w1={w1:.6g} pressure_exponent={pressure_exponent:.6g} "
                  f"pressure_reference_mpa={pressure_reference:.6g} "
                  f"latent_heat={latent_heat:.6g} "
                  f"residual_norm={residual_norm:.6g}", flush=True)
            cache[key] = cost
        return cost

    return objective, unpack


def particle_swarm_optimize(objective, bounds, init_positions: np.ndarray,
                             maxiter: int, tol: float, seed: int | None,
                             pool: concurrent.futures.Executor,
                             inertia: float = 0.7, cognitive: float = 1.5,
                             social: float = 1.5,
                             velocity_clamp_fraction: float = 0.2,
                             stall_limit: int = 5):
    """Hand-rolled particle swarm optimization (no external PSO dependency).

    Each particle has a position and velocity in the bounded parameter box;
    every generation, velocities are updated by a weighted combination of
    inertia (the particle's own momentum), a cognitive term (pull toward
    that particle's own best-seen position), and a social term (pull toward
    the swarm's best-seen position), then positions are advanced and
    clamped back into bounds (with the corresponding velocity component
    damped/reversed on clamping, so a particle that hits a wall doesn't
    just get stuck repeatedly overshooting it). All particles in a
    generation are evaluated concurrently via `pool.map`, mirroring how
    stage 1's random search and the old differential_evolution stage were
    parallelized. Stops early if the global best hasn't improved by more
    than `tol` for `stall_limit` consecutive generations.
    """
    rng = np.random.default_rng(seed)
    lo = np.array([b[0] for b in bounds])
    hi = np.array([b[1] for b in bounds])
    positions = init_positions.copy()
    n_particles = positions.shape[0]
    vmax = velocity_clamp_fraction * (hi - lo)
    velocities = rng.uniform(-1.0, 1.0, size=positions.shape) * vmax

    costs = np.array(list(pool.map(objective, positions)))
    personal_best_pos = positions.copy()
    personal_best_cost = costs.copy()
    g_idx = int(np.argmin(personal_best_cost))
    global_best_pos = personal_best_pos[g_idx].copy()
    global_best_cost = float(personal_best_cost[g_idx])
    print(f"[pso gen 0] global_best_cost={global_best_cost:.6g}", flush=True)

    stall = 0
    for gen in range(1, maxiter + 1):
        r1 = rng.random(positions.shape)
        r2 = rng.random(positions.shape)
        velocities = (inertia * velocities
                      + cognitive * r1 * (personal_best_pos - positions)
                      + social * r2 * (global_best_pos - positions))
        velocities = np.clip(velocities, -vmax, vmax)
        positions = positions + velocities

        below = positions < lo
        above = positions > hi
        positions = np.clip(positions, lo, hi)
        # Damp-and-reverse velocity on any component that hit a wall, so
        # particles don't keep pinning themselves against the boundary.
        velocities[below] *= -0.5
        velocities[above] *= -0.5

        costs = np.array(list(pool.map(objective, positions)))
        improved = costs < personal_best_cost
        personal_best_pos[improved] = positions[improved]
        personal_best_cost[improved] = costs[improved]

        g_idx = int(np.argmin(personal_best_cost))
        candidate_cost = float(personal_best_cost[g_idx])
        if candidate_cost < global_best_cost - tol:
            stall = 0
        else:
            stall += 1
        if candidate_cost < global_best_cost:
            global_best_cost = candidate_cost
            global_best_pos = personal_best_pos[g_idx].copy()

        print(f"[pso gen {gen}] global_best_cost={global_best_cost:.6g} "
              f"stall={stall}/{stall_limit}", flush=True)
        if stall >= stall_limit:
            print(f"PSO converged: no improvement > tol for {stall_limit} "
                  f"generations", flush=True)
            break

    return global_best_pos, global_best_cost, n_particles


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
                         help="value (if fixed) or initial guess (if fit) for "
                              "the optional (P/1MPa)^pressure_exponent "
                              "kinetics correction (dimensionless; 0 == no "
                              "pressure effect, the default)")
    parser.add_argument("--fit-pressure-exponent", dest="fit_pressure_exponent",
                         action="store_true", default=True,
                         help="include pressure_exponent as a free fit "
                              "parameter (default)")
    parser.add_argument("--no-fit-pressure-exponent",
                         dest="fit_pressure_exponent", action="store_false",
                         help="hold pressure_exponent fixed at "
                              "--pressure-exponent0")
    parser.add_argument("--pressure-reference0", type=float,
                         default=DEFAULT_PRESSURE_REFERENCE0_MPA,
                         help="value (if fixed) or initial guess (if fit) for "
                              "pressure_reference [MPa], the P0 in "
                              "(P/P0)^pressure_exponent (default: 1.0 MPa, the "
                              "code default)")
    parser.add_argument("--fit-pressure-reference",
                         dest="fit_pressure_reference", action="store_true",
                         default=True,
                         help="include pressure_reference as a free fit "
                              "parameter (default), strictly bounded above 0 "
                              "-- see --pressure-reference-min/-max")
    parser.add_argument("--no-fit-pressure-reference",
                         dest="fit_pressure_reference", action="store_false",
                         help="hold pressure_reference fixed at "
                              "--pressure-reference0")
    parser.add_argument("--pressure-reference-min", type=float, default=0.1,
                         help="lower search bound for pressure_reference [MPa] "
                              "(default: 0.1; must be > 0, a non-positive "
                              "reference pressure is not physical)")
    parser.add_argument("--pressure-reference-max", type=float, default=10.0,
                         help="upper search bound for pressure_reference [MPa] "
                              "(default: 10.0)")
    parser.add_argument("--latent-heat0", type=float,
                         default=DEFAULT_LATENT_HEAT0,
                         help="value (if fixed) or initial guess (if fit) for "
                              "latent_heat [cal/g] (default: 100, the test's "
                              "literature value)")
    parser.add_argument("--fit-latent-heat", dest="fit_latent_heat",
                         action="store_true", default=True,
                         help="include latent_heat as a free fit parameter "
                              "(default), bounded strictly above 0 -- see "
                              "--latent-heat-min/--latent-heat-max")
    parser.add_argument("--no-fit-latent-heat", dest="fit_latent_heat",
                         action="store_false",
                         help="hold latent_heat fixed at --latent-heat0")
    parser.add_argument("--latent-heat-min", type=float, default=0.1,
                         help="lower search bound for latent_heat [cal/g] "
                              "(default: 0.1; must be > 0, latent_heat <= 0 "
                              "is not physical)")
    parser.add_argument("--latent-heat-max", type=float, default=200.0,
                         help="upper search bound for latent_heat [cal/g] "
                              "(default: 200)")
    parser.add_argument("--lowmach-bin", type=Path, default=DEFAULT_LOWMACH_BIN,
                         help="lowmach binary to run")
    parser.add_argument("--random-samples", type=int, default=24,
                         help="number of Latin-hypercube samples drawn over the "
                              "full parameter box in the stage-1 random search "
                              "(default: 24)")
    parser.add_argument("--pso-seed-count", type=int, default=10,
                         help="number of best random-search points used to seed "
                              "the stage-2 particle-swarm population (default: "
                              "10; must be >= 5)")
    parser.add_argument("--pso-maxiter", type=int, default=15,
                         help="max generations for stage-2 particle swarm "
                              "(default: 15)")
    parser.add_argument("--pso-tol", type=float, default=1.0e-3,
                         help="minimum global-best improvement per generation "
                              "before it counts toward the stall/convergence "
                              "counter (default: 1e-3)")
    parser.add_argument("--pso-stall-limit", type=int, default=5,
                         help="stop stage 2 after this many consecutive "
                              "generations without a > --pso-tol improvement "
                              "in the global best (default: 5)")
    parser.add_argument("--pso-inertia", type=float, default=0.7,
                         help="PSO inertia weight (default: 0.7)")
    parser.add_argument("--pso-cognitive", type=float, default=1.5,
                         help="PSO cognitive (personal-best pull) coefficient "
                              "(default: 1.5)")
    parser.add_argument("--pso-social", type=float, default=1.5,
                         help="PSO social (global-best pull) coefficient "
                              "(default: 1.5)")
    parser.add_argument("--pso-seed", type=int, default=None,
                         help="random seed for reproducibility (default: unseeded)")
    parser.add_argument("--no-polish", dest="polish", action="store_false",
                         default=True,
                         help="skip the final local (L-BFGS-B) polish step after "
                              "the particle swarm converges")
    parser.add_argument("--eval-workers", type=int, default=3,
                         help="number of parameter sets evaluated concurrently "
                              "in both stages (each evaluation itself launches "
                              "its 5 pressures concurrently, so total sim "
                              "processes in flight = eval-workers * 5; default "
                              "3, i.e. up to 15 concurrent lowmach processes)")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.workdir.mkdir(parents=True, exist_ok=True)
    log_path = args.workdir / "iterations.jsonl"

    if args.latent_heat_min <= 0.0:
        raise SystemExit("--latent-heat-min must be > 0 (latent_heat <= 0 is "
                          "not physical)")
    if args.pressure_reference_min <= 0.0:
        raise SystemExit("--pressure-reference-min must be > 0 (a "
                          "non-positive reference pressure is not physical)")

    data = load_experimental_data(args.data)
    targets = nearest_experimental(data, args.fit_pressures)
    print(f"Fitting to {len(args.fit_pressures)} pressures: "
          f"{list(zip(args.fit_pressures, targets))}")

    fixed_w1 = None if args.fit_w1 else args.w1_0
    fixed_pressure_exponent = None if args.fit_pressure_exponent else args.pressure_exponent0
    fixed_pressure_reference = None if args.fit_pressure_reference else args.pressure_reference0
    fixed_latent_heat = None if args.fit_latent_heat else args.latent_heat0

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
    if args.fit_pressure_reference:
        x0_list.append(args.pressure_reference0)
        lo.append(args.pressure_reference_min)
        hi.append(args.pressure_reference_max)
    if args.fit_latent_heat:
        x0_list.append(args.latent_heat0)
        lo.append(args.latent_heat_min)
        hi.append(args.latent_heat_max)
    x0 = np.array(x0_list)
    lo, hi = np.array(lo), np.array(hi)
    bounds = list(zip(lo, hi))
    ndim = len(x0)

    objective, unpack = make_objective(
        args.fit_pressures, targets, args.workdir, args.lowmach_bin, log_path,
        stage="search", fixed_w1=fixed_w1,
        fixed_pressure_exponent=fixed_pressure_exponent,
        fixed_pressure_reference=fixed_pressure_reference,
        fixed_latent_heat=fixed_latent_heat)

    print(f"\n=== Stage 1: random search over the full parameter box "
          f"({args.random_samples} samples, fitting all {args.fit_pressures} "
          f"MPa every evaluation, {args.eval_workers} evaluated concurrently) ===")
    sampler = qmc.LatinHypercube(d=ndim, seed=args.pso_seed)
    unit_samples = sampler.random(n=args.random_samples)
    samples = qmc.scale(unit_samples, lo, hi)
    samples[0] = x0  # always include the current best-known point
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.eval_workers) as pool:
        costs = np.array(list(pool.map(objective, samples)))
    order = np.argsort(costs)
    best_random = samples[order[0]]
    param_names = ["rate_multiplier", "activation_temperature", "w1",
                   "pressure_exponent", "pressure_reference", "latent_heat"]
    best_random_params = dict(zip(param_names, unpack(best_random)))
    print(f"Stage 1 best: cost={costs[order[0]]:.6g} at {best_random_params}")

    seed_count = max(args.pso_seed_count, 5)
    seed_pop = samples[order[:seed_count]]
    if len(seed_pop) < seed_count:
        extra = qmc.scale(sampler.random(n=seed_count - len(seed_pop)), lo, hi)
        seed_pop = np.vstack([seed_pop, extra])

    print(f"\n=== Stage 2: particle swarm seeded from the {seed_count} best "
          f"random-search points, fitting all {args.fit_pressures} MPa every "
          f"evaluation, {args.eval_workers} particles evaluated concurrently "
          f"===")
    pso_objective, _ = make_objective(
        args.fit_pressures, targets, args.workdir, args.lowmach_bin, log_path,
        stage="pso", fixed_w1=fixed_w1,
        fixed_pressure_exponent=fixed_pressure_exponent,
        fixed_pressure_reference=fixed_pressure_reference,
        fixed_latent_heat=fixed_latent_heat)

    with concurrent.futures.ThreadPoolExecutor(max_workers=args.eval_workers) as pool:
        best_x, best_cost, n_particles = particle_swarm_optimize(
            pso_objective, bounds, seed_pop, maxiter=args.pso_maxiter,
            tol=args.pso_tol, seed=args.pso_seed, pool=pool,
            inertia=args.pso_inertia, cognitive=args.pso_cognitive,
            social=args.pso_social, stall_limit=args.pso_stall_limit)

        success = True
        message = "particle swarm stall/maxiter"
        if args.polish:
            print(f"\n=== Polish: local L-BFGS-B from the swarm's global "
                  f"best (cost={best_cost:.6g}) ===")
            polish_result = minimize(
                pso_objective, best_x, method="L-BFGS-B", bounds=bounds)
            print(f"Polish result: cost={polish_result.fun:.6g} "
                  f"success={polish_result.success}: {polish_result.message}")
            if polish_result.fun <= best_cost:
                best_x = polish_result.x
                best_cost = float(polish_result.fun)
            success = bool(polish_result.success)
            message = str(polish_result.message)

    (rate_multiplier, activation_temperature, w1, pressure_exponent,
     pressure_reference, latent_heat) = unpack(best_x)
    print("\n=== Converged (or hit maxiter/stall limit) ===")
    print(f"rate_multiplier         = {rate_multiplier:.6g}")
    print(f"activation_temperature  = {activation_temperature:.6g} K")
    print(f"w1                      = {w1:.6g}")
    print(f"pressure_exponent       = {pressure_exponent:.6g}")
    print(f"pressure_reference      = {pressure_reference:.6g} MPa")
    print(f"latent_heat             = {latent_heat:.6g} cal/g")
    print(f"final cost              = {best_cost:.6g}")
    print(f"success={success}: {message}")

    best = {
        "rate_multiplier": rate_multiplier,
        "activation_temperature": activation_temperature,
        "w1": w1,
        "pressure_exponent": pressure_exponent,
        "pressure_reference_mpa": pressure_reference,
        "latent_heat": latent_heat,
    }
    with open(args.workdir / "best_fit.json", "w") as fh:
        json.dump(best, fh, indent=2)
    print(f"Wrote {args.workdir / 'best_fit.json'}")

    print("\nRunning full validation sweep over all experimental pressures...")
    all_pressures = sorted(data)
    val_dir = args.workdir / "validation"
    val_dir.mkdir(parents=True, exist_ok=True)
    val_rates = run_sweep(rate_multiplier, activation_temperature, w1,
                          pressure_exponent, pressure_reference, latent_heat,
                          all_pressures, val_dir, args.lowmach_bin)

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
                  f"w1={w1:.4g}, pressure_exponent={pressure_exponent:.4g}, "
                  f"pressure_reference={pressure_reference:.4g} MPa, "
                  f"latent_heat={latent_heat:.4g} cal/g")
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
