#!/usr/bin/env python3
"""Broad grid sweep of (rate_multiplier, activation_temperature) against
AP_reg_rate.csv, using tests/LMRFMonoAP/input directly and the same
rate-extraction method as tests/LMRFMonoAP/test (rigid_eta = 0.5 interface
tracking via yt, linear fit over the back half of the run).

Runs here are cheap (tens of seconds to a couple minutes each at these
parameter magnitudes), so this screens the whole grid concurrently instead
of following a gradient from one starting point -- useful for checking
whether the experimental r(P) *slope* (monotonically increasing 2->6 MPa)
is achievable anywhere in the (rate_multiplier, activation_temperature)
plane, not just whether some point matches the magnitude at one pressure.

Writes <workdir>/sweep_results.csv (one row per combo, one column per
pressure) as it goes, so partial progress can be inspected before the
whole grid finishes.
"""

from __future__ import annotations

import csv
import itertools
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
TEST_INPUT = REPO_ROOT / "tests" / "LMRFMonoAP" / "input"
LOWMACH_BIN = REPO_ROOT / "bin" / "lowmach-2d-clang++"

sys.path.insert(0, str(REPO_ROOT / "scripts"))
import testlib  # noqa: E402

PRESSURES = [2.0, 3.0, 4.0, 5.0, 6.0]
RATE_MULTIPLIERS = [1.0e5, 3.0e5, 8.0e5, 2.0e6, 5.0e6, 1.5e7]
ACTIVATION_TEMPERATURES = [1500.0, 2000.0, 2500.0, 3000.0, 3145.0, 3500.0, 4000.0]
RUN_TIMEOUT_S = 240  # kill any single run that gets too stiff to be worth waiting on


def run_and_measure(rate_multiplier: float, activation_temperature: float,
                     pressure_mpa: float, workdir: Path) -> float | None:
    outdir = workdir / f"rm{rate_multiplier:.3g}_Ea{activation_temperature:.0f}_P{pressure_mpa:g}"
    outdir.mkdir(parents=True, exist_ok=True)
    pressure_pa = pressure_mpa * 1.0e6
    args = [
        str(LOWMACH_BIN), str(TEST_INPUT),
        f"Final.density.ic.expression.constant.P={pressure_pa!r}",
        f"component_density.bc.expression.constant.P={pressure_pa!r}",
        f"pressure.ic.constant.value={pressure_pa!r}",
        f"pressure.bc.constant.val.yhi={pressure_pa!r}",
        f"AP_decomposition.phase_change.rate_multiplier={rate_multiplier!r}",
        f"AP_decomposition.phase_change.activation_temperature={activation_temperature!r}_K",
        f"plot_file={outdir}/output",
    ]
    log_path = outdir / "run.log"
    try:
        with open(log_path, "w") as log_fh:
            subprocess.run(args, stdout=log_fh, stderr=subprocess.STDOUT,
                            cwd=str(REPO_ROOT), timeout=RUN_TIMEOUT_S, check=False)
    except subprocess.TimeoutExpired:
        return None

    import glob
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

    times = np.asarray(times)
    positions = np.asarray(positions)
    fit = times >= 0.5 * times[-1]
    if fit.sum() < 2:
        return None
    slope = np.polyfit(times[fit], positions[fit], 1)[0]
    return float(-1000.0 * slope)


def main() -> None:
    workdir = Path(sys.argv[1]) if len(sys.argv) > 1 else REPO_ROOT / "sweep_ap_2026-07-29"
    workdir.mkdir(parents=True, exist_ok=True)
    csv_path = workdir / "sweep_results.csv"

    combos = list(itertools.product(RATE_MULTIPLIERS, ACTIVATION_TEMPERATURES))
    print(f"Sweeping {len(combos)} (rate_multiplier, activation_temperature) combos "
          f"x {len(PRESSURES)} pressures = {len(combos) * len(PRESSURES)} runs", flush=True)

    with open(csv_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["rate_multiplier", "activation_temperature",
                          *[f"rate_{p:g}MPa_mm_s" for p in PRESSURES]])

    pending: dict[tuple[float, float], dict[float, float | None]] = {
        combo: {} for combo in combos
    }
    with ThreadPoolExecutor(max_workers=15) as ex:
        future_to_key = {
            ex.submit(run_and_measure, rm, ea, p, workdir): (rm, ea, p)
            for rm, ea in combos for p in PRESSURES
        }
        for fut in as_completed(future_to_key):
            rm, ea, p = future_to_key[fut]
            pending[(rm, ea)][p] = fut.result()
            if len(pending[(rm, ea)]) == len(PRESSURES):
                row = [rm, ea] + [pending[(rm, ea)][p] if pending[(rm, ea)][p] is not None
                                   else "" for p in PRESSURES]
                with open(csv_path, "a", newline="") as fh:
                    csv.writer(fh).writerow(row)
                rates_str = " ".join(
                    f"{p:g}MPa={pending[(rm, ea)][p]:.4g}"
                    if pending[(rm, ea)][p] is not None else f"{p:g}MPa=FAIL"
                    for p in PRESSURES)
                print(f"rm={rm:.3g} Ea={ea:.0f}  {rates_str}", flush=True)

    print(f"\nWrote {csv_path}")


if __name__ == "__main__":
    main()
