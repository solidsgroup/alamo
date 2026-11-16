#!/usr/bin/env python3
"""Compare matched 1D AMReX regression histories exported to CSV."""

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def load_history(path):
    with path.open(newline="") as source:
        rows = list(csv.DictReader(source))
    if not rows:
        raise ValueError(f"No data rows in {path}")
    required = ("time_s", "front_x_m", "front_temperature_K", "pressure_min_Pa", "pressure_max_Pa")
    for name in required:
        if name not in rows[0]:
            raise ValueError(f"{path} lacks required column {name}")
    time = np.asarray([float(row["time_s"]) for row in rows])
    front = np.asarray([float(row["front_x_m"]) for row in rows])
    temperature = np.asarray([float(row["front_temperature_K"]) for row in rows])
    pressure_min = np.asarray([float(row["pressure_min_Pa"]) for row in rows])
    pressure_max = np.asarray([float(row["pressure_max_Pa"]) for row in rows])
    if not all(np.all(np.isfinite(values)) for values in (time, front, temperature, pressure_min, pressure_max)):
        raise ValueError(f"{path} has nonfinite required data")
    if np.any(np.diff(time) <= 0.0):
        raise ValueError(f"{path} times are not strictly increasing")
    return {
        "time": time,
        "front": front,
        "temperature": temperature,
        "pressure_min": pressure_min,
        "pressure_max": pressure_max,
        "source": path,
    }


def local_burn_rate(time, front, window):
    rate = np.full(time.shape, np.nan)
    half_window = 0.5 * window
    for i, center in enumerate(time):
        mask = np.abs(time - center) <= half_window
        if mask.sum() >= 3:
            rate[i] = -np.polyfit(time[mask] - center, front[mask], 1)[0]
    return rate


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("driven_csv", type=Path)
    parser.add_argument("control_csv", type=Path)
    parser.add_argument("--output", type=Path, required=True, help="PNG output path")
    parser.add_argument("--trim", type=float, default=0.1, help="analysis start time in seconds")
    parser.add_argument("--rate-window", type=float, default=0.01, help="centered fit width in seconds")
    parser.add_argument("--finest-dx", type=float, default=3.125e-6, help="finest cell width in meters")
    args = parser.parse_args()

    driven = load_history(args.driven_csv)
    control = load_history(args.control_csv)
    if driven["time"].shape != control["time"].shape or not np.allclose(
        driven["time"], control["time"], rtol=0.0, atol=1.0e-10
    ):
        raise ValueError("Driven/control snapshots must have matching times")

    time = driven["time"]
    front_diff = driven["front"] - control["front"]
    driven_rate = local_burn_rate(time, driven["front"], args.rate_window)
    control_rate = local_burn_rate(time, control["front"], args.rate_window)
    use = time >= args.trim
    if not np.any(use):
        raise ValueError("Trim starts after the end of the available records")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    ax = axes[0, 0]
    ax.plot(time * 1.0e3, driven["front"] * 1.0e6, label="Driven")
    ax.plot(time * 1.0e3, control["front"] * 1.0e6, label="Zero-drive control")
    ax.axvline(args.trim * 1.0e3, color="0.35", linestyle="--", label="Trim")
    ax.set(title="AP/gas front position", xlabel="Time (ms)", ylabel="eta = 0.5 front x (um)")
    ax.legend()

    ax = axes[0, 1]
    ax.plot(time * 1.0e3, front_diff * 1.0e9)
    ax.axvline(args.trim * 1.0e3, color="0.5", linestyle=":")
    ax.set(title="Driven minus control front", xlabel="Time (ms)", ylabel="Position difference (nm)")
    difference_limit = max(float(np.max(np.abs(front_diff))) * 1.25, 1.0)
    ax.set_ylim(-difference_limit, difference_limit)
    ax.text(
        0.03,
        0.95,
        f"max |difference| = {np.max(np.abs(front_diff)) * 1e9:.3g} nm\n"
        f"finest cell = {args.finest_dx * 1e6:.3g} um",
        transform=ax.transAxes,
        va="top",
    )

    ax = axes[1, 0]
    ax.plot(time * 1.0e3, driven_rate * 1.0e3, label="Driven")
    ax.plot(time * 1.0e3, control_rate * 1.0e3, label="Zero-drive control")
    ax.axvline(args.trim * 1.0e3, color="0.35", linestyle="--")
    ax.set(title=f"Local regression rate ({args.rate_window * 1e3:g} ms linear fit)",
           xlabel="Time (ms)", ylabel="Burn rate (mm/s)")
    ax.legend()

    ax = axes[1, 1]
    p0 = 0.5 * (driven["pressure_min"] + driven["pressure_max"])
    ax.plot(time * 1.0e3, p0 / 1.0e6, label="Driven pressure envelope midpoint")
    ax.fill_between(
        time * 1.0e3,
        driven["pressure_min"] / 1.0e6,
        driven["pressure_max"] / 1.0e6,
        alpha=0.2,
        label="Domain pressure range",
    )
    ax.plot(time * 1.0e3, control["temperature"] * 0.0 + 0.5 * (
        control["pressure_min"] + control["pressure_max"]
    ) / 1.0e6, color="tab:green", alpha=0.6, label="Control pressure midpoint")
    ax.axvline(args.trim * 1.0e3, color="0.35", linestyle="--")
    ax.set(title="Thermodynamic pressure", xlabel="Time (ms)", ylabel="Pressure (MPa)")
    ax.legend(fontsize="small")

    fig.suptitle("1D Beckstead 130 Hz pilot: exploratory forced/control comparison")
    fig.savefig(args.output, dpi=160)
    plt.close(fig)

    comparison_csv = args.output.with_name(f"{args.output.stem}_comparison.csv")
    columns = (
        "time_s",
        "driven_front_x_m",
        "control_front_x_m",
        "front_difference_m",
        "driven_burn_rate_m_s",
        "control_burn_rate_m_s",
        "driven_front_temperature_K",
        "control_front_temperature_K",
        "driven_pressure_min_Pa",
        "driven_pressure_max_Pa",
        "control_pressure_min_Pa",
        "control_pressure_max_Pa",
    )
    with comparison_csv.open("w", newline="") as output:
        writer = csv.DictWriter(output, fieldnames=columns)
        writer.writeheader()
        for i in range(len(time)):
            writer.writerow(
                {
                    "time_s": time[i],
                    "driven_front_x_m": driven["front"][i],
                    "control_front_x_m": control["front"][i],
                    "front_difference_m": front_diff[i],
                    "driven_burn_rate_m_s": driven_rate[i],
                    "control_burn_rate_m_s": control_rate[i],
                    "driven_front_temperature_K": driven["temperature"][i],
                    "control_front_temperature_K": control["temperature"][i],
                    "driven_pressure_min_Pa": driven["pressure_min"][i],
                    "driven_pressure_max_Pa": driven["pressure_max"][i],
                    "control_pressure_min_Pa": control["pressure_min"][i],
                    "control_pressure_max_Pa": control["pressure_max"][i],
                }
            )

    burn_rate_mean = float(np.nanmean(driven_rate[use]))
    control_rate_mean = float(np.nanmean(control_rate[use]))
    resolved_difference = float(np.max(np.abs(front_diff[use])))
    print(f"Snapshots: {len(time)}; time: {time[0]:.6g} to {time[-1]:.6g} s")
    print(f"Trim: {args.trim:.6g} s; post-trim span: {time[-1] - args.trim:.6g} s")
    print(f"Mean post-trim fitted burn rate (driven/control): {burn_rate_mean:.6g} / {control_rate_mean:.6g} m/s")
    print(f"Max post-trim |front difference|: {resolved_difference:.6g} m ({resolved_difference / args.finest_dx:.4g} finest cells)")
    print(f"Plot: {args.output}")
    print(f"Comparison CSV: {comparison_csv}")


if __name__ == "__main__":
    main()
