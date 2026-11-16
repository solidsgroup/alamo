#!/usr/bin/env python3
"""Fit the complex, control-corrected harmonic response from matched 1D histories."""

import argparse
import csv
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


SUMMARY_FIELDS = (
    "frequency_Hz",
    "drive_case",
    "control_case",
    "stop_time_s",
    "analyzed_cycles",
    "control_analyzed_cycles",
    "drive_mean_rate_mps",
    "control_mean_rate_mps",
    "drive_Rp_mag",
    "drive_Rp_phase_deg",
    "control_Rp_mag",
    "control_Rp_phase_deg",
    "control_corrected_Rp_real",
    "control_corrected_Rp_imag",
    "control_corrected_Rp_mag",
    "control_corrected_Rp_phase_deg",
    "control_floor_fractional_rate",
    "control_floor_percent",
    "control_phase_bound_deg",
    "drive_fit_front_residual_rms_m",
    "control_fit_front_residual_rms_m",
    "first_block_Rp_mag",
    "last_block_Rp_mag",
    "relative_magnitude_drift_pct",
    "phase_drift_deg",
    "analysis_trim_start_s",
    "analysis_trim_end_s",
    "drive_input_phase_rad",
    "control_input_phase_rad",
)

BLOCK_FIELDS = (
    "block_start_s",
    "block_end_s",
    "cycles",
    "drive_mean_rate_mps",
    "control_mean_rate_mps",
    "drive_Rp_mag",
    "drive_Rp_phase_deg",
    "control_Rp_mag",
    "control_Rp_phase_deg",
    "control_corrected_Rp_real",
    "control_corrected_Rp_imag",
    "control_corrected_Rp_mag",
    "control_corrected_Rp_phase_deg",
)


def load_front_history(path: Path) -> tuple[np.ndarray, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as source:
        reader = csv.DictReader(source)
        required = {"time_s", "front_x_m"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError(f"{path} must contain time_s and front_x_m columns")
        rows = list(reader)
    if not rows:
        raise ValueError(f"No data rows in {path}")

    try:
        time = np.asarray([float(row["time_s"]) for row in rows])
        front = np.asarray([float(row["front_x_m"]) for row in rows])
    except (TypeError, ValueError) as error:
        raise ValueError(f"{path} contains a missing or invalid front-history value") from error
    if not np.all(np.isfinite(time)) or not np.all(np.isfinite(front)):
        raise ValueError(f"{path} contains nonfinite time or front values")
    if np.any(np.diff(time) <= 0.0):
        raise ValueError(f"{path} times are not strictly increasing")
    return time, front


def fit_front_response(
    time: np.ndarray,
    front: np.ndarray,
    frequency_hz: float,
    pressure_amplitude_pa: float,
    mean_pressure_pa: float,
    forcing_start_s: float,
    forcing_phase_rad: float,
    minimum_cycles: float = 5.0,
    minimum_samples_per_cycle: float = 8.0,
) -> dict[str, float | complex]:
    duration = float(time[-1] - time[0])
    cycles = duration * frequency_hz
    if cycles < minimum_cycles:
        raise ValueError(
            f"Response window has only {cycles:.3g} cycles; "
            f"at least {minimum_cycles:g} are required"
        )
    if len(time) / cycles < minimum_samples_per_cycle:
        raise ValueError(
            f"Response window has {len(time) / cycles:.3g} samples per cycle; "
            f"at least {minimum_samples_per_cycle:g} are required"
        )

    centered_time = time - float(np.mean(time))
    theta = (
        2.0 * math.pi * frequency_hz * (time - forcing_start_s)
        + forcing_phase_rad
    )
    design = np.column_stack(
        (
            np.ones_like(centered_time),
            centered_time,
            centered_time**2,
            np.cos(theta),
            np.sin(theta),
        )
    )
    coefficients, _, rank, _ = np.linalg.lstsq(design, front, rcond=None)
    if rank != design.shape[1]:
        raise ValueError("Harmonic fit is rank deficient")

    mean_rate = -float(coefficients[1])
    if not math.isfinite(mean_rate) or mean_rate <= 0.0:
        raise ValueError(f"Fitted mean burn rate must be positive; got {mean_rate}")

    front_phasor = complex(coefficients[3], -coefficients[4])
    rate_phasor = -1j * (2.0 * math.pi * frequency_hz) * front_phasor
    pressure_fraction = pressure_amplitude_pa / mean_pressure_pa
    response = rate_phasor / (mean_rate * pressure_fraction)
    residual = front - design @ coefficients
    return {
        "mean_rate_mps": mean_rate,
        "response": complex(response),
        "front_residual_rms_m": float(np.sqrt(np.mean(residual**2))),
        "cycles": cycles,
    }


def fit_pair_response(
    time: np.ndarray,
    driven_front: np.ndarray,
    control_front: np.ndarray,
    frequency_hz: float,
    pressure_amplitude_pa: float,
    mean_pressure_pa: float,
    forcing_start_s: float,
    forcing_phase_rad: float,
    minimum_cycles: float = 5.0,
    minimum_samples_per_cycle: float = 8.0,
) -> dict[str, float | complex | dict[str, float | complex]]:
    fit_args = (
        frequency_hz,
        pressure_amplitude_pa,
        mean_pressure_pa,
        forcing_start_s,
        forcing_phase_rad,
        minimum_cycles,
        minimum_samples_per_cycle,
    )
    driven = fit_front_response(time, driven_front, *fit_args)
    control = fit_front_response(time, control_front, *fit_args)
    corrected = driven["response"] - control["response"]
    corrected_magnitude = abs(corrected)
    control_magnitude = abs(control["response"])
    floor_ratio = (
        control_magnitude / corrected_magnitude
        if corrected_magnitude > 0.0
        else math.nan
    )
    return {
        "driven": driven,
        "control": control,
        "corrected": corrected,
        "corrected_magnitude": corrected_magnitude,
        "corrected_phase_deg": math.degrees(np.angle(corrected)),
        "control_floor_fractional_rate": control_magnitude
        * pressure_amplitude_pa
        / mean_pressure_pa,
        "control_floor_percent": 100.0 * floor_ratio,
        "control_phase_bound_deg": math.degrees(
            math.asin(min(floor_ratio, 1.0))
        )
        if math.isfinite(floor_ratio)
        else math.nan,
    }


def case_name(path: Path) -> str:
    if path.parent.name == "csv_export":
        return path.parent.parent.name
    return path.parent.name


def write_rows(path: Path, fieldnames: tuple[str, ...], rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as output:
        writer = csv.DictWriter(output, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("driven_csv", type=Path)
    parser.add_argument("control_csv", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--frequency-hz", type=float, default=130.0)
    parser.add_argument("--pressure-amplitude", type=float, default=68948.0)
    parser.add_argument("--mean-pressure", type=float, default=3447379.0)
    parser.add_argument("--forcing-start", type=float, default=0.00175)
    parser.add_argument("--forcing-phase-rad", type=float, default=-math.pi / 2.0)
    parser.add_argument("--trim", type=float, default=0.1)
    parser.add_argument("--minimum-cycles", type=float, default=5.0)
    parser.add_argument("--block-cycles", type=int, default=5)
    parser.add_argument("--thermal-diffusivity", type=float, default=1.7155e-7)
    args = parser.parse_args()

    positive_values = (
        args.frequency_hz,
        args.pressure_amplitude,
        args.mean_pressure,
        args.minimum_cycles,
        args.thermal_diffusivity,
    )
    if not all(math.isfinite(value) and value > 0.0 for value in positive_values):
        parser.error("frequency, pressure, cycle count, and diffusivity must be positive and finite")
    if args.pressure_amplitude >= args.mean_pressure:
        parser.error("pressure amplitude must be smaller than mean pressure")
    if (
        not math.isfinite(args.forcing_start)
        or not math.isfinite(args.forcing_phase_rad)
        or not math.isfinite(args.trim)
        or args.trim < args.forcing_start
    ):
        parser.error("trim and forcing phase must be finite; trim must follow forcing start")
    if args.block_cycles < 3:
        parser.error("--block-cycles must be at least 3")

    driven_time, driven_front = load_front_history(args.driven_csv)
    control_time, control_front = load_front_history(args.control_csv)
    if driven_time.shape != control_time.shape or not np.allclose(
        driven_time, control_time, rtol=0.0, atol=1.0e-10
    ):
        raise ValueError("Driven/control snapshots must have matching times")

    selected = driven_time >= args.trim
    time = driven_time[selected]
    driven_front = driven_front[selected]
    control_front = control_front[selected]
    if len(time) < 6:
        raise ValueError("Fewer than six matched snapshots remain after the trim")

    fit = fit_pair_response(
        time,
        driven_front,
        control_front,
        args.frequency_hz,
        args.pressure_amplitude,
        args.mean_pressure,
        args.forcing_start,
        args.forcing_phase_rad,
        args.minimum_cycles,
    )

    period = 1.0 / args.frequency_hz
    first_cycle = math.ceil(
        (args.trim - args.forcing_start) / period - 1.0e-12
    )
    final_cycle = math.floor(
        (float(time[-1]) - args.forcing_start) / period + 1.0e-12
    )
    blocks = []
    for cycle_start in range(first_cycle, final_cycle, args.block_cycles):
        cycle_count = min(args.block_cycles, final_cycle - cycle_start)
        if cycle_count < 3:
            break
        block_start = args.forcing_start + cycle_start * period
        block_end = block_start + cycle_count * period
        block_mask = (time >= block_start) & (time < block_end)
        if not np.any(block_mask):
            continue
        block_fit = fit_pair_response(
            time[block_mask],
            driven_front[block_mask],
            control_front[block_mask],
            args.frequency_hz,
            args.pressure_amplitude,
            args.mean_pressure,
            args.forcing_start,
            args.forcing_phase_rad,
            minimum_cycles=2.5,
        )
        driven = block_fit["driven"]
        control = block_fit["control"]
        corrected = block_fit["corrected"]
        blocks.append(
            {
                "block_start_s": block_start,
                "block_end_s": block_end,
                "cycles": cycle_count,
                "drive_mean_rate_mps": driven["mean_rate_mps"],
                "control_mean_rate_mps": control["mean_rate_mps"],
                "drive_Rp_mag": abs(driven["response"]),
                "drive_Rp_phase_deg": math.degrees(np.angle(driven["response"])),
                "control_Rp_mag": abs(control["response"]),
                "control_Rp_phase_deg": math.degrees(np.angle(control["response"])),
                "control_corrected_Rp_real": corrected.real,
                "control_corrected_Rp_imag": corrected.imag,
                "control_corrected_Rp_mag": block_fit["corrected_magnitude"],
                "control_corrected_Rp_phase_deg": block_fit["corrected_phase_deg"],
            }
        )
    if not blocks:
        raise ValueError("No complete response blocks with at least three cycles")

    first_block = blocks[0]
    last_block = blocks[-1]
    first_magnitude = first_block["control_corrected_Rp_mag"]
    magnitude_drift = (
        100.0 * (last_block["control_corrected_Rp_mag"] - first_magnitude)
        / first_magnitude
        if first_magnitude > 0.0
        else math.nan
    )
    phase_drift = (
        (last_block["control_corrected_Rp_phase_deg"]
         - first_block["control_corrected_Rp_phase_deg"]
         + 180.0)
        % 360.0
        - 180.0
    )

    output_dir = args.output_dir.resolve()
    output_files = (
        output_dir / "frequency_response.csv",
        output_dir / "frequency_response_omega.csv",
        output_dir / "response_blocks.csv",
        output_dir / "response_convergence.png",
    )
    existing = [path for path in output_files if path.exists()]
    if existing:
        raise FileExistsError(
            "Refusing to overwrite existing analysis output: "
            + ", ".join(str(path) for path in existing)
        )
    output_dir.mkdir(parents=True, exist_ok=True)

    driven_fit = fit["driven"]
    control_fit = fit["control"]
    corrected = fit["corrected"]
    omega_drive = (
        2.0
        * math.pi
        * args.frequency_hz
        * args.thermal_diffusivity
        / driven_fit["mean_rate_mps"] ** 2
    )
    omega_control = (
        2.0
        * math.pi
        * args.frequency_hz
        * args.thermal_diffusivity
        / control_fit["mean_rate_mps"] ** 2
    )
    summary = {
        "frequency_Hz": args.frequency_hz,
        "drive_case": case_name(args.driven_csv),
        "control_case": case_name(args.control_csv),
        "stop_time_s": float(time[-1]),
        "analyzed_cycles": (float(time[-1]) - float(time[0])) * args.frequency_hz,
        "control_analyzed_cycles": (float(time[-1]) - float(time[0]))
        * args.frequency_hz,
        "drive_mean_rate_mps": driven_fit["mean_rate_mps"],
        "control_mean_rate_mps": control_fit["mean_rate_mps"],
        "drive_Rp_mag": abs(driven_fit["response"]),
        "drive_Rp_phase_deg": math.degrees(np.angle(driven_fit["response"])),
        "control_Rp_mag": abs(control_fit["response"]),
        "control_Rp_phase_deg": math.degrees(np.angle(control_fit["response"])),
        "control_corrected_Rp_real": corrected.real,
        "control_corrected_Rp_imag": corrected.imag,
        "control_corrected_Rp_mag": fit["corrected_magnitude"],
        "control_corrected_Rp_phase_deg": fit["corrected_phase_deg"],
        "control_floor_fractional_rate": fit["control_floor_fractional_rate"],
        "control_floor_percent": fit["control_floor_percent"],
        "control_phase_bound_deg": fit["control_phase_bound_deg"],
        "drive_fit_front_residual_rms_m": driven_fit["front_residual_rms_m"],
        "control_fit_front_residual_rms_m": control_fit["front_residual_rms_m"],
        "first_block_Rp_mag": first_magnitude,
        "last_block_Rp_mag": last_block["control_corrected_Rp_mag"],
        "relative_magnitude_drift_pct": magnitude_drift,
        "phase_drift_deg": phase_drift,
        "analysis_trim_start_s": args.trim,
        "analysis_trim_end_s": float(time[-1]),
        "drive_input_phase_rad": args.forcing_phase_rad,
        "control_input_phase_rad": args.forcing_phase_rad,
    }
    write_rows(output_files[0], SUMMARY_FIELDS, [summary])
    write_rows(
        output_files[1],
        (
            "frequency_Hz",
            "drive_case",
            "control_case",
            "drive_mean_rate_mps",
            "control_mean_rate_mps",
            "thermal_diffusivity_m2_s",
            "Omega_sim_drive",
            "Omega_sim_control",
        ),
        [
            {
                "frequency_Hz": args.frequency_hz,
                "drive_case": summary["drive_case"],
                "control_case": summary["control_case"],
                "drive_mean_rate_mps": driven_fit["mean_rate_mps"],
                "control_mean_rate_mps": control_fit["mean_rate_mps"],
                "thermal_diffusivity_m2_s": args.thermal_diffusivity,
                "Omega_sim_drive": omega_drive,
                "Omega_sim_control": omega_control,
            }
        ],
    )
    write_rows(output_files[2], BLOCK_FIELDS, blocks)

    block_times = np.asarray(
        [0.5 * (row["block_start_s"] + row["block_end_s"]) for row in blocks]
    )
    fig, axes = plt.subplots(2, 1, figsize=(9, 7), sharex=True, constrained_layout=True)
    axes[0].plot(
        block_times,
        [row["control_corrected_Rp_mag"] for row in blocks],
        "o-",
        label="Control-corrected",
    )
    axes[0].plot(
        block_times,
        [row["control_Rp_mag"] for row in blocks],
        "s--",
        label="Control floor",
    )
    axes[0].axhline(
        fit["corrected_magnitude"], color="0.35", linestyle=":", label="Full-window fit"
    )
    axes[0].set_ylabel(r"$|R_p|$")
    axes[0].legend()
    axes[1].plot(
        block_times,
        [row["control_corrected_Rp_phase_deg"] for row in blocks],
        "o-",
    )
    axes[1].axhline(fit["corrected_phase_deg"], color="0.35", linestyle=":")
    axes[1].set(xlabel="Time (s)", ylabel="Phase (degrees)")
    fig.suptitle(f"1D {args.frequency_hz:g} Hz response by {args.block_cycles}-cycle block")
    fig.savefig(output_files[3], dpi=160)
    plt.close(fig)

    print(f"Analyzed cycles: {summary['analyzed_cycles']:.3f}")
    print(
        "Control-corrected response: "
        f"|Rp|={fit['corrected_magnitude']:.6g}, "
        f"phase={fit['corrected_phase_deg']:.3f} deg"
    )
    print(
        f"Control floor: {fit['control_floor_percent']:.2f}% of corrected magnitude; "
        f"first-to-last block magnitude drift={magnitude_drift:.2f}%"
    )
    for path in output_files:
        print(f"Output: {path}")


if __name__ == "__main__":
    main()
