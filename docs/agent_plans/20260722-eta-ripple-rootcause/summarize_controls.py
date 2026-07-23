#!/usr/bin/env python3
"""Aggregate the eta-ripple planar and radial control matrices."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib.pyplot as plt


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
PLANAR = RESULTS / "planar-analysis"
RADIAL = RESULTS / "radial-analysis"


def read_metric(path: Path) -> dict:
    return json.loads((path / "metrics.json").read_text())


def planar_row(name: str, label: str, factor: float | None = None) -> dict:
    metric = read_metric(PLANAR / name)
    interface = metric["interface_metrics"][0]
    return {
        "label": label,
        "factor": factor,
        "dx_m": metric["dx_m"][0],
        "pressure_pa": metric["pressure_fit"]["pressure_pa"],
        "raw_nyquist_pa": interface["raw"]["nyquist_projection_pa"],
        "corrected_nyquist_pa": interface["pressure_corrected"]["nyquist_projection_pa"],
        "corrected_rms_pa": interface["pressure_corrected"]["detrended_rms_pa"],
        "corrected_linf_pa": interface["pressure_corrected"]["detrended_linf_pa"],
        "corrected_constant_linf_pa": metric["corrected_constant"]["linf_about_median_pa"],
        "residual_linf_pa_per_m": metric["profile_residual_linf_pa_per_m"],
        "slope_relative_error": metric["raw_vs_eta_fit"]["relative_slope_error"],
    }


def radial_row(name: str, label: str) -> dict:
    metric = read_metric(RADIAL / name)
    interface = metric["profiles"]["bottom_x"]["interfaces"][0]
    return {
        "label": label,
        "dx_m": metric["dx_m"][0],
        "pressure_pa": metric["pressure_fit"]["pressure_pa"],
        "raw_nyquist_pa": interface["raw"]["nyquist_projection_pa"],
        "corrected_nyquist_pa": interface["pressure_corrected"]["nyquist_projection_pa"],
        "corrected_rms_pa": interface["pressure_corrected"]["detrended_rms_pa"],
        "corrected_linf_pa": interface["pressure_corrected"]["detrended_linf_pa"],
    }


def orders(rows: list[dict], key: str) -> list[float]:
    answer = []
    for coarse, fine in zip(rows, rows[1:]):
        answer.append(
            math.log(coarse[key] / fine[key])
            / math.log(coarse["dx_m"] / fine["dx_m"])
        )
    return answer


def main() -> None:
    pressure = [
        planar_row("x_p0_n128_w0.00113800_s0", "0", 0.0),
        planar_row("x_p431638.6252462386_n128_w0.00113800_s0", "0.1", 0.1),
        planar_row("x_p2158193.126231193_n128_w0.00113800_s0", "0.5", 0.5),
        planar_row("x_p4316386.252462386_n128_w0.00113800_s0", "1", 1.0),
        planar_row("x_p8632772.504924772_n128_w0.00113800_s0", "2", 2.0),
    ]
    phase = {}
    for orientation in ("x", "y"):
        phase[orientation] = [
            planar_row(
                f"{orientation}_p4316386.252462386_n128_w0.00113800_s{shift}",
                shift,
            )
            for shift in ("0", "0.25", "0.5", "0.75")
        ]
    phase["diag"] = [
        planar_row("diag_base", "0"),
        planar_row("diag_s0.25", "0.25"),
        planar_row("diag_s0.5", "0.5"),
        planar_row("diag_s0.75", "0.75"),
    ]
    planar_fixed_width = [
        planar_row("x_p4316386.252462386_n64_w0.00113800_s0", "64"),
        planar_row("x_p4316386.252462386_n128_w0.00113800_s0", "128"),
        planar_row("x_p4316386.252462386_n256_w0.00113800_s0", "256"),
    ]
    diagonal_fixed_width = [
        planar_row("diag_fixedw_n64", "64"),
        planar_row("diag_base", "128"),
        planar_row("diag_fixedw_n256", "256"),
    ]
    radial_fixed_width = [
        radial_row("fixedw_n64", "64"),
        radial_row("fixedw_n128", "128"),
        radial_row("fixedw_n256", "256"),
    ]
    radial_fixed_cells = [
        radial_row("fixedcells_n64", "64"),
        radial_row("fixedw_n128", "128"),
        radial_row("fixedcells_n256", "256"),
    ]
    summary = {
        "pressure_series": pressure,
        "phase_series": phase,
        "fixed_physical_width": {
            "aligned_x": planar_fixed_width,
            "diagonal": diagonal_fixed_width,
            "radial_signed_distance": radial_fixed_width,
            "diagonal_corrected_nyquist_orders": orders(
                diagonal_fixed_width, "corrected_nyquist_pa"
            ),
            "radial_corrected_nyquist_orders": orders(
                radial_fixed_width, "corrected_nyquist_pa"
            ),
        },
        "fixed_width_in_cells": {
            "radial_signed_distance": radial_fixed_cells,
            "radial_corrected_nyquist_orders": orders(
                radial_fixed_cells, "corrected_nyquist_pa"
            ),
        },
    }
    (RESULTS / "control_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    with (RESULTS / "control_summary.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            ("study", "label", "dx_m", "pressure_pa", "raw_nyquist_pa", "corrected_nyquist_pa")
        )
        for study, rows in (
            ("pressure", pressure),
            ("diag_fixed_width", diagonal_fixed_width),
            ("radial_fixed_width", radial_fixed_width),
            ("radial_fixed_cells", radial_fixed_cells),
        ):
            for row in rows:
                writer.writerow(
                    (
                        study,
                        row["label"],
                        row["dx_m"],
                        row["pressure_pa"],
                        row["raw_nyquist_pa"],
                        row["corrected_nyquist_pa"],
                    )
                )

    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2), constrained_layout=True)
    axes[0].plot(
        [row["factor"] for row in pressure],
        [row["raw_nyquist_pa"] / 1e3 for row in pressure],
        "o-",
        label="raw",
    )
    axes[0].plot(
        [row["factor"] for row in pressure],
        [row["corrected_nyquist_pa"] / 1e3 for row in pressure],
        "o-",
        label="corrected",
    )
    axes[0].set(xlabel="pressure / production pressure", ylabel="Nyquist projection [kPa]")
    axes[0].legend()

    for orientation, rows in phase.items():
        axes[1].plot(
            [float(row["label"]) for row in rows],
            [row["corrected_nyquist_pa"] for row in rows],
            "o-",
            label=orientation,
        )
    axes[1].set_yscale("symlog", linthresh=1e-7)
    axes[1].set(xlabel="normal shift [cells]", ylabel="corrected Nyquist [Pa]")
    axes[1].legend()

    for label, rows in (
        ("diagonal", diagonal_fixed_width),
        ("radial SDF", radial_fixed_width),
    ):
        axes[2].loglog(
            [row["dx_m"] for row in rows],
            [row["corrected_nyquist_pa"] for row in rows],
            "o-",
            label=label,
        )
    axes[2].invert_xaxis()
    axes[2].set(xlabel="grid spacing h [m]", ylabel="corrected Nyquist [Pa]")
    axes[2].legend()
    for axis in axes:
        axis.grid(alpha=0.25)
    fig.savefig(RESULTS / "control_summary.png", dpi=180)


if __name__ == "__main__":
    main()
