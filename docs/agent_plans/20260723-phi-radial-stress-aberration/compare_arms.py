#!/usr/bin/env python3
"""Compare two phi-interface radial-stress diagnostic arms."""

from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def read_profiles(path: Path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def numeric(records, axis, name):
    selected = [record for record in records if record["axis"] == axis]
    selected.sort(key=lambda record: float(record["coordinate_m"]))
    return np.array([float(record[name]) for record in selected])


def compare(first_path: Path, second_path: Path, outdir: Path, require_state_match):
    first = json.loads(first_path.read_text())
    second = json.loads(second_path.read_text())
    first_profiles = read_profiles(first_path.parent / "profiles.csv")
    second_profiles = read_profiles(second_path.parent / "profiles.csv")
    outdir.mkdir(parents=True, exist_ok=True)

    first_hashes = first["state_hashes"]
    second_hashes = second["state_hashes"]
    common_hashes = sorted(set(first_hashes) & set(second_hashes))
    mismatched_state = {
        name: [first_hashes[name], second_hashes[name]]
        for name in common_hashes
        if first_hashes[name] != second_hashes[name]
    }
    hierarchy_match = (
        first["finest_level"] == second["finest_level"]
        and first["dx"] == second["dx"]
        and {
            level: value["covered"]
            for level, value in first["level_hashes"].items()
        }
        == {
            level: value["covered"]
            for level, value in second["level_hashes"].items()
        }
    )
    state_match = hierarchy_match and not mismatched_state
    mechanics_changed = sorted(
        name
        for name in set(first["mechanics_hashes"]) & set(second["mechanics_hashes"])
        if first["mechanics_hashes"][name] != second["mechanics_hashes"][name]
    )

    axis_comparison = {}
    for axis in ("bottom", "right"):
        axis_comparison[axis] = {}
        for representation in ("saved", "face_average", "native_face"):
            first_metric = first["axis"][axis][representation]
            second_metric = second["axis"][axis][representation]
            first_nyquist = first_metric["nyquist_amplitude_pa"]
            second_nyquist = second_metric["nyquist_amplitude_pa"]
            reduction = None
            if first_nyquist not in (None, 0.0) and second_nyquist is not None:
                reduction = 1.0 - second_nyquist / first_nyquist
            axis_comparison[axis][representation] = {
                "first_nyquist_pa": first_nyquist,
                "second_nyquist_pa": second_nyquist,
                "fractional_reduction": reduction,
                "first_detrended_linf_pa": first_metric["detrended_linf_pa"],
                "second_detrended_linf_pa": second_metric["detrended_linf_pa"],
                "first_monotone_overshoot_pa": first_metric[
                    "monotone_overshoot_pa"
                ],
                "second_monotone_overshoot_pa": second_metric[
                    "monotone_overshoot_pa"
                ],
            }

    result = {
        "first": str(first_path),
        "second": str(second_path),
        "time": [first["time"], second["time"]],
        "hierarchy_match": hierarchy_match,
        "state_match": state_match,
        "mismatched_state_hashes": mismatched_state,
        "mechanics_fields_changed": mechanics_changed,
        "axis": axis_comparison,
        "finest_phi_band": {
            "first_saved_vs_face_linf_pa": first["finest_phi_band"][
                "saved_vs_face_average_prr_linf_pa"
            ],
            "second_saved_vs_face_linf_pa": second["finest_phi_band"][
                "saved_vs_face_average_prr_linf_pa"
            ],
            "first_residual_linf_pa_per_m": first["finest_phi_band"][
                "face_residual_linf_pa_per_m"
            ],
            "second_residual_linf_pa_per_m": second["finest_phi_band"][
                "face_residual_linf_pa_per_m"
            ],
        },
    }
    (outdir / "comparison.json").write_text(json.dumps(result, indent=2) + "\n")

    figure, axes = plt.subplots(
        2, 2, figsize=(13, 8), constrained_layout=True
    )
    for row, axis in enumerate(("bottom", "right")):
        for records, label, style in (
            (first_profiles, "first", "--"),
            (second_profiles, "second", "-"),
        ):
            radius = numeric(records, axis, "coordinate_m")
            for name, marker in (
                ("saved_prr_pa", "o"),
                ("face_average_prr_pa", "s"),
                ("native_face_normal_pa", "^"),
            ):
                axes[row, 0].plot(
                    radius,
                    numeric(records, axis, name) / 1.0e6,
                    linestyle=style,
                    marker=marker,
                    label=f"{label} {name}",
                )
            axes[row, 1].plot(
                radius,
                numeric(records, axis, "phi_face"),
                linestyle=style,
                label=f"{label} phi",
            )
            axes[row, 1].plot(
                radius,
                numeric(records, axis, "eta_face"),
                linestyle=style,
                label=f"{label} eta",
            )
        axes[row, 0].set_title(f"{axis} axis radial traction")
        axes[row, 0].set_ylabel("$P_{rr}$ [MPa]")
        axes[row, 1].set_title(f"{axis} axis fields")
        axes[row, 1].set_ylim(-0.05, 1.05)
        for column in range(2):
            axes[row, column].set_xlabel("radius [m]")
            axes[row, column].grid(alpha=0.25)
            axes[row, column].legend(fontsize=7)
    figure.savefig(outdir / "comparison.png", dpi=180)
    plt.close(figure)
    print(json.dumps(result, indent=2))
    if require_state_match and not state_match:
        raise SystemExit("paired state or hierarchy does not match")
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("first", type=Path)
    parser.add_argument("second", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--require-state-match", action="store_true")
    args = parser.parse_args()
    if args.out is None:
        common = Path(
            os.path.commonpath(
                [str(args.first.parent), str(args.second.parent)]
            )
        )
        args.out = common / "comparison"
    compare(args.first, args.second, args.out, args.require_state_match)


if __name__ == "__main__":
    main()
