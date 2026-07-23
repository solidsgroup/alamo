#!/usr/bin/env python3
"""Analyze a one-level aligned planar eta/pressure control."""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
ORACLE_PATH = (
    ROOT
    / "docs/agent_plans/20260722-eta-ripple-rootcause/analyze_eta_ripple.py"
)


def load_oracle():
    spec = importlib.util.spec_from_file_location("eta_oracle", ORACLE_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {ORACLE_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


oracle = load_oracle()
prior = oracle.prior


def analyze(case_dir: Path, orientation: str, outdir: Path) -> dict:
    node_dir = case_dir / "00001node"
    cell_dir = case_dir / "00001cell"
    node_header = prior.ca.parse_main_header(node_dir)
    level = int(node_header["finest"])
    if level != 0:
        raise ValueError("The planar oracle requires one AMR level")
    _, node_fields, node_covered, duplicate_spread = prior.load_level(node_dir, level)
    _, cell = oracle.load_eta(cell_dir, level)
    dx = tuple(float(value) for value in node_header["dx"][level])
    frozen = prior.FrozenField(node_fields, node_covered, dx)
    pressure_fit = oracle.fit_pressure(
        cell,
        node_fields["rhs_x"],
        node_fields["rhs_y"],
        node_covered,
        dx,
    )
    pressure = pressure_fit["pressure_pa"]

    records = []
    face_indices = []
    balance = None
    if orientation == "x":
        j = frozen.nj // 2
        for i in range(2, frozen.ni - 3):
            if not frozen.has_stencil(i, j):
                continue
            eta_face = cell.xface(i, j)
            raw = float(frozen.face_stress(i, j, 0)[0, 0])
            coordinate = node_header["prob_lo"][0] + (i + 0.5) * dx[0]
            residual = float(node_fields["res_x"][j, i])
            records.append((coordinate, eta_face, raw, raw + pressure * eta_face, residual))
            face_indices.append(i)
        first, last = face_indices[0], face_indices[-1]
        rhs_integral = float(
            np.sum(node_fields["rhs_x"][j, first + 1 : last + 1]) * dx[0]
        )
        traction_jump = records[-1][2] - records[0][2]
        balance = {
            "elastic_traction_jump_pa": traction_jump,
            "integrated_rhs_pa": rhs_integral,
            "mismatch_pa": traction_jump - rhs_integral,
        }
    elif orientation == "y":
        i = frozen.ni // 2
        for j in range(2, frozen.nj - 3):
            if not frozen.has_stencil(i, j):
                continue
            eta_face = cell.yface(i, j)
            raw = float(frozen.face_stress(i, j, 1)[1, 1])
            coordinate = node_header["prob_lo"][1] + (j + 0.5) * dx[1]
            residual = float(node_fields["res_y"][j, i])
            records.append((coordinate, eta_face, raw, raw + pressure * eta_face, residual))
            face_indices.append(j)
        first, last = face_indices[0], face_indices[-1]
        rhs_integral = float(
            np.sum(node_fields["rhs_y"][first + 1 : last + 1, i]) * dx[1]
        )
        traction_jump = records[-1][2] - records[0][2]
        balance = {
            "elastic_traction_jump_pa": traction_jump,
            "integrated_rhs_pa": rhs_integral,
            "mismatch_pa": traction_jump - rhs_integral,
        }
    elif orientation == "diag":
        normal = np.array([1.0, 1.0]) / np.sqrt(2.0)
        tangent = np.array([-1.0, 1.0]) / np.sqrt(2.0)
        center = np.array([0.02, 0.02])
        bins: dict[int, list[tuple[float, float, float, float, float]]] = {}
        for j in range(3, frozen.nj - 3):
            for i in range(3, frozen.ni - 3):
                if not frozen.has_stencil(i, j):
                    continue
                point = np.array(
                    [
                        node_header["prob_lo"][0] + i * dx[0],
                        node_header["prob_lo"][1] + j * dx[1],
                    ]
                )
                if abs(float(tangent @ (point - center))) > 0.004:
                    continue
                stress = np.column_stack(
                    (
                        0.5
                        * (
                            frozen.face_stress(i - 1, j, 0)[:, 0]
                            + frozen.face_stress(i, j, 0)[:, 0]
                        ),
                        0.5
                        * (
                            frozen.face_stress(i, j - 1, 1)[:, 1]
                            + frozen.face_stress(i, j, 1)[:, 1]
                        ),
                    )
                )
                eta_x = 0.5 * (cell.xface(i - 1, j) + cell.xface(i, j))
                eta_y = 0.5 * (cell.yface(i, j - 1) + cell.yface(i, j))
                eta_normal = normal[0] ** 2 * eta_x + normal[1] ** 2 * eta_y
                raw = float(normal @ stress @ normal)
                corrected = raw + pressure * eta_normal
                residual = float(
                    normal
                    @ np.array(
                        [node_fields["res_x"][j, i], node_fields["res_y"][j, i]]
                    )
                )
                coordinate = float(normal @ (point - center))
                bins.setdefault(i + j, []).append(
                    (coordinate, eta_normal, raw, corrected, residual)
                )
        for key in sorted(bins):
            records.append(tuple(np.mean(np.asarray(bins[key]), axis=0)))
    else:
        raise ValueError(orientation)

    profile = np.asarray(records)
    if len(profile) < 16:
        raise ValueError("Planar profile is too short")
    interfaces, windows = oracle.profile_metrics(
        profile[:, 0],
        profile[:, 1],
        profile[:, 2],
        profile[:, 3],
        dx[0] if orientation == "x" else dx[1] if orientation == "y" else dx[0] / np.sqrt(2.0),
    )
    corrected = profile[:, 3]
    constant = float(np.median(corrected))
    fit = np.polynomial.polynomial.polyfit(profile[:, 1], profile[:, 2], 1)
    slope_scale = abs(pressure)
    relative_slope_error = (
        float(abs(fit[1] + pressure) / slope_scale)
        if slope_scale > 0.0
        else float(abs(fit[1]))
    )
    result = {
        "case_dir": str(case_dir),
        "orientation": orientation,
        "traction_representation": (
            "native aligned face"
            if orientation in ("x", "y")
            else "node-colocated arithmetic average of adjacent native face columns"
        ),
        "level": level,
        "dx_m": list(dx),
        "duplicate_node_max_spread": duplicate_spread,
        "pressure_fit": pressure_fit,
        "raw_vs_eta_fit": {
            "intercept_pa": float(fit[0]),
            "slope_pa": float(fit[1]),
            "expected_slope_pa": -pressure,
            "relative_slope_error": relative_slope_error,
        },
        "corrected_constant": {
            "median_pa": constant,
            "rms_about_median_pa": float(np.sqrt(np.mean((corrected - constant) ** 2))),
            "linf_about_median_pa": float(np.max(np.abs(corrected - constant))),
            "range_pa": float(np.ptp(corrected)),
        },
        "interface_metrics": interfaces,
        "profile_residual_linf_pa_per_m": float(np.max(np.abs(profile[:, 4]))),
        "global_balance_per_unit_tangent": balance,
    }

    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "metrics.json").write_text(json.dumps(result, indent=2) + "\n")
    fig, axes = plt.subplots(3, 1, figsize=(9, 10), sharex=True, constrained_layout=True)
    axes[0].plot(profile[:, 0], profile[:, 2] / 1e6, label="raw elastic")
    axes[0].plot(profile[:, 0], profile[:, 3] / 1e6, label="P + p eta I")
    axes[0].legend()
    axes[0].set_ylabel("normal traction [MPa]")
    axes[1].plot(profile[:, 0], profile[:, 1])
    axes[1].set_ylabel("eta face")
    axes[2].plot(profile[:, 0], (profile[:, 3] - constant) / 1e3)
    axes[2].set_ylabel("corrected - median [kPa]")
    axes[2].set_xlabel("normal coordinate [m]")
    for axis in axes:
        axis.grid(alpha=0.25)
    fig.savefig(outdir / "profile.png", dpi=180)
    plt.close(fig)
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("case_dir", type=Path)
    parser.add_argument("orientation", choices=("x", "y", "diag"))
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(analyze(args.case_dir, args.orientation, args.outdir), indent=2))


if __name__ == "__main__":
    main()
