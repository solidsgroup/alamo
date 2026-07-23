#!/usr/bin/env python3
"""Face-native diagnostics for the remaining eta-interface stress feature.

The pressure load is built from cell-centered eta while displacement, RHS, and
stress are nodal.  This script therefore requires the paired cell and node
plotfiles.  It imports the already-validated face-stress reconstruction from
the preceding diagnostic task and adds the exact pressure-flux accounting.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


ROOT = Path(__file__).resolve().parents[3]
PRIOR_SCRIPT = (
    ROOT
    / "docs/agent_plans/20260722-free-circular-casing/check_stress_ripple.py"
)
DEFAULT_NODE = ROOT / "output_ripple_no_postsolve_regrid/05002node"
DEFAULT_CELL = ROOT / "output_ripple_no_postsolve_regrid/05002cell"
DEFAULT_OUT = (
    ROOT
    / "docs/agent_plans/20260722-eta-ripple-rootcause/results/frozen-baseline"
)


def load_prior_module():
    spec = importlib.util.spec_from_file_location("prior_ripple", PRIOR_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {PRIOR_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


prior = load_prior_module()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def plot_fingerprint(plotdir: Path) -> str:
    digest = hashlib.sha256()
    for path in sorted(p for p in plotdir.rglob("*") if p.is_file()):
        digest.update(str(path.relative_to(plotdir)).encode())
        digest.update(bytes.fromhex(sha256_file(path)))
    return digest.hexdigest()


@dataclass
class CellEta:
    eta: np.ndarray
    covered: np.ndarray

    @property
    def nj(self) -> int:
        return self.eta.shape[0]

    @property
    def ni(self) -> int:
        return self.eta.shape[1]

    def _reflect(self, i: int, j: int) -> tuple[int, int]:
        # The quarter-deck eta BC is even on xhi and ylo.  The identity audit
        # itself uses strict interior rows; reflection is needed only for the
        # two symmetry-axis profiles.
        if i >= self.ni:
            i = 2 * self.ni - 1 - i
        if j < 0:
            j = -j - 1
        if i < 0 or i >= self.ni or j < 0 or j >= self.nj:
            raise IndexError
        if not self.covered[j, i]:
            raise KeyError
        return i, j

    def value(self, i: int, j: int) -> float:
        i, j = self._reflect(i, j)
        return float(self.eta[j, i])

    def xface(self, i: int, j: int) -> float:
        """Eta flux on the nodal x edge from (i,j) to (i+1,j)."""
        return 0.5 * (self.value(i, j) + self.value(i, j - 1))

    def yface(self, i: int, j: int) -> float:
        """Eta flux on the nodal y edge from (i,j) to (i,j+1)."""
        return 0.5 * (self.value(i, j) + self.value(i - 1, j))

    def gradient_on_node(self, i: int, j: int, dx: tuple[float, float]) -> np.ndarray:
        return np.array(
            [
                (self.xface(i, j) - self.xface(i - 1, j)) / dx[0],
                (self.yface(i, j) - self.yface(i, j - 1)) / dx[1],
            ]
        )

    def production_gradient(self, i: int, j: int, dx: tuple[float, float]) -> np.ndarray:
        e00 = self.value(i, j)
        em0 = self.value(i - 1, j)
        e0m = self.value(i, j - 1)
        emm = self.value(i - 1, j - 1)
        return np.array(
            [
                0.5 * (e00 - em0 + e0m - emm) / dx[0],
                0.5 * (e00 - e0m + em0 - emm) / dx[1],
            ]
        )


def load_eta(plotdir: Path, level: int) -> tuple[dict, CellEta]:
    header, fields, covered, _ = prior.load_level(plotdir, level)
    if "eta" not in fields:
        raise ValueError(f"eta is absent from {plotdir}")
    return header, CellEta(fields["eta"], covered)


def common_interior_nodes(cell: CellEta, node_covered: np.ndarray):
    nj, ni = node_covered.shape
    for j in range(1, nj - 1):
        for i in range(1, ni - 1):
            if not node_covered[j, i]:
                continue
            try:
                cell.value(i, j)
                cell.value(i - 1, j)
                cell.value(i, j - 1)
                cell.value(i - 1, j - 1)
            except (IndexError, KeyError):
                continue
            yield i, j


def fit_pressure(
    cell: CellEta,
    rhs_x: np.ndarray,
    rhs_y: np.ndarray,
    node_covered: np.ndarray,
    dx: tuple[float, float],
) -> dict:
    gradients = []
    rhs = []
    identity_error = 0.0
    for i, j in common_interior_nodes(cell, node_covered):
        production = cell.production_gradient(i, j, dx)
        flux_difference = cell.gradient_on_node(i, j, dx)
        identity_error = max(identity_error, float(np.max(np.abs(production - flux_difference))))
        if np.linalg.norm(production) <= 1.0e-12:
            continue
        gradients.append(production)
        rhs.append([rhs_x[j, i], rhs_y[j, i]])
    g = np.asarray(gradients)
    b = np.asarray(rhs)
    if not len(g):
        raise ValueError("No nonzero eta-gradient rows are available for pressure fitting")
    pressure = -float(np.sum(g * b) / np.sum(g * g))
    error = b + pressure * g
    rhs_scale = float(np.max(np.abs(b)))
    relative_linf = (
        float(np.max(np.abs(error)) / rhs_scale)
        if rhs_scale > 0.0
        else float(np.max(np.abs(error)))
    )
    return {
        "pressure_pa": pressure,
        "fit_rows": int(len(g)),
        "pressure_flux_identity_linf_per_m": identity_error,
        "rhs_fit_linf_pa_per_m": float(np.max(np.abs(error))),
        "rhs_fit_rms_pa_per_m": float(np.sqrt(np.mean(error * error))),
        "rhs_fit_relative_linf": relative_linf,
        "rhs_linf_pa_per_m": rhs_scale,
    }


def contiguous_groups(indices: np.ndarray) -> list[np.ndarray]:
    if not len(indices):
        return []
    splits = np.where(np.diff(indices) > 1)[0] + 1
    return [group for group in np.split(indices, splits) if len(group)]


def detrend_metrics(values: np.ndarray, dx: float, degree: int = 3) -> dict:
    if len(values) < degree + 3:
        raise ValueError("Interface window is too short for the locked detrending metric")
    coordinate = (np.arange(len(values)) - 0.5 * (len(values) - 1)) * dx
    coordinate /= max(np.max(np.abs(coordinate)), dx)
    fit = np.polynomial.polynomial.polyval(
        coordinate,
        np.polynomial.polynomial.polyfit(coordinate, values, degree),
    )
    residual = values - fit
    alternating = (-1.0) ** np.arange(len(values))
    return {
        "samples": int(len(values)),
        "fit_degree": degree,
        "nyquist_projection_pa": float(abs(np.dot(residual, alternating)) / len(values)),
        "detrended_rms_pa": float(np.sqrt(np.mean(residual * residual))),
        "detrended_linf_pa": float(np.max(np.abs(residual))),
        "second_difference_linf_pa": float(
            np.max(np.abs(np.diff(values, n=2))) if len(values) >= 3 else 0.0
        ),
    }


def profile_metrics(
    coordinate: np.ndarray,
    eta_face: np.ndarray,
    raw: np.ndarray,
    corrected: np.ndarray,
    dx: float,
    half_window: int = 8,
) -> tuple[list[dict], list[tuple[int, int]]]:
    transition = np.flatnonzero((eta_face > 0.1) & (eta_face < 0.9))
    groups = contiguous_groups(transition)
    metrics = []
    windows = []
    for group in groups:
        lo = max(0, int(group[0]) - half_window)
        hi = min(len(raw), int(group[-1]) + half_window + 1)
        windows.append((lo, hi))
        item = {
            "coordinate_lo_m": float(coordinate[lo]),
            "coordinate_hi_m": float(coordinate[hi - 1]),
            "eta_min": float(np.min(eta_face[lo:hi])),
            "eta_max": float(np.max(eta_face[lo:hi])),
            "raw": detrend_metrics(raw[lo:hi], dx),
            "pressure_corrected": detrend_metrics(corrected[lo:hi], dx),
        }
        metrics.append(item)
    return metrics, windows


def axis_profile(
    frozen,
    cell: CellEta,
    header: dict,
    pressure: float,
    axis: str,
) -> dict:
    dx = tuple(float(x) for x in header["dx"][header["finest"]])
    prob_lo = header["prob_lo"]
    records = []
    if axis == "bottom_x":
        j = 0
        for i in range(frozen.ni - 1):
            try:
                if not frozen.has_stencil(i, j):
                    continue
                eta_face = cell.xface(i, j)
                raw = float(frozen.face_stress(i, j, 0)[0, 0])
            except (IndexError, KeyError, np.linalg.LinAlgError):
                continue
            x = float(prob_lo[0] + (i + 0.5) * dx[0])
            records.append((x, eta_face, raw, raw + pressure * eta_face))
    elif axis == "right_y":
        i = frozen.xhi
        for j in range(frozen.nj - 1):
            try:
                if not frozen.has_stencil(i, j):
                    continue
                eta_face = cell.yface(i, j)
                raw = float(frozen.face_stress(i, j, 1)[1, 1])
            except (IndexError, KeyError, np.linalg.LinAlgError):
                continue
            y = float(prob_lo[1] + (j + 0.5) * dx[1])
            records.append((y, eta_face, raw, raw + pressure * eta_face))
    else:
        raise ValueError(axis)
    if not records:
        raise ValueError(f"No valid {axis} face profile")
    array = np.asarray(records)
    order = np.argsort(array[:, 0])
    array = array[order]
    metrics, windows = profile_metrics(
        array[:, 0], array[:, 1], array[:, 2], array[:, 3], dx[0]
    )
    return {"array": array, "metrics": metrics, "windows": windows}


def analyze(node_dir: Path, cell_dir: Path, outdir: Path, level: int | None) -> dict:
    node_header = prior.ca.parse_main_header(node_dir)
    cell_header = prior.ca.parse_main_header(cell_dir)
    if level is None:
        level = int(node_header["finest"])
    if level != int(cell_header["finest"]):
        raise ValueError("Use the common finest level for the baseline audit")
    _, node_fields, node_covered, duplicate_spread = prior.load_level(node_dir, level)
    _, cell = load_eta(cell_dir, level)
    frozen = prior.FrozenField(node_fields, node_covered, node_header["dx"][level])
    dx = tuple(float(x) for x in node_header["dx"][level])
    pressure = fit_pressure(
        cell,
        node_fields["rhs_x"],
        node_fields["rhs_y"],
        node_covered,
        dx,
    )
    profiles = {
        name: axis_profile(frozen, cell, node_header, pressure["pressure_pa"], name)
        for name in ("bottom_x", "right_y")
    }

    outdir.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 2, figsize=(13, 9), constrained_layout=True)
    for column, name in enumerate(("bottom_x", "right_y")):
        profile = profiles[name]
        a = profile["array"]
        axes[0, column].plot(a[:, 0], a[:, 2] / 1e6, label="raw elastic")
        axes[0, column].plot(a[:, 0], a[:, 3] / 1e6, label="P + p eta I")
        twin = axes[0, column].twinx()
        twin.plot(a[:, 0], a[:, 1], color="black", alpha=0.25, label="eta face")
        twin.set_ylabel("eta face")
        axes[0, column].set_title(name)
        axes[0, column].set_ylabel("native normal coordinate traction [MPa]")
        axes[0, column].legend(loc="best")
        for lo, hi in profile["windows"]:
            values = a[lo:hi, 3]
            coordinate = np.arange(len(values), dtype=float)
            fit = np.polynomial.polynomial.polyval(
                coordinate,
                np.polynomial.polynomial.polyfit(coordinate, values, 3),
            )
            axes[1, column].plot(a[lo:hi, 0], (values - fit) / 1e3)
        axes[1, column].set_title("locked cubic-fit residual")
        axes[1, column].set_ylabel("corrected residual [kPa]")
        axes[1, column].set_xlabel("coordinate [m]")
    for ax in axes.flat:
        ax.grid(alpha=0.25)
    fig.savefig(outdir / "pressure_corrected_profiles.png", dpi=180)
    plt.close(fig)

    serializable_profiles = {
        name: {"interfaces": value["metrics"]} for name, value in profiles.items()
    }
    result = {
        "node_plot": str(node_dir),
        "cell_plot": str(cell_dir),
        "node_plot_sha256": plot_fingerprint(node_dir),
        "cell_plot_sha256": plot_fingerprint(cell_dir),
        "level": level,
        "dx_m": list(dx),
        "duplicate_node_max_spread": duplicate_spread,
        "pressure_fit": pressure,
        "locked_metric": {
            "transition_eta_interval": [0.1, 0.9],
            "half_window_cells": 8,
            "fit_degree": 3,
            "nyquist_definition": "abs(dot(cubic_residual, (-1)^index))/N",
        },
        "profiles": serializable_profiles,
    }
    (outdir / "metrics.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def self_test() -> dict:
    rng = np.random.default_rng(240722)
    eta = rng.random((31, 37))
    cell = CellEta(eta, np.ones_like(eta, dtype=bool))
    dx = (0.0031, 0.0047)
    identity = 0.0
    gradients = []
    for j in range(1, eta.shape[0] - 1):
        for i in range(1, eta.shape[1] - 1):
            a = cell.production_gradient(i, j, dx)
            b = cell.gradient_on_node(i, j, dx)
            identity = max(identity, float(np.max(np.abs(a - b))))
            gradients.append(a)
    g = np.asarray(gradients)
    known_pressure = 4.31639e6
    rhs = -known_pressure * g
    recovered = -float(np.sum(g * rhs) / np.sum(g * g))

    n = 25
    coordinate = np.linspace(-1.0, 1.0, n)
    smooth = 2.0e6 + 3.0e5 * coordinate - 2.0e4 * coordinate**2
    clean = detrend_metrics(smooth, 1.0)
    injected_amplitude = 1234.5
    injected = detrend_metrics(smooth + injected_amplitude * (-1.0) ** np.arange(n), 1.0)
    result = {
        "pressure_flux_identity_linf_per_m": identity,
        "known_pressure_pa": known_pressure,
        "recovered_pressure_pa": recovered,
        "clean_nyquist_projection_pa": clean["nyquist_projection_pa"],
        "injected_nyquist_projection_pa": injected["nyquist_projection_pa"],
        "injected_amplitude_pa": injected_amplitude,
    }
    if identity > 1.0e-10:
        raise AssertionError(result)
    if abs(recovered - known_pressure) > 1.0e-8 * known_pressure:
        raise AssertionError(result)
    if clean["nyquist_projection_pa"] > 1.0e-8:
        raise AssertionError(result)
    if abs(injected["nyquist_projection_pa"] - injected_amplitude) > 0.05 * injected_amplitude:
        raise AssertionError(result)
    return result


def check_manifest(path: Path) -> dict:
    manifest = json.loads(path.read_text())
    incomplete = [case["name"] for case in manifest.get("cases", []) if case.get("status") != "complete"]
    if incomplete:
        raise SystemExit(f"Incomplete cases: {', '.join(incomplete)}")
    return {"cases": len(manifest.get("cases", [])), "complete": True}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--node", type=Path, default=DEFAULT_NODE)
    parser.add_argument("--cell", type=Path, default=DEFAULT_CELL)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--level", type=int)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument("--check-pressure-flux-identity", action="store_true")
    parser.add_argument("--case", choices=("frozen-baseline",))
    parser.add_argument("--check-repeatability", action="store_true")
    parser.add_argument("--repeat-node", type=Path)
    parser.add_argument("--repeat-cell", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--check-complete", action="store_true")
    args = parser.parse_args()

    if args.self_test:
        print(json.dumps(self_test(), indent=2))
        return
    if args.check_complete:
        if args.manifest is None:
            parser.error("--check-complete requires --manifest")
        print(json.dumps(check_manifest(args.manifest), indent=2))
        return

    result = analyze(args.node, args.cell, args.outdir, args.level)
    if args.check_pressure_flux_identity:
        scale = result["pressure_fit"]["rhs_linf_pa_per_m"]
        error = result["pressure_fit"]["pressure_flux_identity_linf_per_m"]
        if error > max(1.0e-12 * scale, 1.0e-10):
            raise SystemExit(f"Pressure-flux identity failed: {error}")
    if args.check_repeatability:
        if args.repeat_node is None or args.repeat_cell is None:
            parser.error("--check-repeatability requires --repeat-node and --repeat-cell")
        repeat_node_hash = plot_fingerprint(args.repeat_node)
        repeat_cell_hash = plot_fingerprint(args.repeat_cell)
        result["repeatability"] = {
            "node_bitwise_identical": repeat_node_hash == result["node_plot_sha256"],
            "cell_bitwise_identical": repeat_cell_hash == result["cell_plot_sha256"],
            "repeat_node_sha256": repeat_node_hash,
            "repeat_cell_sha256": repeat_cell_hash,
        }
        if not all(
            result["repeatability"][key]
            for key in ("node_bitwise_identical", "cell_bitwise_identical")
        ):
            raise SystemExit("Repeated plots are not bitwise identical")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
