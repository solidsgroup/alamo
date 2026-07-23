#!/usr/bin/env python3
"""Diagnose the rod/tube radial-stress feature at the phi interface.

The primary quantities are native coordinate-face normal tractions on the two
symmetry axes.  Nodal face averages and polar projections are secondary
visualization quantities.  This script does not advance the simulation.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


STRESS_VALIDATE = Path("/home/jackplum/Projects/chamberutils/stress_validate")
sys.path.insert(0, str(STRESS_VALIDATE))
import compare_alamo as ca  # noqa: E402


CENTER = np.array([0.0877, 0.0877])
STATE_FIELDS = (
    "phi",
    "eta",
    "temp",
    "model_mu",
    "model_kappa",
    "model_F0xx",
    "model_F0xy",
    "model_F0yx",
    "model_F0yy",
    "rhs_x",
    "rhs_y",
)
MECHANICS_FIELDS = (
    "disp_x",
    "disp_y",
    "stress_xx",
    "stress_xy",
    "stress_yx",
    "stress_yy",
)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_level(plotdir: Path, level: int):
    header = ca.parse_main_header(plotdir)
    names = header["names"]
    boxes = list(ca.read_level_fabs(plotdir, level, header["ncomp"]))
    ni = max(hi[0] for _, hi, _ in boxes) + 1
    nj = max(hi[1] for _, hi, _ in boxes) + 1
    sums = np.zeros((header["ncomp"], nj, ni))
    count = np.zeros((nj, ni), dtype=np.int32)
    first = np.full_like(sums, np.nan)
    duplicate_spread = 0.0
    for (ilo, jlo), (ihi, jhi), data in boxes:
        rows = slice(jlo, jhi + 1)
        cols = slice(ilo, ihi + 1)
        occupied = count[rows, cols] > 0
        if occupied.any():
            previous = first[:, rows, cols]
            duplicate_spread = max(
                duplicate_spread,
                float(np.nanmax(np.abs(previous[:, occupied] - data[:, occupied]))),
            )
        empty = ~occupied
        if empty.any():
            target = first[:, rows, cols]
            target[:, empty] = data[:, empty]
        sums[:, rows, cols] += data
        count[rows, cols] += 1
    fields = np.full_like(sums, np.nan)
    covered = count > 0
    fields[:, covered] = sums[:, covered] / count[covered]
    return (
        header,
        {name: fields[index] for index, name in enumerate(names)},
        covered,
        duplicate_spread,
    )


class FrozenField:
    """Exact 2-D production face-stress reconstruction for one AMR level."""

    def __init__(self, fields, covered, dx):
        self.f = fields
        self.covered = covered
        self.dx = np.asarray(dx)
        self.nj, self.ni = covered.shape
        self.xhi = self.ni - 1
        self.ylo = 0

    def _reflect(self, i, j):
        signs = np.ones(2)
        if i > self.xhi:
            i = 2 * self.xhi - i
            signs[0] *= -1.0
        if j < self.ylo:
            j = 2 * self.ylo - j
            signs[1] *= -1.0
        if i < 0 or i >= self.ni or j < 0 or j >= self.nj:
            raise IndexError
        if not self.covered[j, i]:
            raise KeyError
        return i, j, signs

    def scalar(self, name, i, j):
        i, j, _ = self._reflect(i, j)
        return self.f[name][j, i]

    def vector(self, i, j):
        i, j, signs = self._reflect(i, j)
        value = np.array([self.f["disp_x"][j, i], self.f["disp_y"][j, i]])
        return value * signs

    def model(self, i, j):
        i, j, _ = self._reflect(i, j)
        return (
            self.f["model_mu"][j, i],
            self.f["model_kappa"][j, i],
            np.array(
                [
                    [self.f["model_F0xx"][j, i], self.f["model_F0xy"][j, i]],
                    [self.f["model_F0yx"][j, i], self.f["model_F0yy"][j, i]],
                ]
            ),
        )

    def face_gradient(self, i, j, face):
        ip, jp = i + (face == 0), j + (face == 1)
        gradient = np.zeros((2, 2))
        gradient[:, face] = (
            self.vector(ip, jp) - self.vector(i, j)
        ) / self.dx[face]
        if face != 0:
            gradient[:, 0] = (
                self.vector(i + 1, j)
                - self.vector(i - 1, j)
                + self.vector(ip + 1, jp)
                - self.vector(ip - 1, jp)
            ) / (4.0 * self.dx[0])
        if face != 1:
            gradient[:, 1] = (
                self.vector(i, j + 1)
                - self.vector(i, j - 1)
                + self.vector(ip, jp + 1)
                - self.vector(ip, jp - 1)
            ) / (4.0 * self.dx[1])
        return gradient

    def face_model(self, i, j, face):
        ip, jp = i + (face == 0), j + (face == 1)
        low = self.model(i, j)
        high = self.model(ip, jp)
        return tuple(0.5 * (a + b) for a, b in zip(low, high))

    def face_stress(self, i, j, face):
        deformation = np.eye(2) + self.face_gradient(i, j, face)
        return neo_predeformed_dw(deformation, *self.face_model(i, j, face))

    def reconstructed_stress(self, i, j):
        output = np.zeros((2, 2))
        for face in range(2):
            im, jm = i - (face == 0), j - (face == 1)
            output[:, face] = 0.5 * (
                self.face_stress(im, jm, face)[:, face]
                + self.face_stress(i, j, face)[:, face]
            )
        return output

    def face_divergence(self, i, j):
        output = np.zeros(2)
        for face in range(2):
            im, jm = i - (face == 0), j - (face == 1)
            output += (
                self.face_stress(i, j, face)[:, face]
                - self.face_stress(im, jm, face)[:, face]
            ) / self.dx[face]
        return output

    def has_stencil(self, i, j, radius=2):
        try:
            for di in range(-radius, radius + 1):
                for dj in range(-radius, radius + 1):
                    self._reflect(i + di, j + dj)
            return True
        except (IndexError, KeyError):
            return False


def neo_predeformed_dw(deformation, mu, kappa, f0):
    relative = deformation @ np.linalg.inv(f0)
    determinant = np.linalg.det(relative)
    determinant_23 = abs(determinant) ** (2.0 / 3.0)
    trace = np.sum(relative * relative) + 1.0
    inverse_transpose = np.linalg.inv(relative).T
    return (
        mu
        * (
            relative / determinant_23
            - trace * inverse_transpose / (3.0 * determinant_23)
        )
        + kappa
        * (determinant - 1.0)
        * determinant
        * inverse_transpose
    )


def radial_component(stress, x, y):
    direction = np.array([x - CENTER[0], y - CENTER[1]])
    norm = np.linalg.norm(direction)
    if norm == 0.0:
        return np.nan
    direction /= norm
    return float(direction @ stress @ direction)


def alternating_metrics(coordinate, values):
    coordinate = np.asarray(coordinate, dtype=float)
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(coordinate) & np.isfinite(values)
    coordinate = coordinate[finite]
    values = values[finite]
    if values.size < 5:
        return {
            "samples": int(values.size),
            "nyquist_amplitude_pa": None,
            "detrended_rms_pa": None,
            "detrended_linf_pa": None,
            "monotone_overshoot_pa": None,
        }
    order = np.argsort(coordinate)
    coordinate = coordinate[order]
    values = values[order]
    scaled = coordinate - coordinate.mean()
    scale = np.max(np.abs(scaled))
    if scale > 0.0:
        scaled /= scale
    degree = min(3, values.size - 2)
    design = np.column_stack([scaled**power for power in range(degree + 1)])
    alternating = (-1.0) ** np.arange(values.size)
    design = np.column_stack([design, alternating])
    coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
    smooth = design[:, :-1] @ coefficients[:-1]
    residual = values - smooth
    increasing = values[-1] >= values[0]
    oriented = values if increasing else -values
    monotone_envelope = np.maximum.accumulate(oriented)
    monotone_overshoot = np.max(monotone_envelope - oriented)
    return {
        "samples": int(values.size),
        "nyquist_amplitude_pa": float(abs(coefficients[-1])),
        "detrended_rms_pa": float(np.sqrt(np.mean(residual**2))),
        "detrended_linf_pa": float(np.max(np.abs(residual))),
        "monotone_overshoot_pa": float(monotone_overshoot),
        "coordinate_min_m": float(coordinate.min()),
        "coordinate_max_m": float(coordinate.max()),
    }


def hash_dense_field(values, covered):
    digest = hashlib.sha256()
    digest.update(np.ascontiguousarray(covered).view(np.uint8))
    digest.update(np.ascontiguousarray(values[covered]).view(np.uint8))
    return digest.hexdigest()


def axis_records(header, fields, covered, level, axis):
    ff = FrozenField(fields, covered, header["dx"][level])
    dx = header["dx"][level][0]
    records = []
    if axis == "bottom":
        points = ((i, 0, 0) for i in range(ff.ni - 1))
    elif axis == "right":
        points = ((ff.xhi, j, 1) for j in range(ff.nj - 1))
    else:
        raise ValueError(axis)
    for i, j, face in points:
        if not ff.has_stencil(i, j):
            continue
        ip, jp = i + (face == 0), j + (face == 1)
        phi_face = 0.5 * (
            ff.scalar("phi", i, j) + ff.scalar("phi", ip, jp)
        )
        eta_face = 0.5 * (
            ff.scalar("eta", i, j) + ff.scalar("eta", ip, jp)
        )
        eta_delta = abs(
            ff.scalar("eta", ip, jp) - ff.scalar("eta", i, j)
        )
        if not (0.1 <= phi_face <= 0.9):
            continue
        if eta_delta > 1.0e-8:
            continue
        x = header["prob_lo"][0] + (i + 0.5 * (face == 0)) * dx
        y = header["prob_lo"][1] + (j + 0.5 * (face == 1)) * dx
        native = ff.face_stress(i, j, face)[face, face]
        reconstructed = ff.reconstructed_stress(i, j)
        saved = np.array(
            [
                [fields["stress_xx"][j, i], fields["stress_xy"][j, i]],
                [fields["stress_yx"][j, i], fields["stress_yy"][j, i]],
            ]
        )
        records.append(
            {
                "axis": axis,
                "level": level,
                "i": i,
                "j": j,
                "coordinate_m": float(np.hypot(x - CENTER[0], y - CENTER[1])),
                "x_m": float(x),
                "y_m": float(y),
                "phi_face": float(phi_face),
                "eta_face": float(eta_face),
                "eta_face_delta": float(eta_delta),
                "saved_prr_pa": radial_component(saved, x, y),
                "face_average_prr_pa": radial_component(
                    reconstructed, x, y
                ),
                "native_face_normal_pa": float(native),
            }
        )
    return records


def choose_axis_level(level_data, axis):
    candidates = []
    for level, (header, fields, covered, _) in enumerate(level_data):
        records = axis_records(header, fields, covered, level, axis)
        if records:
            candidates.append((len(records), level, records))
    if not candidates:
        raise RuntimeError(f"No valid {axis} phi-band face records")
    sufficient = [item for item in candidates if item[0] >= 8]
    pool = sufficient if sufficient else candidates
    _, level, records = max(pool, key=lambda item: item[1])
    return level, records


def polar_band_metrics(header, fields, covered, level):
    ff = FrozenField(fields, covered, header["dx"][level])
    dx = header["dx"][level][0]
    saved_values = []
    reconstructed_values = []
    residuals = []
    eta_values = []
    eta_deltas = []
    for j in range(ff.nj):
        for i in range(ff.ni):
            if not covered[j, i] or not ff.has_stencil(i, j):
                continue
            phi = fields["phi"][j, i]
            if not (0.1 <= phi <= 0.9):
                continue
            eta_neighbors = [
                ff.scalar("eta", i + 1, j),
                ff.scalar("eta", i - 1, j),
                ff.scalar("eta", i, j + 1),
                ff.scalar("eta", i, j - 1),
            ]
            eta = fields["eta"][j, i]
            eta_delta = max(abs(value - eta) for value in eta_neighbors)
            if eta_delta > 1.0e-8:
                continue
            x = header["prob_lo"][0] + i * dx
            y = header["prob_lo"][1] + j * dx
            saved = np.array(
                [
                    [fields["stress_xx"][j, i], fields["stress_xy"][j, i]],
                    [fields["stress_yx"][j, i], fields["stress_yy"][j, i]],
                ]
            )
            reconstructed = ff.reconstructed_stress(i, j)
            try:
                divergence = ff.face_divergence(i, j)
            except (IndexError, KeyError, np.linalg.LinAlgError):
                continue
            rhs = np.array([fields["rhs_x"][j, i], fields["rhs_y"][j, i]])
            saved_values.append(radial_component(saved, x, y))
            reconstructed_values.append(
                radial_component(reconstructed, x, y)
            )
            residuals.append(divergence - rhs)
            eta_values.append(eta)
            eta_deltas.append(eta_delta)
    saved_values = np.asarray(saved_values)
    reconstructed_values = np.asarray(reconstructed_values)
    residuals = np.asarray(residuals)
    return {
        "samples": int(saved_values.size),
        "saved_vs_face_average_prr_linf_pa": float(
            np.max(np.abs(saved_values - reconstructed_values))
        ),
        "saved_vs_face_average_prr_rms_pa": float(
            np.sqrt(np.mean((saved_values - reconstructed_values) ** 2))
        ),
        "face_residual_linf_pa_per_m": float(np.max(np.abs(residuals))),
        "face_residual_l2_pa_per_m": float(
            np.sqrt(np.mean(residuals**2))
        ),
        "eta_min": float(np.min(eta_values)),
        "eta_max": float(np.max(eta_values)),
        "eta_neighbor_delta_max": float(np.max(eta_deltas)),
    }


def analyze(plotdir: Path, outdir: Path):
    outdir.mkdir(parents=True, exist_ok=True)
    main_header = ca.parse_main_header(plotdir)
    level_data = [
        load_level(plotdir, level)
        for level in range(main_header["finest"] + 1)
    ]

    chosen = {}
    all_records = []
    for axis in ("bottom", "right"):
        level, records = choose_axis_level(level_data, axis)
        chosen[axis] = level
        all_records.extend(records)

    columns = [
        "axis",
        "level",
        "i",
        "j",
        "coordinate_m",
        "x_m",
        "y_m",
        "phi_face",
        "eta_face",
        "eta_face_delta",
        "saved_prr_pa",
        "face_average_prr_pa",
        "native_face_normal_pa",
    ]
    with (outdir / "profiles.csv").open("w") as stream:
        stream.write(",".join(columns) + "\n")
        for record in all_records:
            stream.write(",".join(str(record[name]) for name in columns) + "\n")

    axis_metrics = {}
    for axis in ("bottom", "right"):
        records = [record for record in all_records if record["axis"] == axis]
        coordinate = [record["coordinate_m"] for record in records]
        axis_metrics[axis] = {
            "level": chosen[axis],
            "saved": alternating_metrics(
                coordinate, [record["saved_prr_pa"] for record in records]
            ),
            "face_average": alternating_metrics(
                coordinate,
                [record["face_average_prr_pa"] for record in records],
            ),
            "native_face": alternating_metrics(
                coordinate,
                [record["native_face_normal_pa"] for record in records],
            ),
            "saved_vs_face_average_linf_pa": float(
                max(
                    abs(
                        record["saved_prr_pa"]
                        - record["face_average_prr_pa"]
                    )
                    for record in records
                )
            ),
            "face_average_vs_native_linf_pa": float(
                max(
                    abs(
                        record["face_average_prr_pa"]
                        - record["native_face_normal_pa"]
                    )
                    for record in records
                )
            ),
            "phi_face_min": float(min(record["phi_face"] for record in records)),
            "phi_face_max": float(max(record["phi_face"] for record in records)),
            "eta_face_min": float(min(record["eta_face"] for record in records)),
            "eta_face_max": float(max(record["eta_face"] for record in records)),
            "eta_face_delta_max": float(
                max(record["eta_face_delta"] for record in records)
            ),
        }

    state_hashes = {}
    mechanics_hashes = {}
    level_hashes = {}
    for level, (_, fields, covered, duplicate_spread) in enumerate(level_data):
        level_hashes[str(level)] = {
            "covered": hashlib.sha256(
                np.ascontiguousarray(covered).view(np.uint8)
            ).hexdigest(),
            "duplicate_spread": duplicate_spread,
        }
        for name in STATE_FIELDS:
            state_hashes[f"L{level}:{name}"] = hash_dense_field(
                fields[name], covered
            )
        for name in MECHANICS_FIELDS:
            mechanics_hashes[f"L{level}:{name}"] = hash_dense_field(
                fields[name], covered
            )

    finest_header, finest_fields, finest_covered, _ = level_data[-1]
    metrics = {
        "plotdir": str(plotdir.resolve()),
        "header_sha256": sha256_file(plotdir / "Header"),
        "time": main_header["time"],
        "finest_level": main_header["finest"],
        "dx": main_header["dx"],
        "selected_axis_levels": chosen,
        "level_hashes": level_hashes,
        "state_hashes": state_hashes,
        "mechanics_hashes": mechanics_hashes,
        "axis": axis_metrics,
        "finest_phi_band": polar_band_metrics(
            finest_header, finest_fields, finest_covered, main_header["finest"]
        ),
    }
    (outdir / "metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")

    figure, axes = plt.subplots(
        2, 2, figsize=(12, 8), constrained_layout=True
    )
    for row, axis in enumerate(("bottom", "right")):
        records = [record for record in all_records if record["axis"] == axis]
        records.sort(key=lambda record: record["coordinate_m"])
        radius = np.array([record["coordinate_m"] for record in records])
        axes[row, 0].plot(
            radius,
            np.array([record["saved_prr_pa"] for record in records]) / 1.0e6,
            "o-",
            label="saved polar",
        )
        axes[row, 0].plot(
            radius,
            np.array(
                [record["face_average_prr_pa"] for record in records]
            )
            / 1.0e6,
            "s-",
            label="face-average node",
        )
        axes[row, 0].plot(
            radius,
            np.array(
                [record["native_face_normal_pa"] for record in records]
            )
            / 1.0e6,
            "^-",
            label="native face",
        )
        axes[row, 0].set_ylabel("$P_{rr}$ [MPa]")
        axes[row, 0].set_title(
            f"{axis} axis, level {chosen[axis]}"
        )
        axes[row, 0].legend()
        axes[row, 1].plot(
            radius,
            [record["phi_face"] for record in records],
            "k-",
            label="phi",
        )
        axes[row, 1].plot(
            radius,
            [record["eta_face"] for record in records],
            "r--",
            label="eta",
        )
        axes[row, 1].set_ylim(-0.05, 1.05)
        axes[row, 1].set_title("target-band fields")
        axes[row, 1].legend()
        for column in range(2):
            axes[row, column].set_xlabel("radius [m]")
            axes[row, column].grid(alpha=0.25)
    figure.savefig(outdir / "profiles.png", dpi=180)
    plt.close(figure)
    print(json.dumps(metrics, indent=2))
    return metrics


def self_test():
    coordinate = np.linspace(0.075, 0.085, 32)
    scaled = (coordinate - coordinate.mean()) / np.ptp(coordinate)
    smooth = 2.0e6 + 4.0e5 * scaled - 1.0e5 * scaled**2 + 5.0e4 * scaled**3
    alternating = (-1.0) ** np.arange(coordinate.size)
    native_amplitude = 1234.5
    nodal_amplitude = 4321.0
    clean = alternating_metrics(coordinate, smooth)
    native = alternating_metrics(
        coordinate, smooth + native_amplitude * alternating
    )
    nodal = alternating_metrics(
        coordinate, smooth + nodal_amplitude * alternating
    )
    result = {
        "clean": clean,
        "native_injection": native,
        "nodal_injection": nodal,
    }
    clean_limit = 1.0e-8 * min(native_amplitude, nodal_amplitude)
    if clean["nyquist_amplitude_pa"] > clean_limit:
        raise AssertionError(result)
    if abs(native["nyquist_amplitude_pa"] - native_amplitude) > 0.02 * native_amplitude:
        raise AssertionError(result)
    if abs(nodal["nyquist_amplitude_pa"] - nodal_amplitude) > 0.02 * nodal_amplitude:
        raise AssertionError(result)
    print(json.dumps(result, indent=2))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--node", type=Path)
    parser.add_argument("--out", type=Path)
    parser.add_argument("--self-test", action="store_true")
    parser.add_argument(
        "--write-visit-comparison",
        action="store_true",
        help="Accepted for the plan oracle; profiles.csv is always written.",
    )
    args = parser.parse_args()
    if args.self_test:
        self_test()
        return
    if args.node is None or args.out is None:
        parser.error("--node and --out are required unless --self-test is used")
    analyze(args.node, args.out)


if __name__ == "__main__":
    main()
