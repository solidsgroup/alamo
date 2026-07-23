#!/usr/bin/env python3
"""Spatially composite AMR nodal stress and reconstruct conservative face stress."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

STRESS_VALIDATE = Path("/home/jackplum/Projects/chamberutils/stress_validate")
sys.path.insert(0, str(STRESS_VALIDATE))
import compare_alamo as ca  # noqa: E402

CENTER = np.array(ca.CENTER)
C_PHI = 0
C_DISP_X = 1
C_DISP_Y = 2
C_STRESS_XX = 5
C_STRESS_XY = 6
C_STRESS_YX = 7
C_STRESS_YY = 8
C_MODEL_MU = 13
C_MODEL_KAPPA = 14
C_MODEL_F0XX = 15
C_MODEL_F0XY = 16
C_MODEL_F0YX = 17
C_MODEL_F0YY = 18
C_ETA = 19


def dense_level(plotdir: Path, level: int, hdr: dict) -> tuple[np.ndarray, list[tuple[tuple[int, int], tuple[int, int]]]]:
    """Load one sparse AMR level into a global-index dense array with NaN gaps."""
    dx, dy = hdr["dx"][level]
    nx = int(round((hdr["prob_hi"][0] - hdr["prob_lo"][0]) / dx)) + 1
    ny = int(round((hdr["prob_hi"][1] - hdr["prob_lo"][1]) / dy)) + 1
    out = np.full((hdr["ncomp"], ny, nx), np.nan)
    boxes = []
    for lo, hi, data in ca.read_level_fabs(plotdir, level, hdr["ncomp"]):
        out[:, lo[1]:hi[1] + 1, lo[0]:hi[0] + 1] = data
        boxes.append((lo, hi))
    return out, boxes


def parse_header(plotdir: Path) -> dict:
    hdr = ca.parse_main_header(plotdir)
    lines = (plotdir / "Header").read_text().splitlines()
    p = 2 + hdr["ncomp"] + 4
    hdr["prob_hi"] = [float(v) for v in lines[p].split()]
    return hdr


def neo_hookean_predeformed_dw(mu, kappa, f0, deformation):
    """Vectorized 2-D implementation of NeoHookeanPredeformed::DW."""
    f0det = f0[..., 0, 0] * f0[..., 1, 1] - f0[..., 0, 1] * f0[..., 1, 0]
    f0inv = np.empty_like(f0)
    f0inv[..., 0, 0] = f0[..., 1, 1] / f0det
    f0inv[..., 0, 1] = -f0[..., 0, 1] / f0det
    f0inv[..., 1, 0] = -f0[..., 1, 0] / f0det
    f0inv[..., 1, 1] = f0[..., 0, 0] / f0det
    f = np.einsum("...ik,...kj->...ij", deformation, f0inv)
    det = f[..., 0, 0] * f[..., 1, 1] - f[..., 0, 1] * f[..., 1, 0]
    j23 = np.abs(det) ** (2.0 / 3.0)
    trace3 = np.sum(f * f, axis=(-2, -1)) + 1.0
    finvt = np.empty_like(f)
    finvt[..., 0, 0] = f[..., 1, 1] / det
    finvt[..., 0, 1] = -f[..., 1, 0] / det
    finvt[..., 1, 0] = -f[..., 0, 1] / det
    finvt[..., 1, 1] = f[..., 0, 0] / det
    return (mu / j23)[..., None, None] * (f - trace3[..., None, None] * finvt / 3.0) \
        + (kappa * (det - 1.0) * det)[..., None, None] * finvt


def model_arrays(data):
    mu = data[C_MODEL_MU]
    kappa = data[C_MODEL_KAPPA]
    f0 = np.empty(mu.shape + (2, 2))
    f0[..., 0, 0] = data[C_MODEL_F0XX]
    f0[..., 0, 1] = data[C_MODEL_F0XY]
    f0[..., 1, 0] = data[C_MODEL_F0YX]
    f0[..., 1, 1] = data[C_MODEL_F0YY]
    return mu, kappa, f0


def face_recovered_stress(data: np.ndarray, dx: float, dy: float) -> np.ndarray:
    """Recover nodal tensor columns from the conservative x/y face tractions."""
    u = np.moveaxis(data[[C_DISP_X, C_DISP_Y]], 0, -1)
    mu, kappa, f0 = model_arrays(data)
    ny, nx = mu.shape

    # Positive x faces: exact normal difference and averaged endpoint y derivative.
    gx = np.full((ny, nx - 1, 2, 2), np.nan)
    gx[..., :, 0] = (u[:, 1:, :] - u[:, :-1, :]) / dx
    gx[1:-1, :, :, 1] = (
        u[2:, :-1, :] - u[:-2, :-1, :] + u[2:, 1:, :] - u[:-2, 1:, :]
    ) / (4.0 * dy)
    mux = 0.5 * (mu[:, :-1] + mu[:, 1:])
    kx = 0.5 * (kappa[:, :-1] + kappa[:, 1:])
    f0x = 0.5 * (f0[:, :-1] + f0[:, 1:])
    fx = gx.copy()
    fx[..., 0, 0] += 1.0
    fx[..., 1, 1] += 1.0
    sx = neo_hookean_predeformed_dw(mux, kx, f0x, fx)

    # Positive y faces.
    gy = np.full((ny - 1, nx, 2, 2), np.nan)
    gy[..., :, 1] = (u[1:, :, :] - u[:-1, :, :]) / dy
    gy[:, 1:-1, :, 0] = (
        u[:-1, 2:, :] - u[:-1, :-2, :] + u[1:, 2:, :] - u[1:, :-2, :]
    ) / (4.0 * dx)
    muy = 0.5 * (mu[:-1, :] + mu[1:, :])
    ky = 0.5 * (kappa[:-1, :] + kappa[1:, :])
    f0y = 0.5 * (f0[:-1, :] + f0[1:, :])
    fy = gy.copy()
    fy[..., 0, 0] += 1.0
    fy[..., 1, 1] += 1.0
    sy = neo_hookean_predeformed_dw(muy, ky, f0y, fy)

    # At a node, average the two adjacent face values for each tensor column.
    recovered = np.full((ny, nx, 2, 2), np.nan)
    recovered[:, 1:-1, :, 0] = 0.5 * (sx[:, :-1, :, 0] + sx[:, 1:, :, 0])
    recovered[1:-1, :, :, 1] = 0.5 * (sy[:-1, :, :, 1] + sy[1:, :, :, 1])
    return recovered


def is_covered(level: int, shape: tuple[int, int], boxes_by_level, hdr) -> np.ndarray:
    covered = np.zeros(shape, dtype=bool)
    for fine in range(level + 1, hdr["finest"] + 1):
        ratio = int(round(hdr["dx"][level][0] / hdr["dx"][fine][0]))
        for lo, hi in boxes_by_level[fine]:
            i0 = max(0, math.ceil(lo[0] / ratio))
            i1 = min(shape[1] - 1, math.floor(hi[0] / ratio))
            j0 = max(0, math.ceil(lo[1] / ratio))
            j1 = min(shape[0] - 1, math.floor(hi[1] / ratio))
            if i0 <= i1 and j0 <= j1:
                covered[j0:j1 + 1, i0:i1 + 1] = True
    return covered


def composite_records(plotdir: Path) -> tuple[pd.DataFrame, dict]:
    hdr = parse_header(plotdir)
    levels, boxes = [], []
    for lev in range(hdr["finest"] + 1):
        data, level_boxes = dense_level(plotdir, lev, hdr)
        levels.append(data)
        boxes.append(level_boxes)

    records = []
    for lev, data in enumerate(levels):
        dx, dy = hdr["dx"][lev]
        ny, nx = data.shape[1:]
        ii, jj = np.meshgrid(np.arange(nx), np.arange(ny))
        x = hdr["prob_lo"][0] + ii * dx
        y = hdr["prob_lo"][1] + jj * dy
        xr, yr = x - CENTER[0], y - CENTER[1]
        r = np.hypot(xr, yr)
        c = np.divide(xr, r, out=np.ones_like(r), where=r > 0)
        s = np.divide(yr, r, out=np.zeros_like(r), where=r > 0)

        node = np.moveaxis(data[[C_STRESS_XX, C_STRESS_XY, C_STRESS_YX, C_STRESS_YY]], 0, -1)
        node = node.reshape(ny, nx, 2, 2)
        face = face_recovered_stress(data, dx, dy)
        srr_node = c * c * node[..., 0, 0] + c * s * node[..., 0, 1] \
            + s * c * node[..., 1, 0] + s * s * node[..., 1, 1]
        srr_face = c * c * face[..., 0, 0] + c * s * face[..., 0, 1] \
            + s * c * face[..., 1, 0] + s * s * face[..., 1, 1]

        valid = np.isfinite(data[C_PHI]) & ~is_covered(lev, (ny, nx), boxes, hdr)
        records.append(pd.DataFrame({
            "r": r[valid], "phi": data[C_PHI][valid], "eta": data[C_ETA][valid],
            "mu": data[C_MODEL_MU][valid], "srr_node": srr_node[valid],
            "srr_face": srr_face[valid], "level": lev,
            "weight": np.full(np.count_nonzero(valid), dx * dy),
        }))
    return pd.concat(records, ignore_index=True), hdr


def radial_profile(df: pd.DataFrame, rmax=0.1, dr=0.00025) -> pd.DataFrame:
    d = df[(df.r >= 0) & (df.r < rmax)].copy()
    d["bin"] = np.floor(d.r / dr).astype(int)
    rows = []
    for b, g in d.groupby("bin"):
        row = {"r": (b + 0.5) * dr, "n": len(g)}
        for col in ("phi", "eta", "mu", "srr_node", "srr_face"):
            finite = np.isfinite(g[col])
            if finite.any():
                row[col] = np.average(g.loc[finite, col], weights=g.loc[finite, "weight"])
            else:
                row[col] = np.nan
        rows.append(row)
    return (pd.DataFrame(rows).dropna(subset=["srr_node", "srr_face"])
            .sort_values("r").reset_index(drop=True))


def case_metrics(profile: pd.DataFrame, hdr: dict) -> dict:
    inner = profile[(profile.r >= 0.065) & (profile.phi >= 0.5) \
                    & profile.srr_node.notna() & profile.srr_face.notna()].copy()
    inner["node_minus_face"] = inner.srr_node - inner.srr_face
    peak = inner.loc[inner.node_minus_face.idxmax()]
    transition = profile[(profile.phi > 0.01) & (profile.phi < 0.99)]
    band = profile[(profile.r >= transition.r.min() - 0.002) \
                   & (profile.r <= transition.r.max() + 0.002) \
                   & profile.srr_node.notna() & profile.srr_face.notna()].copy()
    band["node_minus_face"] = band.srr_node - band.srr_face
    low = band.loc[band.node_minus_face.idxmin()]
    absolute = band.loc[band.node_minus_face.abs().idxmax()]
    return {
        "time": hdr["time"],
        "finest_level": hdr["finest"],
        "finest_dx_m": hdr["dx"][-1][0],
        "inner_peak_r_m": float(peak.r),
        "inner_peak_phi": float(peak.phi),
        "inner_peak_node_pa": float(peak.srr_node),
        "inner_peak_face_pa": float(peak.srr_face),
        "inner_peak_node_minus_face_pa": float(peak.node_minus_face),
        "interface_min_node_minus_face_pa": float(low.node_minus_face),
        "interface_min_r_m": float(low.r),
        "interface_max_abs_node_minus_face_pa": float(abs(absolute.node_minus_face)),
        "interface_max_abs_r_m": float(absolute.r),
        "phi_01_099_width_m": float(transition.r.max() - transition.r.min()) if len(transition) else None,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--case", action="append", required=True,
                        help="LABEL=PLOTDIR (repeatable)")
    parser.add_argument("--outdir", type=Path, required=True)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    profiles, metrics = {}, {}
    for spec in args.case:
        label, path = spec.split("=", 1)
        records, hdr = composite_records(Path(path))
        profile = radial_profile(records)
        profile.to_csv(args.outdir / f"{label}_composite.csv", index=False)
        profiles[label] = profile
        metrics[label] = case_metrics(profile, hdr)

    (args.outdir / "interface_metrics.json").write_text(json.dumps(metrics, indent=2) + "\n")

    fig, axes = plt.subplots(3, 1, figsize=(9, 11), sharex=True)
    for label, p in profiles.items():
        axes[0].plot(p.r, p.srr_node / 1e6, "-", lw=1.2, label=f"{label} nodal")
        axes[0].plot(p.r, p.srr_face / 1e6, "--", lw=1.2, label=f"{label} face")
        axes[1].plot(p.r, (p.srr_node - p.srr_face) / 1e6, lw=1.4, label=label)
        axes[2].plot(p.r, p.phi, lw=1.4, label=label)
    axes[0].set_ylabel(r"$\sigma_{rr}$ [MPa]")
    axes[1].set_ylabel("nodal - face [MPa]")
    axes[2].set_ylabel(r"$\phi$")
    axes[2].set_xlabel("radius [m]")
    axes[2].set_xlim(0.068, 0.090)
    for ax in axes:
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, ncol=2)
    fig.suptitle("Rod/tube casing-interface stress: spatial AMR composite")
    fig.tight_layout()
    fig.savefig(args.outdir / "interface_stress_composite.png", dpi=250, facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    main()
