#!/usr/bin/env python3
"""Frozen-field diagnostics for the rod-and-tube stress ripple.

This deliberately reimplements the branch's 2-D FaceGradient and
NeoHookeanPredeformed::DW formulas.  It never advances the simulation.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

STRESS_VALIDATE = Path("/home/jackplum/Projects/chamberutils/stress_validate")
sys.path.insert(0, str(STRESS_VALIDATE))
import compare_alamo as ca  # noqa: E402


CENTER = np.array(ca.CENTER)


def load_level(plotdir: Path, level: int):
    header = ca.parse_main_header(plotdir)
    names = header["names"]
    boxes = list(ca.read_level_fabs(plotdir, level, header["ncomp"]))
    ni = max(hi[0] for _, hi, _ in boxes) + 1
    nj = max(hi[1] for _, hi, _ in boxes) + 1
    sums = np.zeros((header["ncomp"], nj, ni))
    count = np.zeros((nj, ni), dtype=int)
    duplicate_spread = 0.0
    first = np.full_like(sums, np.nan)
    for (ilo, jlo), (ihi, jhi), data in boxes:
        sl = np.s_[jlo:jhi + 1, ilo:ihi + 1]
        occupied = count[sl] > 0
        if occupied.any():
            duplicate_spread = max(
                duplicate_spread,
                float(np.nanmax(np.abs(first[:, sl[0], sl[1]][:, occupied]
                                       - data[:, occupied]))),
            )
        empty = ~occupied
        if empty.any():
            target = first[:, sl[0], sl[1]]
            target[:, empty] = data[:, empty]
        sums[:, sl[0], sl[1]] += data
        count[sl] += 1
    fields = np.full_like(sums, np.nan)
    covered = count > 0
    fields[:, covered] = sums[:, covered] / count[covered]
    return header, {name: fields[n] for n, name in enumerate(names)}, covered, duplicate_spread


class FrozenField:
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
        return np.array([self.f["disp_x"][j, i], self.f["disp_y"][j, i]]) * signs

    def model(self, i, j):
        i, j, _ = self._reflect(i, j)
        return (
            self.f["model_mu"][j, i],
            self.f["model_kappa"][j, i],
            np.array([[self.f["model_F0xx"][j, i], self.f["model_F0xy"][j, i]],
                      [self.f["model_F0yx"][j, i], self.f["model_F0yy"][j, i]]]),
        )

    def face_gradient(self, i, j, face):
        ip, jp = i + (face == 0), j + (face == 1)
        grad = np.zeros((2, 2))
        grad[:, face] = (self.vector(ip, jp) - self.vector(i, j)) / self.dx[face]
        if face != 0:
            grad[:, 0] = (
                self.vector(i + 1, j) - self.vector(i - 1, j)
                + self.vector(ip + 1, jp) - self.vector(ip - 1, jp)
            ) / (4.0 * self.dx[0])
        if face != 1:
            grad[:, 1] = (
                self.vector(i, j + 1) - self.vector(i, j - 1)
                + self.vector(ip, jp + 1) - self.vector(ip, jp - 1)
            ) / (4.0 * self.dx[1])
        return grad

    def face_model(self, i, j, face):
        ip, jp = i + (face == 0), j + (face == 1)
        a, b = self.model(i, j), self.model(ip, jp)
        return tuple(0.5 * (x + y) for x, y in zip(a, b))

    def face_stress(self, i, j, face):
        return neo_predeformed_dw(np.eye(2) + self.face_gradient(i, j, face),
                                  *self.face_model(i, j, face))

    def reconstructed_stress(self, i, j):
        out = np.zeros((2, 2))
        for face in range(2):
            im, jm = i - (face == 0), j - (face == 1)
            out[:, face] = 0.5 * (
                self.face_stress(im, jm, face)[:, face]
                + self.face_stress(i, j, face)[:, face]
            )
        return out

    def nodal_gradient(self, i, j):
        return np.column_stack((
            (self.vector(i + 1, j) - self.vector(i - 1, j)) / (2.0 * self.dx[0]),
            (self.vector(i, j + 1) - self.vector(i, j - 1)) / (2.0 * self.dx[1]),
        ))

    def nodal_stress(self, i, j):
        return neo_predeformed_dw(np.eye(2) + self.nodal_gradient(i, j), *self.model(i, j))

    def face_divergence(self, i, j):
        out = np.zeros(2)
        for face in range(2):
            im, jm = i - (face == 0), j - (face == 1)
            out += (self.face_stress(i, j, face)[:, face]
                    - self.face_stress(im, jm, face)[:, face]) / self.dx[face]
        return out

    def has_stencil(self, i, j, radius=2):
        try:
            for di in range(-radius, radius + 1):
                for dj in range(-radius, radius + 1):
                    self._reflect(i + di, j + dj)
            return True
        except (IndexError, KeyError):
            return False


def neo_predeformed_dw(F, mu, kappa, F0):
    relative = F @ np.linalg.inv(F0)
    J = np.linalg.det(relative)
    J23 = abs(J) ** (2.0 / 3.0)
    tr = np.sum(relative * relative) + 1.0
    finvt = np.linalg.inv(relative).T
    return mu * (relative / J23 - tr * finvt / (3.0 * J23)) \
        + kappa * (J - 1.0) * J * finvt


def stress_components(stress):
    return np.array([stress[0, 0], stress[0, 1], stress[1, 0], stress[1, 1]])


def roughness(values):
    values = np.asarray(values)
    if values.size < 3:
        return np.nan
    return float(np.max(np.abs(values[2:] - 2.0 * values[1:-1] + values[:-2])))


def analyze_frozen(plotdir, outdir, level):
    header, fields, covered, duplicate_spread = load_level(plotdir, level)
    ff = FrozenField(fields, covered, header["dx"][level])
    dx = header["dx"][level][0]
    records = []
    for j in range(ff.nj):
        for i in range(ff.ni):
            if not covered[j, i] or not ff.has_stencil(i, j):
                continue
            x = header["prob_lo"][0] + i * dx
            y = header["prob_lo"][1] + j * dx
            qx, qy = CENTER[0] - x, y - CENTER[1]
            r = np.hypot(qx, qy)
            angle = np.arctan2(qy, qx)
            try:
                recon = ff.reconstructed_stress(i, j)
                nodal = ff.nodal_stress(i, j)
                saved = np.array([[fields["stress_xx"][j, i], fields["stress_xy"][j, i]],
                                  [fields["stress_yx"][j, i], fields["stress_yy"][j, i]]])
                div = ff.face_divergence(i, j)
            except (IndexError, KeyError, np.linalg.LinAlgError):
                continue
            rhs = np.array([fields["rhs_x"][j, i], fields["rhs_y"][j, i]])
            records.append((i, j, r, angle, *stress_components(saved),
                            *stress_components(recon), *stress_components(nodal),
                            *div, *rhs, fields["phi"][j, i], fields["model_mu"][j, i]))
    columns = ["i", "j", "r", "angle",
               "saved_xx", "saved_xy", "saved_yx", "saved_yy",
               "recon_xx", "recon_xy", "recon_yx", "recon_yy",
               "nodal_xx", "nodal_xy", "nodal_yx", "nodal_yy",
               "div_x", "div_y", "rhs_x", "rhs_y", "phi", "mu"]
    a = np.asarray(records)
    np.savetxt(outdir / "frozen_level2.csv", a, delimiter=",", header=",".join(columns), comments="")
    col = {name: n for n, name in enumerate(columns)}

    interface = (a[:, col["r"]] >= 0.075) & (a[:, col["r"]] <= 0.085)
    solid = a[:, col["mu"]] > 1.0e7
    region = interface & solid
    # The elastic operator imposes its BC directly on physical-boundary rows.
    # A face-divergence balance is therefore only the matching equation in the
    # strict domain interior.  Keep the all-row values for debugging, but do
    # not mistake their boundary contribution for an equilibrium residual.
    interior = (
        (a[:, col["i"]] > 0) & (a[:, col["i"]] < ff.xhi)
        & (a[:, col["j"]] > 0) & (a[:, col["j"]] < ff.nj - 1)
    )
    saved = a[:, [col[x] for x in ("saved_xx", "saved_xy", "saved_yx", "saved_yy")]]
    recon = a[:, [col[x] for x in ("recon_xx", "recon_xy", "recon_yx", "recon_yy")]]
    nodal = a[:, [col[x] for x in ("nodal_xx", "nodal_xy", "nodal_yx", "nodal_yy")]]
    div = a[:, [col["div_x"], col["div_y"]]]
    rhs = a[:, [col["rhs_x"], col["rhs_y"]]]
    resid = div - rhs

    metrics = {
        "level": level,
        "dx_m": dx,
        "samples": int(len(a)),
        "duplicate_node_max_spread": duplicate_spread,
        "saved_vs_reconstructed_linf_pa": float(np.max(np.abs(saved - recon))),
        "saved_vs_reconstructed_interface_linf_pa": float(np.max(np.abs((saved - recon)[region]))),
        "saved_vs_reconstructed_interior_linf_pa": float(np.max(np.abs((saved - recon)[interior]))),
        "saved_vs_reconstructed_interface_interior_linf_pa": float(
            np.max(np.abs((saved - recon)[region & interior]))),
        "nodal_vs_reconstructed_interface_linf_pa": float(np.max(np.abs((nodal - recon)[region]))),
        "face_equilibrium_residual_linf_pa_per_m": float(np.max(np.abs(resid[solid]))),
        "face_equilibrium_residual_l2_pa_per_m": float(np.sqrt(np.mean(resid[solid] ** 2))),
        "face_equilibrium_interior_linf_pa_per_m": float(
            np.max(np.abs(resid[solid & interior]))),
        "face_equilibrium_interior_l2_pa_per_m": float(
            np.sqrt(np.mean(resid[solid & interior] ** 2))),
        "rhs_linf_pa_per_m": float(np.max(np.abs(rhs[solid]))),
        "equilibrium_linf_relative_to_rhs_linf": float(np.max(np.abs(resid[solid])) / np.max(np.abs(rhs[solid]))),
    }

    # Axis-layer profiles and raw matching face tractions in the casing/interface band.
    layer_metrics = {}
    for layer in range(5):
        bottom_j = layer
        right_i = ff.xhi - layer
        bottom = []
        right = []
        xfaces = []
        yfaces = []
        for i in range(ff.ni):
            r = (ff.xhi - i) * dx
            if 0.075 <= r <= 0.085 and ff.has_stencil(i, bottom_j):
                bottom.append(ff.reconstructed_stress(i, bottom_j)[0, 0])
                xfaces.append(ff.face_stress(i, bottom_j, 0)[0, 0])
        for j in range(ff.nj):
            r = j * dx
            if 0.075 <= r <= 0.085 and ff.has_stencil(right_i, j):
                right.append(ff.reconstructed_stress(right_i, j)[1, 1])
                yfaces.append(ff.face_stress(right_i, j, 1)[1, 1])
        layer_metrics[str(layer)] = {
            "bottom_reconstructed_xx_roughness_pa": roughness(bottom),
            "bottom_raw_xface_xx_roughness_pa": roughness(xfaces),
            "right_reconstructed_yy_roughness_pa": roughness(right),
            "right_raw_yface_yy_roughness_pa": roughness(yfaces),
        }
    metrics["axis_layers"] = layer_metrics

    # Transpose comparison on the square finest-level domain.
    transpose = []
    material = []
    for j in range(ff.nj):
        for i in range(ff.ni):
            ip, jp = ff.xhi - j, ff.xhi - i
            if not (0 <= ip < ff.ni and 0 <= jp < ff.nj):
                continue
            if not (covered[j, i] and covered[jp, ip]):
                continue
            material.append(fields["phi"][j, i] - fields["phi"][jp, ip])
            transpose.append(fields["stress_xx"][j, i] - fields["stress_yy"][jp, ip])
    metrics["phi_transpose_rms"] = float(np.sqrt(np.mean(np.asarray(material) ** 2)))
    metrics["stress_xx_yy_transpose_rms_pa"] = float(np.sqrt(np.mean(np.asarray(transpose) ** 2)))

    # Profiles used for direct visual inspection.
    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    for layer, style in ((0, "-"), (1, "--"), (2, ":")):
        mask = (a[:, col["j"]] == layer) & (a[:, col["r"]] >= 0.045) & (a[:, col["r"]] <= 0.088)
        order = np.argsort(a[mask, col["r"]])
        axes[0, 0].plot(a[mask][order, col["r"]], a[mask][order, col["recon_xx"]] / 1e6,
                        style, label=f"layer {layer}")
        mask = (a[:, col["i"]] == ff.xhi - layer) & (a[:, col["r"]] >= 0.045) & (a[:, col["r"]] <= 0.088)
        order = np.argsort(a[mask, col["r"]])
        axes[0, 1].plot(a[mask][order, col["r"]], a[mask][order, col["recon_yy"]] / 1e6,
                        style, label=f"layer {layer}")
    axes[0, 0].set_title("bottom/x-axis reconstructed $P_{xx}$")
    axes[0, 1].set_title("right/y-axis reconstructed $P_{yy}$")
    mask = region
    axes[1, 0].scatter(a[mask, col["angle"]] * 180 / np.pi,
                       (saved - nodal)[mask, 0] / 1e6, s=5)
    axes[1, 0].set_title("saved face reconstruction - nodal $P_{xx}$")
    axes[1, 1].scatter(a[mask, col["angle"]] * 180 / np.pi,
                       np.linalg.norm(resid[mask], axis=1) / 1e6, s=5)
    axes[1, 1].set_title("native face residual magnitude [MPa/m]")
    for ax in axes.flat:
        ax.grid(alpha=0.25)
        ax.set_xlabel("r [m]" if ax in axes[0] else "angle [deg]")
    axes[0, 0].legend()
    axes[0, 1].legend()
    fig.savefig(outdir / "frozen_diagnostics.png", dpi=180)
    plt.close(fig)
    return metrics


def make_manufactured(n=257):
    dx = 0.0877 / (n - 1)
    x = np.arange(n) * dx
    y = 0.0877 + np.arange(n) * dx
    X, Y = np.meshgrid(x, y)
    px, py = X - CENTER[0], Y - CENTER[1]
    r = np.hypot(px, py)
    fields = {}
    covered = np.ones((n, n), dtype=bool)
    fields["disp_x"] = 2.0e-3 * px
    fields["disp_y"] = 2.0e-3 * py
    fields["model_mu"] = np.full_like(X, 26e9)
    fields["model_kappa"] = np.full_like(X, 70e9)
    fields["model_F0xx"] = np.ones_like(X)
    fields["model_F0xy"] = np.zeros_like(X)
    fields["model_F0yx"] = np.zeros_like(X)
    fields["model_F0yy"] = np.ones_like(X)
    ff = FrozenField(fields, covered, (dx, dx))
    affine = np.array([ff.reconstructed_stress(i, j)
                       for j in range(1, n - 1) for i in range(1, n - 1)])

    scale = 1.0e-3 * (1.0 + 20.0 * r ** 2)
    fields["disp_x"] = scale * px
    fields["disp_y"] = scale * py
    transition = 0.5 * (1.0 - np.tanh((r - 0.080) / (4.0 * dx)))
    fields["model_mu"] = 0.5e6 + (26e9 - 0.5e6) * transition
    fields["model_kappa"] = 0.5e6 + (70e9 - 0.5e6) * transition
    ff = FrozenField(fields, covered, (dx, dx))
    transpose_error = []
    shear_asymmetry = []
    for j in range(2, n - 2):
        for i in range(2, n - 2):
            ip, jp = n - 1 - j, n - 1 - i
            if not (2 <= ip < n - 2 and 2 <= jp < n - 2):
                continue
            p = ff.reconstructed_stress(i, j)
            pt = ff.reconstructed_stress(ip, jp)
            transpose_error.append(p[0, 0] - pt[1, 1])
            shear_asymmetry.append(p[0, 1] - p[1, 0])
    return {
        "affine_reconstructed_component_std_pa": [float(v) for v in affine.std(axis=0).ravel()],
        "affine_reconstructed_linf_from_mean_pa": float(np.max(np.abs(affine - affine.mean(axis=0)))),
        "radial_diffuse_xx_yy_transpose_linf_pa": float(np.max(np.abs(transpose_error))),
        "radial_diffuse_xx_yy_transpose_rms_pa": float(np.sqrt(np.mean(np.asarray(transpose_error) ** 2))),
        "radial_diffuse_reconstructed_nonsymmetry_linf_pa": float(np.max(np.abs(shear_asymmetry))),
    }


def compare_output_toggle(off_dir, on_dir):
    results = {}
    for suffix in ("00000node", "00050node", "00100node"):
        off_header = ca.parse_main_header(off_dir / suffix)
        on_header = ca.parse_main_header(on_dir / suffix)
        common = sorted(set(off_header["names"]) & set(on_header["names"]))
        step = {}
        for level in range(off_header["finest"] + 1):
            _, off, off_cov, _ = load_level(off_dir / suffix, level)
            _, on, on_cov, _ = load_level(on_dir / suffix, level)
            mask = off_cov & on_cov
            for name in common:
                diff = np.abs(off[name][mask] - on[name][mask])
                step[f"L{level}:{name}"] = float(np.max(diff)) if diff.size else 0.0
        results[suffix] = {
            "max_common_field_abs_difference": max(step.values(), default=0.0),
            "nonzero_fields": {k: v for k, v in step.items() if v != 0.0},
        }
    return results


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--node", type=Path,
                        default=Path("output_ideal_rod_and_tube_free_circular_t1/05000node"))
    parser.add_argument("--outdir", type=Path,
                        default=Path("docs/agent_plans/20260722-free-circular-casing/results/ripple_checks"))
    parser.add_argument("--level", type=int, default=2)
    args = parser.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    result = {
        "frozen_field": analyze_frozen(args.node, args.outdir, args.level),
        "manufactured": make_manufactured(),
        "output_toggle": compare_output_toggle(
            Path("output_face_stress_off_short"), Path("output_face_stress_on_short")),
    }
    (args.outdir / "metrics.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
