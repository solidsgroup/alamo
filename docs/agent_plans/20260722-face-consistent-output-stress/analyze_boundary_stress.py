#!/usr/bin/env python3
"""Compare legacy and face-consistent stress by physical-boundary distance."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

STRESS_VALIDATE = Path("/home/jackplum/Projects/chamberutils/stress_validate")
sys.path.insert(0, str(STRESS_VALIDATE))
import compare_alamo as ca  # noqa: E402


def fab_map(plotdir: Path, level: int, ncomp: int):
    return {(lo, hi): data for lo, hi, data in ca.read_level_fabs(plotdir, level, ncomp)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("legacy", type=Path)
    parser.add_argument("face_consistent", type=Path)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    old_header = ca.parse_main_header(args.legacy)
    new_header = ca.parse_main_header(args.face_consistent)
    stress_names = ["stress_xx", "stress_xy", "stress_yx", "stress_yy"]
    old_indices = [old_header["names"].index(name) for name in stress_names]
    new_indices = [new_header["names"].index(name) for name in stress_names]
    result = {"levels": {}}

    for level in range(new_header["finest"] + 1):
        old_fabs = fab_map(args.legacy, level, old_header["ncomp"])
        new_fabs = fab_map(args.face_consistent, level, new_header["ncomp"])
        if old_fabs.keys() != new_fabs.keys():
            raise ValueError(f"level {level} has different FAB boxes")
        nx = int(round(0.0877 / new_header["dx"][level][0]))
        ny = int(round(0.0877 / new_header["dx"][level][1]))
        maxima = {"boundary": 0.0, "one_node_in": 0.0, "interior": 0.0}
        counts = {name: 0 for name in maxima}
        old_dense = np.full((old_header["ncomp"], ny + 1, nx + 1), np.nan)
        new_dense = np.full((new_header["ncomp"], ny + 1, nx + 1), np.nan)
        for (lo, hi), new_data in new_fabs.items():
            old_data = old_fabs[(lo, hi)]
            old_dense[:, lo[1]:hi[1] + 1, lo[0]:hi[0] + 1] = old_data
            new_dense[:, lo[1]:hi[1] + 1, lo[0]:hi[0] + 1] = new_data
            error = np.max(
                np.abs(old_data[old_indices] - new_data[new_indices]), axis=0
            )
            ii, jj = np.meshgrid(
                np.arange(lo[0], hi[0] + 1), np.arange(lo[1], hi[1] + 1)
            )
            distance = np.minimum.reduce((ii, jj, nx - ii, ny - jj))
            masks = {
                "boundary": distance == 0,
                "one_node_in": distance == 1,
                "interior": distance >= 2,
            }
            for name, mask in masks.items():
                if np.any(mask):
                    maxima[name] = max(maxima[name], float(np.max(error[mask])))
                    counts[name] += int(np.count_nonzero(mask))
        result["levels"][str(level)] = {
            name: {"max_abs_new_minus_legacy_pa": maxima[name], "nodes": counts[name]}
            for name in maxima
        }
        old_phi = old_dense[old_header["names"].index("phi")]
        new_phi = new_dense[new_header["names"].index("phi")]
        old_sxx = old_dense[old_header["names"].index("stress_xx")]
        new_sxx = new_dense[new_header["names"].index("stress_xx")]
        old_syy = old_dense[old_header["names"].index("stress_yy")]
        new_syy = new_dense[new_header["names"].index("stress_yy")]

        def max_jump(boundary, adjacent, phi):
            mask = np.isfinite(boundary) & np.isfinite(adjacent) \
                & (phi > 0.01) & (phi < 0.99)
            return float(np.max(np.abs(boundary[mask] - adjacent[mask]))) \
                if np.any(mask) else None

        result["levels"][str(level)]["phi_band_boundary_jump_pa"] = {
            "bottom_legacy_sxx": max_jump(old_sxx[0], old_sxx[1], old_phi[0]),
            "bottom_face_consistent_sxx": max_jump(new_sxx[0], new_sxx[1], new_phi[0]),
            "right_legacy_syy": max_jump(old_syy[:, -1], old_syy[:, -2], old_phi[:, -1]),
            "right_face_consistent_syy": max_jump(new_syy[:, -1], new_syy[:, -2], new_phi[:, -1]),
        }

    text = json.dumps(result, indent=2) + "\n"
    print(text, end="")
    if args.output:
        args.output.write_text(text)


if __name__ == "__main__":
    main()
