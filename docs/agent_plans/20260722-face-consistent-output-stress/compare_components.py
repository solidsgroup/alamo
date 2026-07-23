#!/usr/bin/env python3
"""Compare named components in two Alamo nodal plotfiles."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

STRESS_VALIDATE = Path("/home/jackplum/Projects/chamberutils/stress_validate")
sys.path.insert(0, str(STRESS_VALIDATE))
import compare_alamo as ca  # noqa: E402


def load_components(plotdir: Path, names: list[str]):
    header = ca.parse_main_header(plotdir)
    indices = {name: header["names"].index(name) for name in names}
    levels = []
    for level in range(header["finest"] + 1):
        fabs = {}
        for lo, hi, data in ca.read_level_fabs(plotdir, level, header["ncomp"]):
            fabs[(lo, hi)] = {name: data[index] for name, index in indices.items()}
        levels.append(fabs)
    return header, levels


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("left", type=Path)
    parser.add_argument("right", type=Path)
    parser.add_argument("components", nargs="+")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    left_header, left = load_components(args.left, args.components)
    right_header, right = load_components(args.right, args.components)
    if left_header["finest"] != right_header["finest"]:
        raise ValueError("plotfiles have different finest levels")

    result = {}
    for name in args.components:
        max_abs = 0.0
        for level, (left_fabs, right_fabs) in enumerate(zip(left, right)):
            if left_fabs.keys() != right_fabs.keys():
                raise ValueError(f"level {level} has different FAB boxes")
            for key in left_fabs:
                max_abs = max(
                    max_abs,
                    float(np.max(np.abs(left_fabs[key][name] - right_fabs[key][name]))),
                )
        result[name] = {"max_abs_difference": max_abs}

    text = json.dumps(result, indent=2) + "\n"
    print(text, end="")
    if args.output:
        args.output.write_text(text)


if __name__ == "__main__":
    main()
