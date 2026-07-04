#!/usr/bin/env python3
"""Compare two ALAMO thermo.dat files column-by-column."""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path


def read_thermo(path: Path) -> tuple[list[str], list[list[float]]]:
    with path.open("r", encoding="utf-8") as stream:
        header = stream.readline().split()
        rows: list[list[float]] = []
        for line_number, line in enumerate(stream, start=2):
            if not line.strip():
                continue
            try:
                row = [float(value) for value in line.split()]
            except ValueError as exc:
                raise ValueError(f"{path}:{line_number}: non-numeric thermo row") from exc
            if len(row) != len(header):
                raise ValueError(
                    f"{path}:{line_number}: expected {len(header)} columns, got {len(row)}"
                )
            rows.append(row)

    if not header:
        raise ValueError(f"{path}: empty header")
    if not rows:
        raise ValueError(f"{path}: no data rows")
    return header, rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("reference", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--abs-tol", type=float, default=1.0e-8)
    parser.add_argument("--rel-tol", type=float, default=1.0e-6)
    args = parser.parse_args()

    ref_header, ref_rows = read_thermo(args.reference)
    cand_header, cand_rows = read_thermo(args.candidate)

    if ref_header != cand_header:
        print("header mismatch", file=sys.stderr)
        print(f"reference: {ref_header}", file=sys.stderr)
        print(f"candidate: {cand_header}", file=sys.stderr)
        return 2
    if len(ref_rows) != len(cand_rows):
        print(
            f"row-count mismatch: reference={len(ref_rows)} candidate={len(cand_rows)}",
            file=sys.stderr,
        )
        return 2

    failed = False
    print(f"{'column':<24} {'max_abs':>14} {'max_rel':>14} status")
    for column, name in enumerate(ref_header):
        max_abs = 0.0
        max_rel = 0.0
        for row, (ref, cand) in enumerate(zip(ref_rows, cand_rows), start=1):
            ref_value = ref[column]
            cand_value = cand[column]
            if not (math.isfinite(ref_value) and math.isfinite(cand_value)):
                print(f"{name:<24} non-finite at row {row}", file=sys.stderr)
                failed = True
                continue
            abs_err = abs(cand_value - ref_value)
            scale = max(abs(ref_value), abs(cand_value), 1.0)
            rel_err = abs_err / scale
            max_abs = max(max_abs, abs_err)
            max_rel = max(max_rel, rel_err)

        ok = max_abs <= args.abs_tol or max_rel <= args.rel_tol
        failed = failed or not ok
        print(f"{name:<24} {max_abs:14.6e} {max_rel:14.6e} {'ok' if ok else 'FAIL'}")

    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
