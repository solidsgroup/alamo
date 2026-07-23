#!/usr/bin/env python3
"""Task-local oracle for the revision-bound Hydro/Fracture audit scan."""

import csv
import subprocess
import sys
from pathlib import Path


REPO = Path(__file__).resolve().parents[3]
FIELDS = [
    "schema_version", "port_id", "source_revision", "site_id", "file",
    "line", "pattern_id", "state", "evidence",
]
STATES = {"candidate", "converted", "not-applicable", "false-positive"}
REQUIRED_SITES = {
    ("src/Integrator/Hydro.cpp", "466", "GPU-013"),
    ("src/Integrator/Hydro.cpp", "467", "GPU-013"),
    ("src/Integrator/Hydro.cpp", "468", "GPU-013"),
    ("src/Integrator/Hydro.cpp", "469", "GPU-013"),
    ("src/Integrator/Fracture.H", "696", "GPU-013"),
}


def main():
    default = Path(__file__).parent / "results/current-source-coverage.csv"
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else default
    head = subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=REPO, text=True
    ).strip()
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        if reader.fieldnames != FIELDS:
            raise SystemExit(f"bad coverage header: {reader.fieldnames!r}")
    if not rows:
        raise SystemExit("coverage is empty")
    if {row["schema_version"] for row in rows} != {"3"}:
        raise SystemExit("coverage is not schema v3")
    if {row["port_id"] for row in rows} != {"hostile-current"}:
        raise SystemExit("coverage has wrong or mixed port IDs")
    if {row["source_revision"] for row in rows} != {head}:
        raise SystemExit("coverage revision does not equal current HEAD")
    if any(row["state"] not in STATES or not row["evidence"] for row in rows):
        raise SystemExit("coverage has invalid state or missing evidence")
    if len({row["site_id"] for row in rows}) != len(rows):
        raise SystemExit("coverage has duplicate site IDs")
    order = [
        (row["file"], row["pattern_id"], int(row["line"]), row["state"], row["site_id"])
        for row in rows
    ]
    if order != sorted(order):
        raise SystemExit("coverage order is not deterministic")
    sites = {(row["file"], row["line"], row["pattern_id"]) for row in rows}
    missing = REQUIRED_SITES - sites
    if missing:
        raise SystemExit(f"coverage misses reproduced aggregate sites: {sorted(missing)!r}")
    required_ids = {"GPU-002", "GPU-004", "GPU-013"}
    for target in ("src/Integrator/Hydro.cpp", "src/Integrator/Fracture.H"):
        ids = {row["pattern_id"] for row in rows if row["file"] == target}
        if required_ids - ids:
            raise SystemExit(f"{target} misses pattern classes {sorted(required_ids - ids)!r}")
    print(f"COVERAGE_ROWS={len(rows)}; REVISION={head}; REQUIRED_SITES={len(REQUIRED_SITES)}")


if __name__ == "__main__":
    main()
