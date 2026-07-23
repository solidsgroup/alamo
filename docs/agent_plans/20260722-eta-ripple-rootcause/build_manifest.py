#!/usr/bin/env python3
"""Build the machine-readable manifest for valid eta-ripple runs."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
import re
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
RESULTS = HERE / "results"
EXPERIMENTAL_EXECUTABLE_SHA256 = (
    "ea8bd70daf511fab48d2bc1dbcc60359a0d0d30a3f60b7f080a7f5206be10708"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_plot_reader():
    path = ROOT / "docs/agent_plans/20260722-free-circular-casing/check_stress_ripple.py"
    spec = importlib.util.spec_from_file_location("ripple_plot_reader", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module.ca


ca = load_plot_reader()


def plot_metadata(case_dir: Path) -> dict:
    node = case_dir / "00001node"
    cell = case_dir / "00001cell"
    header = ca.parse_main_header(node)
    level = int(header["finest"])
    return {
        "levels": level + 1,
        "dx_m": [float(value) for value in header["dx"][level]],
        "node_header_sha256": sha256(node / "Header"),
        "cell_header_sha256": sha256(cell / "Header"),
    }


def main() -> None:
    head = subprocess.check_output(
        ["git", "rev-parse", "--short=10", "HEAD"], cwd=ROOT, text=True
    ).strip()
    input_hash = sha256(HERE / "input_planar")
    cases = []
    planar_pattern = re.compile(
        r"^(x|y|diag)_p(.+)_n([0-9]+)_w(.+)_s(.+)$"
    )
    for case_dir in sorted((RESULTS / "planar").iterdir()):
        match = planar_pattern.match(case_dir.name)
        if not match or not (case_dir / "00001node/Header").is_file():
            continue
        orientation, pressure_text, ncell_text, width_text, shift_text = match.groups()
        metadata = plot_metadata(case_dir)
        width = float(width_text)
        metadata.update(
            {
                "name": f"planar-{case_dir.name}",
                "status": "complete",
                "geometry": f"planar-{orientation}",
                "command": (
                    "bash docs/agent_plans/20260722-eta-ripple-rootcause/"
                    f"run_planar.sh {orientation} {pressure_text} {ncell_text} "
                    f"{width_text} {shift_text}"
                ),
                "input_sha256": input_hash,
                "executable_sha256": EXPERIMENTAL_EXECUTABLE_SHA256,
                "mpi_ranks": 1,
                "pressure_pa": float(pressure_text),
                "eta_width_parameter_m": width,
                "eta_10_90_width_cells": (
                    2.0 * math.atanh(0.8) * width / metadata["dx_m"][0]
                ),
                "normal_shift_cells": float(shift_text),
            }
        )
        cases.append(metadata)

    radial_pattern = re.compile(r"^p(.+)_n([0-9]+)_w(.+)_s(.+)$")
    for case_dir in sorted((RESULTS / "radial").iterdir()):
        match = radial_pattern.match(case_dir.name)
        if not match or not (case_dir / "00001node/Header").is_file():
            continue
        pressure_text, ncell_text, width_text, shift_text = match.groups()
        metadata = plot_metadata(case_dir)
        width = float(width_text)
        metadata.update(
            {
                "name": f"radial-{case_dir.name}",
                "status": "complete",
                "geometry": "quarter-circle-signed-distance",
                "command": (
                    "bash docs/agent_plans/20260722-eta-ripple-rootcause/"
                    f"run_radial.sh {pressure_text} {ncell_text} {width_text} {shift_text}"
                ),
                "input_sha256": input_hash,
                "executable_sha256": EXPERIMENTAL_EXECUTABLE_SHA256,
                "mpi_ranks": 1,
                "pressure_pa": float(pressure_text),
                "eta_width_parameter_m": width,
                "eta_10_90_width_cells": (
                    2.0 * math.atanh(0.8) * width / metadata["dx_m"][0]
                ),
                "normal_shift_cells": float(shift_text),
            }
        )
        cases.append(metadata)

    for label in ("frozen-final-a", "frozen-final-b"):
        case_dir = RESULTS / "runs" / label
        node = case_dir / "05003node"
        cell = case_dir / "05003cell"
        header = ca.parse_main_header(node)
        level = int(header["finest"])
        cases.append(
            {
                "name": label,
                "status": "complete",
                "geometry": "frozen-production-double-circle",
                "command": (
                    "bash docs/agent_plans/20260722-eta-ripple-rootcause/"
                    f"run_frozen_baseline.sh {label}"
                ),
                "input_sha256": sha256(ROOT / "input_rt1s_ideal"),
                "executable_sha256": EXPERIMENTAL_EXECUTABLE_SHA256,
                "mpi_ranks": 1,
                "levels": level + 1,
                "dx_m": [float(value) for value in header["dx"][level]],
                "eta_10_90_width_cells": 8.0,
                "pressure_pa": 4316386.252462386,
                "node_header_sha256": sha256(node / "Header"),
                "cell_header_sha256": sha256(cell / "Header"),
                "repeatability": "bitwise-identical fields and logs",
            }
        )

    manifest = {
        "task": "eta-ripple-rootcause",
        "repository_head": head,
        "repository_was_dirty": True,
        "experimental_executable": {
            "path": "bin/alamo-2d-g++",
            "sha256_at_run_time": EXPERIMENTAL_EXECUTABLE_SHA256,
            "temporary_frozen_restart_hook": "reverted after evidence capture",
        },
        "post_revert_executable": {
            "path": "bin/alamo-2d-g++",
            "sha256": sha256(ROOT / "bin/alamo-2d-g++"),
        },
        "cases": cases,
        "invalid_attempts": [
            {
                "name": "early-frozen-restart-attempts",
                "location": "results/runs",
                "reason": "retained as failed diagnostics; not evidence cases",
            },
            {
                "name": "planar-input-unit-syntax-attempt",
                "location": "results/planar/x_p4316386.252462386_n128_w0.001138_s0",
                "reason": "temperature syntax rejected before a valid solve",
            },
            {
                "name": "radial-expression-whitespace-attempt",
                "location": "results/radial/invalid_expression_p4316386.252462386_n64_w0.00113800_s0",
                "reason": "CLI expression tokenized as constant eta=0.5",
            },
        ],
    }
    (RESULTS / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
