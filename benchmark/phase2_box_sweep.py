#!/usr/bin/env python3
"""Phase 2.1 box/grid sweep for the Flame phase-field GPU path."""

from __future__ import annotations

import csv
import json
import math
import os
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "benchmark" / "phase2_box_sweep"


@dataclass(frozen=True)
class SweepCase:
    case_id: str
    n_cell: str
    max_level: int
    blocking_factor: int
    max_grid_size: int
    grid_eff: float
    base_regrid_int: int = 1000000
    nsubsteps: int = 2


CASES = (
    SweepCase("base_512_bf8_mgs32", "64 64 64", 3, 8, 32, 0.7),
    SweepCase("wide_512_bf16_mgs64", "64 64 64", 3, 16, 64, 0.9),
    SweepCase("wide_512_bf32_mgs128", "64 64 64", 3, 32, 128, 0.9),
    SweepCase("wide_1024_bf32_mgs128", "128 128 128", 3, 32, 128, 0.9),
    SweepCase("shallow_1024_bf32_mgs128", "256 256 256", 2, 32, 128, 0.9),
)


PROFILES = {
    "cpu_np8": ("bin/alamo-2d-clang++", 8),
    "gpu_fast": ("bin/alamo_gpu-2d-cuda86-g++", 1),
}


COMMON_OVERRIDES = (
    "allow_unused=True",
    "max_step=5",
    "stop_time=1e99_s",
    "elastic.tstart=1000000000.0",
    "elastic.solver.verbose=0",
    "elastic.print_model=0",
    "amr.plot_int=-1",
    "amr.thermo.plot_int=1",
    "amr.thermo.int=1",
)


def read_thermo(path: Path) -> tuple[list[str], list[list[float]]]:
    with path.open("r", encoding="utf-8") as stream:
        header = stream.readline().split()
        rows = [[float(value) for value in line.split()] for line in stream if line.strip()]
    if not header or not rows:
        raise RuntimeError(f"{path} is missing thermo data")
    return header, rows


def compare_thermo(reference: Path, candidate: Path, abs_tol: float = 1.0e-8, rel_tol: float = 1.0e-5) -> tuple[bool, float, float]:
    ref_header, ref_rows = read_thermo(reference)
    cand_header, cand_rows = read_thermo(candidate)
    if ref_header != cand_header or len(ref_rows) != len(cand_rows):
        return False, math.inf, math.inf
    max_abs = 0.0
    max_rel = 0.0
    ok = True
    for ref_row, cand_row in zip(ref_rows, cand_rows):
        for ref, cand in zip(ref_row, cand_row):
            if not (math.isfinite(ref) and math.isfinite(cand)):
                return False, math.inf, math.inf
            abs_err = abs(cand - ref)
            rel_err = abs_err / max(abs(ref), abs(cand), 1.0)
            max_abs = max(max_abs, abs_err)
            max_rel = max(max_rel, rel_err)
            if abs_err > abs_tol and rel_err > rel_tol:
                ok = False
    return ok, max_abs, max_rel


def case_overrides(case: SweepCase) -> list[str]:
    return [
        f"amr.n_cell={case.n_cell}",
        f"amr.max_level={case.max_level}",
        f"amr.blocking_factor={case.blocking_factor}",
        f"amr.max_grid_size={case.max_grid_size}",
        f"amr.grid_eff={case.grid_eff}",
        f"amr.base_regrid_int={case.base_regrid_int}",
        f"amr.nsubsteps={case.nsubsteps}",
    ]


def run_one(case: SweepCase, profile: str, binary: str, np: int) -> dict[str, object]:
    run_dir = OUT / case.case_id / profile
    run_dir.mkdir(parents=True, exist_ok=True)
    plot_dir = run_dir / "plot"
    log_file = run_dir / "run.log"
    cmd = [
        "mpiexec",
        "-np",
        str(np),
        str(ROOT / binary),
        "input",
        *COMMON_OVERRIDES,
        *case_overrides(case),
        f"plot_file={plot_dir}",
    ]
    start = time.monotonic()
    with log_file.open("w", encoding="utf-8") as log:
        completed = subprocess.run(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, text=True)
    elapsed = time.monotonic() - start
    result = {
        "case": case.case_id,
        "profile": profile,
        "binary": binary,
        "np": np,
        "elapsed_s": elapsed,
        "returncode": completed.returncode,
        "log": str(log_file.relative_to(ROOT)),
        "plot": str(plot_dir.relative_to(ROOT)),
        "command": cmd,
    }
    (run_dir / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    if completed.returncode != 0:
        raise RuntimeError(f"{case.case_id}/{profile} failed; see {log_file}")
    return result


def write_reports(rows: list[dict[str, object]]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    csv_path = OUT / "summary.csv"
    fields = [
        "case",
        "n_cell",
        "max_level",
        "blocking_factor",
        "max_grid_size",
        "grid_eff",
        "cpu_elapsed_s",
        "gpu_elapsed_s",
        "gpu_vs_cpu",
        "correctness",
        "max_abs",
        "max_rel",
    ]
    with csv_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)

    md = [
        "# Phase 2.1 Box/Grid Sweep",
        "",
        "Elastic solve is parsed but skipped with `elastic.tstart=1000000000.0` so this isolates the phase-field/thermal path.",
        "",
        "| Case | n_cell | max_level | blocking | max_grid | grid_eff | CPU np8 s | GPU s | GPU/CPU | Correctness |",
        "| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
    ]
    for row in rows:
        md.append(
            f"| {row['case']} | {row['n_cell']} | {row['max_level']} | "
            f"{row['blocking_factor']} | {row['max_grid_size']} | {row['grid_eff']} | "
            f"{float(row['cpu_elapsed_s']):.3f} | {float(row['gpu_elapsed_s']):.3f} | "
            f"{float(row['gpu_vs_cpu']):.3f} | {row['correctness']} |"
        )
    md.append("")
    md.append(f"CSV: `{csv_path.relative_to(ROOT)}`")
    (OUT / "README.md").write_text("\n".join(md) + "\n", encoding="utf-8")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for case in CASES:
        print(f"=== {case.case_id} ===", flush=True)
        results = {}
        for profile, (binary, np) in PROFILES.items():
            print(f"running {profile}", flush=True)
            results[profile] = run_one(case, profile, binary, np)
        cpu_thermo = OUT / case.case_id / "cpu_np8" / "plot" / "thermo.dat"
        gpu_thermo = OUT / case.case_id / "gpu_fast" / "plot" / "thermo.dat"
        ok, max_abs, max_rel = compare_thermo(cpu_thermo, gpu_thermo)
        cpu_elapsed = float(results["cpu_np8"]["elapsed_s"])
        gpu_elapsed = float(results["gpu_fast"]["elapsed_s"])
        rows.append(
            {
                "case": case.case_id,
                "n_cell": case.n_cell,
                "max_level": case.max_level,
                "blocking_factor": case.blocking_factor,
                "max_grid_size": case.max_grid_size,
                "grid_eff": case.grid_eff,
                "cpu_elapsed_s": cpu_elapsed,
                "gpu_elapsed_s": gpu_elapsed,
                "gpu_vs_cpu": gpu_elapsed / cpu_elapsed,
                "correctness": "ok" if ok else "FAIL",
                "max_abs": max_abs,
                "max_rel": max_rel,
            }
        )
        write_reports(rows)
    return 0 if all(row["correctness"] == "ok" for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
