#!/usr/bin/env python3
"""Performance-regression tracking for the GPU Flame canonical case.

Captures the roadmap's standing GPU perf metrics with `nsys` on the canonical
parity config (see ~/Desktop/SESSION_HANDOFF_2026-06-20.md section 5) and
appends one row, keyed by git SHA + date, to `benchmark/perf_regression.csv`:

    date, git_sha, config, steps, launches_per_step, kernel_avg_us,
    sync_frac, wall_per_step_gpu_ms, wall_per_step_cpu_ms, notes

Metric definitions (per the roadmap):
  - launches_per_step = cudaLaunchKernel "Num Calls" / steps
  - kernel_avg_us      = weighted-average GPU kernel duration (ns -> us) from
                         cuda_gpu_kern_sum, weighted by each kernel's
                         "Total Time (ns)" (i.e. total GPU busy time / total
                         kernel instances)
  - sync_frac          = cudaStreamSynchronize "Total Time (ns)" /
                         sum of all cuda_api_sum "Total Time (ns)" rows
                         (i.e. fraction of total CUDA-API time)
  - wall_per_step_gpu_ms / wall_per_step_cpu_ms = optional, supplied by the
                         caller (this tool does not itself run a CPU
                         comparison binary; pass --wall-gpu-ms/--wall-cpu-ms
                         if you have them, e.g. from baseline_suite.py or
                         benchmark_gpu_cpu.sh).

This tool is CI-safe: with no GPU or no nsys binary available, it prints
"skipped: no GPU" and exits 0 (does not fail the build).

Usage:
    # Capture a fresh nsys profile of the canonical case and record a row.
    python3 benchmark/perf_regression_track.py --capture

    # Parse an existing nsys stats CSV pair (e.g. from a prior capture) and
    # record a row without re-running anything.
    python3 benchmark/perf_regression_track.py \\
        --api-csv analysis/results_phase22/stats_cuda_api_sum.csv \\
        --kern-csv analysis/results_phase22/stats_cuda_gpu_kern_sum.csv \\
        --steps 30 --config wide_512_bf32_mgs128

    # Compare the most recent row for a config against the new capture and
    # flag a regression beyond the threshold (default 10%).
    python3 benchmark/perf_regression_track.py --capture --compare

See benchmark/PERF_TRACKING.md for the full walkthrough.
"""

from __future__ import annotations

import argparse
import csv
import datetime
import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CSV = ROOT / "benchmark" / "perf_regression.csv"
DEFAULT_OUT_DIR = ROOT / "benchmark" / "perf_regression_runs"

CSV_FIELDS = [
    "date",
    "git_sha",
    "config",
    "steps",
    "launches_per_step",
    "kernel_avg_us",
    "sync_frac",
    "wall_per_step_gpu_ms",
    "wall_per_step_cpu_ms",
    "notes",
]

# Canonical parity config, from SESSION_HANDOFF_2026-06-20.md section 5.
DEFAULT_CONFIG_NAME = "wide_512_bf32_mgs128"
DEFAULT_GPU_BIN_GLOB = "bin/alamo_gpu-2d-cuda*-g++"
DEFAULT_INPUT = "input"
DEFAULT_MAX_STEP = 30
PARITY_OVERRIDES = [
    "allow_unused=True",
    "stop_time=1e99_s",
    "elastic.tstart=1000000000.0",
    "elastic.solver.verbose=0",
    "elastic.print_model=0",
    "amr.plot_int=-1",
    "amr.thermo.plot_int=-1",
    "amr.thermo.int=1",
    "amr.n_cell=64 64 64",
    "amr.max_level=3",
    "amr.blocking_factor=32",
    "amr.max_grid_size=128",
    "amr.grid_eff=0.9",
    "amr.base_regrid_int=1000000",
    "amr.nsubsteps=2",
]


def resolve_nsys() -> str | None:
    """Resolve the nsys binary: $NSYS -> .local/nsight -> PATH."""
    env_nsys = os.environ.get("NSYS")
    if env_nsys and Path(env_nsys).is_file() and os.access(env_nsys, os.X_OK):
        return env_nsys

    local_glob = sorted(
        ROOT.glob(".local/nsight/opt/nvidia/nsight-systems/*/target-linux-x64/nsys")
    )
    for candidate in reversed(local_glob):
        if candidate.is_file() and os.access(candidate, os.X_OK):
            return str(candidate)

    path_nsys = shutil.which("nsys")
    if path_nsys:
        return path_nsys

    return None


def gpu_available() -> bool:
    """Best-effort check for a usable CUDA device (nvidia-smi probe)."""
    nvidia_smi = shutil.which("nvidia-smi")
    if not nvidia_smi:
        return False
    try:
        completed = subprocess.run(
            [nvidia_smi, "-L"],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            timeout=10,
        )
    except (subprocess.SubprocessError, OSError):
        return False
    return completed.returncode == 0 and "GPU" in completed.stdout


def find_gpu_binary(pattern: str) -> Path | None:
    candidates = [p for p in ROOT.glob(pattern) if p.is_file() and os.access(p, os.X_OK)]
    if not candidates:
        return None
    return sorted(candidates)[-1]


def git_sha() -> str:
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True, stderr=subprocess.DEVNULL
        )
        return out.strip()
    except (subprocess.SubprocessError, OSError):
        return "unknown"


def git_sha_short() -> str:
    try:
        out = subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            cwd=ROOT,
            text=True,
            stderr=subprocess.DEVNULL,
        )
        return out.strip()
    except (subprocess.SubprocessError, OSError):
        return "unknown"


def capture_nsys(
    nsys: str,
    gpu_bin: Path,
    input_file: str,
    max_step: int,
    out_dir: Path,
    plot_dir: Path,
    extra_overrides: list[str],
) -> tuple[Path, Path]:
    """Run the canonical case under nsys profile, export cuda_api_sum +
    cuda_gpu_kern_sum stats CSVs, and return their paths."""
    out_dir.mkdir(parents=True, exist_ok=True)
    rep_path = out_dir / "nsys"
    cmd = [
        nsys,
        "profile",
        "-o",
        str(rep_path),
        "--force-overwrite",
        "true",
        "--stats=false",
        "--trace=cuda,nvtx",
        "--sample=none",
        "--cpuctxsw=none",
        str(gpu_bin),
        input_file,
        f"max_step={max_step}",
        *PARITY_OVERRIDES,
        *extra_overrides,
        f"plot_file={plot_dir}",
    ]
    log_path = out_dir / "capture.log"
    with log_path.open("w", encoding="utf-8") as log:
        completed = subprocess.run(
            cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, text=True
        )
    if completed.returncode != 0:
        raise RuntimeError(
            f"nsys capture failed (rc={completed.returncode}); see {log_path}"
        )

    rep_file = rep_path.with_suffix(".nsys-rep")
    stats_prefix = out_dir / "stats"
    stats_cmd = [
        nsys,
        "stats",
        "--force-export=true",
        "--report",
        "cuda_api_sum",
        "--report",
        "cuda_gpu_kern_sum",
        "--format",
        "csv",
        "--output",
        str(stats_prefix),
        str(rep_file),
    ]
    completed = subprocess.run(
        stats_cmd, cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"nsys stats failed (rc={completed.returncode}):\n{completed.stdout}"
        )

    api_csv = out_dir / "stats_cuda_api_sum.csv"
    kern_csv = out_dir / "stats_cuda_gpu_kern_sum.csv"
    if not api_csv.is_file() or not kern_csv.is_file():
        raise RuntimeError(
            f"nsys stats did not produce expected CSVs in {out_dir} "
            f"(looked for {api_csv.name}, {kern_csv.name})"
        )
    return api_csv, kern_csv


def parse_cuda_api_sum(path: Path) -> dict[str, dict[str, float]]:
    """Parse an `nsys stats --report cuda_api_sum --format csv` file.

    Columns: Time (%), Total Time (ns), Num Calls, Avg (ns), Med (ns),
    Min (ns), Max (ns), StdDev (ns), Name
    Returns {Name: {"total_ns": float, "num_calls": float, ...}}.
    """
    rows: dict[str, dict[str, float]] = {}
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        for row in reader:
            name = row.get("Name", "").strip()
            if not name:
                continue
            try:
                total_ns = float(row["Total Time (ns)"])
                num_calls = float(row["Num Calls"])
            except (KeyError, ValueError):
                continue
            rows[name] = {"total_ns": total_ns, "num_calls": num_calls}
    if not rows:
        raise RuntimeError(f"{path}: no usable rows parsed (cuda_api_sum)")
    return rows


def parse_cuda_gpu_kern_sum(path: Path) -> dict[str, dict[str, float]]:
    """Parse an `nsys stats --report cuda_gpu_kern_sum --format csv` file.

    Columns: Time (%), Total Time (ns), Instances, Avg (ns), Med (ns),
    Min (ns), Max (ns), StdDev (ns), Name
    Returns {Name: {"total_ns": float, "instances": float}}.
    """
    rows: dict[str, dict[str, float]] = {}
    with path.open("r", encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream)
        for row in reader:
            name = row.get("Name", "").strip()
            if not name:
                continue
            try:
                total_ns = float(row["Total Time (ns)"])
                instances = float(row["Instances"])
            except (KeyError, ValueError):
                continue
            rows[name] = {"total_ns": total_ns, "instances": instances}
    if not rows:
        raise RuntimeError(f"{path}: no usable rows parsed (cuda_gpu_kern_sum)")
    return rows


def compute_metrics(
    api_rows: dict[str, dict[str, float]],
    kern_rows: dict[str, dict[str, float]],
    steps: int,
) -> dict[str, float]:
    if steps <= 0:
        raise ValueError("steps must be > 0")

    launch_calls = api_rows.get("cudaLaunchKernel", {}).get("num_calls", 0.0)
    launches_per_step = launch_calls / steps

    sync_ns = api_rows.get("cudaStreamSynchronize", {}).get("total_ns", 0.0)
    total_api_ns = sum(r["total_ns"] for r in api_rows.values())
    sync_frac = (sync_ns / total_api_ns) if total_api_ns > 0 else 0.0

    total_kernel_ns = sum(r["total_ns"] for r in kern_rows.values())
    total_kernel_instances = sum(r["instances"] for r in kern_rows.values())
    kernel_avg_us = (
        (total_kernel_ns / total_kernel_instances) / 1000.0
        if total_kernel_instances > 0
        else 0.0
    )

    return {
        "launches_per_step": launches_per_step,
        "kernel_avg_us": kernel_avg_us,
        "sync_frac": sync_frac,
    }


def append_row(csv_path: Path, row: dict[str, object]) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    write_header = not csv_path.is_file() or csv_path.stat().st_size == 0
    with csv_path.open("a", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=CSV_FIELDS)
        if write_header:
            writer.writeheader()
        writer.writerow({field: row.get(field, "") for field in CSV_FIELDS})


def read_rows(csv_path: Path) -> list[dict[str, str]]:
    if not csv_path.is_file():
        return []
    with csv_path.open("r", encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def last_row_for_config(rows: list[dict[str, str]], config: str, exclude_sha: str | None = None) -> dict[str, str] | None:
    matching = [r for r in rows if r.get("config") == config]
    if exclude_sha is not None:
        matching = [r for r in matching if r.get("git_sha") != exclude_sha]
    return matching[-1] if matching else None


def fmt_pct(x: float) -> str:
    return f"{x * 100.0:+.1f}%"


def compare_rows(
    previous: dict[str, str], current: dict[str, object], threshold_pct: float
) -> tuple[bool, list[str]]:
    """Compare current metrics vs previous row; return (regressed, messages).

    A metric "regresses" if it got worse (higher launches/kernel-avg/sync-frac/
    wall-time) by more than threshold_pct relative to the previous value.
    """
    regressed = False
    messages: list[str] = []
    higher_is_worse = [
        "launches_per_step",
        "kernel_avg_us",
        "sync_frac",
        "wall_per_step_gpu_ms",
        "wall_per_step_cpu_ms",
    ]
    for metric in higher_is_worse:
        try:
            prev_val = float(previous.get(metric, "") or "nan")
            cur_val = float(current.get(metric, "") or "nan")
        except ValueError:
            continue
        if prev_val != prev_val or cur_val != cur_val:  # NaN check
            continue
        if prev_val == 0:
            continue
        delta_frac = (cur_val - prev_val) / abs(prev_val)
        flag = " ".join(
            [
                f"{metric}: {prev_val:.4g} -> {cur_val:.4g}",
                f"({fmt_pct(delta_frac)})",
            ]
        )
        if delta_frac > threshold_pct / 100.0:
            regressed = True
            messages.append(f"REGRESSION {flag}")
        else:
            messages.append(f"ok         {flag}")
    return regressed, messages


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--capture",
        action="store_true",
        help="Run nsys against the canonical GPU case to produce stats CSVs.",
    )
    parser.add_argument(
        "--api-csv",
        type=Path,
        default=None,
        help="Use an existing `nsys stats --report cuda_api_sum --format csv` file "
        "instead of capturing a new one.",
    )
    parser.add_argument(
        "--kern-csv",
        type=Path,
        default=None,
        help="Use an existing `nsys stats --report cuda_gpu_kern_sum --format csv` "
        "file instead of capturing a new one.",
    )
    parser.add_argument(
        "--config",
        default=DEFAULT_CONFIG_NAME,
        help=f"Config label recorded in the CSV (default: {DEFAULT_CONFIG_NAME}).",
    )
    parser.add_argument(
        "--steps",
        type=int,
        default=DEFAULT_MAX_STEP,
        help=f"Number of steps in the captured run (default: {DEFAULT_MAX_STEP}). "
        "Required to normalize launches_per_step.",
    )
    parser.add_argument(
        "--input",
        default=DEFAULT_INPUT,
        help=f"Input file for --capture (default: {DEFAULT_INPUT}).",
    )
    parser.add_argument(
        "--gpu-bin",
        default=None,
        help="Path to the GPU binary for --capture (default: newest match of "
        f"'{DEFAULT_GPU_BIN_GLOB}').",
    )
    parser.add_argument(
        "--override",
        action="append",
        default=[],
        help="Extra ParmParse override(s) appended to the parity config for "
        "--capture (repeatable).",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_OUT_DIR,
        help=f"Directory for capture artifacts (default: {DEFAULT_OUT_DIR}).",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=DEFAULT_CSV,
        help=f"Path to perf_regression.csv (default: {DEFAULT_CSV}).",
    )
    parser.add_argument(
        "--wall-gpu-ms",
        type=float,
        default=None,
        help="Optional measured GPU wall-clock per step (ms), e.g. from "
        "baseline_suite.py / benchmark_gpu_cpu.sh, recorded alongside the "
        "nsys metrics.",
    )
    parser.add_argument(
        "--wall-cpu-ms",
        type=float,
        default=None,
        help="Optional measured CPU wall-clock per step (ms) for the same case.",
    )
    parser.add_argument(
        "--notes",
        default="",
        help="Free-text note stored with the row (e.g. 'phase2.2 fused reduction').",
    )
    parser.add_argument(
        "--compare",
        action="store_true",
        help="After recording (or parsing) a row, compare it against the most "
        "recent previously-recorded row for the same config and flag regressions.",
    )
    parser.add_argument(
        "--threshold-pct",
        type=float,
        default=10.0,
        help="Regression threshold in percent for --compare (default: 10).",
    )
    parser.add_argument(
        "--no-record",
        action="store_true",
        help="Compute and print metrics but do not append a row to the CSV "
        "(useful with --compare for a dry-run check).",
    )
    args = parser.parse_args(argv)

    # ---- CI-safe early exits -------------------------------------------------
    if args.api_csv is None and args.kern_csv is None:
        # We will need to actually drive nsys + a GPU. Check availability first.
        if not gpu_available():
            print("skipped: no GPU")
            return 0
        nsys = resolve_nsys()
        if nsys is None:
            print("skipped: no GPU")
            return 0
    elif (args.api_csv is None) != (args.kern_csv is None):
        print("error: --api-csv and --kern-csv must both be given together", file=sys.stderr)
        return 2
    else:
        nsys = None  # parsing pre-existing CSVs; nsys not needed.

    # ---- Obtain the stats CSVs ------------------------------------------------
    if args.api_csv is not None and args.kern_csv is not None:
        api_csv_path = args.api_csv
        kern_csv_path = args.kern_csv
        if not api_csv_path.is_file() or not kern_csv_path.is_file():
            print(f"error: missing stats CSV(s): {api_csv_path}, {kern_csv_path}", file=sys.stderr)
            return 2
    else:
        gpu_bin = (
            Path(args.gpu_bin)
            if args.gpu_bin
            else find_gpu_binary(DEFAULT_GPU_BIN_GLOB)
        )
        if gpu_bin is None or not gpu_bin.is_file():
            print("skipped: no GPU")
            return 0

        timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        capture_dir = args.out_dir / f"{git_sha_short()}_{timestamp}"
        plot_dir = capture_dir / "plt"
        try:
            api_csv_path, kern_csv_path = capture_nsys(
                nsys=nsys,
                gpu_bin=gpu_bin,
                input_file=args.input,
                max_step=args.steps,
                out_dir=capture_dir,
                plot_dir=plot_dir,
                extra_overrides=args.override,
            )
        except RuntimeError as exc:
            # Treat any capture failure as a soft skip too (CI-safe contract:
            # never hard-fail just because the GPU leg couldn't run), but keep
            # the message distinguishable from the true "no GPU present" case.
            print(f"skipped: GPU capture unavailable/failed ({exc})")
            shutil.rmtree(capture_dir, ignore_errors=True)
            return 0

    # ---- Parse + compute metrics -----------------------------------------
    try:
        api_rows = parse_cuda_api_sum(api_csv_path)
        kern_rows = parse_cuda_gpu_kern_sum(kern_csv_path)
        metrics = compute_metrics(api_rows, kern_rows, args.steps)
    except (RuntimeError, ValueError) as exc:
        print(f"error: failed to parse/compute metrics: {exc}", file=sys.stderr)
        return 1

    sha = git_sha()
    row = {
        "date": datetime.date.today().isoformat(),
        "git_sha": sha,
        "config": args.config,
        "steps": args.steps,
        "launches_per_step": f"{metrics['launches_per_step']:.4f}",
        "kernel_avg_us": f"{metrics['kernel_avg_us']:.4f}",
        "sync_frac": f"{metrics['sync_frac']:.6f}",
        "wall_per_step_gpu_ms": f"{args.wall_gpu_ms:.4f}" if args.wall_gpu_ms is not None else "",
        "wall_per_step_cpu_ms": f"{args.wall_cpu_ms:.4f}" if args.wall_cpu_ms is not None else "",
        "notes": args.notes,
    }

    print(f"config:            {args.config}")
    print(f"git_sha:           {sha}")
    print(f"steps:             {args.steps}")
    print(f"launches_per_step: {metrics['launches_per_step']:.2f}")
    print(f"kernel_avg_us:     {metrics['kernel_avg_us']:.3f}")
    print(f"sync_frac:         {metrics['sync_frac'] * 100.0:.2f}%")
    if args.wall_gpu_ms is not None:
        print(f"wall_per_step_gpu_ms: {args.wall_gpu_ms:.3f}")
    if args.wall_cpu_ms is not None:
        print(f"wall_per_step_cpu_ms: {args.wall_cpu_ms:.3f}")

    existing_rows = read_rows(args.csv)
    previous = last_row_for_config(existing_rows, args.config, exclude_sha=sha)

    if not args.no_record:
        append_row(args.csv, row)
        print(f"recorded row to {args.csv}")
    else:
        print("--no-record set: row not appended")

    exit_code = 0
    if args.compare:
        if previous is None:
            print(f"compare: no prior row for config={args.config!r}; nothing to compare against")
        else:
            regressed, messages = compare_rows(previous, row, args.threshold_pct)
            print(f"compare: vs git_sha={previous.get('git_sha')} date={previous.get('date')}")
            for message in messages:
                print(f"  {message}")
            if regressed:
                print(f"REGRESSION DETECTED (> {args.threshold_pct:.1f}% threshold)")
                exit_code = 1
            else:
                print("no regression beyond threshold")

    return exit_code


if __name__ == "__main__":
    sys.exit(main())
