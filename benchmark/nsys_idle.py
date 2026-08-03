#!/usr/bin/env python3
"""Compute GPU idle fractions inside selected projected NVTX ranges."""

import argparse
import bisect
import csv
import statistics
import sys
import tempfile
from pathlib import Path


DEFAULT_RANGES = (
    ":Integrator::Evolve",
    ":Integrator::Flame::TimeStepBegin",
    ":Integrator::Base::Mechanics::TimeStepBegin",
    ":MLMG::solve()",
    ":Integrator::TimeStep",
)


def numeric(row, key):
    try:
        return int(row[key])
    except (KeyError, TypeError, ValueError):
        return None


def read_gpu_intervals(path):
    intervals = []
    with path.open(newline="", errors="replace") as handle:
        for row in csv.DictReader(handle):
            start = numeric(row, "Start (ns)")
            duration = numeric(row, "Duration (ns)")
            if start is None or duration is None or duration <= 0:
                continue
            intervals.append((start, start + duration))
    intervals.sort()
    merged = []
    for start, end in intervals:
        if merged and start <= merged[-1][1]:
            merged[-1] = (merged[-1][0], max(merged[-1][1], end))
        else:
            merged.append((start, end))
    return merged


def read_ranges(path, selected):
    ranges = {name: [] for name in selected}
    with path.open(newline="", errors="replace") as handle:
        for row in csv.DictReader(handle):
            name = row.get("Name", "")
            if name not in ranges:
                continue
            start = numeric(row, "Projected Start (ns)")
            duration = numeric(row, "Projected Duration (ns)")
            if start is None or duration is None or duration <= 0:
                continue
            ranges[name].append((start, start + duration))
    return ranges


def busy_time(start, end, intervals, interval_ends):
    index = bisect.bisect_right(interval_ends, start)
    busy = 0
    while index < len(intervals):
        op_start, op_end = intervals[index]
        if op_start >= end:
            break
        busy += max(0, min(end, op_end) - max(start, op_start))
        index += 1
    return busy


def summarize(gpu_path, nvtx_path, selected):
    intervals = read_gpu_intervals(gpu_path)
    if not intervals:
        raise ValueError(f"no GPU intervals in {gpu_path}")
    interval_ends = [end for _, end in intervals]
    ranges = read_ranges(nvtx_path, selected)
    result = []
    for name in selected:
        samples = []
        projected_total = 0
        busy_total = 0
        for start, end in ranges[name]:
            projected = end - start
            busy = busy_time(start, end, intervals, interval_ends)
            projected_total += projected
            busy_total += busy
            samples.append(1.0 - busy / projected)
        if not samples:
            continue
        result.append(
            {
                "range": name,
                "instances": len(samples),
                "projected_ms": projected_total / 1.0e6,
                "gpu_busy_ms": busy_total / 1.0e6,
                "idle_ms": (projected_total - busy_total) / 1.0e6,
                "idle_fraction": 1.0 - busy_total / projected_total,
                "median_instance_idle_fraction": statistics.median(samples),
            }
        )
    return result


def render(rows):
    header = (
        "range\tinstances\tprojected_ms\tgpu_busy_ms\tidle_ms\t"
        "idle_fraction\tmedian_instance_idle_fraction"
    )
    lines = [header]
    for row in rows:
        lines.append(
            "{range}\t{instances}\t{projected_ms:.3f}\t{gpu_busy_ms:.3f}\t"
            "{idle_ms:.3f}\t{idle_fraction:.6f}\t"
            "{median_instance_idle_fraction:.6f}".format(**row)
        )
    return "\n".join(lines)


def unit():
    with tempfile.TemporaryDirectory() as temp:
        root = Path(temp)
        gpu = root / "gpu.csv"
        nvtx = root / "nvtx.csv"
        gpu.write_text(
            "Start (ns),Duration (ns),Name\n"
            "0,10,a\n"
            "5,10,b\n"
            "20,10,c\n"
        )
        nvtx.write_text(
            "Name,Projected Start (ns),Projected Duration (ns)\n"
            ":test,0,40\n"
        )
        rows = summarize(gpu, nvtx, (":test",))
        assert len(rows) == 1
        assert rows[0]["gpu_busy_ms"] == 25 / 1.0e6
        assert rows[0]["idle_fraction"] == 0.375
    print("nsys_idle: unit OK")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gpu_trace", nargs="?", type=Path)
    parser.add_argument("nvtx_trace", nargs="?", type=Path)
    parser.add_argument(
        "--range",
        dest="ranges",
        action="append",
        help="exact NVTX range name; repeat for multiple ranges",
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--unit", action="store_true")
    args = parser.parse_args()
    if args.unit:
        unit()
        return
    if args.gpu_trace is None or args.nvtx_trace is None:
        parser.error("provide CUDA GPU and NVTX projected-trace CSVs")
    selected = tuple(args.ranges) if args.ranges else DEFAULT_RANGES
    try:
        rows = summarize(args.gpu_trace, args.nvtx_trace, selected)
    except (OSError, csv.Error, ValueError) as error:
        print(f"nsys_idle: {error}", file=sys.stderr)
        raise SystemExit(2) from error
    if not rows:
        print("nsys_idle: none of the selected NVTX ranges were found", file=sys.stderr)
        raise SystemExit(3)
    text = render(rows)
    if args.output:
        args.output.write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
