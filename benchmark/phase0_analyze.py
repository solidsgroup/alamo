#!/usr/bin/env python3
"""Summarise the stdlib-only artifacts emitted by benchmark/phase0_capture.sh."""
import argparse
import csv
import io
import re
import sys
import tempfile
from pathlib import Path


def read(path):
    try:
        return path.read_text(errors="replace")
    except OSError:
        return ""


def fmt(x):
    try:
        return f"{float(x):,.3f}"
    except (TypeError, ValueError):
        return x or "—"


def timing(d):
    out = []
    for mode in ("managed", "device"):
        p = d / "timing" / mode / "timing.txt"
        text = read(p)
        if not text:
            out.append((mode, "MISSING"))
            continue
        vals = dict(re.findall(r"([A-Za-z_]+)=([^\s]+)", text))
        out.append((mode, f"wall median {vals.get('wall_median_s', '—')} s; "
                           f"steady/step median {vals.get('steady_per_step_median_s', '—')} s; "
                           f"failed reps {vals.get('failed_reps', '—')}"))
    return out


def inventory_fields(text):
    fields = {}
    for key in (
        "host",
        "date",
        "pushed_from_host",
        "pushed_at",
        "local_branch",
        "local_head",
        "tree_hash",
        "src_hash",
        "local_dirty_files",
        "manifest_files",
    ):
        match = re.search(rf"(?:^|\s){re.escape(key)}=([^\s]+)", text)
        if match:
            fields[key] = match.group(1)
    return fields


def arena_table(d):
    rows = []
    ad = d / "arena"
    if not ad.exists():
        return rows
    steps = sorted((p for p in ad.iterdir() if p.is_dir()), key=lambda p: int(re.search(r"\d+", p.name).group()) if re.search(r"\d+", p.name) else 0)
    for step in steps:
        text = read(step / "run.log")
        if not text:
            rows.append((step.name, "MISSING", "", [], {}))
            continue
        if "CAPTURE_FAILED" in text or (step / "CAPTURE_FAILED").exists():
            rows.append((step.name, "CAPTURE_FAILED", "", [], {}))
            continue
        # Final AMReX summary is less ambiguous than request tables.
        m = re.search(r"\[The\s+Arena\] max space \(MB\) used\s+spread across MPI: \[([^]]+)]", text)
        used = m.group(1).strip() if m else "—"
        m = re.search(r"\[The\s+Arena\] max space \(MB\) allocated\s+spread across MPI: \[([^]]+)]", text)
        alloc = m.group(1).strip() if m else "—"
        requests = []
        totals = {"device": 0, "managed": 0, "pinned": 0}
        section = None
        for line in text.splitlines():
            heading = line.strip()
            if heading == "Device Memory Usage:":
                section = "device"
                continue
            if heading == "Managed Memory Usage:":
                section = "managed"
                continue
            if heading == "Pinned Memory Usage:":
                section = "pinned"
                continue
            if section is None or "Nalloc" in line or "MaxMem" in line:
                continue
            m = re.match(r"^\s*(.*?)\s+(\d+)\s+.*?(\d+(?:\.\d+)?\s+(?:B|KiB|MiB|GiB))\s*$", line)
            if not m:
                continue
            name = m.group(1).strip()
            if "Arena::Initialize()" in name:
                continue
            nalloc = int(m.group(2))
            totals[section] += nalloc
            if section == "device":
                requests.append((str(nalloc), name, m.group(3).strip()))
        rows.append(
            (
                step.name,
                f"used {used} MB",
                f"allocated {alloc} MB",
                requests[:3],
                totals,
            )
        )
    return rows


def csv_rows(path):
    try:
        with path.open(newline="", errors="replace") as fh:
            return list(csv.DictReader(fh))
    except (OSError, csv.Error):
        return []


def nsys(d):
    nd = d / "nsys"
    kern = next(iter(nd.glob("discovered_top10.tsv")), None)
    if kern is None:
        kern = next(iter(d.rglob("discovered_top10.tsv")), None)
    kr = []
    if kern:
        for line in read(kern).splitlines()[1:11]:
            if line.strip():
                kr.append(line)
    else:
        files = list(nd.glob("*cuda_gpu_kern_sum*.csv"))
        if files:
            for row in csv_rows(files[0])[:10]:
                kr.append(f"{row.get('Time (%)', '—')}% | {row.get('Total Time (ns)', '—')} ns | {row.get('Name', '—')}")
    mem = []
    for path in nd.glob("*cuda_gpu_mem_size_sum*.csv"):
        for row in csv_rows(path):
            mem.append(f"{row.get('Operation', '—')}: {row.get('Total (MB)', '—')} MB ({row.get('Count', '—')} calls)")
    api = []
    for path in nd.glob("*cuda_api_sum*.csv"):
        for row in csv_rows(path)[:10]:
            api.append(f"{row.get('Name', '—')}: {row.get('Total Time (ns)', '—')} ns ({row.get('Num Calls', '—')} calls)")
    sync = [x for x in api if re.search(r"Synchronize|DeviceSynchronize|StreamSynchronize", x, re.I)]
    idle = []
    idle_path = nd / "gpu_idle_summary.tsv"
    if idle_path.exists():
        with idle_path.open(newline="", errors="replace") as handle:
            idle = list(csv.DictReader(handle, delimiter="\t"))
    return kr, mem, api, sync, idle


def ncu(d):
    result = []
    for path in sorted((d / "ncu").glob("*.csv")):
        rows = csv_rows(path)
        metrics = {}
        for row in rows:
            name = row.get("Metric Name", "").strip()
            val = row.get("Metric Value", "").strip()
            unit = row.get("Metric Unit", "").strip()
            if name and val and re.search(
                r"achieved|throughput|occupancy|dram|sm__|speed of light|"
                r"register|block limit|duration",
                name,
                re.I,
            ):
                metrics.setdefault(name, f"{val} {unit}".strip())
        if metrics:
            result.append((path.stem, metrics))
    return result


def capture(d):
    lines = [f"## Capture: `{d.name}`", "", f"Source: `{d}`"]
    inv = read(d / "env" / "inventory.txt")
    provenance = [
        f"**{key}** `{value}`" for key, value in inventory_fields(inv).items()
    ]
    lines += ["### Provenance", "", "; ".join(provenance) if provenance else "MISSING inventory.txt", "", "### Timing", ""]
    lines += [f"- **{mode}**: {summary}" for mode, summary in timing(d)]
    lines += ["", "### Arena step 1/full", ""]
    ar = arena_table(d)
    lines += (
        [
            "| Capture | High-water used | Allocation | Device requests | "
            "Managed requests | Pinned requests |",
            "|---|---:|---:|---:|---:|---:|",
        ]
        + [
            f"| {step} | {used} | {allocated} | "
            f"{totals.get('device', '—')} | {totals.get('managed', '—')} | "
            f"{totals.get('pinned', '—')} |"
            for step, used, allocated, _, totals in ar
        ]
        if ar
        else ["MISSING arena artifacts"]
    )
    for step, _, _, req, _ in ar:
        if req:
            lines += ["", f"Top request rows ({step}):", "", "| Nalloc | Region | MaxMem |", "|---:|---|---:|"]
            lines += [f"| {n} | {name} | {mx} |" for n, name, mx in req[:10]]
    failures = [p for p in d.rglob("CAPTURE_FAILED")]
    # failures.txt is a report, not a marker: only surface it when non-empty.
    failures += [p for p in d.rglob("failures.txt") if read(p).strip()]
    flip = read(d / "flip" / "failures.txt")
    lines += ["", f"Flip failures: **{len([x for x in flip.splitlines() if x.strip()]) if flip else 0}**"]
    kr, mem, api, sync, idle = nsys(d)
    lines += ["", "### Nsight Systems", "", "Top kernels (top 10):"]
    lines += [f"- {x}" for x in kr] or ["- MISSING kernel summary"]
    lines += ["", "CUDA transfers:"] + ([f"- {x}" for x in mem] or ["- MISSING CUDA memory summary"])
    lines += ["", "CUDA API top rows:"] + ([f"- {x}" for x in api] or ["- MISSING CUDA API summary"])
    lines += [f"", f"Synchronization rows: {len(sync)}"]
    lines += ["", "GPU idle fractions:"]
    if idle:
        lines += [
            "",
            "| NVTX range | Instances | Idle fraction | Median instance idle |",
            "|---|---:|---:|---:|",
        ]
        lines += [
            f"| {row.get('range', '—')} | {row.get('instances', '—')} | "
            f"{row.get('idle_fraction', '—')} | "
            f"{row.get('median_instance_idle_fraction', '—')} |"
            for row in idle
        ]
    else:
        lines += ["- MISSING idle summary"]
    metrics = ncu(d)
    lines += ["", "### NCU dynamic metrics", ""]
    lines += [f"- **{name}**: " + "; ".join(f"{k}={v}" for k, v in vals.items()) for name, vals in metrics] or ["MISSING NCU CSV metrics"]
    if failures:
        lines += ["", "Markers: " + ", ".join(str(x.relative_to(d)) for x in failures)]
    return "\n".join(lines)


def unit():
    with tempfile.TemporaryDirectory() as td:
        d = Path(td) / "_phase0_fixture"
        (d / "env").mkdir(parents=True)
        (d / "timing" / "managed").mkdir(parents=True)
        (d / "nsys").mkdir(parents=True)
        (d / "env" / "inventory.txt").write_text(
            "host=test date=2026-07-31T00:00:00-05:00\n"
            "local_head=abc tree_hash=def src_hash=123\n"
        )
        (d / "timing" / "managed" / "timing.txt").write_text("mode=managed wall_median_s=1.2 failed_reps=0\n")
        (d / "nsys" / "gpu_idle_summary.tsv").write_text(
            "range\tinstances\tidle_fraction\tmedian_instance_idle_fraction\n"
            ":test\t1\t0.125000\t0.125000\n"
        )
        out = capture(d)
        assert all(
            expected in out
            for expected in (
                "**host** `test`",
                "**date** `2026-07-31T00:00:00-05:00`",
                "**local_head** `abc`",
                "wall median 1.2",
                "device",
                "MISSING",
                "| :test | 1 | 0.125000 | 0.125000 |",
            )
        )
    print("phase0_analyze: unit OK")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("captures", nargs="*", type=Path)
    ap.add_argument("--output", type=Path)
    ap.add_argument("--unit", action="store_true")
    args = ap.parse_args()
    if args.unit:
        unit(); return
    if not args.captures:
        ap.error("provide one or more capture directories (or --unit)")
    text = "\n\n".join(capture(p) for p in args.captures)
    if args.output:
        args.output.write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
