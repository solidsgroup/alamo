#!/usr/bin/env python3
"""Summarise the stdlib-only artifacts emitted by benchmark/phase0_capture.sh."""
import argparse
import csv
import io
import re
import statistics
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
                           f"MAD {vals.get('steady_per_step_mad_s', '—')} s; "
                           f"sample sd {vals.get('steady_per_step_sd_s', '—')} s; "
                           f"startup median {vals.get('startup_wall_median_s', '—')} s; "
                           f"failed reps {vals.get('failed_reps', '—')}"))
    return out


def paired_timing(d):
    def reps(mode):
        result = {}
        for line in read(d / "timing" / mode / "reps.txt").splitlines():
            fields = line.split()
            if len(fields) < 4:
                continue
            try:
                rep = int(fields[0])
                rc = int(fields[1])
                per_step = float(fields[3])
                short_rc = int(fields[4]) if len(fields) >= 6 else 0
            except ValueError:
                continue
            if rc == 0 and short_rc == 0:
                result[rep] = per_step
        return result

    managed = reps("managed")
    device = reps("device")
    common = sorted(managed.keys() & device.keys())
    if not common:
        return None
    deltas = [device[rep] - managed[rep] for rep in common]
    percents = [
        100.0 * (device[rep] - managed[rep]) / managed[rep]
        for rep in common
        if managed[rep] != 0.0
    ]
    median_delta = statistics.median(deltas)
    mad_delta = statistics.median(
        abs(delta - median_delta) for delta in deltas
    )
    median_percent = statistics.median(percents) if percents else float("nan")
    lower = median_delta - 2.0 * mad_delta
    upper = median_delta + 2.0 * mad_delta
    if upper < 0.0:
        verdict = "device faster (2-MAD band excludes zero)"
    elif lower > 0.0:
        verdict = "device slower (2-MAD band excludes zero)"
    else:
        verdict = "inconclusive (2-MAD band overlaps zero)"
    paired_protocol = "startup_wall_median_s" in read(
        d / "timing" / "managed" / "timing.txt"
    )
    return {
        "count": len(common),
        "median_delta": median_delta,
        "median_percent": median_percent,
        "mad_delta": mad_delta,
        "verdict": verdict,
        "protocol": (
            "paired short/long"
            if paired_protocol
            else "legacy single-startup calibration; supporting only"
        ),
    }


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


def profile_memory_rows(text, heading):
    if heading not in text:
        return []
    section = text.rsplit(heading, 1)[1]
    section = re.split(
        r"\n(?:Device|Managed|Pinned) Memory Usage:|"
        r"\nTotal GPU global memory",
        section,
        maxsplit=1,
    )[0]
    rows = []
    for line in section.splitlines():
        if "Nalloc" in line or set(line.strip()) <= {"-"}:
            continue
        match = re.match(
            r"^\s*(.*?)\s+(\d+)\s+.*?"
            r"(\d+(?:\.\d+)?\s+(?:B|KiB|MiB|GiB))\s*$",
            line,
        )
        if not match:
            continue
        name = match.group(1).strip()
        if "Arena::Initialize()" in name:
            continue
        rows.append((name, int(match.group(2)), match.group(3).strip()))
    return rows


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
    mem_rows = []
    for path in nd.glob("*cuda_gpu_mem_size_sum*.csv"):
        mem_rows = csv_rows(path)
        for row in mem_rows:
            mem.append(f"{row.get('Operation', '—')}: {row.get('Total (MB)', '—')} MB ({row.get('Count', '—')} calls)")
    api = []
    api_rows = []
    for path in nd.glob("*cuda_api_sum*.csv"):
        api_rows = csv_rows(path)
        for row in api_rows[:10]:
            api.append(f"{row.get('Name', '—')}: {row.get('Total Time (ns)', '—')} ns ({row.get('Num Calls', '—')} calls)")
    sync = [x for x in api if re.search(r"Synchronize|DeviceSynchronize|StreamSynchronize", x, re.I)]
    idle = []
    idle_path = nd / "gpu_idle_summary.tsv"
    if idle_path.exists():
        with idle_path.open(newline="", errors="replace") as handle:
            idle = list(csv.DictReader(handle, delimiter="\t"))
    um_paths = sorted(nd.glob("um*.csv"))
    if any(path.stat().st_size for path in um_paths):
        um_status = "rows present"
    elif um_paths:
        um_status = "EMPTY/UNAVAILABLE (inspect stats.log before claiming zero)"
    else:
        um_status = "MISSING"
    run_text = read(nd / "run.log")
    managed_rows = profile_memory_rows(run_text, "Managed Memory Usage:")
    trace_steps = len(re.findall(
        r"^STEP\s+\d+\s+starts", run_text, re.MULTILINE
    ))
    normalized = {
        "trace_steps": trace_steps,
        "transfers": [],
        "blocking_transfer_calls_per_step": None,
        "sync_calls_per_step": None,
        "sync_ms_per_step": None,
    }
    if trace_steps:
        for row in mem_rows:
            try:
                total_mb = float(row.get("Total (MB)", ""))
                count = int(row.get("Count", ""))
            except ValueError:
                continue
            normalized["transfers"].append(
                (
                    row.get("Operation", "—"),
                    total_mb / trace_steps,
                    count / trace_steps,
                )
            )
        blocking_calls = 0
        sync_calls = 0
        sync_ns = 0
        for row in api_rows:
            name = row.get("Name", "")
            try:
                count = int(row.get("Num Calls", ""))
                total_ns = int(row.get("Total Time (ns)", ""))
            except ValueError:
                continue
            if re.match(r"^cudaMemcpy(?!.*Async)", name):
                blocking_calls += count
            if re.search(r"(?:Device|Stream|Event|Thread)Synchronize", name):
                sync_calls += count
                sync_ns += total_ns
        normalized["blocking_transfer_calls_per_step"] = (
            blocking_calls / trace_steps
        )
        normalized["sync_calls_per_step"] = sync_calls / trace_steps
        normalized["sync_ms_per_step"] = sync_ns / trace_steps / 1.0e6
    return kr, mem, api, sync, idle, um_status, normalized, managed_rows


def ncu_target_matches(target, kernel_name):
    selector = target.get("selector", "")
    if not selector.startswith("regex:"):
        return False
    if "ReduceOps" in selector and "::value" in selector:
        expected = tuple(re.findall(r"ReduceOp[A-Z][A-Za-z0-9_]*", selector))
        start = kernel_name.find("ReduceOps<")
        if start < 0 or "::value" not in kernel_name:
            return False
        start += len("ReduceOps<")
        depth = 1
        cursor = start
        while cursor < len(kernel_name) and depth:
            if kernel_name[cursor] == "<":
                depth += 1
            elif kernel_name[cursor] == ">":
                depth -= 1
            cursor += 1
        actual = tuple(
            re.findall(
                r"ReduceOp[A-Z][A-Za-z0-9_]*",
                kernel_name[start : cursor - 1],
            )
        )
        return bool(expected) and actual == expected
    pattern = selector.removeprefix("regex:")
    try:
        if re.search(pattern, kernel_name):
            return True
    except re.error:
        return False
    # NCU's demangled-name export sometimes removes the owning namespace/type
    # even though its kernel filter matched the full symbol (for example it
    # emits `FillBoundary(...)` for `BC::Constant::FillBoundary`).  Accept that
    # lossy form only when the callable itself still matches.  If NCU preserved
    # an explicit owner, the full selector above must match it.
    label = target.get("label", "")
    if not label or not re.search(
        rf"(?<![A-Za-z0-9_]){re.escape(label)}(?:<|\()", kernel_name
    ):
        return False
    explicit_owner = re.search(
        rf"[A-Za-z_][A-Za-z0-9_]*(?:<[^>]*>)?::{re.escape(label)}(?:<|\()",
        kernel_name,
    )
    return explicit_owner is None


def ncu(d):
    result = []
    targets = d / "nsys" / "discovered_ncu_targets.tsv"
    if not targets.exists():
        targets = d / "ncu" / "discovered_ncu_targets.tsv"
    target_rows = []
    if targets.exists():
        target_rows = list(
            csv.DictReader(io.StringIO(read(targets)), delimiter="\t")
        )
    targets_by_rank = {row.get("rank", ""): row for row in target_rows}
    semantic_issues = []
    report_paths = sorted((d / "ncu").glob("rank*.csv"))
    for path in report_paths:
        rows = csv_rows(path)
        kernel_name = rows[0].get("Kernel Name", "").strip() if rows else ""
        display_kernel = (
            kernel_name
            if len(kernel_name) <= 240
            else kernel_name[:237] + "..."
        )
        metrics = {"Kernel": display_kernel or "MISSING"}
        rank_match = re.match(r"rank(\d+)_", path.stem)
        target = targets_by_rank.get(rank_match.group(1), {}) if rank_match else {}
        if target and not ncu_target_matches(target, kernel_name):
            semantic_issues.append(
                f"{path.stem} expected {target.get('label', '?')}"
            )
        for row in rows:
            name = row.get("Metric Name", "").strip()
            val = row.get("Metric Value", "").strip()
            unit = row.get("Metric Unit", "").strip()
            if name and val and re.search(
                r"achieved|throughput|occupancy|dram|sm__|speed of light|"
                r"register|block limit|block size|grid size|waves per sm|"
                r"active warps|duration",
                name,
                re.I,
            ):
                metrics.setdefault(name, f"{val} {unit}".strip())
        if metrics:
            result.append((path.stem, metrics))
    expected = len(target_rows) if targets.exists() else None
    produced = len(report_paths)
    semantic_markers = list((d / "ncu").glob("*VALIDATION_FAILED*"))
    if semantic_markers or semantic_issues:
        detail = "; ".join(semantic_issues) or "validation marker present"
        coverage = (
            f"INVALID ({produced}/{expected if expected is not None else '?'} "
            f"files; semantic validation failed: {detail})"
        )
    elif expected is None:
        coverage = "UNAVAILABLE (discovery table missing)"
    elif produced == expected:
        coverage = f"COMPLETE ({produced}/{expected})"
    else:
        coverage = f"INCOMPLETE ({produced}/{expected})"
    return result, coverage


def synchronization_inventory(d):
    path = d / "env" / "sync_inventory.tsv"
    if not path.exists():
        return []
    with path.open(newline="", errors="replace") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def capture(d):
    lines = [f"## Capture: `{d.name}`", "", f"Source: `{d}`"]
    inv = read(d / "env" / "inventory.txt")
    if not inv:
        inv = read(d / "env" / "provenance.txt")
    provenance = [
        f"**{key}** `{value}`" for key, value in inventory_fields(inv).items()
    ]
    scheduler_out = read(d / "env" / "scheduler.out")
    scheduler_err_path = d / "env" / "scheduler.err"
    scheduler_err = read(scheduler_err_path)
    scheduler_err_status = (
        "MISSING"
        if not scheduler_err_path.exists()
        else ("NONEMPTY" if scheduler_err.strip() else "empty")
    )
    required_matches = re.findall(
        r"^required_leg_failures=(\d+)\s*$", scheduler_out, re.MULTILINE
    )
    required_verdict = required_matches[-1] if required_matches else "MISSING"
    lines += [
        "### Provenance",
        "",
        "; ".join(provenance) if provenance else "MISSING inventory.txt",
        "",
        f"Scheduler required-leg failures: **{required_verdict}**; "
        f"scheduler stderr: **{scheduler_err_status}**.",
        "",
        "### Timing",
        "",
    ]
    lines += [f"- **{mode}**: {summary}" for mode, summary in timing(d)]
    paired = paired_timing(d)
    if paired:
        lines += [
            "",
            "Paired arena comparison (device − managed): "
            f"**{paired['median_delta']:+.5f} s/step** "
            f"(**{paired['median_percent']:+.2f}%**), "
            f"MAD **{paired['mad_delta']:.5f} s/step**, "
            f"n=**{paired['count']}**; **{paired['verdict']}**. "
            f"Protocol: {paired['protocol']}.",
        ]
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
    if len(ar) >= 2:
        first = ar[0]
        last = ar[-1]
        first_step = re.search(r"\d+", first[0])
        last_step = re.search(r"\d+", last[0])
        if first_step and last_step:
            added_steps = int(last_step.group()) - int(first_step.group())
            if added_steps > 0:
                deltas = {
                    arena: (
                        last[4].get(arena, 0) - first[4].get(arena, 0)
                    ) / added_steps
                    for arena in ("device", "managed", "pinned")
                }
                lines += [
                    "",
                    "Request delta per added step: "
                    f"device **{deltas['device']:.3f}**, "
                    f"managed **{deltas['managed']:.3f}**, "
                    f"pinned **{deltas['pinned']:.3f}**.",
                ]
    if ar:
        lines += [
            "",
            "T5b live allocation count: **UNAVAILABLE**. The log prints the "
            "memory table after balanced teardown, not at the requested "
            "in-evolution endpoint, and contains Nalloc/AvgMem/MaxMem only.",
        ]
    failures = [p for p in d.rglob("*FAILED*") if p.is_file()]
    # failures.txt is a report, not a marker: only surface it when non-empty.
    failures += [p for p in d.rglob("failures.txt") if read(p).strip()]
    flip = read(d / "flip" / "failures.txt")
    lines += ["", f"Flip failures: **{len([x for x in flip.splitlines() if x.strip()]) if flip else 0}**"]
    (
        kr,
        mem,
        api,
        sync,
        idle,
        um_status,
        normalized,
        managed_rows,
    ) = nsys(d)
    lines += ["", "### Nsight Systems", "", "Top kernels (top 10):"]
    lines += [f"- {x}" for x in kr] or ["- MISSING kernel summary"]
    lines += ["", "CUDA transfers:"] + ([f"- {x}" for x in mem] or ["- MISSING CUDA memory summary"])
    if normalized["trace_steps"]:
        lines += [
            "",
            f"Normalized over **{normalized['trace_steps']}** coarse steps:",
            "",
            "| Direction | MB/step | Transfers/step |",
            "|---|---:|---:|",
        ]
        lines += [
            f"| {operation} | {mb_per_step:.6f} | {count_per_step:.3f} |"
            for operation, mb_per_step, count_per_step
            in normalized["transfers"]
        ]
        lines += [
            "",
            "Blocking CUDA transfer calls per step "
            f"(non-Async `cudaMemcpy*`): "
            f"**{normalized['blocking_transfer_calls_per_step']:.3f}**.",
            "CUDA synchronization calls per step: "
            f"**{normalized['sync_calls_per_step']:.3f}** "
            f"(**{normalized['sync_ms_per_step']:.3f} ms/step** in API time).",
        ]
    lines += ["", f"Unified-memory page-fault reports: **{um_status}**"]
    lines += ["", "Managed-pool application requests (profile run):"]
    if managed_rows:
        lines += [
            "",
            "| Region | Nalloc | MaxMem |",
            "|---|---:|---:|",
        ]
        lines += [
            f"| {name} | {nalloc} | {maxmem} |"
            for name, nalloc, maxmem in managed_rows[:10]
        ]
    else:
        lines += [
            "- No non-initialization managed-pool rows found "
            "(interpret only after pinning the T1 configuration)."
        ]
    lines += ["", "CUDA API top rows:"] + ([f"- {x}" for x in api] or ["- MISSING CUDA API summary"])
    lines += [f"", f"Synchronization rows: {len(sync)}"]
    lines += ["", "GPU idle fractions (fine-NVTX diagnostic; not the coarse T6b gate):"]
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
    metrics, ncu_coverage = ncu(d)
    lines += ["", "### NCU dynamic metrics", ""]
    lines += [f"Target coverage: **{ncu_coverage}**", ""]
    lines += [f"- **{name}**: " + "; ".join(f"{k}={v}" for k, v in vals.items()) for name, vals in metrics] or ["MISSING NCU CSV metrics"]
    inventory = synchronization_inventory(d)
    lines += ["", "### Mechanical synchronization inventory", ""]
    if inventory:
        counts = {}
        for row in inventory:
            kind = row.get("kind", "unknown")
            counts[kind] = counts.get(kind, 0) + 1
        lines += [
            "; ".join(f"**{kind}** {count}" for kind, count in counts.items()),
            "",
            "| Kind | Site | Code |",
            "|---|---|---|",
        ]
        for row in inventory:
            code = row.get("code", "—").replace("|", "&#124;").replace("`", "'")
            lines.append(
                f"| {row.get('kind', '—')} | "
                f"{row.get('path', '—')}:{row.get('line', '—')} | "
                f"`{code}` |"
            )
    else:
        lines += ["MISSING sync_inventory.tsv"]
    if failures:
        lines += ["", "Markers: " + ", ".join(str(x.relative_to(d)) for x in failures)]
    return "\n".join(lines)


def unit():
    valid_reduce = {
        "label": "value",
        "selector": (
            r"regex:.*::ReduceOps<[^,>]*ReduceOpSum[^,>]*>::value.*"
        ),
    }
    assert ncu_target_matches(
        valid_reduce,
        "void amrex::ReduceOps<amrex::ReduceOpSum>::value<T>()",
    )
    assert not ncu_target_matches(
        valid_reduce,
        "void amrex::MaybeDeviceRunnable<T, void>::value ResizeRandomSeed()",
    )
    assert not ncu_target_matches(
        valid_reduce,
        "void amrex::ReduceOps<amrex::ReduceOpLogicalOr>::value<T>()",
    )
    assert not ncu_target_matches(
        valid_reduce,
        "void amrex::ReduceOps<amrex::ReduceOpSum, "
        "amrex::ReduceOpSum>::value<T>()",
    )
    flame_advance = {
        "label": "Advance",
        "selector": r"regex:.*::Flame.*::Advance.*",
    }
    assert ncu_target_matches(
        flame_advance,
        "void Integrator::Flame::Advance(int)::[lambda()]()",
    )
    assert not ncu_target_matches(
        flame_advance,
        "void Integrator::Base::Mechanics<Model>::Advance(int)::[lambda()]()",
    )
    assert ncu_target_matches(
        {
            "label": "FillBoundary",
            "selector": r"regex:.*::Constant.*::FillBoundary.*",
        },
        "void launch_global<FillBoundary(BaseFab<double>&)::[lambda()]>()",
    )
    with tempfile.TemporaryDirectory() as td:
        d = Path(td) / "_phase0_fixture"
        (d / "env").mkdir(parents=True)
        (d / "timing" / "managed").mkdir(parents=True)
        (d / "timing" / "device").mkdir(parents=True)
        (d / "nsys").mkdir(parents=True)
        (d / "env" / "provenance.txt").write_text(
            "host=test date=2026-07-31T00:00:00-05:00\n"
            "local_head=abc tree_hash=def src_hash=123\n"
        )
        (d / "env" / "scheduler.out").write_text(
            "phase0 capture fixture\nrequired_leg_failures=0\n"
        )
        (d / "env" / "scheduler.err").write_text("")
        (d / "timing" / "managed" / "timing.txt").write_text(
            "mode=managed wall_median_s=1.2 "
            "steady_per_step_median_s=0.10 startup_wall_median_s=0.4 "
            "failed_reps=0\n"
        )
        (d / "timing" / "device" / "timing.txt").write_text(
            "mode=device wall_median_s=1.1 "
            "steady_per_step_median_s=0.09 startup_wall_median_s=0.4 "
            "failed_reps=0\n"
        )
        (d / "timing" / "managed" / "reps.txt").write_text(
            "1 0 1.2 0.10 0 0.4\n"
            "2 0 1.3 0.11 0 0.4\n"
            "3 0 1.1 0.09 0 0.4\n"
        )
        (d / "timing" / "device" / "reps.txt").write_text(
            "1 0 1.1 0.09 0 0.4\n"
            "2 0 1.2 0.10 0 0.4\n"
            "3 0 1.0 0.08 0 0.4\n"
        )
        (d / "nsys" / "gpu_idle_summary.tsv").write_text(
            "range\tinstances\tidle_fraction\tmedian_instance_idle_fraction\n"
            ":test\t1\t0.125000\t0.125000\n"
        )
        (d / "nsys" / "run.log").write_text(
            "STEP 1 starts ...\nSTEP 1 ends.\n"
            "STEP 2 starts ...\nSTEP 2 ends.\n"
            "Managed Memory Usage:\n"
            "------------------------------------------------------------\n"
            "Name                              Nalloc    AvgMem    MaxMem\n"
            "------------------------------------------------------------\n"
            "The_Arena::Initialize()                1      8 B      1 GiB\n"
            "Integrator::Flame::Regrid              2     16 B     64 MiB\n"
        )
        (d / "nsys" / "cuda_gpu_mem_size_sum.csv").write_text(
            "Total (MB),Count,Operation\n"
            "90,180,[CUDA memcpy Host-to-Device]\n"
        )
        (d / "nsys" / "cuda_api_sum.csv").write_text(
            "Total Time (ns),Num Calls,Name\n"
            "90000000,90,cudaStreamSynchronize\n"
            "1000,1,cudaMemcpy\n"
            "2000,2,cudaMemcpyAsync\n"
        )
        (d / "env" / "sync_inventory.tsv").write_text(
            "kind\tpath\tline\tcode\n"
            "explicit_stream_sync\tsrc/probe.H\t9\tamrex::Gpu::streamSynchronizeAll();\n"
        )
        out = capture(d)
        assert all(
            expected in out
            for expected in (
                "**host** `test`",
                "**date** `2026-07-31T00:00:00-05:00`",
                "**local_head** `abc`",
                "Scheduler required-leg failures: **0**",
                "wall median 1.2",
                "Paired arena comparison (device − managed): "
                "**-0.01000 s/step**",
                "device faster (2-MAD band excludes zero)",
                "| :test | 1 | 0.125000 | 0.125000 |",
                "Normalized over **2** coarse steps",
                "| [CUDA memcpy Host-to-Device] | 45.000000 | 90.000 |",
                "non-Async `cudaMemcpy*`): **0.500**",
                "CUDA synchronization calls per step: **45.000**",
                "| Integrator::Flame::Regrid | 2 | 64 MiB |",
                "| explicit_stream_sync | src/probe.H:9 |",
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
