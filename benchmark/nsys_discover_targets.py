#!/usr/bin/env python3
"""Turn an Nsight Systems kernel summary into ranked NCU selectors."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Kernel:
    rank: int
    time_percent: float
    total_time_ms: float
    instances: int
    label: str
    selector: str
    name: str


def scoped_call_label(name: str) -> str:
    """Return the innermost scoped callable enclosing the generated lambda."""
    prefix = name.split("::[lambda", 1)[0]
    labels: list[str] = []

    for position, character in enumerate(prefix):
        if character != "(":
            continue

        cursor = position - 1
        while cursor >= 0 and prefix[cursor].isspace():
            cursor -= 1

        if cursor >= 0 and prefix[cursor] == ">":
            depth = 0
            while cursor >= 0:
                if prefix[cursor] == ">":
                    depth += 1
                elif prefix[cursor] == "<":
                    depth -= 1
                    if depth == 0:
                        cursor -= 1
                        break
                cursor -= 1
            while cursor >= 0 and prefix[cursor].isspace():
                cursor -= 1

        end = cursor + 1
        while cursor >= 0 and (prefix[cursor].isalnum() or prefix[cursor] == "_"):
            cursor -= 1
        start = cursor + 1
        if start == end or prefix[max(0, start - 2) : start] != "::":
            continue
        labels.append(prefix[start:end])

    if not labels:
        raise ValueError(f"could not derive a scoped callable from kernel name: {name}")
    return labels[-1]


def selector_for(name: str, label: str) -> str:
    """Build a demangled-name regex broad enough to survive tool formatting."""
    selector = rf".*::{re.escape(label)}.*"
    # A bare "::value" selector is poisonous: virtually every AMReX launch
    # contains the trait MaybeDeviceRunnable<...>::value, so NCU can match an
    # unrelated early kernel before it reaches the ranked ReduceOps::value
    # launch.  Anchor this callable to its owning reduction type.
    if label == "value" and "::ReduceOps<" in name:
        selector = r".*::ReduceOps.*::value.*"
    if label == "placementNew":
        match = re.search(r"::placementNew<([A-Za-z_]\w*(?:::[A-Za-z_]\w*)+)", name)
        if match:
            selector = rf".*::placementNew.*{re.escape(match.group(1))}.*"
    return f"regex:{selector}"


def read_kernels(path: Path, top: int) -> list[Kernel]:
    kernels: list[Kernel] = []
    with path.open(newline="", encoding="utf-8") as stream:
        for rank, row in enumerate(csv.DictReader(stream), 1):
            if rank > top:
                break
            label = scoped_call_label(row["Name"])
            kernels.append(
                Kernel(
                    rank=rank,
                    time_percent=float(row["Time (%)"]),
                    total_time_ms=float(row["Total Time (ns)"]) / 1.0e6,
                    instances=int(row["Instances"]),
                    label=label,
                    selector=selector_for(row["Name"], label),
                    name=row["Name"],
                )
            )
    if not kernels:
        raise ValueError(f"{path} contains no kernel rows")
    return kernels


def write_top(path: Path, kernels: list[Kernel]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(
            ("rank", "time_percent", "total_time_ms", "instances", "label", "kernel_name")
        )
        for kernel in kernels:
            writer.writerow(
                (
                    kernel.rank,
                    f"{kernel.time_percent:.6f}",
                    f"{kernel.total_time_ms:.6f}",
                    kernel.instances,
                    kernel.label,
                    kernel.name,
                )
            )


def write_targets(path: Path, kernels: list[Kernel], count: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    selected: list[Kernel] = []
    labels: set[str] = set()
    for kernel in kernels:
        if kernel.label in labels:
            continue
        selected.append(kernel)
        labels.add(kernel.label)
        if len(selected) == count:
            break

    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(("rank", "time_percent", "total_time_ms", "instances", "label", "selector"))
        for kernel in selected:
            writer.writerow(
                (
                    kernel.rank,
                    f"{kernel.time_percent:.6f}",
                    f"{kernel.total_time_ms:.6f}",
                    kernel.instances,
                    kernel.label,
                    kernel.selector,
                )
            )


def unit_test() -> None:
    names = (
        "void amrex::launch_global<Operator::Elastic<(int)1>::Fapply(int)::[lambda(int)]>()",
        "void amrex::launch_global<amrex::placementNew<Set::Matrix4<(int)3, (int)1>>(T1 *, long)::[lambda(long)]>()",
        "void amrex::launch_global<Operator::Elastic<(int)1>::SetModel(int)::[lambda(int)]>()",
        "void amrex::launch_global<Operator::Operator<(Grid)1>::Fsmooth(int)::[lambda(int)]>()",
        "void amrex::launch_global<T1::Type amrex::ReduceOps<amrex::ReduceOpSum>::value<T2>(T1&)::[lambda()]>()",
    )
    expected = ("Fapply", "placementNew", "SetModel", "Fsmooth", "value")
    actual = tuple(scoped_call_label(name) for name in names)
    if actual != expected:
        raise AssertionError(f"labels: expected {expected}, got {actual}")
    if selector_for(names[1], actual[1]) != r"regex:.*::placementNew.*Set::Matrix4.*":
        raise AssertionError("placementNew selector lost its discovered value type")
    value_selector = selector_for(names[4], actual[4])
    if value_selector != r"regex:.*::ReduceOps.*::value.*":
        raise AssertionError("ReduceOps::value selector lost its owning type")
    unrelated = (
        "void amrex::launch_global<std::enable_if<"
        "amrex::MaybeDeviceRunnable<T2, void>::value, void>::type "
        "amrex::ParallelFor<ResizeRandomSeed::[lambda()]>>()"
    )
    if re.search(value_selector.removeprefix("regex:"), unrelated):
        raise AssertionError("ReduceOps::value selector matches an unrelated trait")
    print("nsys target discovery unit tests passed")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", nargs="?", type=Path)
    parser.add_argument("--top", type=int, default=10)
    parser.add_argument("--ncu-targets", type=int, default=4)
    parser.add_argument("--top-output", type=Path)
    parser.add_argument("--target-output", type=Path)
    parser.add_argument("--unit", action="store_true")
    args = parser.parse_args()

    if args.unit:
        unit_test()
        return 0
    if args.summary is None or args.top_output is None or args.target_output is None:
        parser.error("summary, --top-output, and --target-output are required")
    if args.top < 1 or args.ncu_targets < 1:
        parser.error("--top and --ncu-targets must be positive")

    kernels = read_kernels(args.summary, args.top)
    write_top(args.top_output, kernels)
    write_targets(args.target_output, kernels, args.ncu_targets)
    for kernel in kernels:
        print(
            f"{kernel.rank:2d} {kernel.time_percent:7.3f}% "
            f"{kernel.total_time_ms:10.3f} ms x{kernel.instances:<7d} {kernel.label}"
        )
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (OSError, ValueError) as error:
        print(f"nsys target discovery failed: {error}", file=sys.stderr)
        raise SystemExit(2)
