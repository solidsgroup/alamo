#!/usr/bin/env python3
"""Mechanically enumerate host-landing and explicit synchronization sites."""

import argparse
import re
import sys
import tempfile
from pathlib import Path


DEFAULT_PATHS = (
    Path("src/Integrator/Integrator.H"),
    Path("src/Integrator/Integrator.cpp"),
    Path("src/Integrator/Flame.cpp"),
    Path("src/Integrator/Base/Mechanics.H"),
    Path("src/Operator/Elastic.cpp"),
    Path("src/Solver/Nonlocal/Newton.H"),
    Path("src/Util/MPI.H"),
    Path("src/Util/Util.H"),
)

PATTERNS = (
    (
        "explicit_stream_sync",
        re.compile(r"\b(?:streamSynchronize(?:All)?|Device::synchronize|Gpu::synchronize)\s*\("),
    ),
    (
        "device_result_landing",
        re.compile(r"(?:reduce_data|flag)\.value\s*\(|\.dataValue\s*\("),
    ),
    (
        "multifab_host_norm",
        re.compile(r"\.norm0\s*\("),
    ),
    (
        "blocking_collective",
        re.compile(r"\b(?:ParallelDescriptor::Reduce\w+|ParallelAllReduce::\w+)\s*\("),
    ),
    (
        "blocking_copy",
        re.compile(r"\b(?:dtoh_memcpy|htod_memcpy|copyDtoH|copyHtoD)\s*\("),
    ),
)


def source_files(paths):
    for path in paths:
        if path.is_dir():
            yield from sorted(
                child
                for child in path.rglob("*")
                if child.suffix in {".H", ".cpp", ".cc"}
            )
        elif path.is_file():
            yield path


def inventory(paths):
    rows = []
    seen = set()
    for path in source_files(paths):
        if path in seen:
            continue
        seen.add(path)
        try:
            lines = path.read_text(errors="replace").splitlines()
        except OSError as error:
            raise ValueError(f"cannot read {path}: {error}") from error
        for lineno, raw in enumerate(lines, 1):
            code = raw.split("//", 1)[0]
            if not code.strip():
                continue
            for kind, pattern in PATTERNS:
                if pattern.search(code):
                    rows.append((kind, str(path), lineno, code.strip()))
    return rows


def render(rows):
    lines = ["kind\tpath\tline\tcode"]
    lines.extend(
        f"{kind}\t{path}\t{lineno}\t{code}"
        for kind, path, lineno, code in rows
    )
    return "\n".join(lines)


def unit():
    with tempfile.TemporaryDirectory() as temp:
        path = Path(temp) / "probe.cpp"
        path.write_text(
            "amrex::Gpu::streamSynchronizeAll();\n"
            "auto x = reduce_data.value(reduce_op);\n"
            "amrex::ParallelDescriptor::ReduceRealSum(x);\n"
            "// flag.value() in a comment is not a site\n"
        )
        rows = inventory((path,))
        assert [row[0] for row in rows] == [
            "explicit_stream_sync",
            "device_result_landing",
            "blocking_collective",
        ]
    print("sync_inventory: unit OK")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "paths",
        nargs="*",
        type=Path,
        help="source files/directories (default: Flame/Elastic/Newton closure)",
    )
    parser.add_argument("--output", type=Path)
    parser.add_argument("--unit", action="store_true")
    args = parser.parse_args()
    if args.unit:
        unit()
        return
    try:
        rows = inventory(tuple(args.paths) if args.paths else DEFAULT_PATHS)
    except ValueError as error:
        print(f"sync_inventory: {error}", file=sys.stderr)
        raise SystemExit(2) from error
    text = render(rows)
    if args.output:
        args.output.write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
