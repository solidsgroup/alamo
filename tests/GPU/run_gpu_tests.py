#!/usr/bin/env python3
"""Driver for the ALAMO GPU test suite.

Usage:
    python3 tests/GPU/run_gpu_tests.py [--test NAME] [--dry-run] [--timeout SEC]

Each sub-directory under tests/GPU/ that contains both an 'input' file and a
'test.py' script is treated as a test case. test.py receives the output directory
as argv[1] and exits:
    0  -> PASS
    1  -> FAIL
    77 -> SKIP (binary not found or yt not available)
"""

from __future__ import annotations

import argparse
import os
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
GPU_DIR = Path(__file__).resolve().parent


# ---------------------------------------------------------------------------
# Binary discovery helpers
# ---------------------------------------------------------------------------

def _find_bin(pattern: str) -> Path | None:
    candidates = [
        p for p in ROOT.glob(pattern)
        if p.is_file() and os.access(p, os.X_OK)
    ]
    return sorted(candidates)[-1] if candidates else None


def resolve_binaries() -> dict[str, Path | None]:
    return {
        "ALAMO_GPU_BIN":        Path(os.environ["ALAMO_GPU_BIN"])        if "ALAMO_GPU_BIN"        in os.environ else _find_bin("bin/alamo_gpu-2d-cuda86-g++"),
        "ALAMO_GPU_STRICT_BIN": Path(os.environ["ALAMO_GPU_STRICT_BIN"]) if "ALAMO_GPU_STRICT_BIN" in os.environ else _find_bin("bin/alamo_gpu-2d-nofast-cuda86-g++"),
        "ALAMO_GPU_3D_BIN":     Path(os.environ["ALAMO_GPU_3D_BIN"])     if "ALAMO_GPU_3D_BIN"     in os.environ else _find_bin("bin/alamo_gpu-3d-cuda86-g++"),
        "ALAMO_CPU_BIN":        Path(os.environ["ALAMO_CPU_BIN"])        if "ALAMO_CPU_BIN"        in os.environ else _find_bin("bin/alamo-2d-g++"),
    }


# ---------------------------------------------------------------------------
# Test discovery
# ---------------------------------------------------------------------------

def discover_tests() -> list[Path]:
    tests = []
    for d in sorted(GPU_DIR.iterdir()):
        if d.is_dir() and (d / "input").exists() and (d / "test.py").exists():
            tests.append(d)
    return tests


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description="Run ALAMO GPU tests")
    parser.add_argument("--test", metavar="NAME", help="Run only the test with this directory name")
    parser.add_argument("--dry-run", action="store_true", help="List tests and binaries, then exit")
    parser.add_argument("--timeout", type=int, default=660, help="Per-test timeout in seconds (default: 660)")
    args = parser.parse_args()

    tests = discover_tests()
    if args.test:
        tests = [t for t in tests if t.name == args.test]
        if not tests:
            print(f"ERROR: no test named '{args.test}' found under {GPU_DIR}", file=sys.stderr)
            return 1

    bins = resolve_binaries()

    if args.dry_run:
        print("=== GPU Test Suite — Dry Run ===")
        print(f"Root:     {ROOT}")
        print(f"Tests dir:{GPU_DIR}")
        print()
        print("Discovered tests:")
        for t in tests:
            print(f"  {t.name}")
        print()
        print("Binary paths:")
        for k, v in bins.items():
            status = str(v) if (v and v.exists()) else "NOT FOUND"
            print(f"  {k:28s} {status}")
        return 0

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    env = {**os.environ}
    for k, v in bins.items():
        if v is not None:
            env[k] = str(v)

    results: list[tuple[str, str, float]] = []

    for test_dir in tests:
        test_py = test_dir / "test.py"
        out_dir = test_dir / f"output_{timestamp}"
        out_dir.mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}")
        print(f"Running: {test_dir.name}")
        print(f"Output:  {out_dir}")
        print(f"{'='*60}")

        t0 = time.monotonic()
        try:
            proc = subprocess.run(
                [sys.executable, str(test_py), str(out_dir)],
                env=env,
                timeout=args.timeout,
                cwd=str(ROOT),
            )
            rc = proc.returncode
        except subprocess.TimeoutExpired:
            rc = -999
        elapsed = time.monotonic() - t0

        if rc == 0:
            status = "PASS"
        elif rc == 77:
            status = "SKIP"
        elif rc == -999:
            status = f"TIMEOUT (>{args.timeout}s)"
        else:
            status = f"FAIL (exit {rc})"

        results.append((test_dir.name, status, elapsed))
        print(f"\n{status}  [{elapsed:.1f}s]  {test_dir.name}")

    # Summary table
    print(f"\n{'='*60}")
    print("GPU Test Summary")
    print(f"{'='*60}")
    print(f"{'Test':<40} {'Result':<20} {'Time':>8}")
    print("-" * 60)
    passes = fails = skips = 0
    for name, status, elapsed in results:
        print(f"{name:<40} {status:<20} {elapsed:>7.1f}s")
        if status == "PASS":
            passes += 1
        elif status.startswith("SKIP"):
            skips += 1
        else:
            fails += 1
    print("-" * 60)
    print(f"{passes} passed, {fails} failed, {skips} skipped  ({len(results)} total)")

    return 0 if fails == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
