#!/usr/bin/env python3
"""Phase 1 task 1.B -- local-A1000 validation runner.

Local-only: no hardware auto-detect across clusters, no Slurm. See
benchmark/validate/README.md ("Two entry points, not one auto-detecting
script") for why this is deliberately separate from the NOVA path
(run_validation_nova.slurm). This script assumes it is running on the box
where the binaries already live in `bin/`; the case-running/extraction core
is shared with the NOVA driver via validation_common.py.

For each case in cases.manifest.yaml whose `hardware` list includes `local`,
and for each requested profile (cpu / gpu_strict / gpu_fast), this:
  1. runs the matching alamo binary via mpiexec,
  2. lands run.log + thermo.dat + the final <N>node/<N>cell plotfile pair into
     a bundle directory per the schema in README.md,
  3. calls extract_metrics.py to populate metrics.json + field_norms.json.

Each profile gets its own bundle directory (device differs: cpu vs
a1000_sm86) -- two sibling bundles from one invocation is what task 1.E's two
named references (cpu_strict / gpu_alpha1_local) are built from.

Usage:
    python3 benchmark/validate/run_validation_local.py
    python3 benchmark/validate/run_validation_local.py --profiles gpu_strict --case canonical_2d_elastic
    python3 benchmark/validate/run_validation_local.py --profiles cpu,gpu_strict,gpu_fast  # gpu_fast = perf smoke row, not a CORRECTNESS claim
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import validation_common as vc


def device_tag(profile: str, arch: str) -> str:
    if profile == "cpu":
        return "cpu"
    if profile == "gpu_strict":
        return f"a1000_sm{arch}_strict"
    if profile == "gpu_fast":
        return f"a1000_sm{arch}_fast"
    raise ValueError(f"unknown profile: {profile}")


def mpiexec_launch(profile: str, np: int) -> list[str]:
    return ["mpiexec", "-np", str(np)]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--profiles", default="cpu,gpu_strict")
    parser.add_argument("--case", default=None, help="restrict to a single case id")
    parser.add_argument("--manifest", type=Path, default=vc.DEFAULT_MANIFEST)
    parser.add_argument("--budget", type=Path, default=vc.DEFAULT_BUDGET)
    parser.add_argument("--runs-dir", type=Path, default=vc.DEFAULT_RUNS_DIR)
    parser.add_argument("--binary", type=Path, help="exact absolute executable (required for campaign evidence)")
    parser.add_argument("--bundle-dir", type=Path, help="exact absolute new bundle directory")
    parser.add_argument("--build-command")
    parser.add_argument("--build-flags")
    args = parser.parse_args()

    if args.binary is not None:
        if not args.binary.is_absolute() or not args.binary.is_file() or not os.access(args.binary, os.X_OK):
            parser.error("--binary must be an existing absolute executable")
    if args.bundle_dir is not None:
        if not args.bundle_dir.is_absolute():
            parser.error("--bundle-dir must be absolute")
        if args.bundle_dir.exists():
            parser.error(f"--bundle-dir already exists: {args.bundle_dir}")

    profiles = [p.strip() for p in args.profiles.split(",") if p.strip()]
    if (args.binary is None) != (args.bundle_dir is None):
        parser.error("--binary and --bundle-dir must be supplied together")
    if args.binary is not None and (len(profiles) != 1 or args.case is None):
        parser.error("exact mode requires exactly one profile and --case")
    if args.binary is not None and (not args.build_command or not args.build_flags):
        parser.error("exact mode requires --build-command and --build-flags")
    bad = [p for p in profiles if p not in vc.PROFILES]
    if bad:
        parser.error(f"unknown profile(s): {bad}; choose from {vc.PROFILES}")

    cases = vc.load_cases(args.manifest, "local", args.case)
    if not cases:
        print(f"ERROR: no local-hardware cases matched (filter={args.case!r})", file=sys.stderr)
        return 2

    return vc.run_all(
        profiles=profiles, cases=cases, arch=vc.cuda_arch(), host=vc.hostname(),
        runs_dir=args.runs_dir, budget_path=args.budget, device_tag_fn=device_tag,
        launch_builder=mpiexec_launch, label="local validation",
        exact_binary=args.binary, exact_bundle_dir=args.bundle_dir,
        manifest_extra={"build_command": args.build_command, "build_flags": args.build_flags},
        manifest_path=args.manifest,
    )


if __name__ == "__main__":
    raise SystemExit(main())
