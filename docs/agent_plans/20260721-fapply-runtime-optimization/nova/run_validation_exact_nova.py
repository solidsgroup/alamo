#!/usr/bin/env python3
"""Run one exact validation binary inside an existing NOVA GPU allocation."""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path


ROOT = Path.cwd()
sys.path.insert(0, str(ROOT / "benchmark" / "validate"))
import validation_common as vc  # noqa: E402


def srun_launch(profile: str, np: int) -> list[str]:
    cmd = ["srun", "--mpi=pmix", "-n", str(np)]
    if profile != "cpu":
        cmd += ["--gpus-per-task=1"]
    return cmd


def device_tag(profile: str, arch: str) -> str:
    suffix = "strict" if profile == "gpu_strict" else "fast"
    return f"a100_sm{arch}_{suffix}"


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile", default="gpu_fast", choices=vc.PROFILES)
    parser.add_argument("--case", required=True)
    parser.add_argument("--manifest", type=Path, default=vc.DEFAULT_MANIFEST)
    parser.add_argument("--budget", type=Path, default=vc.DEFAULT_BUDGET)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--bundle-dir", type=Path, required=True)
    parser.add_argument("--build-command", required=True)
    parser.add_argument("--build-flags", required=True)
    args = parser.parse_args()

    binary = args.binary.resolve()
    bundle_dir = args.bundle_dir.resolve()
    if not binary.is_file() or not os.access(binary, os.X_OK):
        parser.error("--binary must be an existing executable")
    if args.bundle_dir.exists():
        parser.error("--bundle-dir must not already exist")

    cases = vc.load_cases(args.manifest, "nova", args.case)
    if len(cases) != 1:
        parser.error(f"expected exactly one NOVA case, found {len(cases)}")

    arch = vc.cuda_arch()
    if arch != "80":
        parser.error(f"expected A100 compute capability 8.0, found {arch!r}")

    np = int(os.environ.get("SLURM_NTASKS", "1"))
    case = cases[0]
    case_dir = bundle_dir / case["id"]
    argv = [
        *srun_launch(args.profile, np), str(binary), case["input"],
        f"max_step={case['max_step']}", *case["overrides"],
        f"plot_file={case_dir / '_plot'}",
    ]
    bundle_dir.mkdir(parents=True)
    case_dir.mkdir(parents=True)
    (case_dir / "command.json").write_text(
        json.dumps({"argv": argv, "overrides": case["overrides"]}, indent=2),
        encoding="utf-8",
    )
    if not vc.run_case(case, args.profile, binary, case_dir, np, srun_launch):
        return 1
    vc.extract(case_dir, args.budget)

    protected = [
        "src/Operator/Elastic.cpp",
        "src/Operator/Elastic.H",
        "src/Set/Matrix4_Major.H",
    ]
    source_diff = subprocess.check_output(
        ["git", "diff", "HEAD", "--binary", "--", *protected], cwd=ROOT
    )
    manifest_path = args.manifest.resolve()
    budget_path = args.budget.resolve()
    oracle_paths = [
        manifest_path,
        budget_path,
        ROOT / "benchmark/validate/extract_metrics.py",
        ROOT / "benchmark/validate/compare_validation.py",
        ROOT / "benchmark/validate/validation_common.py",
        Path(__file__).resolve(),
    ]
    manifest = {
        "timestamp_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ"),
        "git_sha": subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        "host": vc.hostname(),
        "device": device_tag(args.profile, arch),
        "profile": args.profile,
        "gpu_name": vc.nvidia_smi_field("name") or None,
        "cuda_compute_cap": vc.nvidia_smi_field("compute_cap") or None,
        "driver_version": vc.nvidia_smi_field("driver_version") or None,
        "gpu_uuid": vc.nvidia_smi_field("uuid") or None,
        "gpu_type": "a100",
        "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
        "binary_path": str(binary),
        "binary_sha256": sha256_file(binary),
        "scoped_source_diff_sha256": hashlib.sha256(source_diff).hexdigest(),
        "build_command": args.build_command,
        "build_flags": args.build_flags,
        "oracle_scripts": {
            str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path): sha256_file(path)
            for path in oracle_paths
        },
        "cases": [{
            "id": case["id"],
            "input": case["input"],
            "dim": case["dim"],
            "max_step": case["max_step"],
            "np": np,
            "input_sha256": sha256_file(ROOT / case["input"]),
            "overrides": case["overrides"],
            "overrides_sha256": hashlib.sha256(
                "\n".join(case["overrides"]).encode()
            ).hexdigest(),
            "executed_argv": argv,
        }],
    }
    (bundle_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8"
    )
    print(f"wrote exact NOVA bundle {bundle_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
