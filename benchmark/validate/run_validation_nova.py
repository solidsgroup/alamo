#!/usr/bin/env python3
"""Phase 1 task 1.B -- NOVA validation driver.

Invoked from inside a Slurm allocation by run_validation_nova.slurm (NOT run
directly -- it needs SLURM_NTASKS and a GPU gres allocation already in place
for GPU profiles). See benchmark/validate/README.md ("Two entry points, not
one auto-detecting script") for why this is a separate path from
run_validation_local.py; the two share their case-running/extraction core via
validation_common.py.

KNOWN GAP, read before using --profiles gpu_strict here: NOVA's existing
build scripts (build_alamo_nova.sh / build_alamo_nova_3d.sh) only ever build
`--profile` targets, which configure's own help text says are "fully
optimized (--use_fast_math etc. on GPU)" -- i.e. FAST-MATH builds. There is
currently no NOVA build step that passes `--cuda-fp strict` (the flag that
produces the `-nofast-` binaries this script looks for). Until that's added
to the NOVA build scripts, `--profiles gpu_strict` here will fail with a
clear "binary not found" error telling you the exact configure/make command
that would produce it (see binary_for_profile_nova below) -- it does NOT
silently fall back to the fast-math binary, because v3 guiding principle 3 is
explicit that correctness claims use the strict binary only.

Usage (normally via the .slurm wrapper, not directly):
    python3 benchmark/validate/run_validation_nova.py --profile gpu_strict --gpu-type a100
    python3 benchmark/validate/run_validation_nova.py --profile cpu
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import validation_common as vc

GPU_ARCH = {"v100": "70", "a100": "80", "h200": "90"}


def device_tag(profile: str, gpu_type: str, arch: str) -> str:
    if profile == "cpu":
        return "cpu"
    suffix = "strict" if profile == "gpu_strict" else "fast"
    return f"{gpu_type}_sm{arch}_{suffix}"


def srun_launch(profile: str, np: int) -> list[str]:
    cmd = ["srun", "--mpi=pmix", "-n", str(np)]
    if profile != "cpu":
        cmd += ["--gpus-per-task=1"]
    return cmd


def binary_for_profile_nova(profile: str, dim: int, arch: str) -> Path | None:
    if profile == "cpu":
        # build_alamo_nova[_3d].sh also builds a CPU baseline tagged -profile-
        return vc.find_binary(f"bin/alamo-{dim}d-profile-g++") or vc.find_binary(f"bin/alamo-{dim}d-g++")
    if profile == "gpu_fast":
        return vc.find_binary(f"bin/alamo_gpu-{dim}d-profile-cuda{arch}-g++")
    if profile == "gpu_strict":
        # See module docstring KNOWN GAP. Try both plausible postfix orderings
        # from configure's naming logic (profile flag is appended before the
        # cuda_fp-strict "-nofast" suffix, so "profile-nofast" is the form a
        # `--profile --cuda-fp strict` build would actually produce).
        return (vc.find_binary(f"bin/alamo_gpu-{dim}d-profile-nofast-cuda{arch}-g++")
                or vc.find_binary(f"bin/alamo_gpu-{dim}d-nofast-cuda{arch}-g++"))
    raise ValueError(f"unknown profile: {profile}")


def slurm_np_resolver(case, profile: str) -> int:  # noqa: ANN001 (dict-shaped case)
    env_np = os.environ.get("SLURM_NTASKS")
    if env_np:
        return int(env_np)
    return vc.default_np_resolver(case, profile)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--profile", required=True, choices=vc.PROFILES,
                         help="single profile per Slurm submission (resource request differs per profile -- "
                              "see run_validation_nova.slurm)")
    parser.add_argument("--gpu-type", default=os.environ.get("GPU_TYPE", "a100"), choices=sorted(GPU_ARCH))
    parser.add_argument("--case", default=None, help="restrict to a single case id")
    parser.add_argument("--manifest", type=Path, default=vc.DEFAULT_MANIFEST)
    parser.add_argument("--budget", type=Path, default=vc.DEFAULT_BUDGET)
    parser.add_argument("--runs-dir", type=Path, default=vc.DEFAULT_RUNS_DIR)
    args = parser.parse_args()

    arch = GPU_ARCH[args.gpu_type]
    cases = vc.load_cases(args.manifest, "nova", args.case)
    if not cases:
        print(f"ERROR: no nova-hardware cases matched (filter={args.case!r})", file=sys.stderr)
        return 2

    if args.profile == "gpu_strict":
        missing = [c for c in cases if binary_for_profile_nova("gpu_strict", c["dim"], arch) is None]
        if missing:
            dims = sorted({c["dim"] for c in missing})
            print("ERROR: no strict/no-fast-math NOVA binary found for dim(s) "
                  f"{dims} at sm_{arch}. NOVA's build scripts do not currently produce one "
                  "(see this script's module docstring). To build it:\n"
                  f"    ./configure --comp=g++ --dim <dim> --cuda {arch} --cuda-fp strict --profile --get-eigen\n"
                  "    make -j64 bin/alamo_gpu\n"
                  f"  -> expected: bin/alamo_gpu-<dim>d-profile-nofast-cuda{arch}-g++",
                  file=sys.stderr)
            return 2

    return vc.run_all(
        profiles=[args.profile], cases=cases, arch=arch, host=vc.hostname(),
        runs_dir=args.runs_dir, budget_path=args.budget,
        device_tag_fn=lambda profile, a: device_tag(profile, args.gpu_type, a),
        launch_builder=srun_launch, label=f"nova validation ({args.gpu_type})",
        manifest_extra={"gpu_type": args.gpu_type, "slurm_job_id": os.environ.get("SLURM_JOB_ID")},
        binary_finder=binary_for_profile_nova, np_resolver=slurm_np_resolver,
    )


if __name__ == "__main__":
    raise SystemExit(main())
