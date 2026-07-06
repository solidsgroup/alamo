#!/usr/bin/env python3
"""Standardized Flame CPU/GPU baseline and regression suite."""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_REF_DIR = ROOT / "benchmark" / "baseline_references"
DEFAULT_RUN_DIR = ROOT / "benchmark" / "baseline_runs"


@dataclass(frozen=True)
class Case:
    case_id: str
    input_file: str
    max_step: int
    abs_tol: float = 1.0e-8
    rel_tol: float = 1.0e-6
    overrides: tuple[str, ...] = field(default_factory=tuple)


CASES: tuple[Case, ...] = (
    Case(
        case_id="canonical_step1",
        input_file="input",
        max_step=1,
    ),
    Case(
        case_id="canonical_step2",
        input_file="input",
        max_step=2,
        rel_tol=1.0e-5,
    ),
    Case(
        case_id="eta_expression_step1",
        input_file="input",
        max_step=1,
        rel_tol=1.0e-5,
        overrides=(
            "allow_unused=1",
            "pf.eta.ic.type=expression",
            "pf.eta.ic.expression.region0=0.5 + 0.5*tanh((x-0.0877_m)/0.005_m)",
        ),
    ),
    # Worst-case high-contrast elastic regression: full 2D rod-and-tube (a
    # near-floating stiff rod coupled to a stiff tube only through a thin soft
    # void seam) exercises the persistent-MLMG resync path (SyncCoefficients on
    # Newton relinearization). interval=1 fires the elastic solve at step 2;
    # thermo.dat's disp_*/trac_* columns capture the solve. A resync/MLMG-recipe
    # regression trips this via divergence->abort or changed boundary tractions.
    Case(
        case_id="rod_and_tube_step2",
        input_file="input_rod_and_tube_2d",
        max_step=2,
        rel_tol=1.0e-5,
        overrides=("elastic.interval=1",),
    ),
)


PROFILES: tuple[str, ...] = ("cpu", "gpu_fast", "gpu_strict")


def find_binary(patterns: list[str]) -> str | None:
    candidates: list[Path] = []
    for pattern in patterns:
        candidates.extend(ROOT.glob(pattern))
    executable = [path for path in candidates if path.is_file() and os.access(path, os.X_OK)]
    if not executable:
        return None
    # Newest mtime, not lexicographic max: a name-sorted pick can prefer a
    # stale binary (e.g. alamo-2d-perf-clang++) over the one just built.
    newest = max(executable, key=lambda path: path.stat().st_mtime)
    return str(newest.relative_to(ROOT))


def local_cuda_arch() -> str | None:
    if "CUDA_ARCH" in os.environ:
        return os.environ["CUDA_ARCH"]
    if shutil.which("nvidia-smi") is None:
        return None
    try:
        output = subprocess.check_output(
            ["nvidia-smi", "--query-gpu=compute_cap", "--format=csv,noheader"],
            cwd=ROOT,
            text=True,
            stderr=subprocess.DEVNULL,
        )
    except subprocess.SubprocessError:
        return None
    first = output.splitlines()[0] if output.splitlines() else ""
    return first.replace(".", "").replace(" ", "") or None


def default_binaries() -> dict[str, str | None]:
    arch = local_cuda_arch()
    gpu_fast_patterns = [f"bin/alamo_gpu-2d-cuda{arch}-*"] if arch else []
    gpu_fast_patterns.append("bin/alamo_gpu-2d-cuda*-*")
    gpu_strict_patterns = [f"bin/alamo_gpu-2d-nofast-cuda{arch}-*"] if arch else []
    gpu_strict_patterns.append("bin/alamo_gpu-2d-nofast-cuda*-*")
    return {
        "cpu": os.environ.get("CPU_BIN") or find_binary(["bin/alamo-2d-*"]),
        "gpu_fast": os.environ.get("GPU_FAST_BIN") or find_binary(gpu_fast_patterns),
        "gpu_strict": os.environ.get("GPU_STRICT_BIN") or find_binary(gpu_strict_patterns),
    }


def read_thermo(path: Path) -> tuple[list[str], list[list[float]]]:
    with path.open("r", encoding="utf-8") as stream:
        header = stream.readline().split()
        rows = [[float(value) for value in line.split()] for line in stream if line.strip()]
    if not header or not rows:
        raise RuntimeError(f"{path} is missing thermo data")
    return header, rows


def compare_rows(
    header: list[str],
    reference: list[list[float]],
    candidate: list[list[float]],
    abs_tol: float,
    rel_tol: float,
) -> tuple[bool, list[dict[str, float | str]]]:
    if len(reference) != len(candidate):
        return False, [{"column": "<rows>", "max_abs": math.inf, "max_rel": math.inf, "status": "FAIL"}]

    results: list[dict[str, float | str]] = []
    ok_all = True
    for column, name in enumerate(header):
        max_abs = 0.0
        max_rel = 0.0
        ok = True
        for ref_row, cand_row in zip(reference, candidate):
            ref = ref_row[column]
            cand = cand_row[column]
            if not (math.isfinite(ref) and math.isfinite(cand)):
                ok = False
                max_abs = math.inf
                max_rel = math.inf
                break
            abs_err = abs(cand - ref)
            rel_err = abs_err / max(abs(ref), abs(cand), 1.0)
            max_abs = max(max_abs, abs_err)
            max_rel = max(max_rel, rel_err)
        ok = ok and (max_abs <= abs_tol or max_rel <= rel_tol)
        ok_all = ok_all and ok
        results.append(
            {
                "column": name,
                "max_abs": max_abs,
                "max_rel": max_rel,
                "status": "ok" if ok else "FAIL",
            }
        )
    return ok_all, results


def profile_np(profile: str) -> str:
    env_name = {
        "cpu": "CPU_NP",
        "gpu_fast": "GPU_FAST_NP",
        "gpu_strict": "GPU_STRICT_NP",
    }[profile]
    return os.environ.get(env_name) or os.environ.get("NP", "1")


def run_case(case: Case, profile: str, binary: str, run_dir: Path) -> dict[str, object]:
    profile_dir = run_dir / case.case_id / profile
    profile_dir.mkdir(parents=True, exist_ok=True)
    plot_dir = profile_dir / "plot"
    log_file = profile_dir / "run.log"

    cmd = [
        "mpiexec",
        "-np",
        profile_np(profile),
        str(ROOT / binary),
        case.input_file,
        f"max_step={case.max_step}",
        "stop_time=1e99_s",
        "amr.plot_int=-1",
        "amr.thermo.plot_int=1",
        "amr.thermo.int=1",
        "elastic.solver.verbose=0",
        "elastic.print_model=0",
        *case.overrides,
        f"plot_file={plot_dir}",
    ]

    started = time.monotonic()
    with log_file.open("w", encoding="utf-8") as log:
        completed = subprocess.run(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, text=True)
    elapsed = time.monotonic() - started
    if completed.returncode != 0:
        raise RuntimeError(f"{case.case_id}/{profile} failed with {completed.returncode}; see {log_file}")

    header, rows = read_thermo(plot_dir / "thermo.dat")
    return {
        "case": case.case_id,
        "profile": profile,
        "input": case.input_file,
        "max_step": case.max_step,
        "binary": binary,
        "mpi_np": int(profile_np(profile)),
        "command": cmd,
        "elapsed_s": elapsed,
        "abs_tol": case.abs_tol,
        "rel_tol": case.rel_tol,
        "overrides": list(case.overrides),
        "header": header,
        "rows": rows,
        "final": dict(zip(header, rows[-1])),
    }


def write_json(path: Path, data: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def load_json(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def print_compare(case: Case, left_name: str, left: dict[str, object], right_name: str, right: dict[str, object]) -> bool:
    header = left["header"]
    if header != right["header"]:
        print(f"{case.case_id}: header mismatch {left_name} vs {right_name}")
        return False
    ok, columns = compare_rows(
        header,
        left["rows"],
        right["rows"],
        case.abs_tol,
        case.rel_tol,
    )
    status = "ok" if ok else "FAIL"
    print(f"{case.case_id}: {right_name} vs {left_name}: {status}")
    for item in columns:
        if item["status"] != "ok":
            print(
                f"  {item['column']}: abs={item['max_abs']:.6e} "
                f"rel={item['max_rel']:.6e} {item['status']}"
            )
    return ok


def compare_records(case: Case, reference: dict[str, object], candidate: dict[str, object]) -> tuple[bool, list[dict[str, float | str]]]:
    if reference["header"] != candidate["header"]:
        return False, [{"column": "<header>", "max_abs": math.inf, "max_rel": math.inf, "status": "FAIL"}]
    return compare_rows(
        reference["header"],
        reference["rows"],
        candidate["rows"],
        case.abs_tol,
        case.rel_tol,
    )


def report(cases: list[Case], profiles: list[str], ref_dir: Path) -> int:
    ok_all = True
    for case in cases:
        records: dict[str, dict[str, object]] = {}
        for profile in profiles:
            path = ref_dir / case.case_id / f"{profile}.json"
            if not path.exists():
                raise SystemExit(f"Missing reference: {path}")
            records[profile] = load_json(path)

        cpu = records.get("cpu")
        print(f"\n{case.case_id}")
        print("profile      np  elapsed_s  vs_cpu   correctness")
        print("-----------  --  ---------  -------  -----------")
        for profile in profiles:
            record = records[profile]
            elapsed = float(record["elapsed_s"])
            if profile == "cpu" or cpu is None:
                speed = "1.000x"
                status = "reference"
            else:
                ok, _ = compare_records(case, cpu, record)
                ok_all = ok_all and ok
                speed = f"{elapsed / float(cpu['elapsed_s']):.3f}x"
                status = "ok" if ok else "FAIL"
            print(
                f"{profile:<11}  {int(record.get('mpi_np', 1)):>2}  "
                f"{elapsed:>9.3f}  {speed:>7}  {status}"
            )
    return 0 if ok_all else 1


def selected_cases(names: list[str]) -> list[Case]:
    if not names:
        return list(CASES)
    by_id = {case.case_id: case for case in CASES}
    missing = [name for name in names if name not in by_id]
    if missing:
        raise SystemExit(f"Unknown case(s): {', '.join(missing)}")
    return [by_id[name] for name in names]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=["list", "record", "check", "report", "unit"])
    parser.add_argument("--case", action="append", default=[])
    parser.add_argument("--profiles", default=",".join(PROFILES))
    parser.add_argument("--ref-dir", type=Path, default=DEFAULT_REF_DIR)
    parser.add_argument("--run-dir", type=Path, default=DEFAULT_RUN_DIR)
    args = parser.parse_args()

    if args.mode == "unit":
        sample_header = ["time", "value"]
        reference = [[0.0, 1.0], [1.0, 2.0]]
        candidate = [[0.0, 1.0 + 1.0e-12], [1.0, 2.0 - 1.0e-12]]
        ok, _ = compare_rows(sample_header, reference, candidate, 1.0e-8, 1.0e-8)
        if not ok:
            print("compare_rows unit test failed", file=sys.stderr)
            return 1
        print("baseline suite unit tests passed")
        return 0

    cases = selected_cases(args.case)
    profiles = [profile.strip() for profile in args.profiles.split(",") if profile.strip()]
    bad_profiles = [profile for profile in profiles if profile not in PROFILES]
    if bad_profiles:
        raise SystemExit(f"Unknown profile(s): {', '.join(bad_profiles)}")

    if args.mode == "list":
        print("cases:")
        for case in cases:
            print(f"  {case.case_id}: input={case.input_file} max_step={case.max_step}")
        print("profiles:")
        for profile, binary in default_binaries().items():
            print(f"  {profile}: np={profile_np(profile)} binary={binary or '<missing>'}")
        return 0

    if args.mode == "report":
        return report(cases, profiles, args.ref_dir)

    binaries = default_binaries()
    missing = [profile for profile in profiles if not binaries.get(profile)]
    if missing:
        raise SystemExit(f"Missing binaries for profiles: {', '.join(missing)}")

    ok_all = True
    if args.mode == "record":
        for case in cases:
            results: dict[str, dict[str, object]] = {}
            for profile in profiles:
                print(f"running {case.case_id}/{profile}")
                result = run_case(case, profile, binaries[profile] or "", args.run_dir)
                results[profile] = result
                write_json(args.ref_dir / case.case_id / f"{profile}.json", result)

            if "cpu" in results:
                for profile, result in results.items():
                    if profile == "cpu":
                        continue
                    ok_all = print_compare(case, "cpu", results["cpu"], profile, result) and ok_all
        return 0 if ok_all else 1

    if args.mode == "check":
        for case in cases:
            for profile in profiles:
                ref_path = args.ref_dir / case.case_id / f"{profile}.json"
                if not ref_path.exists():
                    raise SystemExit(f"Missing reference: {ref_path}")
                print(f"checking {case.case_id}/{profile}")
                fresh = run_case(case, profile, binaries[profile] or "", args.run_dir)
                reference = load_json(ref_path)
                ok_all = print_compare(case, "reference", reference, "fresh", fresh) and ok_all
        return 0 if ok_all else 1

    return 2


if __name__ == "__main__":
    raise SystemExit(main())
