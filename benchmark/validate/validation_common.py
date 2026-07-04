"""Shared orchestration for the Phase 1 validation runners (roadmap task 1.B).

Both `run_validation_local.py` (plain `mpiexec`, this machine's `bin/`) and
the NOVA-side driver invoked by `run_validation_nova.slurm` (`srun --mpi=pmix`,
NOVA's Slurm/module environment) need the same case-running and
metrics-extraction logic -- only the launch command and binary/hardware
discovery differ. This module is that shared core; the two entry points stay
thin and environment-specific. See README.md ("Two entry points, not one
auto-detecting script") for why the entry points themselves are not merged.
"""

from __future__ import annotations

import datetime as dt
import json
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Callable

import yaml

ROOT = Path(__file__).resolve().parents[2]
VALIDATE_DIR = ROOT / "benchmark" / "validate"
DEFAULT_MANIFEST = VALIDATE_DIR / "cases.manifest.yaml"
DEFAULT_BUDGET = VALIDATE_DIR / "physics_budget.yaml"
DEFAULT_RUNS_DIR = VALIDATE_DIR / "runs"

PROFILES = ("cpu", "gpu_strict", "gpu_fast")

LaunchBuilder = Callable[[str, int], list[str]]  # (profile, np) -> command prefix, binary/args appended after


def sh(cmd: list[str]) -> str:
    return subprocess.check_output(cmd, text=True, stderr=subprocess.DEVNULL).strip()


def git_sha() -> str:
    try:
        return sh(["git", "rev-parse", "--short", "HEAD"])
    except subprocess.SubprocessError:
        return "nogit"


def hostname() -> str:
    try:
        return sh(["hostname", "-s"])
    except (subprocess.SubprocessError, FileNotFoundError):
        return "unknown"


def nvidia_smi_field(query: str) -> str:
    try:
        out = sh(["nvidia-smi", f"--query-gpu={query}", "--format=csv,noheader"])
        return out.splitlines()[0].strip() if out else ""
    except (subprocess.SubprocessError, FileNotFoundError):
        return ""


def cuda_arch() -> str:
    cc = nvidia_smi_field("compute_cap")
    return cc.replace(".", "").replace(" ", "")


def find_binary(pattern: str) -> Path | None:
    candidates = sorted(ROOT.glob(pattern))
    executables = [p for p in candidates if p.is_file()]
    return executables[-1] if executables else None


def binary_for_profile(profile: str, dim: int, arch: str) -> Path | None:
    if profile == "cpu":
        return find_binary(f"bin/alamo-{dim}d-g++")
    if profile == "gpu_strict":
        return find_binary(f"bin/alamo_gpu-{dim}d-nofast-cuda{arch}-g++")
    if profile == "gpu_fast":
        return find_binary(f"bin/alamo_gpu-{dim}d-cuda{arch}-g++")
    raise ValueError(f"unknown profile: {profile}")


def load_cases(manifest_path: Path, hardware_tag: str, case_filter: str | None) -> list[dict[str, Any]]:
    data = yaml.safe_load(manifest_path.read_text(encoding="utf-8"))
    default_overrides = data.get("defaults", {}).get("overrides", [])
    cases = []
    for c in data["cases"]:
        if hardware_tag not in c.get("hardware", []):
            continue
        if case_filter and c["id"] != case_filter:
            continue
        c = dict(c)
        c["overrides"] = list(default_overrides) + list(c.get("overrides", []))
        cases.append(c)
    return cases


def run_case(case: dict[str, Any], profile: str, binary: Path, case_dir: Path, np: int,
             launch_builder: LaunchBuilder) -> bool:
    case_dir.mkdir(parents=True, exist_ok=True)
    plot_dir = case_dir / "_plot"
    log_path = case_dir / "run.log"

    cmd = [
        *launch_builder(profile, np), str(binary), case["input"],
        f"max_step={case['max_step']}",
        *case["overrides"],
        f"plot_file={plot_dir}",
    ]
    with log_path.open("w", encoding="utf-8") as log:
        log.write(f"# command: {' '.join(cmd)}\n")
        log.flush()
        completed = subprocess.run(cmd, cwd=ROOT, stdout=log, stderr=subprocess.STDOUT, text=True)

    if completed.returncode != 0:
        print(f"ERROR: {case['id']}/{profile} failed (exit {completed.returncode}); see {log_path}",
              file=sys.stderr)
        return False

    if (plot_dir / "thermo.dat").exists():
        shutil.copy(plot_dir / "thermo.dat", case_dir / "thermo.dat")
    # ALAMO writes <step>node/<step>cell pairs directly under plot_dir; copy
    # every pair (extract_metrics.py picks the highest-numbered "node" dir as
    # the final plotfile -- see its find_final_plotfile()).
    for sub in sorted(plot_dir.glob("[0-9]*")):
        if sub.is_dir() and (sub.name.endswith("node") or sub.name.endswith("cell")):
            shutil.copytree(sub, case_dir / sub.name, dirs_exist_ok=True)
    shutil.rmtree(plot_dir, ignore_errors=True)

    region_times = case_dir / "region_times.txt"
    log_text = log_path.read_text(encoding="utf-8", errors="replace")
    if "TinyProfiler" in log_text:
        idx = log_text.index("TinyProfiler")
        region_times.write_text(log_text[idx:], encoding="utf-8")
    else:
        region_times.write_text(
            "# no TinyProfiler region summary found in run.log "
            "(profiling not enabled for this binary/run)\n",
            encoding="utf-8",
        )

    return True


def extract(case_dir: Path, budget_path: Path) -> None:
    subprocess.run(
        [sys.executable, str(VALIDATE_DIR / "extract_metrics.py"), str(case_dir), "--budget", str(budget_path)],
        check=True,
    )


def write_manifest(bundle_dir: Path, profile: str, dev_tag: str, host: str,
                    cases_run: list[dict[str, Any]], extra: dict[str, Any] | None = None) -> None:
    manifest = {
        "timestamp_utc": dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ"),
        "git_sha": git_sha(),
        "host": host,
        "device": dev_tag,
        "profile": profile,
        "gpu_name": nvidia_smi_field("name") or None,
        "cuda_compute_cap": nvidia_smi_field("compute_cap") or None,
        "driver_version": nvidia_smi_field("driver_version") or None,
        "cases": cases_run,
    }
    if extra:
        manifest.update(extra)
    (bundle_dir / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")


def default_np_resolver(case: dict[str, Any], profile: str) -> int:
    return case.get("np", {}).get("cpu" if profile == "cpu" else "gpu", 1)


def run_all(*, profiles: list[str], cases: list[dict[str, Any]], arch: str, host: str,
            runs_dir: Path, budget_path: Path, device_tag_fn: Callable[[str, str], str],
            launch_builder: LaunchBuilder, label: str, manifest_extra: dict[str, Any] | None = None,
            binary_finder: Callable[[str, int, str], Path | None] = binary_for_profile,
            np_resolver: Callable[[dict[str, Any], str], int] = default_np_resolver) -> int:
    """Shared per-profile/per-case loop. Returns process exit code (0 ok, 1 any failure).

    `binary_finder`/`np_resolver` default to the local-machine conventions
    (binary_for_profile's plain glob, manifest's `np.cpu`/`np.gpu`) but the
    NOVA driver overrides both: NOVA's build scripts use a different binary
    naming convention (see run_validation_nova.py), and rank count there
    comes from the Slurm allocation (SLURM_NTASKS), not the manifest.
    """
    timestamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    sha = git_sha()
    print(f"=== {label}: {len(cases)} case(s), profiles={profiles} ===")

    any_failure = False
    for profile in profiles:
        dev_tag = device_tag_fn(profile, arch)
        bundle_dir = runs_dir / f"{timestamp}_{sha}_{host}_{dev_tag}"
        bundle_dir.mkdir(parents=True, exist_ok=True)
        print(f"\n--- profile={profile} -> {bundle_dir} ---")

        cases_run = []
        for case in cases:
            binary = binary_finder(profile, case["dim"], arch)
            if binary is None:
                print(f"WARNING: no binary for profile={profile} dim={case['dim']}d, "
                      f"skipping case {case['id']}", file=sys.stderr)
                any_failure = True
                continue
            np = np_resolver(case, profile)
            case_dir = bundle_dir / case["id"]
            print(f"  case={case['id']} bin={binary.relative_to(ROOT)} np={np}")
            ok = run_case(case, profile, binary, case_dir, np, launch_builder)
            if not ok:
                any_failure = True
                continue
            extract(case_dir, budget_path)
            cases_run.append({
                "id": case["id"], "input": case["input"], "dim": case["dim"],
                "max_step": case["max_step"], "np": np,
            })

        write_manifest(bundle_dir, profile, dev_tag, host, cases_run, manifest_extra)
        print(f"  wrote bundle {bundle_dir}")

    print(f"\ndone. bundles under {runs_dir}/{timestamp}_{sha}_{host}_*")
    return 1 if any_failure else 0
