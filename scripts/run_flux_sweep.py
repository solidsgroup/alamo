#!/usr/bin/env python3
"""Run and compare Hydro flux schemes for the diffuse flame input."""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
from pathlib import Path
import shlex
import shutil
import subprocess
import sys
import time


RUNS: list[tuple[str, list[str]]] = [
    ("riemann", ["hydro.flux_scheme=riemann"]),
    ("advect_upwind", ["timestep=2e-8", "hydro.flux_scheme=advect", "hydro.advection.type=upwind"]),
    #("advect_centered", ["timestep=2e-8", "hydro.flux_scheme=advect", "hydro.advection.type=centered"]),
    #("advect_muscl_koren", ["timestep=2e-8", "hydro.flux_scheme=advect", "hydro.advection.type=muscl", "hydro.advection.muscl.limiter.type=koren"]),
    #("advect_muscl_mc", ["timestep=2e-8", "hydro.flux_scheme=advect", "hydro.advection.type=muscl", "hydro.advection.muscl.limiter.type=mc"]),
    ("advect_muscl_minmod", ["timestep=2e-8", "hydro.flux_scheme=advect", "hydro.advection.type=muscl", "hydro.advection.muscl.limiter.type=minmod"]),
    ("advect_muscl_superbee", ["timestep=2e-8", "hydro.flux_scheme=advect", "hydro.advection.type=muscl", "hydro.advection.muscl.limiter.type=superbee"]),
    ("advect_muscl_umist", ["timestep=2e-8", "hydro.flux_scheme=advect", "hydro.advection.type=muscl", "hydro.advection.muscl.limiter.type=umist"]),
    ("advect_muscl_vanalbada", ["timestep=2e-8", "hydro.flux_scheme=advect", "hydro.advection.type=muscl", "hydro.advection.muscl.limiter.type=vanalbada"]),
    #("advect_muscl_vanleer", ["timestep=2e-8", "hydro.flux_scheme=advect", "hydro.advection.type=muscl", "hydro.advection.muscl.limiter.type=vanleer"]),
    #("advect_quick", ["timestep=2e-8", "hydro.flux_scheme=advect", "hydro.advection.type=quick"]),
    ("advect_weno5", ["timestep=2e-8", "hydro.flux_scheme=advect", "hydro.advection.type=weno5"]),
]

RUN_LABELS = {
    "riemann": "Riemann",
    "advect_upwind": "Advect upwind",
    #"advect_centered": "Advect centered",
    #"advect_muscl_koren": "Advect MUSCL-Koren",
    #"advect_muscl_mc": "Advect MUSCL-MC",
    "advect_muscl_minmod": "Advect MUSCL-MinMod",
    "advect_muscl_superbee": "Advect MUSCL-Superbee",
    "advect_muscl_umist": "Advect MUSCL-UMIST",
    "advect_muscl_vanalbada": "Advect MUSCL-VanAlbada",
    #"advect_muscl_vanleer": "Advect MUSCL-VanLeer",
    #"advect_quick": "Advect QUICK",
    "advect_weno5": "Advect WENO5",
}

PROFILE_VARIABLES = [
    ("temperature", "Temperature"),
    ("velocity", "Speed |u|"),
    ("pressure", "Pressure"),
    ("eta", "Eta"),
]


def repo_dir() -> Path:
    return Path(__file__).resolve().parents[1]


def default_exe(root: Path) -> str:
    if (root / "bin/alamo").is_file():
        return "bin/alamo"
    return "bin/alamo-2d-6species-perf-clang++"


def timestamp() -> str:
    return time.strftime("%Y%m%d_%H%M%S")


def parse_args() -> argparse.Namespace:
    root = repo_dir()
    parser = argparse.ArgumentParser(
        description="Run ALAMO flux-scheme sweep and compare centerline profiles."
    )
    parser.add_argument("--exe", default=os.environ.get("EXE", default_exe(root)))
    parser.add_argument(
        "--input",
        default=os.environ.get("INPUT", "input.sandwich_rocfire_flame_diffuse"),
    )
    parser.add_argument(
        "--out-root",
        default=os.environ.get("OUT_ROOT", f"/tmp/alamo_flux_sweep_{timestamp()}"),
    )
    #parser.add_argument("--max-step", type=int, default=int(os.environ.get("MAX_STEP", "50")))
    #parser.add_argument("--plot-int", type=int, default=int(os.environ.get("PLOT_INT", "50")))
    parser.add_argument(
        "--no-build",
        action="store_true",
        default=os.environ.get("BUILD", "1") == "0",
        help="Skip make before running the sweep.",
    )
    return parser.parse_args()


def require_safe_out_root(out_root: Path, root: Path) -> None:
    resolved = out_root.resolve()
    unsafe = {Path("/"), Path("/tmp"), root.resolve()}
    if resolved in unsafe:
        raise SystemExit(f"Refusing unsafe OUT_ROOT={out_root}")


def maybe_build(root: Path, exe: str, no_build: bool) -> None:
    if no_build:
        return
    if exe.startswith("bin/"):
        print(f"Building {exe}")
        subprocess.run(["make", exe], cwd=root, check=True)
    else:
        print(f"Skipping make for non-repo executable {exe}")


def run_sweep(root: Path, args: argparse.Namespace) -> list[dict[str, str]]:
    exe = root / args.exe if not Path(args.exe).is_absolute() else Path(args.exe)
    input_file = root / args.input if not Path(args.input).is_absolute() else Path(args.input)
    out_root = Path(args.out_root)
    logs_dir = out_root / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    if not input_file.is_file():
        raise SystemExit(f"Input file not found: {input_file}")
    if not exe.is_file() or not os.access(exe, os.X_OK):
        raise SystemExit(f"Executable not found or not executable: {exe}")

    summary_rows: list[dict[str, str]] = []
    summary_csv = out_root / "run_summary.csv"

    print(f"Writing sweep outputs to {out_root}")
    print()

    for name, extra_args in RUNS:
        plot_file = out_root / name
        log_file = logs_dir / f"{name}.log"
        if plot_file.exists():
            shutil.rmtree(plot_file)

        cmd = [
            str(exe),
            str(input_file),
            f"plot_file={plot_file}",
            *extra_args,
        ]

        print(f"==> {name}")
        print(f"    {shlex.join(cmd)}")

        start = time.perf_counter()
        with log_file.open("w", encoding="utf-8") as log:
            proc = subprocess.run(
                cmd,
                cwd=root,
                stdout=log,
                stderr=subprocess.STDOUT,
                check=False,
            )
        runtime = time.perf_counter() - start

        row = {
            "name": name,
            "status": str(proc.returncode),
            "runtime_s": f"{runtime:.6f}",
            "plot_file": str(plot_file),
            "log": str(log_file),
            "args": " ".join(extra_args),
        }
        summary_rows.append(row)
        print(f"    status={proc.returncode} runtime={runtime:.2f}s log={log_file}")

    with summary_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=["name", "status", "runtime_s", "plot_file", "log", "args"]
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    return summary_rows


def latest_plotfile(plot_root: Path) -> Path:
    visit = plot_root / "celloutput.visit"
    if visit.exists():
        entries: list[Path] = []
        for raw in visit.read_text(encoding="utf-8").splitlines():
            line = raw.strip()
            if not line:
                continue
            if line.endswith("/Header"):
                line = line[: -len("/Header")]
            path = plot_root / line
            if path.exists():
                entries.append(path)
        if entries:
            return entries[-1]

    candidates = sorted(Path(p) for p in glob.glob(str(plot_root / "*cell")))
    candidates += sorted(Path(p) for p in glob.glob(str(plot_root / "*cell.h5")))
    if not candidates:
        raise FileNotFoundError(f"No cell plotfile found under {plot_root}")
    return candidates[-1]


def field_ref(ds, *names: str):
    candidates = []
    for name in names:
        candidates.extend([("boxlib", name), ("gas", name), ("stream", name), name])
    for candidate in candidates:
        try:
            ds._get_field_info(candidate)
            return candidate
        except Exception:
            pass
    available = ", ".join(str(field) for field in ds.field_list)
    raise KeyError(f"Fields {names!r} not found. Available fields: {available}")


def dataframe_column_name(field) -> str:
    if isinstance(field, tuple):
        return field[-1]
    return field


def extract_profile(name: str, plot_root: Path, profiles_dir: Path):
    import numpy as np
    import yt

    yt.set_log_level(50)

    plotfile = latest_plotfile(plot_root)
    ds = yt.load(str(plotfile))
    dims = [int(v) for v in ds.domain_dimensions]
    dim = sum(1 for value in dims if value > 1)
    geom_lo = [float(value) for value in ds.domain_left_edge.d]
    geom_hi = [float(value) for value in ds.domain_right_edge.d]

    x_center = 0.5 * (geom_lo[0] + geom_hi[0])
    z_center = 0.0
    if len(geom_lo) > 2:
        z_center = 0.5 * (geom_lo[2] + geom_hi[2])

    fields = {
        "temperature": field_ref(ds, "temperature", "temp"),
        "velocityx": field_ref(ds, "velocityx"),
        "velocityy": field_ref(ds, "velocityy"),
        "pressure": field_ref(ds, "pressure"),
        "eta": field_ref(ds, "eta"),
    }
    if dim == 3:
        fields["velocityz"] = field_ref(ds, "velocityz")

    start = [x_center, geom_lo[1], z_center]
    end = [x_center, geom_hi[1], z_center]
    ray = ds.ray(start, end)
    requested = [("gas", "x"), ("gas", "y"), ("gas", "z"), *fields.values()]
    df = ray.to_dataframe(requested).reset_index(drop=True)

    rename = {
        dataframe_column_name(ref): logical
        for logical, ref in fields.items()
    }
    df = df.rename(columns=rename)

    cols = ["x", "y", "z", "temperature", "velocityx", "velocityy", "pressure", "eta"]
    if "velocityz" in df:
        cols.insert(cols.index("pressure"), "velocityz")
        df["velocity"] = np.sqrt(
            df["velocityx"] ** 2 + df["velocityy"] ** 2 + df["velocityz"] ** 2
        )
    else:
        df["velocity"] = np.sqrt(df["velocityx"] ** 2 + df["velocityy"] ** 2)
    cols.insert(cols.index("pressure"), "velocity")

    df = df[cols].sort_values("y").reset_index(drop=True)
    df = df.groupby("y", as_index=False).mean(numeric_only=True)

    csv_path = profiles_dir / f"{name}.csv"
    df.to_csv(csv_path, index=False)
    return df, plotfile, csv_path


def compare_profile(base, prof, var: str) -> tuple[float, float, float]:
    import numpy as np

    y_base = base["y"].to_numpy()
    base_v = base[var].to_numpy()
    y = prof["y"].to_numpy()
    values = prof[var].to_numpy()

    y_min = max(np.min(y_base), np.min(y))
    y_max = min(np.max(y_base), np.max(y))
    mask = (y_base >= y_min) & (y_base <= y_max)
    y_common = y_base[mask]
    if y_common.size == 0:
        return math.nan, math.nan, math.nan

    base_common = base_v[mask]
    values_common = np.interp(y_common, y, values)
    diff = values_common - base_common
    max_abs = float(np.max(np.abs(diff)))
    rms = float(np.sqrt(np.mean(diff * diff)))
    scale = float(np.max(np.abs(base_common)))
    rel_linf = max_abs / scale if scale > 0.0 else math.nan
    return max_abs, rms, rel_linf


def fmt_float(value, precision: int = 4) -> str:
    if value == "" or value is None:
        return ""
    try:
        float_value = float(value)
    except Exception:
        return str(value)
    if not math.isfinite(float_value):
        return str(float_value)
    return f"{float_value:.{precision}g}"


def plot_centerline_profiles(out_root: Path, profiles: dict[str, object]) -> list[Path]:
    plot_dir = out_root / "plots"
    plot_dir.mkdir(exist_ok=True)
    mpl_config = out_root / "matplotlib"
    mpl_config.mkdir(exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_config))

    import matplotlib

    matplotlib.use("Agg", force=True)
    import matplotlib.pyplot as plt

    created: list[Path] = []
    for var, ylabel in PROFILE_VARIABLES:
        fig, ax = plt.subplots(figsize=(8.0, 5.0), constrained_layout=True)
        plotted = 0
        for name, _ in RUNS:
            profile = profiles.get(name)
            if profile is None or var not in profile.columns:
                continue
            x = profile["y"].to_numpy() * 1.0e6
            y = profile[var].to_numpy()
            if name == "riemann":
                ax.plot(x, y, label=RUN_LABELS.get(name, name), color="black", linewidth=2.4)
                ax.set_ylim([None,max(y)*1.1])
            else:
                ax.plot(x, y, label=RUN_LABELS.get(name, name), linewidth=1.6)
            plotted += 1

        if plotted == 0:
            plt.close(fig)
            continue

        ax.set_xlabel("y (micron)")
        ax.set_ylabel(ylabel)
        ax.set_title(f"Centerline {ylabel} after sweep final step")
        ax.grid(True, alpha=0.3, linewidth=0.7)
        ax.legend(loc="best", fontsize=9, frameon=True)

        png = plot_dir / f"centerline_{var}.png"
        pdf = plot_dir / f"centerline_{var}.pdf"
        fig.savefig(png, dpi=200)
        fig.savefig(pdf)
        plt.close(fig)
        created.extend([png, pdf])

    return created


def analyze(out_root: Path, rows: list[dict[str, str]]) -> None:
    import pandas as pd

    profiles_dir = out_root / "profiles"
    profiles_dir.mkdir(exist_ok=True)

    profiles = {}
    for row in rows:
        if row["status"] != "0":
            continue
        name = row["name"]
        try:
            profile, plotfile, profile_csv = extract_profile(
                name, Path(row["plot_file"]), profiles_dir
            )
            row["plotfile_used"] = str(plotfile)
            row["profile_csv"] = str(profile_csv)
            profiles[name] = profile
        except Exception as exc:
            row["status"] = f"profile_error: {exc}"

    if "riemann" not in profiles:
        raise RuntimeError("Riemann baseline profile was not available")

    baseline = profiles["riemann"]
    base_runtime = float(next(row["runtime_s"] for row in rows if row["name"] == "riemann"))
    vars_to_compare = [name for name, _ in PROFILE_VARIABLES]
    comparison_rows = []

    for row in rows:
        name = row["name"]
        comp = {
            "name": name,
            "status": row["status"],
            "runtime_s": row["runtime_s"],
            "runtime_ratio": "",
            "profile_csv": row.get("profile_csv", ""),
        }
        if row["status"] == "0" and name in profiles:
            runtime = float(row["runtime_s"])
            comp["runtime_ratio"] = runtime / base_runtime if base_runtime > 0.0 else math.nan
            for var in vars_to_compare:
                max_abs, rms, rel_linf = compare_profile(baseline, profiles[name], var)
                comp[f"{var}_max_abs"] = max_abs
                comp[f"{var}_rms"] = rms
                comp[f"{var}_rel_linf"] = rel_linf
        comparison_rows.append(comp)

    comparison_csv = out_root / "comparison.csv"
    pd.DataFrame(comparison_rows).to_csv(comparison_csv, index=False)
    plot_paths = plot_centerline_profiles(out_root, profiles)

    summary_md = out_root / "summary.md"
    with summary_md.open("w", encoding="utf-8") as f:
        f.write("# Flux Sweep Summary\n\n")
        f.write(f"- Output root: `{out_root}`\n")
        f.write("- Baseline: `riemann` (`hydro.flux_scheme=riemann`)\n")
        f.write("- Centerline: constant domain-center `x`, sampled along `y`\n")
        f.write("- Velocity comparison uses `sqrt(velocityx^2 + velocityy^2)` for 2D runs.\n\n")
        f.write("| Flux model | Status | Runtime (s) | Runtime / Riemann | max abs dT | max abs d_speed | max abs dP | max abs deta |\n")
        f.write("|---|---:|---:|---:|---:|---:|---:|---:|\n")
        for row in comparison_rows:
            f.write(
                "| {name} | {status} | {runtime} | {ratio} | {dt} | {du} | {dp} | {de} |\n".format(
                    name=row["name"],
                    status=row["status"],
                    runtime=fmt_float(row["runtime_s"]),
                    ratio=fmt_float(row.get("runtime_ratio", "")),
                    dt=fmt_float(row.get("temperature_max_abs", "")),
                    du=fmt_float(row.get("velocity_max_abs", "")),
                    dp=fmt_float(row.get("pressure_max_abs", "")),
                    de=fmt_float(row.get("eta_max_abs", "")),
                )
            )
        f.write("\n")
        f.write(f"Detailed metrics: `{comparison_csv}`\n\n")
        f.write("Profile CSVs:\n")
        for row in comparison_rows:
            if row.get("profile_csv"):
                f.write(f"- `{row['name']}`: `{row['profile_csv']}`\n")
        if plot_paths:
            f.write("\nPlots:\n")
            for path in plot_paths:
                f.write(f"- `{path}`\n")

    print(f"Wrote {comparison_csv}")
    if plot_paths:
        print(f"Wrote plots to {out_root / 'plots'}")
    print(f"Wrote {summary_md}")
    print()
    print(summary_md.read_text(encoding="utf-8"))


def main() -> int:
    root = repo_dir()
    args = parse_args()
    out_root = Path(args.out_root)
    require_safe_out_root(out_root, root)
    maybe_build(root, args.exe, args.no_build)
    rows = run_sweep(root, args)
    print()
    print("Extracting centerline profiles and computing baseline comparisons")
    analyze(out_root, rows)
    print()
    print(f"Done. Summary: {out_root / 'summary.md'}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise SystemExit(130)
