#!/usr/bin/env python3
"""Summarize the paired Step 2 smoothing measurements."""

from __future__ import annotations

import json
import statistics
import sys
from pathlib import Path


REPO = Path("/home/jackplum/Projects/alamo")
TASK = REPO / "docs/agent_plans/20260721-fapply-runtime-optimization"
ROOT = TASK / "artifacts/a1000-sm86-2e6a8f8f-20260721/step2/smoothing_sweep"
sys.path.insert(0, str(REPO / "benchmark"))
from compare_tinyprofiler import parse  # noqa: E402


def distribution(values: list[float]) -> dict[str, object]:
    med = statistics.median(values)
    return {
        "values": values,
        "median": med,
        "mad": statistics.median(abs(value - med) for value in values),
    }


def measured_dirs(root: Path) -> list[Path]:
    return [root / f"rep{i}" for i in range(1, 6)]


def summarize_setting(root: Path) -> dict[str, object]:
    result: dict[str, object] = {}
    for mode in ("unprofiled", "profiled"):
        dirs = measured_dirs(root / mode)
        if not all((directory / "wall_seconds.txt").is_file() for directory in dirs):
            raise RuntimeError(f"incomplete measured arm: {root / mode}")
        result[f"{mode}_external_wall_seconds"] = distribution(
            [float((directory / "wall_seconds.txt").read_text()) for directory in dirs]
        )
        if mode == "profiled":
            parsed = [parse(directory / "stdout") for directory in dirs]
            for label, region in (
                ("mlmg_solve", "MLMG::solve()"),
                ("fapply", "Operator::Elastic::Fapply()"),
                ("mlmg_one_iter", "MLMG::oneIter()"),
                ("mlmg_vcycle", "MLMG::mgVcycle()"),
            ):
                if not all(region in item for item in parsed):
                    raise RuntimeError(f"missing {region} in {root}")
                result[f"{label}_inclusive_wall_seconds"] = distribution(
                    [item[region]["inclusive_wall_seconds"] for item in parsed]
                )
                result[f"{label}_calls"] = [item[region]["ncalls"] for item in parsed]
    return result


def pct_gain(reference: float, candidate: float) -> float:
    return 100.0 * (reference - candidate) / reference


def main() -> int:
    summary: dict[str, object] = {}
    for regime in ("2d_conservative", "3d_psi"):
        candidate = summarize_setting(ROOT / regime / "2x2")
        drift = summarize_setting(ROOT / regime / "drift_4x4")
        candidate["matched_gain_pct"] = {
            "unprofiled_external_wall": pct_gain(
                drift["unprofiled_external_wall_seconds"]["median"],
                candidate["unprofiled_external_wall_seconds"]["median"],
            ),
            "mlmg_solve": pct_gain(
                drift["mlmg_solve_inclusive_wall_seconds"]["median"],
                candidate["mlmg_solve_inclusive_wall_seconds"]["median"],
            ),
            "fapply_inclusive_wall": pct_gain(
                drift["fapply_inclusive_wall_seconds"]["median"],
                candidate["fapply_inclusive_wall_seconds"]["median"],
            ),
            "fapply_calls": pct_gain(
                statistics.median(drift["fapply_calls"]),
                statistics.median(candidate["fapply_calls"]),
            ),
        }
        summary[regime] = {"2x2": candidate, "drift_4x4": drift}

    out_json = ROOT.parent / "smoothing_summary.json"
    out_md = ROOT.parent / "smoothing_summary.md"
    out_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    lines = [
        "# Step 2 paired smoothing summary",
        "",
        "All values below exclude the warmup and use five measured repetitions.",
        "",
        "| regime | setting | external wall median ± MAD (s) | MLMG solve median ± MAD (s) | FApply calls | FApply wall median ± MAD (s) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for regime in ("2d_conservative", "3d_psi"):
        for setting in ("drift_4x4", "2x2"):
            arm = summary[regime][setting]
            ext = arm["unprofiled_external_wall_seconds"]
            solve = arm["mlmg_solve_inclusive_wall_seconds"]
            fapply = arm["fapply_inclusive_wall_seconds"]
            calls = statistics.median(arm["fapply_calls"])
            lines.append(
                f"| {regime} | {setting} | {ext['median']:.6g} ± {ext['mad']:.3g} | "
                f"{solve['median']:.6g} ± {solve['mad']:.3g} | {calls:g} | "
                f"{fapply['median']:.6g} ± {fapply['mad']:.3g} |"
            )
        gains = summary[regime]["2x2"]["matched_gain_pct"]
        lines.extend([
            "",
            f"- {regime} matched 2x2 gains: external {gains['unprofiled_external_wall']:.3f}%, "
            f"MLMG solve {gains['mlmg_solve']:.3f}%, FApply wall "
            f"{gains['fapply_inclusive_wall']:.3f}%, calls {gains['fapply_calls']:.3f}%.",
            "",
        ])
    out_md.write_text("\n".join(lines) + "\n")
    print(out_md.read_text(), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
