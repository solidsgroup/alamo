#!/usr/bin/env python3
"""Summarize the isolated Step 3 baseline/candidate timing matrix."""

from __future__ import annotations

import json
import statistics
import sys
from pathlib import Path


REPO = Path("/home/jackplum/Projects/alamo")
TASK = REPO / "docs/agent_plans/20260721-fapply-runtime-optimization"
ROOT = TASK / "artifacts/a1000-sm86-2e6a8f8f-20260721/candidate-step3/timing"
sys.path.insert(0, str(REPO / "benchmark"))
from compare_tinyprofiler import parse  # noqa: E402


def distribution(values: list[float]) -> dict[str, object]:
    median = statistics.median(values)
    return {
        "values": values,
        "median": median,
        "mad": statistics.median(abs(value - median) for value in values),
    }


def measured_dirs(root: Path) -> list[Path]:
    return [root / f"rep{i}" for i in range(1, 6)]


def summarize_arm(root: Path) -> dict[str, object]:
    result: dict[str, object] = {}
    for mode in ("unprofiled", "profiled"):
        dirs = measured_dirs(root / mode)
        complete = all((directory / "wall_seconds.txt").is_file() for directory in dirs)
        if mode == "profiled" and not complete:
            result["profiled_complete"] = False
            continue
        if not complete:
            raise RuntimeError(f"incomplete measured arm: {root / mode}")
        result[f"{mode}_external_wall_seconds"] = distribution(
            [float((directory / "wall_seconds.txt").read_text()) for directory in dirs]
        )
        if mode != "profiled":
            continue
        result["profiled_complete"] = True
        parsed = [parse(directory / "stdout") for directory in dirs]
        for label, region in (
            ("mlmg_solve", "MLMG::solve()"),
            ("fapply", "Operator::Elastic::Fapply()"),
        ):
            if not all(region in item for item in parsed):
                raise RuntimeError(f"missing {region} in {root}")
            result[f"{label}_inclusive_wall_seconds"] = distribution(
                [item[region]["inclusive_wall_seconds"] for item in parsed]
            )
            result[f"{label}_calls"] = [item[region]["ncalls"] for item in parsed]
        result["fapply_seconds_per_call"] = distribution([
            item["Operator::Elastic::Fapply()"]["inclusive_wall_seconds"]
            / item["Operator::Elastic::Fapply()"]["ncalls"]
            for item in parsed
        ])
    return result


def gain(reference: float, candidate: float) -> float:
    return 100.0 * (reference - candidate) / reference


def main() -> int:
    summary: dict[str, object] = {}
    for regime in ("2d_conservative_2x2", "3d_psi_2x2"):
        baseline = summarize_arm(ROOT / regime / "baseline")
        candidate = summarize_arm(ROOT / regime / "candidate")
        gains = {
            "unprofiled_external_wall": gain(
                baseline["unprofiled_external_wall_seconds"]["median"],
                candidate["unprofiled_external_wall_seconds"]["median"],
            ),
        }
        if baseline["profiled_complete"] and candidate["profiled_complete"]:
            gains.update({
                "mlmg_solve": gain(
                    baseline["mlmg_solve_inclusive_wall_seconds"]["median"],
                    candidate["mlmg_solve_inclusive_wall_seconds"]["median"],
                ),
                "fapply_inclusive_wall": gain(
                    baseline["fapply_inclusive_wall_seconds"]["median"],
                    candidate["fapply_inclusive_wall_seconds"]["median"],
                ),
                "fapply_seconds_per_call": gain(
                    baseline["fapply_seconds_per_call"]["median"],
                    candidate["fapply_seconds_per_call"]["median"],
                ),
            })
        summary[regime] = {
            "baseline": baseline,
            "candidate": candidate,
            "candidate_gain_pct": gains,
        }

    out_json = ROOT.parent / "timing_summary.json"
    out_md = ROOT.parent / "timing_summary.md"
    out_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    lines = [
        "# Step 3 isolated timing summary",
        "",
        "Warmups are excluded; each entry is the median ± MAD of five runs.",
        "",
        "| regime | arm | external wall (s) | MLMG solve (s) | FApply calls | FApply wall (s) | FApply per call (µs) |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for regime in ("2d_conservative_2x2", "3d_psi_2x2"):
        for arm in ("baseline", "candidate"):
            data = summary[regime][arm]
            external = data["unprofiled_external_wall_seconds"]
            if data["profiled_complete"]:
                solve = data["mlmg_solve_inclusive_wall_seconds"]
                fapply = data["fapply_inclusive_wall_seconds"]
                per_call = data["fapply_seconds_per_call"]
                calls = statistics.median(data["fapply_calls"])
                profiled = (
                    f"{solve['median']:.6g} ± {solve['mad']:.3g} | {calls:g} | "
                    f"{fapply['median']:.6g} ± {fapply['mad']:.3g} | "
                    f"{1e6 * per_call['median']:.6g} ± {1e6 * per_call['mad']:.3g}"
                )
            else:
                profiled = "not run (fail-fast) | — | — | —"
            lines.append(
                f"| {regime} | {arm} | {external['median']:.6g} ± {external['mad']:.3g} | "
                f"{profiled} |"
            )
        gains = summary[regime]["candidate_gain_pct"]
        detail = f"external {gains['unprofiled_external_wall']:.3f}%"
        if "mlmg_solve" in gains:
            detail += (
                f", MLMG solve {gains['mlmg_solve']:.3f}%, FApply wall "
                f"{gains['fapply_inclusive_wall']:.3f}%, per-call "
                f"{gains['fapply_seconds_per_call']:.3f}%"
            )
        lines.extend(["", f"- {regime} candidate gains: {detail}.", ""])
    out_md.write_text("\n".join(lines) + "\n")
    print(out_md.read_text(), end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
