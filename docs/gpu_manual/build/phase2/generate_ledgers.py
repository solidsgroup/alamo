#!/usr/bin/env python3
import csv
import re
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
MANUAL = ROOT / "docs/gpu_manual"
BUILD = MANUAL / "build"
PHASE2 = BUILD / "phase2"

RETIRED = {
    "NUM-001": (
        "Adds the conservative face-flux elastic formulation.",
        "5ae7dffa35ac5ec0bfeb9e2b9ced9200885f8227",
        "device-independent solver formulation retired from NUM-001",
    ),
    "NUM-002": (
        "Adds explicit coarse-fine ghost-row and interpolation policies.",
        "5ae7dffa35ac5ec0bfeb9e2b9ced9200885f8227",
        "device-independent AMR policy retired from NUM-002",
    ),
    "NUM-003": (
        "Carries component layouts through AMR transfers.",
        "5ae7dffa35ac5ec0bfeb9e2b9ced9200885f8227",
        "device-independent AMR component behavior retired from NUM-003",
    ),
    "NUM-004": (
        "Adds bounded line-search damping to Newton updates.",
        "dd4054056aa5ef24bde948ac1b04013b5b106eb6",
        "device-independent nonlinear solver behavior retired from NUM-004",
    ),
}


def read_csv(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def latest_file_commit(path):
    result = subprocess.run(
        [
            "git", "log", "-1", "--format=%H",
            "a6f40a3b4daca4de8dab4979750b359ddce99b35..chamber-gpu",
            "--", path,
        ],
        cwd=ROOT,
        check=True,
        text=True,
        capture_output=True,
    )
    commit = result.stdout.strip()
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise SystemExit(f"no evidence commit for {path}")
    return commit


old_text = (MANUAL / "ONE_OFFS.md").read_text()
commit_by_id = {}
for line in old_text.splitlines():
    match = re.match(
        r"^- (?P<id>\S+)(?: \[NUM\])? \| .* \| commit (?P<commit>[0-9a-f]{40}) \|",
        line,
    )
    if match:
        commit_by_id[match.group("id")] = match.group("commit")

review_by_id = {}
for part in ("A", "B"):
    for row in read_csv(PHASE2 / "feature_review" / f"{part}.csv"):
        review_by_id[row["hunk_id"]] = row

retired_by_id = {
    row["hunk_id"]: row
    for row in read_csv(PHASE2 / "RETIRED_PATTERNS.csv")
}
map_rows = read_csv(BUILD / "HUNK_MAP.csv")


def entry(row):
    hunk_id = row["hunk_id"]
    tag = " [NUM]" if row["num_tag"] == "yes" else ""
    if hunk_id in review_by_id:
        review = review_by_id[hunk_id]
        summary = review["summary"].rstrip(".") + "."
        commit = commit_by_id.get(hunk_id) or latest_file_commit(row["file"])
        if row["classification"] == "ONEOFF":
            reason = "file-specific BASE GPU-port residue with no repeatable transform"
        else:
            reason = "independent chamber-gpu capability; explicit task-level opt-in required"
    elif hunk_id in retired_by_id:
        pattern_id = retired_by_id[hunk_id]["retired_pattern"]
        summary, commit, reason = RETIRED[pattern_id]
    else:
        raise SystemExit(f"no ledger metadata for {hunk_id}")
    summary = summary.replace("|", "/")
    return f"- {hunk_id}{tag} | {summary} | commit {commit} | {reason}"


oneoffs = [row for row in map_rows if row["classification"] == "ONEOFF"]
features = [row for row in map_rows if row["classification"] == "FEATURE"]

oneoff_lines = [
    "# One-off changes",
    "",
    "Policy: A port worker touching a file with a [NUM] entry must surface it to the user and never apply it silently. Every entry is file-specific GPU-port residue, not a reusable transform.",
    "",
    *(entry(row) for row in oneoffs),
    "",
]
feature_lines = [
    "# Chamber-gpu features",
    "",
    "Policy: This is both an upstream-merge roadmap and an explicit do-not-import list. A port worker using chamber-gpu as reference must never apply a FEATURE without task-level user opt-in. [NUM] entries additionally require explicit numerical-behavior opt-in.",
    "",
    *(entry(row) for row in features),
    "",
]

(MANUAL / "ONE_OFFS.md").write_text("\n".join(oneoff_lines))
(MANUAL / "FEATURES.md").write_text("\n".join(feature_lines))

print(f"wrote {len(oneoffs)} ONEOFF and {len(features)} FEATURE entries")
