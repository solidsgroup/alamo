#!/usr/bin/env python3
import re
import subprocess
import sys
import csv
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
PATTERN_DIR = ROOT / "docs/gpu_manual/patterns"
EXPECTED = {
    "GPU-001": "correctness",
    "GPU-002": "correctness",
    "GPU-003": "correctness",
    "GPU-004": "correctness",
    "GPU-005": "correctness",
    "GPU-006": "correctness",
    "GPU-007": "correctness",
    "GPU-008": "correctness",
    "GPU-009": "correctness",
    "GPU-010": "scaffolding",
    "GPU-011": "correctness",
    "GPU-012": "correctness",
    "GPU-013": "correctness",
    "GPU-015": "performance",
    "GPU-016": "correctness",
    "GPU-017": "correctness",
    "GPU-018": "performance",
    "GPU-019": "performance",
    "GPU-020": "performance",
    "GPU-021": "correctness",
    "GPU-022": "correctness",
    "GPU-023": "correctness",
    "GPU-024": "correctness",
    "GPU-025": "scaffolding",
    "GPU-030": "correctness",
}
FIELDS = [
    "Status",
    "Class",
    "Recognizer",
    "Applies",
    "Transform",
    "Constraints",
    "Verify",
    "Failure modes",
    "Evidence",
]


def git_ok(*args: str) -> bool:
    return subprocess.run(
        ["git", *args], cwd=ROOT, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
    ).returncode == 0


def main() -> int:
    errors = []
    files = sorted(PATTERN_DIR.glob("*.md"))
    found = {}
    recognizers = {}
    for path in files:
        match = re.match(r"((?:GPU|NUM)-\d{3})-.*\.md$", path.name)
        if not match:
            errors.append(f"{path.name}: invalid filename")
            continue
        pattern_id = match.group(1)
        found[pattern_id] = path
        text = path.read_text()
        lines = text.splitlines()
        if not lines or not lines[0].startswith(f"# {pattern_id}: "):
            errors.append(f"{path.name}: H1 does not match filename ID")

        top_fields = []
        values = {}
        for line in lines[1:]:
            field = re.match(r"^([A-Za-z][A-Za-z ]*):(?:\s*(.*))?$", line)
            if field:
                top_fields.append(field.group(1))
                values[field.group(1)] = field.group(2) or ""
        if top_fields != FIELDS:
            errors.append(
                f"{path.name}: top-level fields {top_fields!r}, expected frozen order"
            )
        if values.get("Status") != "draft":
            errors.append(f"{path.name}: Phase 2 Status must be draft")
        if values.get("Class") != EXPECTED.get(pattern_id):
            errors.append(
                f"{path.name}: Class {values.get('Class')!r}, expected {EXPECTED.get(pattern_id)!r}"
            )
        recognizer = values.get("Recognizer", "").lower()
        recognizers[pattern_id] = values.get("Recognizer", "")
        if not recognizer.startswith(("regex:", "build-error:", "manual:")):
            errors.append(f"{path.name}: invalid Recognizer type")
        elif recognizer.startswith("regex:"):
            expression = values["Recognizer"][len("regex:"):].strip().strip("`")
            try:
                re.compile(expression)
            except re.error as error:
                errors.append(f"{path.name}: invalid Python regex: {error}")
        if not lines or not lines[-1].startswith("Evidence:"):
            errors.append(f"{path.name}: Evidence must be the final line")

        words = len(re.findall(r"\b[\w'-]+\b", text))
        tokens = words * 1.3
        if not 200 <= tokens <= 400:
            errors.append(
                f"{path.name}: approximate token count {tokens:.1f} ({words} words) outside 200-400"
            )

        evidence = values.get("Evidence", "")
        hashes = re.findall(r"\b[0-9a-f]{40}\b", evidence)
        if not hashes:
            errors.append(f"{path.name}: Evidence has no full commit hash")
        for commit in hashes:
            if not git_ok("cat-file", "-e", f"{commit}^{{commit}}"):
                errors.append(f"{path.name}: missing evidence commit {commit}")
            elif not git_ok("merge-base", "--is-ancestor", commit, "chamber-gpu"):
                errors.append(f"{path.name}: evidence {commit} is not on chamber-gpu")

        constraint = values.get("Constraints", "").lower().replace("-", " ")
        if pattern_id.startswith("GPU-") and values.get("Class") == "correctness":
            if "mandatory" not in constraint or "recognizer" not in constraint:
                errors.append(
                    f"{path.name}: GPU correctness constraint must make matched recognizers mandatory"
                )
        if values.get("Class") == "performance":
            if "profiling justification" not in constraint:
                errors.append(f"{path.name}: performance requires profiling justification")
            if not any(word in constraint for word in ("golden", "correctness")):
                errors.append(f"{path.name}: performance must preserve golden correctness")
        if pattern_id in {"GPU-010", "GPU-025"}:
            required = ("temporary", "incremental port", "shrink", "never grow")
            for phrase in required:
                if phrase not in constraint:
                    errors.append(f"{path.name}: scaffolding constraint lacks {phrase!r}")
        if pattern_id == "GPU-010" and "user approval" not in constraint:
            errors.append(f"{path.name}: hard-conversion quarantine requires user approval")
        if pattern_id == "GPU-001":
            if not recognizer.startswith("build-error:"):
                errors.append(f"{path.name}: GPU-001 must use build-error recognizer")
            if "calling a __host__ function from a __device__ function" not in values.get("Recognizer", ""):
                errors.append(f"{path.name}: GPU-001 lacks stable nvcc diagnostic fragment")

    missing = sorted(set(EXPECTED) - set(found))
    extra = sorted(set(found) - set(EXPECTED))
    if missing:
        errors.append(f"missing patterns: {', '.join(missing)}")
    if extra:
        errors.append(f"unexpected patterns: {', '.join(extra)}")
    if any(pattern_id in found for pattern_id in (
        "GPU-014", "GPU-026", "GPU-027", "GPU-028", "GPU-029",
        "NUM-001", "NUM-002", "NUM-003", "NUM-004",
    )):
        errors.append("reserved/renamed IDs must not have pattern files")

    table_path = ROOT / "docs/gpu_manual/recognizers/table.csv"
    with table_path.open(newline="") as handle:
        table_reader = csv.DictReader(handle)
        table_rows = list(table_reader)
        if table_reader.fieldnames != ["pattern_id", "type", "expression", "notes"]:
            errors.append("recognizers/table.csv: invalid header")
    table_ids = [row["pattern_id"] for row in table_rows]
    if table_ids != sorted(EXPECTED):
        errors.append("recognizers/table.csv: IDs must exactly match active patterns in order")
    for row in table_rows:
        pattern_id = row["pattern_id"]
        if pattern_id not in recognizers:
            continue
        kind, payload = recognizers[pattern_id].split(":", 1)
        kind = kind.lower()
        payload = payload.strip().strip("`")
        if row["type"] != kind:
            errors.append(f"recognizers/table.csv: {pattern_id} type differs from Tier 1")
        if kind == "manual":
            if row["expression"]:
                errors.append(f"recognizers/table.csv: {pattern_id} manual row has expression")
        elif row["expression"] != payload:
            errors.append(f"recognizers/table.csv: {pattern_id} expression differs from Tier 1")
        if not row["notes"]:
            errors.append(f"recognizers/table.csv: {pattern_id} lacks notes")

    regex_ids = {row["pattern_id"] for row in table_rows if row["type"] == "regex"}

    def validate_coverage(path, require_all_regex):
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            rows = list(reader)
            if reader.fieldnames != ["file", "pattern_id", "hits"]:
                errors.append(f"{path.name}: invalid header")
                return
        seen_regex = set()
        order = []
        for number, row in enumerate(rows, 2):
            pattern_id = row["pattern_id"]
            if pattern_id not in regex_ids:
                errors.append(f"{path.name}:{number}: non-regex or inactive ID {pattern_id}")
            else:
                seen_regex.add(pattern_id)
            try:
                if int(row["hits"]) <= 0:
                    raise ValueError
            except ValueError:
                errors.append(f"{path.name}:{number}: hits must be a positive integer")
            if Path(row["file"]).is_absolute() or Path(row["file"]).suffix not in {".H", ".cpp", ".cu", ".cc"}:
                errors.append(f"{path.name}:{number}: invalid source path {row['file']}")
            order.append((row["file"], pattern_id))
        if order != sorted(order):
            errors.append(f"{path.name}: rows are not deterministically sorted")
        if require_all_regex and seen_regex != regex_ids:
            errors.append(
                f"{path.name}: missing BASE true positives for {', '.join(sorted(regex_ids - seen_regex))}"
            )

    validate_coverage(ROOT / "docs/gpu_manual/COVERAGE.csv", False)
    base_coverage = ROOT / "docs/gpu_manual/build/phase4/BASE_COVERAGE.csv"
    if base_coverage.exists():
        validate_coverage(base_coverage, True)

    one_off_path = ROOT / "docs/gpu_manual/ONE_OFFS.md"
    one_off_text = one_off_path.read_text()
    if "[NUM] entry must surface it to the user and never apply it silently" not in one_off_text:
        errors.append("ONE_OFFS.md: missing [NUM] user-surfacing policy")
    feature_path = ROOT / "docs/gpu_manual/FEATURES.md"
    feature_text = feature_path.read_text()
    if "explicit do-not-import list" not in feature_text or "never apply a FEATURE without task-level user opt-in" not in feature_text:
        errors.append("FEATURES.md: missing do-not-import/task-level opt-in policy")
    with (ROOT / "docs/gpu_manual/build/HUNK_MAP.csv").open(newline="") as handle:
        hunk_rows = list(csv.DictReader(handle))
    expected_one_offs = {
        row["hunk_id"]: row["num_tag"] for row in hunk_rows
        if row["classification"] == "ONEOFF"
    }
    expected_features = {
        row["hunk_id"]: row["num_tag"] for row in hunk_rows
        if row["classification"] == "FEATURE"
    }
    entry_pattern = re.compile(
        r"^- (?P<id>\S+)(?P<num> \[NUM\])? \| .+ \| commit (?P<commit>[0-9a-f]{40}) \| .+$"
    )

    def validate_ledger(name, text, expected):
        found = {}
        for number, line in enumerate(text.splitlines(), 1):
            if not line.startswith("- "):
                continue
            match = entry_pattern.match(line)
            if not match:
                errors.append(f"{name}:{number}: malformed entry or non-full commit hash")
                continue
            hunk_id = match.group("id")
            if hunk_id in found:
                errors.append(f"{name}:{number}: duplicate {hunk_id}")
            found[hunk_id] = "yes" if match.group("num") else "no"
            commit = match.group("commit")
            if not git_ok("cat-file", "-e", f"{commit}^{{commit}}"):
                errors.append(f"{name}:{number}: missing commit {commit}")
            elif not git_ok("merge-base", "--is-ancestor", commit, "chamber-gpu"):
                errors.append(f"{name}:{number}: commit {commit} is not on chamber-gpu")
        missing = sorted(set(expected) - set(found))
        extra = sorted(set(found) - set(expected))
        if missing:
            errors.append(f"{name}: missing ledger IDs: {', '.join(missing)}")
        if extra:
            errors.append(f"{name}: extra ledger IDs: {', '.join(extra)}")
        for hunk_id in sorted(set(expected) & set(found)):
            if expected[hunk_id] != found[hunk_id]:
                errors.append(f"{name}: [NUM] tag mismatch for {hunk_id}")

    validate_ledger("ONE_OFFS.md", one_off_text, expected_one_offs)
    validate_ledger("FEATURES.md", feature_text, expected_features)

    md_map_path = ROOT / "docs/gpu_manual/MD_MAP.csv"
    if md_map_path.exists():
        with md_map_path.open(newline="") as handle:
            md_rows = list(csv.DictReader(handle))
        if len(md_rows) != 155:
            errors.append(f"MD_MAP.csv: expected 155 rows, found {len(md_rows)}")
        allowed_dispositions = {
            "mapped", "anti-pattern", "evidence-only", "stale-dropped", "unmined"
        }
        for number, row in enumerate(md_rows, 2):
            if row["disposition"] not in allowed_dispositions:
                errors.append(f"MD_MAP.csv:{number}: bad disposition {row['disposition']}")
            for pattern_id in filter(None, row["pattern_ids"].split(";")):
                if pattern_id not in EXPECTED:
                    errors.append(f"MD_MAP.csv:{number}: inactive pattern {pattern_id}")
        inventory_lines = (
            ROOT / "docs/gpu_manual/build/MD_INVENTORY.txt"
        ).read_text().splitlines()
        inventory_paths = [line.split("\t", 2)[2] for line in inventory_lines[5:]]
        if Counter(row["path"] for row in md_rows) != Counter(inventory_paths):
            errors.append("MD_MAP.csv: path multiset does not match MD_INVENTORY.txt")

    if errors:
        print("pattern validation failed:")
        for error in errors:
            print(f"- {error}")
        return 1
    print(f"validated {len(files)} Tier 1 patterns")
    return 0


if __name__ == "__main__":
    sys.exit(main())
