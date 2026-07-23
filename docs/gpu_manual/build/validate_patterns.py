#!/usr/bin/env python3
"""Structural oracle for the generalized GPU manual."""

import csv
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[3]
MANUAL = ROOT / "docs/gpu_manual"
PATTERNS = MANUAL / "patterns"
EXPECTED_CLASSES = {
    **{f"GPU-{number:03d}": "correctness" for number in range(1, 14)},
    "GPU-015": "optimization",
    "GPU-016": "correctness",
    "GPU-017": "correctness",
    "GPU-018": "optimization",
    "GPU-019": "optimization",
    "GPU-020": "optimization",
    "GPU-021": "correctness",
    "GPU-022": "correctness",
    "GPU-023": "correctness",
    "GPU-024": "correctness",
    "GPU-025": "scaffolding",
    "GPU-030": "correctness",
    "GPU-031": "correctness",
}
EXPECTED_CLASSES["GPU-010"] = "scaffolding"
FIELDS = [
    "Transform status",
    "Class",
    "Detection",
    "Invariant",
    "Port contract",
    "Transform",
    "Corpus example",
    "Constraints",
    "Verify",
    "Failure modes",
    "Evidence",
]
TABLE_FIELDS = [
    "schema_version", "pattern_id", "type", "candidate_expression",
    "converted_expression", "exclude_expression", "confirmation", "notes",
]
COVERAGE_FIELDS = [
    "schema_version", "port_id", "source_revision", "site_id", "file",
    "line", "pattern_id", "state", "evidence",
]
STATUSES = {"draft", "file-verified", "transfer-verified", "cross-family"}
SITE_STATES = {"candidate", "converted", "not-applicable", "false-positive"}
IDENTIFIER_FITTED_CANDIDATES = {
    "trac_hi", "disp_hi", "massflux", "mdot", "volume", "DDW", "ximg",
    "m_bc", "GetBC", "static_polymorphism_parser", "Set::Garbage",
}


def top_fields(text):
    names = []
    values = {}
    for line in text.splitlines()[1:]:
        match = re.match(r"^([A-Za-z][A-Za-z ]*):(?:\s*(.*))?$", line)
        if match:
            names.append(match.group(1))
            values[match.group(1)] = (match.group(2) or "").strip()
    return names, values


def require_text(errors, path, snippets):
    text = path.read_text()
    for snippet in snippets:
        if snippet not in text:
            errors.append(f"{path.relative_to(ROOT)}: missing {snippet!r}")


def validate_csv_header(errors, relative, expected):
    path = ROOT / relative
    with path.open(newline="") as handle:
        reader = csv.reader(handle)
        header = next(reader, [])
    if header != expected:
        errors.append(f"{relative}: header {header!r}, expected {expected!r}")


def main():
    errors = []
    found = {}
    forbidden_verify = re.compile(
        r"(?:Flame|PFC|CahnHilliard|HeatConduction|Elastic|benchmark/|make\s)",
        re.IGNORECASE,
    )

    for path in sorted(PATTERNS.glob("GPU-*.md")):
        match = re.match(r"(GPU-\d{3})-.*\.md$", path.name)
        if not match:
            errors.append(f"{path.name}: invalid pattern filename")
            continue
        pattern_id = match.group(1)
        text = path.read_text()
        names, values = top_fields(text)
        found[pattern_id] = values

        if not text.startswith(f"# {pattern_id}: "):
            errors.append(f"{path.name}: heading does not match ID")
        if names != FIELDS:
            errors.append(f"{path.name}: top-level fields {names!r}, expected {FIELDS!r}")
        expected_status = "file-verified" if pattern_id in {"GPU-007", "GPU-016"} else "draft"
        if values.get("Transform status") != expected_status:
            errors.append(f"{path.name}: status must be {expected_status}")
        if values.get("Transform status") not in STATUSES:
            errors.append(f"{path.name}: invalid status")
        if values.get("Class") != EXPECTED_CLASSES.get(pattern_id):
            errors.append(f"{path.name}: wrong class {values.get('Class')!r}")
        for field in ("Detection", "Invariant", "Port contract", "Corpus example", "Constraints", "Verify", "Failure modes", "Evidence"):
            if not values.get(field):
                errors.append(f"{path.name}: empty {field}")
        if "advisory" not in values.get("Detection", "").lower():
            errors.append(f"{path.name}: Detection must say advisory")
        if not re.search(r"corpus|chamber-gpu", values.get("Corpus example", ""), re.I):
            errors.append(f"{path.name}: Corpus example is not explicitly labeled")
        verify = values.get("Verify", "")
        if forbidden_verify.search(verify):
            errors.append(f"{path.name}: port-specific command/physics in Verify")
        if values.get("Class") == "correctness":
            if "VALIDATION.md" not in verify or "tolerance rationale" not in verify:
                errors.append(f"{path.name}: correctness Verify must use validation contract and tolerance rationale")
            if not any(category in verify for category in (
                "strict-build", "analytic-exact", "conservation",
                "golden-regression", "restart-parity", "multi-box", "sanitizer",
            )):
                errors.append(f"{path.name}: correctness Verify lacks a contract category")
        if values.get("Class") == "optimization":
            if "PERFORMANCE.md" not in verify or "baseline" not in verify.lower():
                errors.append(f"{path.name}: optimization must follow the baseline contract")
        if pattern_id in {"GPU-010", "GPU-012", "GPU-025"}:
            if "ARCHITECTURE_POLICIES.md" not in text:
                errors.append(f"{path.name}: missing single policy-home reference")
        evidence = values.get("Evidence", "")
        if "Primary:" not in evidence or "primary-sources.md" not in evidence:
            errors.append(f"{path.name}: missing primary evidence")
        if "Corpus:" not in evidence or not re.search(r"\b[0-9a-f]{40}\b", evidence):
            errors.append(f"{path.name}: missing full-hash corpus evidence")
        if not text.rstrip().splitlines()[-1].startswith("Evidence:"):
            errors.append(f"{path.name}: Evidence must be final")
        tokens = len(re.findall(r"\b[\w'-]+\b", text)) * 1.3
        if not 200 <= tokens <= 400:
            errors.append(f"{path.name}: approximate tokens {tokens:.1f} outside 200-400")

    if set(found) != set(EXPECTED_CLASSES):
        errors.append("pattern inventory does not match the active 26 IDs")

    table_path = MANUAL / "recognizers/table.csv"
    with table_path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        table_rows = list(reader)
        if reader.fieldnames != TABLE_FIELDS:
            errors.append("recognizers/table.csv: invalid v3 header")
    if [row.get("pattern_id") for row in table_rows] != sorted(EXPECTED_CLASSES):
        errors.append("recognizers/table.csv: IDs must be sorted and complete")
    table_types = {}
    for row in table_rows:
        pattern_id = row.get("pattern_id", "")
        table_types[pattern_id] = row.get("type")
        if row.get("schema_version") != "3":
            errors.append(f"{pattern_id}: recognizer schema_version must be 3")
        if row.get("type") not in {"regex", "manual", "build-error"}:
            errors.append(f"{pattern_id}: invalid recognizer type")
        if not row.get("confirmation") or not row.get("notes"):
            errors.append(f"{pattern_id}: recognizer lacks confirmation/notes")
        expression_fields = ("candidate_expression", "converted_expression", "exclude_expression")
        if row.get("type") == "regex" and not row.get("candidate_expression"):
            errors.append(f"{pattern_id}: regex lacks candidate expression")
        candidate = row.get("candidate_expression", "")
        for literal in IDENTIFIER_FITTED_CANDIDATES:
            if literal in candidate:
                errors.append(f"{pattern_id}: identifier-fitted candidate contains {literal!r}")
        if row.get("type") != "regex" and any(row.get(field) for field in expression_fields):
            errors.append(f"{pattern_id}: manual/build-error row has lexical expression")
        for field in expression_fields:
            if row.get(field):
                try:
                    re.compile(row[field], re.MULTILINE)
                except re.error as error:
                    errors.append(f"{pattern_id}: invalid {field}: {error}")

    if (MANUAL / "COVERAGE.csv").exists():
        errors.append("COVERAGE.csv: static root coverage must not exist")

    required_artifacts = [
        "BRIEF.md", "INDEX.md", "VALIDATION.md", "ONBOARDING.md",
        "ARCHITECTURE_POLICIES.md", "GPU_NATIVE_SHAPE.md", "PERFORMANCE.md", "RECOGNIZERS.md",
        "STATUS.md", "BLIND_SPOTS.md", "evidence/primary-sources.md",
        "evidence/GPU-002-value-dispatch-example.md",
        "evidence/host-only-numerical-kernels.md",
        "recognizers/scan.py", "templates/SCOPE.md", "templates/CLOSURE.csv",
        "templates/INSPECTION_LEDGER.csv", "templates/VALIDATION.csv",
        "templates/EFFICIENCY.md", "templates/FEATURE_DECISIONS.csv",
        "templates/HARVEST.md", "templates/PORT_STATUS.csv",
        "templates/COVERAGE.csv", "templates/FIELD_LAYOUT.csv",
        "templates/KERNEL_GRAPH.csv", "templates/SHAPE_PROFILE.md",
    ]
    for relative in required_artifacts:
        if not (MANUAL / relative).is_file():
            errors.append(f"missing artifact: {relative}")

    validate_csv_header(errors, "docs/gpu_manual/templates/CLOSURE.csv", [
        "schema_version", "port_id", "node", "kind", "required_by",
        "discovery", "evidence", "status", "owner", "retirement_trigger",
    ])
    validate_csv_header(errors, "docs/gpu_manual/templates/INSPECTION_LEDGER.csv", [
        "schema_version", "port_id", "file", "anchor_line", "taxonomy",
        "question", "source", "pattern_id", "disposition", "evidence", "owner",
    ])
    require_text(errors, MANUAL / "templates/INSPECTION_LEDGER.csv", [
        "numerical-kernel", "field-layout", "kernel-graph",
    ])
    validate_csv_header(errors, "docs/gpu_manual/templates/VALIDATION.csv", [
        "schema_version", "port_id", "id", "category", "applicability",
        "oracle", "command", "reference", "result", "tolerance_rationale",
        "evidence", "owner",
    ])
    validate_csv_header(errors, "docs/gpu_manual/templates/PORT_STATUS.csv", [
        "schema_version", "port_id", "physics_family", "scope", "closure",
        "inspection", "closed_book_onboarding", "closed_book_transform",
        "validation", "gpu_safety", "gpu_native_shape", "baseline_efficiency", "harvest",
        "optimization", "evidence",
    ])
    validate_csv_header(errors, "docs/gpu_manual/templates/COVERAGE.csv", COVERAGE_FIELDS)
    validate_csv_header(errors, "docs/gpu_manual/templates/FIELD_LAYOUT.csv", [
        "schema_version", "port_id", "field", "storage", "ncomp",
        "component_order", "ghost_region", "kernel_consumers", "access_pattern",
        "residency", "transfer_bytes_per_step", "coalescing_locality_rationale",
        "disposition", "evidence", "owner",
    ])
    validate_csv_header(errors, "docs/gpu_manual/templates/KERNEL_GRAPH.csv", [
        "schema_version", "port_id", "kernel", "phase", "inputs", "outputs",
        "iteration_space", "mfiter_granularity", "launches_per_step",
        "components_per_launch", "dependencies", "synchronization",
        "host_device_bytes", "branch_iteration_notes", "dispatch_state",
        "reduction_atomic", "allocation_lifetime", "disposition", "evidence", "owner",
    ])

    index = (MANUAL / "INDEX.md").read_text()
    if len(re.findall(r"\b[\w'-]+\b", index)) * 1.3 > 1500:
        errors.append("INDEX.md: exceeds approximate 1500-token budget")
    for marker in ("`[I]`", "`[P]`", "`[C]`"):
        if marker not in index:
            errors.append(f"INDEX.md: missing legend marker {marker}")
    pattern_section = index.split("## Patterns", 1)[-1].split("## Contracts", 1)[0]
    for pattern_id in sorted(EXPECTED_CLASSES):
        if pattern_section.count(f"- {pattern_id} ") != 1:
            errors.append(f"INDEX.md: missing/duplicate one-liner for {pattern_id}")

    require_text(errors, MANUAL / "VALIDATION.md", [
        "analytic-exact", "conservation", "golden-regression", "restart-parity",
        "multi-box", "sanitizer", "Tolerance rationale", "pending-pilot",
    ])
    require_text(errors, MANUAL / "ONBOARDING.md", [
        "compiler-first", "host-loop", "launch", "receiver-type", "diagnostic",
        "capture", "lifetime", "reduction", "dispatch", "numerical-kernel",
        "field-layout", "kernel-graph", "zero `open`", "--port-id",
        "--source-revision",
    ])
    require_text(errors, MANUAL / "ARCHITECTURE_POLICIES.md", [
        "## Device-side error propagation", "## Scaffolding lifecycle",
        "## FEATURE and [NUM] surfacing", "Owner decision:", "Worked example:",
    ])
    require_text(errors, MANUAL / "PERFORMANCE.md", [
        "No hot host loop", "no per-cell or per-tile host", "per-tile device-wide",
        "runtime virtual/plugin", "once per component", "dominate the steady-state",
        "GPU_NATIVE_SHAPE=pass", "Field layout", "registers/spills",
    ])
    require_text(errors, MANUAL / "GPU_NATIVE_SHAPE.md", [
        "Field layout and component ordering", "Kernel graph and iteration granularity",
        "Residency and transfer accounting", "Numerical call-chain complexity",
        "registers per thread", "Shared memory, atomics, and allocation",
        "GPU_NATIVE_SHAPE=pass",
    ])
    require_text(errors, MANUAL / "STATUS.md", [
        "file-verified", "transfer-verified", "cross-family", "completed-port",
        "Frozen unseen-target gate",
        "Harvest obligation", "A port with no harvested\nmanual change is incomplete",
    ])
    require_text(errors, MANUAL / "RECOGNIZERS.md", [
        "Compiler diagnostics and the\nrecorded inspection ledger outrank regex",
        "zero `candidate` rows", "--previous", "--port-id", "--source-revision",
    ])
    phase6_counts = (MANUAL / "build/phase6/COUNTS.txt").read_text()
    if re.search(r"^GATE=pass$", phase6_counts, re.MULTILINE):
        errors.append("phase6/COUNTS.txt: unqualified global GATE=pass")

    if errors:
        print("manual validation failed:", file=sys.stderr)
        for error in errors:
            print(f"- {error}", file=sys.stderr)
        return 1
    print(f"VALIDATED_PATTERNS={len(found)}")
    print("FILE_VERIFIED=2/26; TRANSFER_VERIFIED=0/26; CROSS_FAMILY=0/26")
    print("RECOGNIZER_ORACLE=advisory-lifecycle; formal precision/recall out-of-scope")
    print("PILOT=required")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
