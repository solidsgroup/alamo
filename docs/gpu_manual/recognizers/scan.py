#!/usr/bin/env python3
"""Site-level advisory recognizer scanner (stdlib only).

The scanner records candidates; compiler diagnostics and inspection ledgers are
authoritative.  Reviewed dispositions can be carried forward with --previous.
"""
import argparse
import csv
import hashlib
import re
import sys
from pathlib import Path

SCHEMA = "3"
FIELDS = ["schema_version", "port_id", "source_revision", "site_id", "file", "line", "pattern_id", "state", "evidence"]
TABLE_FIELDS = ["schema_version", "pattern_id", "type", "candidate_expression",
                "converted_expression", "exclude_expression", "confirmation", "notes"]
STATES = {"candidate", "converted", "not-applicable", "false-positive"}
EXTENSIONS = {".H", ".cpp", ".cu", ".cc"}
SKIP_COMPONENTS = {".git", "ext", "bin", "build", "docs"}


def read_table(path):
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames != TABLE_FIELDS:
            raise ValueError(f"recognizer table schema must be {','.join(TABLE_FIELDS)}")
        rows = list(reader)
    for row in rows:
        if row["schema_version"] != SCHEMA or not row["pattern_id"]:
            raise ValueError("recognizer table has missing/unsupported schema_version or pattern_id")
        if row["type"] not in {"regex", "manual", "build-error"}:
            raise ValueError(f"{row['pattern_id']}: invalid type {row['type']}")
        if row["type"] != "regex" and (row["candidate_expression"] or row["converted_expression"] or row["exclude_expression"]):
            raise ValueError(f"{row['pattern_id']}: non-regex expressions must be empty")
        for field in ("candidate_expression", "converted_expression", "exclude_expression"):
            if row[field]:
                try:
                    re.compile(row[field], re.MULTILINE)
                except re.error as exc:
                    raise ValueError(f"{row['pattern_id']}: invalid {field}: {exc}") from exc
    return rows


def site_id(file_name, pattern, matched, occurrence):
    """Stable identity: path/content, deliberately independent of line number."""
    digest = hashlib.sha1((file_name + "\0" + pattern + "\0" + matched).encode()).hexdigest()[:16]
    return f"{pattern}:{digest}:{occurrence}"


def previous_rows(path, port_id, source_revision):
    if not path:
        return {}
    with path.open(newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames != FIELDS:
            raise ValueError(f"previous scan schema must be {','.join(FIELDS)}")
        rows = list(reader)
    kept = {}
    for row in rows:
        if (row["schema_version"] != SCHEMA or row.get("port_id") != port_id
                or row.get("source_revision") != source_revision
                or row["state"] not in {"not-applicable", "false-positive"}):
            continue
        kept[row["site_id"]] = row
    return kept


def scan(root, table, prior, port_id, source_revision):
    output = []
    paths = []
    for p in root.rglob("*"):
        if not p.is_file() or p.suffix not in EXTENSIONS:
            continue
        relative_parts = p.relative_to(root).parts
        # Skip repository metadata/generated trees by default. If --root points
        # directly at one of those trees, its own files remain scannable.
        if any(part in SKIP_COMPONENTS for part in relative_parts[:-1]):
            continue
        paths.append(p)
    for path in sorted(paths):
        source = path.read_text(errors="ignore")
        rel = path.relative_to(root).as_posix()
        for row in table:
            if row["type"] != "regex":
                continue
            candidate = re.compile(row["candidate_expression"], re.MULTILINE)
            converted = re.compile(row["converted_expression"], re.MULTILINE) if row["converted_expression"] else None
            excluded = re.compile(row["exclude_expression"], re.MULTILINE) if row["exclude_expression"] else None
            seen = {}
            candidate_spans = set()
            for match in candidate.finditer(source):
                candidate_spans.add(match.span())
                text = match.group(0)
                occurrence = seen.get(text, 0)
                seen[text] = occurrence + 1
                sid = site_id(rel, row["pattern_id"], text, occurrence)
                line = source.count("\n", 0, match.start()) + 1
                state = "candidate"
                evidence = "candidate regex"
                line_text = source.splitlines()[line - 1] if source.splitlines() else text
                if excluded and excluded.search(line_text):
                    state, evidence = "not-applicable", "table exclusion"
                elif sid in prior:
                    state, evidence = prior[sid]["state"], prior[sid]["evidence"]
                output.append({"schema_version": SCHEMA, "port_id": port_id, "source_revision": source_revision,
                               "site_id": sid, "file": rel,
                               "line": str(line), "pattern_id": row["pattern_id"],
                               "state": state, "evidence": evidence})
            if converted:
                seen = {}
                for match in converted.finditer(source):
                    # An overlapping rule is unresolved, not converted. Keep the
                    # candidate visible and avoid two states for one source site.
                    if match.span() in candidate_spans:
                        continue
                    text = match.group(0)
                    occurrence = seen.get(text, 0)
                    seen[text] = occurrence + 1
                    sid = site_id(rel, row["pattern_id"], text, occurrence)
                    line = source.count("\n", 0, match.start()) + 1
                    output.append({"schema_version": SCHEMA, "port_id": port_id, "source_revision": source_revision,
                                   "site_id": sid, "file": rel,
                                   "line": str(line), "pattern_id": row["pattern_id"],
                                   "state": "converted", "evidence": "converted regex"})
    return sorted(output, key=lambda r: (r["file"], r["pattern_id"], int(r["line"]), r["state"], r["site_id"]))


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--table", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--previous", type=Path)
    parser.add_argument("--port-id", required=True)
    parser.add_argument("--source-revision", required=True)
    args = parser.parse_args(argv)
    try:
        rows = scan(args.root, read_table(args.table),
                    previous_rows(args.previous, args.port_id, args.source_revision),
                    args.port_id, args.source_revision)
    except (OSError, ValueError) as exc:
        parser.error(str(exc))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
