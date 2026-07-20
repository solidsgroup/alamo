#!/usr/bin/env python3
import csv
from pathlib import Path


ROOT = Path(__file__).resolve().parents[4]
MANUAL = ROOT / "docs/gpu_manual"
PHASE3 = MANUAL / "build/phase3"
ACTIVE = {
    "GPU-001", "GPU-002", "GPU-003", "GPU-004", "GPU-005", "GPU-006",
    "GPU-007", "GPU-008", "GPU-009", "GPU-010", "GPU-011", "GPU-012",
    "GPU-013", "GPU-015", "GPU-016", "GPU-017", "GPU-018", "GPU-019",
    "GPU-020", "GPU-021", "GPU-022", "GPU-023", "GPU-024", "GPU-025",
    "GPU-030",
}
DISPOSITIONS = {
    "evidence-only": "evidence-only",
    "EVIDENCE_ONLY": "evidence-only",
    "ROUTE_PATTERN": "mapped",
    "ANTI_PATTERN_CANDIDATE": "mapped",
    "STALE_DROP": "stale-dropped",
    "UNMINED": "unmined",
}


rows = []
for part in ("A", "B"):
    with (PHASE3 / f"{part}_MD_MAP.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            disposition = DISPOSITIONS[row["disposition"]]
            raw_ids = row["pattern_ids"].replace(",", ";").split(";")
            pattern_ids = []
            for pattern_id in raw_ids:
                pattern_id = pattern_id.strip()
                if not pattern_id or pattern_id.startswith("NUM-"):
                    continue
                if pattern_id not in ACTIVE:
                    raise SystemExit(f"unknown active pattern {pattern_id}: {row['path']}")
                if pattern_id not in pattern_ids:
                    pattern_ids.append(pattern_id)
            rows.append((row["path"], disposition, ";".join(pattern_ids)))

if len(rows) != 155:
    raise SystemExit(f"expected 155 inventory rows, found {len(rows)}")

with (MANUAL / "MD_MAP.csv").open("w", newline="") as handle:
    writer = csv.writer(handle, lineterminator="\n")
    writer.writerow(["path", "disposition", "pattern_ids"])
    writer.writerows(rows)

counts = {}
for _, disposition, _ in rows:
    counts[disposition] = counts.get(disposition, 0) + 1
print(" ".join(f"{key}={counts[key]}" for key in sorted(counts)))
