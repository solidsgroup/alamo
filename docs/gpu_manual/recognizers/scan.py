import argparse
import csv
import re
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(usage="%(prog)s --root DIR --table CSV --out CSV")
    parser.add_argument("--root", required=True, type=Path)
    parser.add_argument("--table", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()

    with args.table.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    regexes = [
        (row["pattern_id"], re.compile(row["expression"], re.MULTILINE))
        for row in rows
        if row["type"] == "regex"
    ]

    hits = []
    extensions = {".H", ".cpp", ".cu", ".cc"}
    for path in sorted(p for p in args.root.rglob("*") if p.suffix in extensions):
        source = path.read_text(errors="ignore")
        relative = path.relative_to(args.root).as_posix()
        for pattern_id, expression in regexes:
            count = sum(1 for _ in expression.finditer(source))
            if count:
                hits.append((relative, pattern_id, count))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.writer(handle, lineterminator="\n")
        writer.writerow(("file", "pattern_id", "hits"))
        writer.writerows(hits)


if __name__ == "__main__":
    main()
