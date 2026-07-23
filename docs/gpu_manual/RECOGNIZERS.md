# Recognizer lifecycle

Recognizer output is an advisory worklist, never a correctness oracle. The
version-3 table in `recognizers/table.csv` records candidate and post-state
expressions, exclusions, confirmation questions, and notes. Manual and
build-error rows intentionally have no lexical expression.

For each port and pinned source revision, run:

```bash
python3 docs/gpu_manual/recognizers/scan.py \
  --root <closure-root> --table docs/gpu_manual/recognizers/table.csv \
  --port-id <stable-port-id> --source-revision <git-hash-or-content-id> \
  --out <port-evidence>/COVERAGE.csv [--previous <prior-port-coverage.csv>]
```

Copy `templates/COVERAGE.csv` only for its header; never reuse its example row.
Every candidate hit becomes a row in the inspection/port ledger, or is reviewed
and retained as `false-positive`/`not-applicable`. Compiler diagnostics and the
recorded inspection ledger outrank regex hits. A compiler or inspection finding
that the scanner missed is a candidate rule to add to the table; a scan hit with
no ledger row is unexplained and blocks closure. Rerun after each disposition;
completed ports end with zero `candidate` rows, while unresolved candidates are
explicitly visible.

When `--root` is the repository root, traversal skips `.git`, `ext`, `bin`,
`build`, and `docs` components to avoid vendored/generated material. Pointing
`--root` directly at an isolated closure still scans that directory's own files.

Use `--previous` with a prior site CSV only when both `port_id` and
`source_revision` are unchanged. A different port or revision cannot inherit
those dispositions: every source revision reopens candidates so a newly
device-reachable site cannot retain a stale `not-applicable` decision. Candidate
and converted rows are always recomputed.
Site IDs hash the pattern and matched text (plus stable duplicate occurrence),
so line movement alone does not discard a reviewed disposition. Exclusions emit
`not-applicable` rows rather than hiding source sites, and converted matches do
not erase candidate matches in the same file.

There is no root `COVERAGE.csv`: a static repository snapshot is stale by
construction and must not select a port's files or patterns. A task may retain a
revision-bound scan as evidence, but only the consuming port's report is its
status worklist.
