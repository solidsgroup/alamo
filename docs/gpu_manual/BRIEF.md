# Agent Brief: GPU Pattern Manual Builder

Mission: mine the `chamber-gpu` branch of ALAMO (diffs + markdown failure/success library) and produce a token-efficient, three-tier pattern manual that a fresh LLM session can use to convert files on the `gpu` branch without access to `chamber-gpu`. This brief is the agent's complete instruction set. Runtime-agnostic (Claude Code, Codex CLI, or similar).

---

## 0. Kickoff prompt (paste into the agent session)

```
You are the manual-builder agent for the ALAMO repository. Read
PATTERN_MANUAL_AGENT.md in full before acting. Execute Phases 0 through 6 in
order; do not proceed past a phase until its exit criteria are met and
reported. All writes go under docs/gpu_manual/ on a new branch named
manual-build. The branches chamber-gpu, gpu, and the main branch are
read-only. Begin with Phase 0 and report its inventory before continuing.
```

---

## 1. Ground rules

1. Evidence hierarchy, strongest first:
   a. Code diffs on `chamber-gpu` (ground truth)
   b. Build/runtime results the agent reproduces
   c. The markdown library (explains why; may be stale)
   d. Inference
   Never write a pattern supported only by (d). Prose that contradicts a diff is stale: drop it and log the drop.
2. Schemas in Section 3 are frozen. Do not restructure them. No narrative documents anywhere in the manual.
3. Write boundary: `docs/gpu_manual/**` on branch `manual-build` only. Source branches read-only. Never run destructive git (`reset --hard`, `clean -fdx`, force push, branch deletion).
4. Token budgets: `INDEX.md` <= 1500 tokens. Each pattern file 200-400 tokens. Evidence tier unbounded but never loaded by default.
5. Every pattern claim carries evidence: commit hash and/or markdown path.
6. Ask the user when: branch names differ from those assumed here, the markdown library location is ambiguous, or a stop condition in Section 6 fires.

---

## 2. Deliverables and layout

```
docs/gpu_manual/
  INDEX.md                     Tier 0: invariants, triage tree, pattern one-liners
  patterns/GPU-NNN-<slug>.md   Tier 1: one file per pattern
  evidence/<topic>.md          Tier 2: full narratives, benchmarks, links
  recognizers/table.csv        pattern_id,type,expression,notes
  recognizers/scan.py          recognizer runner (spec in Phase 4)
  COVERAGE.csv                 file,pattern_id,hits  (scan of gpu branch)
  ONE_OFFS.md                  meaningful changes that are not patterns
  BUILD_LOG.md                 decisions, drops, conflicts
  build/                       scratch: diffs, inventories, ledgers
```

`COVERAGE.csv` doubles as the port plan for the `gpu` branch: task list, effort estimate, parallelization boundaries.

---

## 3. Frozen schemas

### 3.1 Tier 1 pattern file (`patterns/GPU-NNN-<slug>.md`)

```
# GPU-NNN: <imperative name>
Status: draft | verified
Recognizer: <grep -E / regex for the bad shape; MANUAL if not greppable>
Applies: <code shape, one line>
Transform:
  Before:
    <minimal snippet>
  After:
    <minimal snippet>
Constraints: <when NOT to apply>
Verify: <build/test command + expected result>
Failure modes: <compile/runtime signature if misapplied or if pattern absent>
Evidence: <commit hash(es); markdown path(s)>
```

Anti-patterns (documented dead ends with no surviving diff) use the same schema with `Transform: avoid` and the failed approach shown under `Before`.

### 3.2 Tier 0 (`INDEX.md`)

```
# GPU Pattern Index
## Invariants
<rules with zero exceptions, one line each>
## Triage
<code symptom> -> GPU-NNN[, GPU-MMM]
## Patterns
GPU-001: <one line>
...
## One-offs
See ONE_OFFS.md
```

### 3.3 `ONE_OFFS.md` entry

```
- <file>:<hunk range> | <what changed, one line> | commit <hash> | <why not a pattern>
```

### 3.4 Ledgers (in `build/`)

```
HUNK_MAP.csv:  hunk_id,file,classification   classification in {GPU-NNN, ONEOFF, NOISE}
MD_MAP.csv:    path,disposition,pattern_ids  disposition in {mapped, anti-pattern, evidence-only, stale-dropped}
```

### 3.5 `BUILD_LOG.md` entry

```
<date> | <phase> | <decision> | <evidence> | <artifacts affected>
```

---

## 4. Phase plan

### Phase 0: Recon (read-only)

1. Verify branch names. Assumed: `chamber-gpu`, `gpu`, and a main branch (`main` or `master` or `development`). If assumptions fail, stop and ask.
2. `BASE=$(git merge-base <main> chamber-gpu)`
3. `git diff --stat $BASE..chamber-gpu > docs/gpu_manual/build/DIFF_STAT.txt`
4. One diff file per changed source file into `build/diffs/` (`git diff $BASE..chamber-gpu -- <path>`).
5. Locate the markdown library. List every file with path and size into `build/MD_INVENTORY.txt`. If location ambiguous, ask.

Exit: DIFF_STAT, diffs/, MD_INVENTORY exist. Report counts (files changed, hunks approx, md files) to user.

### Phase 1: Scaffold

Create the Section 2 layout on branch `manual-build`. Empty scaffolds plus this brief copied to `docs/gpu_manual/BRIEF.md` for provenance.

Exit: layout committed.

### Phase 2: Diff mining

1. Split each per-file diff into hunks. `hunk_id = <file>:<old-start>-<old-end>`.
2. Classify every hunk:
   - Transformation class (candidate pattern): cluster similar hunks across files, assign `GPU-NNN`.
   - `NOISE`: formatting, include reordering, comment-only.
   - `ONEOFF`: meaningful but unique. Record per Section 3.3.
3. Draft one Tier 1 file per cluster. Use the smallest clean hunk in the cluster for Before/After. `Status: draft`. Populate Evidence with commit hashes (`git log -L` or `git blame` on chamber-gpu as needed).

Exit: `HUNK_MAP.csv` complete, zero unassigned hunks. Report cluster count and ONEOFF/NOISE counts.

### Phase 3: Markdown mining

1. For each md file: map content to existing pattern IDs. Extract into those patterns: Failure modes, Constraints, Verify commands. Append path to Evidence.
2. Content matching no diff-derived pattern:
   - Documented dead end -> anti-pattern file (Section 3.1 variant).
   - Useful context, not actionable -> `evidence/<topic>.md`, linked from nearest pattern.
   - Contradicts a diff -> stale. Drop, log in BUILD_LOG.
3. Compress: no md prose is copied verbatim into Tier 1; rewrite to schema fields.

Exit: `MD_MAP.csv` complete, every md file dispositioned.

### Phase 4: Recognizers and coverage

1. Fill `recognizers/table.csv`. Prefer regex; `MANUAL` only where the shape is not greppable.
2. Write `recognizers/scan.py`. Spec: Python 3 stdlib only; inputs `--root <dir> --table <csv> --out <csv>`; walks source files (verify extensions in repo, expect `.H .cpp .cu`); applies regex rows; emits `file,pattern_id,hits`. No comments beyond one usage line.
3. True-positive check: scan the `chamber-gpu` pre-transform tree (`git worktree add build/base-scan $BASE`). Every regex recognizer must hit the code its pattern later transformed. Zero hits here = broken recognizer; fix before proceeding.
4. Scan the target: `git worktree add build/gpu-scan gpu`, run scan, write `COVERAGE.csv`. Remove worktrees after (`git worktree remove`).
5. Patterns with zero hits on `gpu`: keep, mark in BUILD_LOG as `n/a this port` unless recognizer is at fault.

Exit: COVERAGE.csv committed, per-pattern hit counts reported.

### Phase 5: Tier 0 assembly (last, not first)

1. Invariants: distill from patterns whose Constraints show no exceptions.
2. Triage tree: build from Recognizer/Applies lines, symptom -> pattern IDs.
3. One-liners for every pattern.
4. Token check: approx tokens = words x 1.3. Trim until <= 1500.

Exit: INDEX.md within budget, committed.

### Phase 6: Closed-book validation gate

1. From COVERAGE.csv pick the smallest `gpu`-branch file with >= 2 distinct pattern hits.
2. Fresh session receives only: INDEX.md, the hit pattern files, the target file. No repo access, no chamber-gpu.
3. Session produces the converted file.
4. Lead (with chamber-gpu access) reviews: correct transforms applied, constraints respected, nothing invented.
5. Any gap: patch the responsible pattern, log in BUILD_LOG, rerun on a different file.
6. Gate passes on two consecutive clean runs on different files. Flip validated patterns to `Status: verified`.

Exit: gate passed, results logged.

---

## 5. Model routing (optional, single-agent fallback is valid)

- Cheap tier, batchable: Phase 0 inventories, Phase 2 hunk splitting, Phase 3 per-file md mapping.
- Mid tier: Phase 2 clustering and pattern drafting, Phase 4 recognizer authoring.
- Top tier: Phase 5 distillation, Phase 6 review judgment.

---

## 6. Stop conditions

- Diff and markdown conflict -> prefer diff, log, continue.
- Hunk unclassifiable after two attempts -> ONEOFF with note, continue.
- More than 20 percent of non-NOISE hunks landing in ONEOFF -> stop, report; clustering is probably wrong.
- Any write needed outside `docs/gpu_manual/`, or any destructive git -> stop, ask.
- Markdown library or branch names not found as assumed -> stop, ask.

---

## 7. Definition of done

- [ ] Every diff hunk in HUNK_MAP.csv as GPU-NNN, ONEOFF, or NOISE
- [ ] Every md file dispositioned in MD_MAP.csv
- [ ] Every regex recognizer passes the true-positive check on the BASE tree
- [ ] COVERAGE.csv generated from the `gpu` branch
- [ ] Zero-hit patterns dispositioned in BUILD_LOG
- [ ] INDEX.md <= 1500 tokens; every pattern file 200-400 tokens
- [ ] Every pattern has populated Evidence
- [ ] Closed-book gate: two consecutive passes; validated patterns marked verified
