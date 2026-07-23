# Adversarial review

## Scope

Review targeted instance-as-invariant leakage, conflicting policy homes,
scanner overclaiming, schema drift, Flame-shaped verification, misleading pass
language, missing retirement/harvest obligations, and broken manual references.

## Findings fixed during integration

1. Scanner site IDs initially omitted the file path, allowing identical matches
   in different files to share a disposition. IDs now include relative path,
   pattern, match, and duplicate occurrence; a cross-file test covers this.
2. The first root scan included dependency and documentation trees. Traversal
   now skips `.git`, `ext`, `bin`, `build`, and `docs` unless one is itself the
   explicitly selected root.
3. Several converted recognizers were too broad. GPU-017 and GPU-023 were
   tightened; unsound GPU-006/GPU-019/GPU-020 post-state regexes were removed.
4. Candidate and converted regexes could label the same span twice. Candidate
   now wins and remains visible; a regression test covers the overlap.
5. Historical phase-6 output used unqualified `GATE=pass`. It now reports
   `PATTERN_APPLICATION_GATE=pass` with GPU-007/GPU-016 file-conversion scope.
6. The first status wording required completed-port material for `verified`,
   contradicting the historical file-level evidence. `verified` now means a
   successful closed-book target exercise; `cross-family` requires two
   completed-port harvests.
7. Two Tier-0 invariant bullets lacked visible grounding. INDEX now links the
   relevant primary-source anchors and separates CPU-path preservation as a
   port contract.

## Final assessment

The scoped structural oracle and seven scanner tests pass. The delegated
cross-file audit found no remaining high/medium documentation or schema defect,
no port-specific correctness Verify line, and no unqualified global pass.
Policy ownership, retirement, baseline efficiency, harvest, and scanner
precedence are consistent at their points of use.

The following are acceptance dependencies, not review defects:

- no completed non-Flame validation-contract instantiation;
- no pilot run of the expanded onboarding closed-book stage;
- no completed-port recognizer harvest or cross-family status;
- current named-GPU sanitizer evidence blocked by the absent binary;
- current repository golden gate reported failed/contradictory evidence.

