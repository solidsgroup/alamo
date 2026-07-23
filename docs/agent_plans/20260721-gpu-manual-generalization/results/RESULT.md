# Result: GPU manual generalization

## Outcome

The documentation and recognizer-tooling scope is implemented. The manual now
separates invariant semantics, port contracts, and chamber-gpu corpus examples;
ships generic onboarding, validation, efficiency, decision, status, and harvest
schemas; and treats regex as a versioned advisory status instrument.

The larger port-level definition of done is not claimed. A completed non-Flame
pilot, its filled oracle deck, the expanded closed-book onboarding run, and
named-GPU runtime/profile evidence remain explicitly `pending-pilot` or
`blocked`.

## Workstream results

| Workstream | Result | Evidence |
|------------|--------|----------|
| WS-A | implemented | INDEX legend; all 25 Tier 1 files use Invariant/Port contract/Corpus example; primary source map; evidence-ledger labels |
| WS-B | contract implemented; pilot pending | VALIDATION.md and template cover strict build, analytic/exact, conservation, golden, restart, multi-box, sanitizer, and tolerance rationale |
| WS-C | implemented; pilot exercise pending | ONBOARDING.md; scope, closure, inspection, decision, and status templates; two-stage closed-book gate in STATUS.md |
| WS-D | implemented | recognizer/coverage schema v2, stateful line-addressable scanner, durable dispositions, feedback procedure, seven tests |
| WS-E | implemented | one policy home with owner decisions and worked examples for error propagation, scaffolding retirement, and FEATURE/[NUM] surfacing |
| WS-F | implemented | physics-agnostic profiling procedure and six-item baseline checklist; four former performance patterns relabeled optimization |
| WS-G | implemented; cross-port evidence pending | draft/verified/cross-family status, independent port axes, harvest template, and blind-spot register |

## Counts

- Tier 1: 25 valid patterns; 2 verified, 23 draft, 0 cross-family.
- Recognizer lifecycle: schema v2; legacy repository scan has 211 sites: 155
  candidate and 56 converted across 13 hit-producing regex patterns.
- Templates: eight (`SCOPE`, `CLOSURE`, `INSPECTION_LEDGER`, `VALIDATION`,
  `EFFICIENCY`, `FEATURE_DECISIONS`, `HARVEST`, `PORT_STATUS`).
- Tier 0: approximately 800 tokens, below the 1500-token limit.

## Verification

```text
python3 docs/gpu_manual/build/validate_patterns.py
  VALIDATED_PATTERNS=25
  TRANSFORMS_VERIFIED=2/25; CROSS_FAMILY=0/25
  RECOGNIZER_ORACLE=advisory-lifecycle; formal precision/recall out-of-scope
  PILOT=required

python3 -m unittest discover -s docs/gpu_manual/recognizers -p 'test_*.py'
  Ran 7 tests ... OK

python3 -m py_compile docs/gpu_manual/recognizers/scan.py docs/gpu_manual/build/validate_patterns.py
  PASS

git diff --check
  PASS
```

The independent mechanical audit reproduced the validator and tests and found
no high/medium structural blocker. Its Tier-0 primary-source visibility finding
was fixed before closeout.

## Repository gate baseline

`benchmark/status.sh` reported device-lint PASS and golden-compare FAIL. The
golden log contains both a `canonical_step1/cpu` exit 131 and a later PASS
banner, so the harness result needs separate adjudication. The A100 sanitizer is
blocked because `bin/alamo_gpu-3d-cuda86-g++` is absent. This docs-only task did
not reclassify either condition or modify source.

## Deviations and ownership

- Held-out recognizer labels and formal precision/recall from the earlier v2
  draft were deliberately not built; the user's plan excludes team-scale
  metrics and assumes one maintainer.
- No integrator-specific port plan, closure, oracle deck, or optimization was
  invented. That makes the non-Flame acceptance item an explicit external
  dependency, not a silent pass.
- No commits were created in the shared dirty worktree, and the user-owned
  `docs/gpu_manual.zip` and preceding reliability task folder were not changed.

