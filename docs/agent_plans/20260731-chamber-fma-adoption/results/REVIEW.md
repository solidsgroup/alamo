# REVIEW — chamber-fma-adoption

Reviewer: fresh adversarial reviewer with no campaign context

Final disposition: **ACCEPT WITH THE EXPLICITLY AUTHORIZED STALE-REFERENCE
EXCEPTION**

## Source review

No source defect was found.

- `BC.H` performs periodic/inter-box exchange before deciding whether a tile
  needs a physical-boundary fill. It grows by the exact MultiFab ghost width,
  converts the physical domain to the tile index type, and skips only when the
  grown region is wholly inside the physical domain. Physical-edge tiles keep
  the existing BC path.
- `BaseField.H` exchanges internal ghosts first. Its skip is confined to the
  cell-field path, whose allocation and domain are cell-indexed, and uses the
  same ghost width as the subsequent physical fill.
- The source-only patch hash independently matches
  `6369167a6b823b60760d6a258cdbc71bc6d2dd72a29f96418e2614183f634585`.
  No test, reference, tolerance, or benchmark-oracle file changed.

## Findings and adjudication

### 1. Strict target oracle remains red

Initial disposition: blocking under the literal task oracle.

Adjudication: the four `rod_and_tube_step2/gpu_strict` traction mismatches were
independently reproduced with identical deltas at pre-adoption target
`508a8785d`. The prior thermoelastic-chain-rule task already classified that
reference as stale. The user explicitly authorized continuing under that
adjudication. No reference or test was changed, and the exception remains
visible rather than being called a pass.

Final disposition: authorized exception; no evidence of an adoption regression.

### 2. Sanitizer harness accepted aborted applications

Initial disposition: blocking. The first logs reported zero sanitizer findings
but the application later aborted in a deliberately under-iterated elastic
solve. The repository wrapper checks the sanitizer summary without also
requiring the application return code.

Remediation: supersede those logs with exact-candidate runs that stop after five
completed steps, before the under-iterated solve. Tier 1, memcheck, initcheck,
and racecheck each record an actual return code of zero, reach normal AMReX
finalization, and contain no abort. Memcheck/initcheck report zero errors;
racecheck reports zero hazards, errors, and warnings.

Final disposition: closed. The reviewer validated all return-code and log files.
The short run is correctly scoped as a memory-safety smoke; separate strict and
800-step runs cover scientific behavior.

### 3. Physics comparison lacked manifest enforcement

Initial disposition: high. The first physics comparison did not use
`--require-compatible-manifest`, and the extracted bundles lacked manifests.

Remediation: supplement each bundle from captured job, per-arm, source, build,
GPU, runner, input, override, and oracle provenance. Rerun
`compare_validation.py --require-compatible-manifest --gate`.

Final disposition: closed. The reviewer independently reran the command; it
exited zero and reproduced `COMPARE_MANIFEST_GATED.json` byte-for-byte.

### 4. Partial result was stale after resumed work

Initial disposition: medium.

Remediation: replace the paused result with the completed evidence and retain
the invalid initial jobs as explicitly quarantined records.

Final disposition: closed by this closeout.

## Evidence verification

- Valid timing job `11824925`: all 16 rows return zero; 7 measured repetitions
  per arm; `0.060375` to `0.053175 s/step`, an 11.93% median reduction.
- Valid correctness job `11824926`: both rows return zero; all 32 scientific
  plot/checkpoint files and full `thermo.dat` are byte-identical.
- Manifest-gated physics result: PASS, with 21 correctness and 12 engineering
  trajectory rows passing and identical solver counts.
- Invalid jobs `11824918` and `11824919`: correctly quarantined; no result used.
- Later commit `962eb4b3c` does not overlap either adopted file and is outside
  the exact target-bound differential.

## Reviewer recommendation

There is no remaining sanitizer, manifest, correctness, performance, or source
blocker. Retention is acceptable under the user's explicit authorization of the
pre-existing stale strict-GPU reference; that exception should continue to be
tracked separately.
