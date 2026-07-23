# Fresh antagonistic review

Reviewer scope: frozen retained 2/2 configuration, rejected Step 3/4 artifacts,
Step 6 memory gate, Step 7 launch trace, and protected source cleanliness.

## Findings and adjudication

1. **Medium: Step 7 trace did not initially persist command/binary attribution.**
   Confirmed as an evidence-packaging defect. The trace SQLite itself embeds the
   exact frozen baseline executable path and A1000 UUID; the executable is
   read-only and its current SHA-256 matches the frozen baseline checksum. Plot
   metadata and embedded process arguments record the effective 2/2 and two-step
   parameters. These facts and all artifact hashes are now recorded in
   `step7/nsys/2d_conservative_2x2/ATTRIBUTION.md`. The trace is used only for the
   narrow launch-count conclusion, not as a new performance baseline.

2. **Low: 2/2 stability evidence is only two steps.** Confirmed. The accepted
   claim is limited to the frozen `max_step=2` correctness/physics horizon. The
   3D final residual is `9.72181e-09`, within the `1e-08` absolute gate but closer
   to it than 4/4 (`4.45907e-09`). No longer-evolution stability claim is made.

3. **Low: rejected Step 4 timing provenance is incomplete.** Confirmed. The
   isolated binaries/objects and per-run outputs exist, but there is no complete
   timing manifest/checksum file. Because the candidate was rejected on a
   conservative fail-fast screen and fully reverted, its figures are not used to
   support the retained result. No remediation run is warranted.

## Verified non-findings

- Protected Elastic source and `Matrix4_Major.H` match HEAD; the Step 3 test file
  is absent and both rejected source experiments are fully reverted.
- Frozen Step 2 binary hashes match the baseline checksum manifest; timing arms
  were interleaved with only smoothing overrides changed. The 2D/3D 2/2 gains
  are supported by both profiled and unprofiled modes.
- Step 3 rejection is supported by isolated checksums and its measured threshold
  miss/regression.
- Step 6 decisively fails the numerical memory-headroom gate.
- Step 7 contains exactly 8,904 FApply kernel launches for 8,904 FApply calls;
  no per-box launch multiplicity exists to fuse in the target case.

## NOVA/A100 evidence follow-up

The same reviewer independently reconstructed the copied A100 artifacts and
found no blocker to the narrow confirmation or retaining 2/2. It reproduced
the five-run medians/MADs and call counts, verified the exact HEAD/binary/input/
source hashes and actual A100 UUID, confirmed that only the smoothing overrides
differ, and confirmed the two-step physics gate PASS.

The original timing job's nonzero status is not a benchmark failure: all timing
repetitions were complete before a validation-wrapper API `TypeError`. The
separate validation tail completed successfully against those preserved timing
artifacts. The compact local bundle intentionally omits the 3 GB raw validation
plot payloads retained on NOVA; copied metrics, field norms, commands, logs,
manifests, comparison report, and profiler evidence are present and verified.
