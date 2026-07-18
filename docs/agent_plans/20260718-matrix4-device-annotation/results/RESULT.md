# RESULT — matrix4-device-annotation (2026-07-18)

## What changed

- src/Set/Matrix4_Isotropic.H: AMREX_GPU_HOST_DEVICE added to element-access
  `operator ()` (retrieval-only, pure mu/lambda arithmetic).
- src/Set/Matrix4_Diagonal.H: same, on the pure array-read `operator ()`.

## Why

chamber-gpu at dd0f86c55 did not compile for CUDA (any DIM): 5ae7dffa3's
conservative-face-flux branch in Operator::Elastic::Diagonal() calls both
element-access operators inside an AMREX_GPU_DEVICE lambda; nvcc rejects the
host-only members (Elastic.cpp:459-460). Found 2026-07-18 while building the
2D CUDA binary to run the GPU suite before pushing the 14-commit backlog.

Process finding: the newest CUDA binary predated the recovery merge by 9
days (3D, Jul 8) — no CUDA compile had ever exercised the unpushed set, and
the a100-sanitizer PASS gate line was stale evidence from the Jul 8 binary.

## Evidence

- 2D fast, 2D strict (nofast), 3D fast CUDA builds: all exit 0, zero errors
  after the two annotations (no further hidden device-compile breaks).
- GPU suite (tests/GPU/run_gpu_tests.py) on fresh binaries: 9 passed,
  0 failed, 0 skipped — C1 elastic parity, C2 restart roundtrip, C3
  multibox stress, C4 AMR correctness, F1/F2 smoke, P1-P3 perf.
  (P* thresholds informational on the shared 50W A1000.)
- ci_golden_compare.sh after CPU rebuild with edited headers: PASS all 4
  cases (macro is empty on host builds; verified, not assumed).

## Deviations from plan

None. Two-line fix sufficed; Step 1's "iterate on further errors" loop
terminated after the first rebuild.
