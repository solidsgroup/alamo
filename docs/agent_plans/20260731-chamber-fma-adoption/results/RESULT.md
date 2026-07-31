# RESULT — chamber-fma-adoption

Status: **COMPLETE WITH AN AUTHORIZED, PRE-EXISTING ORACLE EXCEPTION**

## Recommendation

Retain and commit the selectively transplanted OPT-10 source pair on
`chamber-gpu`. Fresh target-bound evidence shows an 11.93% median reduction in
800-step wall time with byte-identical scientific output, identical convergence
and physics observables, and clean sanitizer runs. Do not merge the campaign
branch or adopt any rejected configuration.

The only red comparison is the already-known strict-GPU
`rod_and_tube_step2` traction reference. It was reproduced with identical
deltas at the pre-adoption target, so it is not caused by this patch. The user
explicitly authorized continuing under the prior stale-reference adjudication;
no test, tolerance, or reference was changed.

## Adopted source state

The pre-adoption `chamber-gpu` target was
`508a8785d972191235a80394614fc2e33a7ef57c`. The adoption plan was committed as
`c1f49caba6f415dd301d189f973cff993482d3c5`, followed by the reviewed source
pair:

- `8110c5e5c5380230e65e651b6b99f35e9d0698bc` —
  `perf: skip interior physical boundary fills`
- `f73154f7e94e8b45a92d3bf70cffa0c1448b6d0c` —
  `fix: preserve index type in BC interior skip`

The source diff is exactly 13 insertions in `src/BC/BC.H` and
`src/Integrator/BaseField.H`. Its patch SHA-256 is
`6369167a6b823b60760d6a258cdbc71bc6d2dd72a29f96418e2614183f634585`,
identical to the frozen campaign patch. No test, benchmark oracle, tolerance,
or stored reference changed in the adoption differential.

While target-bound validation was running, the branch independently advanced
to `962eb4b3c` (`fix: restore dynamic mechanics parser branches`). That commit
changes `src/Integrator/Base/Mechanics.H`, does not overlap either adopted file,
and is outside the `508a8785d..f73154f7e` A/B differential. The claims below
apply to the exact post-adoption candidate `f73154f7e`; the final repository
status gate separately covers the current composite branch head.

## Local correctness

- Fresh strict CPU gate:
  `GOLDEN_MODE=cpu BUILD_JOBS=8 benchmark/ci_golden_compare.sh` passed all four
  golden cases and the NaN tripwire.
- The strict CUDA build used the self-contained CUDA 12.6.3 redistributable
  because the installed CUDA 12.0 `nvcc` wrapper delegates to a missing system
  compiler.
- Strict GPU cases `canonical_step1`, `canonical_step2`, and
  `eta_expression_step1` passed, and the strict-GPU NaN tripwire exited zero.
- `rod_and_tube_step2/gpu_strict` retains these four traction-reference
  mismatches:

| Observable | Absolute delta | Relative delta |
| --- | ---: | ---: |
| `trac_xhi_x` | `1.384000e+02` | `2.265232e-03` |
| `trac_xhi_y` | `4.585000e-01` | `2.349704e-02` |
| `trac_yhi_x` | `6.512000e-01` | `2.642621e-01` |
| `trac_yhi_y` | `1.256000e+02` | `2.100370e-03` |

An independent detached-worktree run at pre-adoption target `508a8785d`
reproduced the same four failures and deltas. This is the authorized exception;
the strict GPU gate is not represented as fully green.

## Clean-exit sanitizer evidence

The first memcheck/initcheck/racecheck attempt reached zero-error summaries but
then aborted in a deliberately under-iterated elastic solve. Fresh adversarial
review correctly rejected those logs as sanitizer passes.

The superseding runs use the exact adopted local binary
(`3be02c88710df79635b3c7a3d043035dd779efe2e66b0888d7fbb86a4a4c8051`),
`input_3d_centre_bore_128_a2`, MGS64, and `stop_time=5e-4`. This stops after
five completed steps, before that deliberately under-iterated solve. Every
underlying command's actual return code is recorded and required to be zero;
every log reaches normal AMReX finalization with no abort:

| Check | Return code | Result |
| --- | ---: | --- |
| Runtime-strict tier 1 | 0 | PASS |
| Compute Sanitizer memcheck | 0 | `ERROR SUMMARY: 0 errors` |
| Compute Sanitizer initcheck | 0 | `ERROR SUMMARY: 0 errors` |
| Compute Sanitizer racecheck | 0 | `0 hazards`, `0 errors`, `0 warnings` |

The CUDA 12.6.3-built binary/runtime was inspected by the available
Compute Sanitizer 2022.4.1 installation. This is a pre-elastic memory-safety
smoke, not a convergence test; strict correctness and the completed 800-step
NOVA differential provide the scientific and solver coverage.

The sanitizer run directory was initially named with a `20260730` artifact
date, so the recorded `plot_file` arguments retain that path. The complete
directory was moved without changing its command, log, or return-code files to
the correct `20260731` artifact root; their hashes bind the relocation.

## Target-bound NOVA differential

Both A100 binaries were built from isolated, hash-verified source trees with
matched 2-D, sm_80, plain/no-fast-math flags:

| Arm | Commit | `src` tree | Binary SHA-256 |
| --- | --- | --- | --- |
| Baseline | `508a8785d` | `f8e977da18b3d5b6e051f5fabb53f251daf54b37` | `85d2babe6a364604a431ad49c00f5cd15415b4af51634ae67e163d137db96cf0` |
| Candidate | `f73154f7e` | `2d6a6c043f17fa2aab6fdf9a413df2700f839d1e` | `a135e0257d8191238335a74a9ecefc06c6fe21e3eaa565d23ba170a758e14a17` |

The shared `input_copy` SHA-256 is
`afbbd60a4573727eb2df87dee874d23c7ae380011301ae016c7e759410b0e7f2`.
Both arms use MGS32/BF8, 4/4 smoothing, synchronized MLMG, and 800 steps.

NOVA job `11824925` completed all 16 timing rows with return code zero: one
warmup pair and seven interleaved measured repetitions per arm.

| Arm | Median wall/step | MAD | Samples |
| --- | ---: | ---: | --- |
| Baseline | `0.060375 s` | `0.0003875 s` | `0.060350, 0.060800, 0.0613125, 0.0599875, 0.061300, 0.060075, 0.060375` |
| Candidate | `0.053175 s` | `0.0001625 s` | `0.053800, 0.0530125, 0.054975, 0.053175, 0.0532625, 0.053025, 0.0527625` |

The candidate reduces median wall/step by **11.9255%**.

NOVA job `11824926` completed both output-enabled correctness arms with return
code zero. All 32 initial/final scientific plotfile payloads per arm have the
same aggregate hash
`1fff7943d9e3c9f879d67fe5bc37a168739217214d1607e85b425097d2bcdfb3`.
Both complete `thermo.dat` files have SHA-256
`b44d7c4bc927ae570802e86c49e73289547c6a3d54e22cbcb3a21c115576e90c`.

The physics comparator was rerun with
`--require-compatible-manifest --gate`. The supplemented manifests bind the
same host, A100 UUID, CUDA/toolchain/build flags, runner, input, overrides, and
oracle hashes while preserving distinct baseline/candidate source and binary
identities. It passes 21 correctness observables and 12 engineering
trajectories with zero differential. Solver health is also identical: 20
Newton iterations and 448 MLMG V-cycles.

The first timing/correctness submissions, jobs `11824918` and `11824919`, are
quarantined. Every arm exited 6 before simulation because an external bitmap
asset was not staged; their empty summaries are not used.

## Review and evidence

Fresh adversarial review found no defect in the physical-vs-interior BC
semantics, periodic/inter-box ordering, index conversion, ghost width, or
cell-field path. It initially rejected closeout because the aborted sanitizer
logs and missing manifest enforcement were invalid evidence. Both findings were
accepted and remediated. Focused re-review independently reproduced the
manifest-gated PASS and accepted the clean-exit sanitizer evidence. It found no
remaining blocker beyond the user-authorized stale GPU-reference exception.

The external evidence root is:

`/home/jackplum/Projects/alamo-chamber-fma-adoption-artifacts-20260731`

`results/KEY_EVIDENCE_SHA256SUMS` binds the key files under that root.

## Retained configuration and limitations

- Retained: OPT-10 plus the mandatory index-type repair, current MGS32/BF8
  grid, 4/4 smoothing, synchronized MLMG.
- Not adopted: MGS64, MGS128, 2/2 smoothing, or no-sync.
- The 256-cubed campaign case remains outside A100-80GB capacity.
- The sanitizer smoke stops before an elastic solve; the 800-step differential
  and strict golden runs cover solver/scientific behavior.
- The unrelated later `962eb4b3c` mechanics-parser change was not part of the
  performance differential.
