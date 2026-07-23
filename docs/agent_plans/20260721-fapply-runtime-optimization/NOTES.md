# Execution notes

## 2026-07-21 — Step 0 workspace freeze

- Human approval: the user instructed the orchestrator to execute the adversarially reviewed plan on 2026-07-21. This satisfies the Step 0 written-approval gate.
- Repository: `/home/jackplum/Projects/alamo`
- Branch: `chamber-gpu`
- HEAD: `2e6a8f8f5430a58dd6a85b07b08d472e2bd1d5cd`
- Hooks: `.githooks`
- Scoped status before campaign edits:
  - `?? benchmark/validate/runs/`
  - `?? docs/agent_plans/20260721-fapply-runtime-optimization/`
- Protected source diff: clean for `src/Operator/Elastic.cpp`,
  `src/Operator/Elastic.H`, and `src/Set/Matrix4_Major.H`.
- `benchmark/status.sh`: device-lint PASS; the script reported 101 dirty files in
  the broader pre-existing worktree. Those artifacts are outside this task and
  must remain untouched.
- Host compiler: `/usr/bin/g++`, GCC 13.3.0
- CUDA compiler: `.local/cuda/bin/nvcc`, CUDA 12.6
  (`Build cuda_12.6.r12.6/compiler.35059454_0`)
- MPI: Open MPI 4.1.6
- Python: 3.12.3
- AMReX: `ext/AMReX-Codes/amrex` at
  `bf15fdce52093c715a1dd9a9da64ba66f0233be1`; the dependency worktree was
  already dirty (`AMReX_MLCGSolver.H`, `AMReX_MLMG.H`, and generated build
  directories). Tracked binary-diff SHA-256 at freeze:
  `37a7887d5edc174289045065b3c406f405670359c16cc43e31adb0dbb4004f6a`.
- GPU: NVIDIA RTX A1000, UUID
  `GPU-ff00e057-b36d-c833-9da1-8f71516e7d73`, driver 595.71.05, 8188 MiB,
  compute capability 8.6, 50 W limit. Freeze query observed P5, 427 MHz SM,
  810 MHz memory, and no active compute processes.
- Nsight Systems: 2026.1.3.243
- Nsight Compute: 2022.4.1.0 (build 32308335)
- Compute Sanitizer: 2022.4.1
- The separate launch-bounds helper/campaign is not present on this branch and
  will not be recreated or merged during this A/B series.

## 2026-07-21 — Step 1 frozen artifact paths and cases

All paths below were fixed while the Step 1 harness paths were still clean.
Every output directory named here must be absent before its producing command.

- Artifact root:
  `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721`
- Baseline binary root:
  `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/baseline/bin`
- Baseline object root:
  `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/baseline/obj`
- Baseline raw-log root:
  `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/baseline/raw`
- Frozen baseline binaries under the baseline binary root:
  - `alamo_gpu-2d-strict-sm86-baseline`
  - `alamo_gpu-2d-fast-sm86-baseline`
  - `alamo_gpu-3d-strict-sm86-baseline`
  - `alamo_gpu-3d-fast-sm86-baseline`
  - `alamo_gpu-2d-profile-fast-sm86-baseline`
  - `alamo_gpu-3d-profile-fast-sm86-baseline`
- Frozen baseline validation bundles:
  - strict 2D: `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/baseline/validation/strict_2d`
  - strict 3D: `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/baseline/validation/strict_3d`
  - fast 2D: `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/baseline/validation/fast_2d`
  - fast 3D: `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/baseline/validation/fast_3d`
- Frozen baseline SoftVoid oracle run root:
  `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/baseline/oracle/softvoid`
- Candidate Step 3 validation root:
  `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/candidate-step3/validation`
- Candidate Step 3 SoftVoid oracle root:
  `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/candidate-step3/oracle/softvoid`
- Candidate Step 4 validation root:
  `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/candidate-step4/validation`
- Candidate Step 4 SoftVoid oracle root:
  `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/candidate-step4/oracle/softvoid`
- Combined-finalist validation root:
  `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/combined/validation`
- Combined-finalist SoftVoid oracle root:
  `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/combined/oracle/softvoid`

Each candidate validation root has exact leaf directories named `strict_2d`,
`strict_3d`, `fast_2d`, and `fast_3d`; none may exist before its run.

Frozen regimes:

1. Conservative 2D correctness and performance: `tests/ElasticSoftVoid/input`
   with `max_step=2`, the input's `elastic.use_psi=0` and 4/4 smoothing, and
   the exact `SOFTVOID_COMMON` overrides in the approved plan. Single-box and
   `amr.max_grid_size=32` multi-box forms are correctness arms. Timed arms add
   only output-suppression/verbosity overrides shared by A and B.
2. Psi-enabled 3D correctness and performance:
   `input_3d_centre_bore_128_a2`, case
   `centre_bore_3d_128_a2_converged`, `max_step=2`, the manifest defaults,
   `amr.max_grid_size=64`, `amrex.the_arena_is_managed=1`, and
   `elastic.solver.nriters=2`. The input supplies `elastic.use_psi=1` and 4/4
   smoothing. Timed arms add only output-suppression/verbosity overrides
   shared by A and B.
3. `canonical_2d_elastic` is retained only as the required psi-path validation
   case; it is never labeled as the conservative performance regime.

Timing protocol is one warmup plus five measured repetitions per binary and
mode, interleaved A/B for source candidates. Profiler-enabled and unprofiled
runs have separate output roots and are compared only within the same mode.

## 2026-07-21 — Step 1 oracle stop

The literal approved `run_softvoid_gpu` command is not executable as written.
It creates `${out}`, redirects stdout/stderr into that directory, and then
passes the same directory as `plot_file`. ALAMO preserves an existing plot
directory by renaming it to `${out}.old.<timestamp>` before creating the new
plot output. Consequently, the redirected logs move into the `.old` sibling
and `tests/ElasticSoftVoid/test "${out}"` fails with missing
`${out}/stdout` even though the simulation completed.

Observed failed arm:

- Plot output:
  `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/baseline/oracle/softvoid/strict_single`
- Displaced logs:
  `/home/jackplum/Projects/alamo/docs/agent_plans/20260721-fapply-runtime-optimization/artifacts/a1000-sm86-2e6a8f8f-20260721/baseline/oracle/softvoid/strict_single.old.20260721155829`
- Failure: `FileNotFoundError` for `strict_single/stdout` at
  `tests/ElasticSoftVoid/test:102`.

Execution stopped under Operating rule 2. Proposed narrow amendment, pending
human approval: preserve the failed root, use a new absolute retry root named
`baseline/oracle/softvoid_retry1`, redirect logs to a sibling log directory
that is not passed as `plot_file`, copy the completed logs into each plot
output, and then invoke the unchanged property oracle. Binary, input,
overrides, layout, and physics checks remain unchanged.

## 2026-07-21 — Step 1 oracle amendment authorized and passed

The human authorized the narrow log-placement amendment. The failed
`baseline/oracle/softvoid` root was preserved unchanged and the retry used the
new root `baseline/oracle/softvoid_retry1`. Before the retry, the device matched
the frozen A1000 UUID and no competing compute process was reported.

All four amended arms exited zero and passed the unchanged
`tests/ElasticSoftVoid/test` property oracle:

- `strict_single`
- `strict_multibox` with `amr.max_grid_size=32`
- `fast_single`
- `fast_multibox` with `amr.max_grid_size=32`

Each simulation wrote to a previously absent plot directory. Its stdout and
stderr were captured in the sibling `softvoid_retry1/logs/<label>` directory,
copied into the completed plot directory, and then consumed by the oracle. The
input, overrides, binaries, and physics checks were unchanged. The exact shell
trace, binary/input/oracle hashes, oracle exits, artifact paths, and post-run
GPU state are retained in `baseline/raw/softvoid_gpu_oracle_retry1.log`.

## 2026-07-21 — Step 1 frozen validation bundles

The exact frozen strict/fast 2D and 3D binaries produced all four planned
same-A1000 baseline bundles at the absolute paths recorded above. The 2D
bundles use `canonical_2d_elastic`; the 3D bundles use
`centre_bore_3d_128_a2_converged`. Each run extracted 27 observables.

All four manifests were checked for nonempty binary path/hash, full HEAD,
scoped-source-diff hash, build command/flags, host, device/profile, GPU
name/UUID/driver/architecture, oracle hashes, input/override hashes, command
shape, and full executed argv. Each bundle passed a self-comparison with
`--gate --require-compatible-manifest`, establishing that its frozen physics
output is readable and satisfies the comparison gate at zero delta. Raw runner
logs are retained as `baseline/raw/validation_{strict,fast}_{2d,3d}.log`.

## 2026-07-21 — Step 1 baseline timing distributions

Both frozen A1000 performance regimes completed one warmup plus five measured
repetitions in each mode. Unprofiled and profiler-enabled binaries were never
compared across modes. Exact shell traces, input and binary hashes, commands,
and per-run GPU/process snapshots are retained in
`baseline/raw/timing_2d_conservative_4x4.log` and
`baseline/raw/timing_3d_psi_4x4.log`; per-run logs are under
`baseline/timing/{2d_conservative_4x4,3d_psi_4x4}`.

| regime / metric | measured values | median | MAD |
|---|---:|---:|---:|
| 2D conservative unprofiled external wall (s) | 1.75, 1.74, 1.74, 1.82, 1.73 | 1.74 | 0.01 |
| 2D conservative profiled external wall (s) | 1.77, 1.77, 1.79, 1.78, 1.77 | 1.77 | 0.00 |
| 2D profiled `MLMG::solve()` inclusive wall (s) | 0.9772, 0.9701, 0.9819, 0.9695, 0.9663 | 0.9701 | 0.0038 |
| 2D profiled `Operator::Elastic::Fapply()` inclusive wall (s) | 0.4613, 0.4621, 0.4648, 0.4593, 0.4594 | 0.4613 | 0.0019 |
| 3D psi unprofiled external wall (s) | 287.12, 287.87, 287.20, 287.11, 286.95 | 287.12 | 0.08 |
| 3D psi profiled external wall (s) | 286.94, 287.24, 286.65, 287.02, 287.22 | 287.02 | 0.20 |
| 3D profiled `MLMG::solve()` inclusive wall (s) | 252.0, 252.4, 251.7, 252.3, 252.2 | 252.2 | 0.2 |
| 3D profiled `Operator::Elastic::Fapply()` inclusive wall (s) | 226.2, 226.5, 226.1, 226.3, 226.3 | 226.3 | 0.1 |

The `Fapply` call count was invariant across measured runs: 11,808 in the 2D
conservative case and 15,477 in the 3D psi case. No pre-run process snapshot
reported a concurrent compute workload. Several long 3D runs emitted only the
host desktop warning `Authorization required, but no authorization protocol
specified`; all solver exits were zero, timing remained tightly clustered, and
no CUDA, AMReX, convergence, or sanitizer error was present.

## 2026-07-21 — Step 1 Nsight Compute permission stop

Static `cuobjdump` resource reporting completed for both baseline binaries.
The relevant `Elastic<1>::Fapply` rows are `REG:96 STACK:0` in 2D and
`REG:254 STACK:440` in 3D, confirming that the 3D pressure targeted by Step 4
is material. Raw tables are retained in
`baseline/raw/fapply_resources_{2d,3d}.log`.

Runtime counter collection then stopped before launching a profiled solver:
Nsight Compute 2022.4.1 returned `ERR_NVGPUCTRPERM` even for
`--query-metrics`. The read-only driver state reports
`RmProfilingAdminOnly: 1`; the campaign user is non-root but is a member of
the `sudo` group. Consequently, achieved occupancy and hardware-counter spill
evidence cannot be collected as the current user. No Nsight target run was
started, no driver setting was changed, Step 2 has not begun, and the protected
source diff remains clean. Human direction is required on privileged counter
collection versus a documented static-analysis amendment.

## 2026-07-22 — Step 1 privileged NCU retry required

The human authorized root profiling and manually ran the task-local
`run_root_ncu.sh`. Counter and section enumeration succeeded as root, proving
the permissions hurdle resolved, and the binary/input hashes matched the frozen
baseline. Neither 2D nor 3D produced a report because Nsight Compute treated
the instrumented solver as a child process and required
`--target-processes all`. Both applications exited normally; the failed
artifacts are preserved in `baseline/ncu` and
`baseline/raw/ncu_fapply_root.log`.

The retry was amended to use new paths `baseline/ncu_retry1` and
`baseline/raw/ncu_fapply_root_retry1.log`, add `--target-processes all`, and
retain the exact demangled `Operator::Elastic<1>::Fapply` kernel filter. Static
demangling confirmed that filter exists in both frozen binaries. The retry
paths were absent and the script passed `bash -n` and `git diff --check` before
requesting the second privileged invocation.

The child-process retry still missed both kernels because the demangled regex
did not match Nsight Compute's kernel-name representation. It was canceled
during the otherwise normal 3D run after the 2D miss was observed; artifacts
are preserved under `baseline/ncu_retry1` and its matching raw log.

An exact-mangled-symbol retry then successfully captured the frozen 2D
`Elastic<1>::Fapply` kernel into
`baseline/ncu_retry2/2d_fapply.ncu-rep`. The 3D kernel matched, but ten-pass
kernel replay failed with exit code 9 while backing the managed-memory working
set to disk. No 3D report was produced. The successful 2D report and its raw
CSV are preserved with hashes in `baseline/ncu_retry2`.

The next privileged retry is 3D-only at `baseline/ncu_retry3`, uses the same
exact mangled 3D symbol, switches to Nsight's recommended application replay,
and reduces collection to LaunchStats, Occupancy, duration, local load/store
instruction counts, and local load/store bytes. All requested metrics were
confirmed present in the root-generated GA107 metric catalog, 249 GiB remained
free, retry paths were absent, the GPU was idle, and the protected source diff
was clean before requesting execution.

## 2026-07-22 — Step 1 privileged NCU collection completed

The human ran the amended root profiler. The 3D application-replay arm exited
zero and produced `baseline/ncu_retry3/3d_fapply.ncu-rep` plus its raw CSV. The
successful 2D and 3D captures describe the exact frozen `Elastic<1>::Fapply`
kernels:

| metric | 2D conservative | 3D psi |
|---|---:|---:|
| captured-kernel duration | 339.936 us | 32.575232 ms |
| block / grid | 256 / 52 | 256 / 1,158 |
| registers per thread | 96 | 254 |
| static stack bytes | 0 | 440 |
| register-limited blocks per SM | 2 | 1 |
| waves per SM | 1.44 | 64.33 |
| achieved occupancy | 28.092176% | 16.463341% |
| local-load instructions / bytes | 0 / 0 | 173,261 / 44.352384 MB |
| local-store instructions / bytes | 0 / 0 | 173,261 / 44.352384 MB |

The 3D result independently confirms that register pressure, local traffic,
and per-call cost remain material enough to satisfy the Step 4 precondition.
Report/CSV SHA-256 pairs are
`db5a542e...` / `0dc97296...` for 2D and
`ecd490ab...` / `78a1ad37...` for 3D; the complete hashes and command traces
remain in the retry directories and matching raw logs.

## 2026-07-22 — Step 1 launch trace and final pre-approval audit

Nsight Systems 2026.1.3 traced the exact frozen 2D conservative performance
case at `baseline/nsys/2d_conservative_4x4`. The target simulation completed
and finalized normally. Nsight's integrated `--stats=true` post-processing
then hung after it had written a valid `trace.nsys-rep` and `trace.sqlite`; it
was interrupted after more than three minutes, so the wrapper exit is 130 and
is preserved in `baseline/raw/nsys_2d_conservative_4x4.log`. Running the three
CUDA reports directly against the completed SQLite database exited zero.

The trace contains 76,062 GPU kernel launches totaling 585.803196 ms of
exclusive device time. Exactly 11,808 are the intended `Elastic<1>::Fapply`
kernel, totaling 416.100172 ms exclusive device time (71.0% of all traced
kernel time; 35.2388 us mean and 26.784 us median per launch). The CUDA API
summary independently contains 76,062 `cudaLaunchKernel` calls. Report hashes:

- `trace.nsys-rep`: `cbb77867bc5791c74982e878a7d27717b477ae61e9badda526c182ceada577af`
- `trace.sqlite`: `bef9cecd3fa9231b9ae7446d6e209c99f8ae8902d974181574ef1aa1a0a4d344`
- kernel summary CSV: `373059be407b8a87f8b7cdef70cc3c76eb8815efdba7bfd88fc86aa52d47bcc6`
- kernel trace CSV: `3b9d61c879585dde95bfbfb94481fabe4e1be51ae92eea53a01f6f0fe6bee282`
- CUDA API summary CSV: `55dcab46c24061a51b5a03da826649315b06d51fd7fdf62f27bddd9f249d49ce`

The Step 1 letter-and-spirit audit now passes pending the required human
baseline approval: exact binaries/objects and four manifest-compatible physics
bundles are frozen; strict/fast, 2D/3D, conservative/psi, single-/multi-box,
CPU golden, device lint, compute-sanitizer, and property oracles passed; both
timing regimes have one warmup plus five measurements per mode with stable
noise; register/stack/local-memory/occupancy data and a launch trace exist; and
the protected source paths remain clean. Final harness checks (`py_compile`,
both parser/comparator self-tests, `bash -n`, device lint, `git diff --check`)
also passed.

Known gaps are unchanged from the plan: this does not prove every nodal overlap
ordering or hierarchy invalidation order, it is A1000 rather than A100 evidence,
and five repetitions bound local noise rather than establish broad statistical
certainty. No source experiment, Step 2 sweep, NOVA run, or unrelated-worktree
mutation has begun.

## 2026-07-22 — Step 2 4/4 versus 2/2 smoothing checkpoint

The human explicitly approved the frozen Step 1 baseline and authorized the
serial Step 2 smoothing sweep. All experiments used the frozen Step 1 binaries;
the independent variable was supplied only through command-line
`elastic.solver.pre_smooth` and `elastic.solver.post_smooth` overrides. Protected
source remained clean, the local A1000 was the only compute device, and no NOVA
access or source experiment occurred.

Correctness passed before timing. The conservative 2D 2/2 arm passed the
ElasticSoftVoid property oracle in both the single-box and `max_grid_size=32`
multi-box layouts. The psi-enabled 3D 2/2 arm passed the full physics-budget gate
against the frozen fast 4/4 bundle: every correctness and engineering-trajectory
observable passed. The non-gating solver-health report records 30 MLMG V-cycles
per solve for 2/2 versus 24 for 4/4; Newton iterations remained 2.

This stability/correctness statement is limited to the frozen `max_step=2`
horizon. The 3D 2/2 final residual is `9.72181e-09`, within but close to the
`1e-08` absolute gate, versus `4.45907e-09` for 4/4; no longer-evolution
stability claim is made.

The timing protocol was one warmup plus five measured repetitions for both
unprofiled and profiler-enabled modes, with 2/2 followed by a fresh 4/4 drift arm
in every pair. Warmups are excluded below.

| regime | setting | external wall median ± MAD (s) | MLMG solve median ± MAD (s) | FApply calls | FApply wall median ± MAD (s) |
|---|---:|---:|---:|---:|---:|
| 2D conservative | drift 4/4 | 1.72 ± 0.01 | 0.9444 ± 0.0040 | 11,808 | 0.4630 ± 0.0020 |
| 2D conservative | 2/2 | 1.50 ± 0.00 | 0.7181 ± 0.0073 | 8,904 | 0.3496 ± 0.0026 |
| 3D psi | drift 4/4 | 288.58 ± 0.15 | 230.5 ± 4.1 | 15,477 | 211.9 ± 3.5 |
| 3D psi | 2/2 | 227.84 ± 2.02 | 176.5 ± 1.4 | 12,660 | 161.6 ± 1.1 |

Matched 2/2 gains are 12.791% external / 23.962% solve / 24.492% FApply /
24.593% calls in 2D, and 21.048% external / 23.427% solve / 23.738% FApply /
18.201% calls in 3D. The 2D TinyProfiler iteration counts increased from 82 to
86 `MLMG::oneIter` calls (492 to 516 `mgVcycle` region calls); the 3D counts
increased from 44 to 57 (308 to 399). Thus 2/2 wins by making each cycle cheaper
despite requiring more cycles.

For this configuration A/B, the fresh matched 4/4 drift distribution is the
baseline arm used by the performance contract. The resulting solve-gain
thresholds are 3.0% in 2D and approximately 3.56% in 3D (the larger of 3% and
twice that baseline arm's median absolute deviation as a fraction of its
median). Both solve gains clear those thresholds, and neither the matched
unprofiled external-wall primary nor the matched profiled solve primary regresses.
Together with the correctness passes, 2/2 satisfies the written performance
contract in both regimes. A 3/3 arm is not justified because the result is large
and directionally consistent rather than non-monotonic.

The 4/4 drift arms remained close to the frozen Step 1 unprofiled external wall
(2D 1.74 to 1.72 s; 3D 287.12 to 288.58 s). The 3D profiled 4/4 regions ran
faster than the Step 1 profile (MLMG solve 252.2 to 230.5 s and FApply 226.3 to
211.9 s), outside the earlier profile noise. This is disclosed as a profiler-run
drift, but it does not reverse the result: the first four profiled A/B pairs were
interleaved in the same session and produced 22.4–23.8% solve gains; the later
recovered fifth 4/4 arm yields a 26.1% pairwise gain but is not treated as
same-session evidence. The independent unprofiled external-wall gain was
18.0–22.6% per pair.

The original 3D timing wrapper was interrupted during the final profiled 4/4
`rep5`, leaving an empty wall-time file and no wrapper PASS marker. A bounded
recovery reran only that missing arm with the exact frozen binary and overrides;
it exited zero at 277.65 s external wall, 242.3 s MLMG solve, 222.1 s FApply, and
15,477 FApply calls. Its raw log ends
`STEP2_TIMING_3D_REP5_RECOVERY_PASS`. The attempted preservation name
`rep5_interrupted_20260722T1333` disappeared during the shared-workspace run even
though the move had initially succeeded; the old partial plot survives as
`rep5/plot.old.20260722134034`, while the old partial stdout and empty wall file
do not. No task process reported removing it. The recovered stdout, wall, GPU
snapshots, empty compute-process snapshots, and recovery log remain intact; their
hashes are frozen in `step2/raw/timing_pair_3d_rep5_recovery.sha256`. This
artifact-preservation anomaly does not change a measured value but is retained
here as an evidence-chain caveat.

Canonical summaries are `step2/smoothing_summary.json` and
`step2/smoothing_summary.md`; raw pass evidence is in
`step2/raw/correctness_2x2.log`, `step2/raw/timing_pair_2d.log`, and
`step2/raw/timing_pair_3d_rep5_recovery.log`. Recommendation: retain 2/2 as the
successful Step 2 experimental configuration, but do not change any canonical
input or proceed to Step 3 until the required human keep/reject checkpoint.

## 2026-07-22 — Step 3 pre-edit dispatch and test freeze

The human approved retaining 2/2 subject to stability and authorized Step 3.
Protected source was clean, the frozen 2D conservative kernel remained material,
and device lint passed before any Step 3 edit.

The four raw `Elastic::Fapply` states are object-level representable because
`SetConservativeFaceFlux` and `SetPsi` are independent. Production `Flame` wiring
uses only conservative/no-psi and general/psi. The mixed general/no-psi state is
well-defined. The mixed conservative/psi state is not: Newton's preparation path
requires both conservative enabled and a null psi pointer, whereas current
`Fapply` dispatches only on the conservative flag and would silently omit psi
weighting and its gradient correction. Step 3 therefore treats
conservative/psi as unsupported and must reject it on the host before any kernel
launch rather than dispatching by `use_psi` or silently executing inconsistent
equations.

The exact focused-test path is frozen as `src/Test/Operator/Elastic.H`, registered
through the existing unit-test harness. The focused unit enumerates all four
flag combinations against the exact production host classifier: general/no-psi
and general/psi select the general lambda, conservative/no-psi selects the
specialized lambda, and conservative/psi is classified unsupported. `Fapply`
checks that result before domain construction, `MFIter`, or a device launch.
Boundary/interior numerical behavior is covered by the full strict/fast,
single-/multi-box Oracle rather than duplicated in a synthetic operator fixture.
The required harness-registration files were added to the plan's Context budget
before inspection or implementation.

## 2026-07-22 — Step 3 specialization checkpoint and rejection

The Step 3 dispatch specialization passed device lint, CPU unit/golden checks,
compute-sanitizer, strict/fast 2D and 3D validation, single-/multi-box
ElasticSoftVoid checks, and the combined retained-2/2 physics gate. The initial
full candidate build had incorporated unrelated concurrent dirty-worktree
objects and failed sanitizer in unchanged `Elastic::Diagonal`; replacing only
the candidate `Elastic.cpp.o` in the frozen baseline object trees removed that
attribution error, and the isolated binary passed sanitizer with zero errors.

The isolated candidate did not satisfy the performance contract. In 2D the
external wall median improved from 1.49 to 1.46 s (2.013%, below the required
3%), FApply improved only 0.553%, and the profiled MLMG solve medians were
0.7008 versus 0.7577 s. In 3D the unprofiled wall median regressed from
224.93 +/- 0.64 s to 228.35 +/- 5.20 s (1.520%). The remaining 3D profiler
matrix was omitted fail-fast after the primary regression was established.
Correctness therefore passed, but Step 3 was rejected for performance. The
human authorized continued execution; the classifier, split kernels, and
focused dispatch test were removed while unrelated concurrent Newton tests were
preserved. Retained 2/2 smoothing is unchanged.

## 2026-07-22 — Step 4 sequential coefficient-gradient start

Step 4 changes only the nonuniform general `Fapply` coefficient-gradient block.
It computes and contracts `Cgrad1`, `Cgrad2`, and (in 3D) `Cgrad3` in successive
scopes, accumulating the smaller vector in x/y/z order before the unchanged
`psi_avg` scaling. Device lint and `git diff --check` pass. The isolated 3D fast
object hash is `debf279384ef0d59519389bbef57f53086faf2b1d16016b42ddae8b6afa0e3a3`;
the linked isolated binary hash is
`e9f42b5f1ff8c4f960b9515e9766e6e14a26c9ddf5b544066778d38265fdbda2`.
Static resource usage for the hot `Elastic<1>::Fapply` kernel improved from
254 registers / 440-byte stack to 250 registers / 288-byte stack, clearing the
cheap resource screen and justifying a short retained-2/2 3D runtime A/B.

The shortened runtime screen used one warmup and three alternating measured
runs per arm. Baseline times were 229.83, 229.49, and 229.73 s; candidate times
were 227.39, 230.77, and 227.48 s. The medians were therefore 229.73 versus
227.48 s, only a 0.979% candidate gain. A final matched profiler screen showed
MLMG solve 192.6 versus 191.6 s (0.519% gain) and FApply 172.2 versus 171.2 s
(0.581% gain), with 12,660 calls in both arms. This is far below the 3% keep
threshold, so the five-run matrix and full candidate oracle were omitted
fail-fast. Step 4 was rejected and its only source hunk was reverted; protected
Elastic source is again identical to HEAD.

## 2026-07-22 — Steps 5–7 disposition

Step 5 has no source combination to evaluate: Step 3 and Step 4 were rejected,
so the retained set is exactly the already-measured Step 2 configuration-only
2/2 smoothing result.

Step 6 fails its mandatory memory gate and is skipped. The 3D baseline reports
7,829 MB physical GPU memory but a 26,059 MB managed-arena high-water mark, with
only 43 MB free at finalize. Projected baseline peak plus any nodal psi cache
therefore cannot be below either 85% of physical memory or physical memory minus
1 GiB, and cannot retain 1 GiB observed headroom.

Step 7 used a fresh retained-2/2 Nsight Systems trace at
`step7/nsys/2d_conservative_2x2`. It contains 58,622 CUDA kernel launches and
8,904 `Elastic<1>::Fapply` launches totaling 306.949066 ms (70.13% of all kernel
time). Every FApply launch is immediately followed by `cudaStreamSynchronize`;
the exact FApply launch API calls total 21.698539 ms and their paired stream
synchronizations total 387.391909 ms. However, TinyProfiler also reports exactly
8,904 FApply calls. The performance case already has one FApply kernel launch per
FApply invocation, so a FabArray-wide multi-box launch cannot eliminate any of
these launches. The apparent host/device gap is per solver invocation rather
than extra per-box launch overhead. A structural child plan is therefore not
justified for this workload and Step 7 is closed without implementation.

## 2026-07-22 — NOVA/A100 confirmation staging

The user authorized NOVA use and established a reusable `ssh -MN nova` control
socket. The existing `/work/brunnels/jackplum/alamo` checkout is old and dirty,
so it was left untouched. An isolated checkout was created at
`/work/brunnels/jackplum/alamo-fapply-20260722` from a local git bundle at exact
HEAD `2e6a8f8f5430a58dd6a85b07b08d472e2bd1d5cd`. Its AMReX checkout is exact
commit `bf15fdce52093c715a1dd9a9da64ba66f0233be1`, with the two baseline tracked
AMReX modifications copied byte-for-byte; the input, `Elastic.cpp`, and both
modified AMReX file hashes match the local baseline.

Build job `11751553` requests 8 CPUs / 32 GiB and builds only the 3D sm_80
profile-fast binary offline. A100 A/B job `11751561` has
`afterok:11751553` dependency and will run one warmup plus five alternating
4/4-versus-2/2 measurements per arm, two manifest-driven physics bundles and
their gated comparison, and static FApply resource reporting. The earlier
larger build submissions `11751458` and `11751463` were canceled while still
pending before consuming resources.

Preflight review found that the workstation exact-binary runner would reject
the candidate manifest's original NOVA-only hardware tag, mislabel an A100 as
an A1000, and use `mpiexec` rather than the standard NOVA `srun` launch. The
manifest now permits exact-runner selection, and a campaign-local exact NOVA
wrapper records an `a100_sm80_fast` device tag and launches through
`srun --mpi=pmix --gpus-per-task=1`. NOVA's default Python lacked `yt`, so an
isolated `.nova-validation-venv` was created inside the campaign checkout with
`yt==4.3.0`; its import is checked before any measured run. The first dependent
A/B job `11751524` was canceled while still pending and replaced by `11751542`
so Slurm captured the corrected batch script. The build job was unchanged.

The scavenger scheduler later projected no build start until approximately
23:12 CDT. A dry-run submission identified an immediate one-hour slot in the
`interactive` partition. Still-pending jobs `11751471` and `11751542` were
canceled without resource use, then replaced by build `11751553` and dependent
A100 A/B `11751554`; the replacement build began immediately.

A subsequent dry-run showed that a normal-partition A100 submission would be
queued far beyond the campaign window, while the authorized `interactive`
partition had an immediate A100 slot. Still-pending dependent job `11751554`
was canceled without resource use and replaced by `11751561`, retaining the
same `afterok:11751553` gate and identical A/B script.

## 2026-07-22 — NOVA/A100 confirmation complete

Build job `11751553` completed in 8:24 with exit `0:0` and marker
`NOVA_A100_BUILD_PASS`. The sm_80 fast-profile binary hash is
`036c7ec520c6ae301b31d1e75a107b5dff2e343f02729c880f7be3ceeda98c84`.
Timing job `11751561` completed one warmup plus five alternating measurements
per arm on an NVIDIA A100 80GB PCIe. Medians ± MAD were:

| metric | 4/4 | 2/2 | reduction |
|---|---:|---:|---:|
| external wall | 15.58 ± 0.02 s | 13.00 ± 0.01 s | 16.560% |
| MLMG solve | 11.42 ± 0.01 s | 8.846 ± 0.002 s | 22.539% |
| FApply | 8.388 ± 0.008 s | 6.434 ± 0.001 s | 23.295% |
| FApply calls | 15,477 | 12,660 | 18.201% |

After all timing reps, the job encountered a campaign-wrapper `TypeError`
because the isolated checkout's committed `validation_common.py` predates the
local exact-binary API. The campaign-local wrapper was made independent of that
new API. Recovery job `11751579` ran only the two validation bundles, comparison,
resource report, and checksums against the preserved timing directory. It
completed in 1:53 with exit `0:0` and `NOVA_A100_VALIDATION_TAIL_PASS`. Overall
and gate verdicts are PASS for the frozen `max_step=2` case; the 2/2 final
residual remains `9.72181e-09` versus the `1e-08` absolute gate.

Nsight Compute default kernel-replay attempts `11751583` (32 GiB) and
`11751584` (64 GiB) hit Slurm host-memory cgroups while snapshotting the large
managed arena. They were retained as failed-tooling records. Job `11751588`
instead collected only LaunchStats and Occupancy with application replay and
completed in 24 seconds with `NOVA_A100_NCU_PASS`. The captured FApply launch
uses 256 threads, 254 registers/thread, one register-limited block/SM, 12.50%
theoretical occupancy, and 12.02% achieved occupancy.

The compact local evidence is under
`artifacts/nova-a100-sm80-2e6a8f8f-20260722/`. It contains all timing, logs,
commands, manifests, field norms, metrics, comparison, and the 181 MB NCU
report. The 3 GB raw validation plot payloads remain on NOVA and were omitted
locally; all 103 copied entries from the remote checksum manifest verified,
with 72 plotfile entries intentionally absent. Fresh reviewer follow-up found
no blocker to the narrow A100 confirmation or retaining 2/2.

Closeout checks passed: protected Elastic/Matrix4 source has no diff from HEAD,
the transient Step 3 test header is absent, `git diff --check` is clean, and
`benchmark/lint_device_patterns.sh` reports zero non-allowlisted violations.
`benchmark/status.sh` records branch `chamber-gpu`, the expected HEAD, 123 dirty
workspace files from the known shared worktree, and `device-lint: PASS`.
