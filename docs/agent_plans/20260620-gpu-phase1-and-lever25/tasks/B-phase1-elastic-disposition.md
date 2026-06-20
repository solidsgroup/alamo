# Task B: Phase 1 — elastic divergence + path disposition (D1)

## Goal
Resolve the high-res MLMG divergence on the GPU elastic solve and decide, with
data, whether elastic lives on the GPU, the CPU, or a hybrid. Produce report R1
and the D1 verdict {GPU | CPU | HYBRID}. This is a **decision/measurement** task —
**no source edits, no builds.**

## Context (only what this task needs)
- Repo `/home/jackplum/Projects/alamo`, branch `chamber-gpu` @ `fe3844f31`, main tree.
- Read `../PLAN.md` for ownership rules. You OWN: a **private copy** of the nofast
  binary, `benchmark/phase1_elastic_2048/`, `analysis/results_phase1*/`, and input
  copies prefixed `input_phase1_*`. You must NOT `./configure`, `make`, or build;
  must NOT edit `src/`; and must **NEVER modify the shared `input` file** (Track A's
  gate reads it — copy it instead).
- **First action:** `cp bin/alamo_gpu-2d-nofast-cuda86-g++ bin/alamo_gpu_phase1-nofast`
  and use that private copy for all GPU elastic runs (Track A rebuilds `bin/`).
- GPU: single RTX A1000, **8 GB**, sm_86 — shared with Track A; be robust to OOM.
- `source benchmark/local_cuda_env.sh` before GPU/nsys/ncu runs.
- **Prior work to build on (read first):** `benchmark/PHASE1_ELASTIC_DISPOSITION.md`
  and `benchmark/phase1_elastic_2048/` already contain `cpu_np8`,
  `gpu_strict_512`, and `gpu_strict_1024` logs+plots. Do NOT redo what's done —
  read what was already concluded and extend it. Note the dir is named "2048" but
  prior runs only reached 1024².
- Background (memory `gpu_perf_nsys_findings.md`, `chamber_gpu_port.md`): GPU
  elastic SIGABRTs ("MLMG failing") on the first elastic solve at high res under
  `--use_fast_math`; CPU solves the identical config fine; at coarse 64²-base GPU
  elastic was ~1.65× faster than CPU. Benchmarks normally SKIP elastic via
  `elastic.tstart=1e9`. To EXERCISE elastic set `elastic.on=1`, `elastic.type=static`,
  a sane `elastic.tstart`, and the stability settings from memory (`psi_floor`,
  clamped BCs). Make these edits in an `input_phase1_*` copy only.

## Files to read first
- `benchmark/PHASE1_ELASTIC_DISPOSITION.md`, `benchmark/phase1_elastic_2048/*.log`.
- `../PLAN.md`. The shared `input` (read-only — copy, never edit).
- `benchmark/g0_ncu_capture.sh` (ncu harness for step 1.3).

## Files allowed to modify / create
- `input_phase1_*` (copies of `input`), `benchmark/phase1_elastic_2048/**`,
  `analysis/results_phase1*/**`, the report file `../results/B-RESULT.md` and an
  R1 report doc (e.g. `benchmark/PHASE1_ELASTIC_DISPOSITION.md` — extend it, or a
  new dated file).
- The private binary copy `bin/alamo_gpu_phase1-nofast`.

## Files NOT allowed to modify
- `src/**`, `configure`, `Makefile`, `.alamo`, the shared `input`, `bin/` (other
  than your private copy), and anything Track A owns.

## Implementation steps (roadmap Phase 1)
1. **1.1 — no-fast-math GPU elastic at high res.** Using the private nofast binary,
   run the elastic case that SIGABRTs under fast-math. Target 2048²; **if it OOMs
   on 8 GB, document the memory ceiling and run the largest fitting resolution
   (e.g. 1536²/1024²) — the convergence question is answerable below 2048² too.**
   Record MLMG convergence (does it converge, or still SIGABRT/diverge?) and the
   iteration count, vs the CPU (`bin/alamo-2d-clang++`, np8) on the identical system.
2. **1.2 — branch on the result (D1):**
   - Converges ⇒ fast-math × conditioning hypothesis confirmed; fix = the build-flag
     split (fast vs strict binaries already exist). 
   - Still diverges ⇒ it is NOT the flag; open an operator/conditioning
     investigation (psi formulation, `psi_floor`, masked-operator conditioning at
     thin features). Freeze device-elastic perf work; elastic stays CPU-resident
     until convergence is restored.
3. **1.3 — occupancy/register profile** of `Fapply` / `Diagonal` / `Newton` with
   `maxrregcount` removed, via ncu (`NCU=.local/nsight-compute/usr/bin/ncu`,
   `benchmark/g0_ncu_capture.sh`). Capture achieved occupancy + registers/thread;
   quantify the Eigen + `Matrix4` per-thread state cost. If high-res elastic won't
   run on GPU, profile at the largest resolution that does.
4. **1.4 — fair benchmark:** GPU elastic vs **fully-subscribed CPU** MLMG
   (`bin/alamo-2d-clang++`, np8) at ≥1024² finest, counting only the solve regions.
   State explicitly whether the coarse "1.65× faster" result survives a multi-core
   baseline.
5. **D1 verdict** per the decision tree: occupancy <25% ⇒ register-bound/weak fit,
   bias CPU/hybrid; GPU slower than CPU node ⇒ CPU-resident; within ~1.5× ⇒ hybrid;
   GPU faster ⇒ GPU-resident. Output `elastic disposition = {GPU | CPU | HYBRID}`
   with the numbers behind it.

## Invariants
- No builds, no source edits, never modify the shared `input`. Use the private
  binary copy and `input_phase1_*` copies only.
- Be robust to shared-GPU OOM: retry, lower resolution, and record what actually
  ran vs what was attempted.

## Build and test commands
None — you do not build. Run binaries directly, e.g.:
```bash
source benchmark/local_cuda_env.sh
mpiexec -np 1 bin/alamo_gpu_phase1-nofast input_phase1_elastic <overrides>
mpiexec -np 8 bin/alamo-2d-clang++        input_phase1_elastic <overrides>
```

## Expected result
A defensible D1 verdict with: the high-res convergence result (and the actual
resolution achieved + memory ceiling), occupancy/register table for the three
elastic kernels, the fair GPU-vs-CPU-node benchmark, and an honest statement of
whether the coarse GPU win survives. Report R1 written. **Do NOT git commit.**

## Non-goals
- Editing solver source / "fixing" the operator in code (1.2b is a diagnosis here,
  not a code change — flag it for the lead).
- Phase 2 / 2.5 / Phase 3 work. 3D. Multi-GPU.

## Stop conditions
- Even the lowest sensible resolution OOMs or SIGABRTs unrecoverably on GPU →
  record it as the result (that itself informs D1) and continue with CPU/ncu data.
- ncu unavailable/permission-blocked → record and proceed with what you have.

## Final report: write `../results/B-RESULT.md`
```
# Result: Task B — Phase 1 elastic disposition
## Summary (D1 verdict: GPU | CPU | HYBRID, with confidence)
## 1.1 convergence: resolution achieved, memory ceiling, converge vs diverge, iters vs CPU
## 1.2 branch taken and rationale
## 1.3 occupancy/registers for Fapply / Diagonal / Newton
## 1.4 fair GPU-vs-CPU-node benchmark (solve regions only)
## Does the coarse 1.65x win survive the multi-core baseline?
## Files changed/created
## Issues found
## Deviations from task
## Follow-up needed (e.g. operator root-cause if diverging)
```
