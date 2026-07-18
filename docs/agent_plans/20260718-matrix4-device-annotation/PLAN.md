# TASK: matrix4-device-annotation

---

## Header

| Field        | Value                                                        |
|--------------|--------------------------------------------------------------|
| Risk tier    | 2 (src/Set/ headers; device-annotation only, no logic)       |
| Model        | opus (found during suite run; fix is mechanical)             |
| Verification | full-oracle (CUDA compile + GPU suite + golden gate)         |
| Est. scope   | 2 files, 2 lines (src/Set/Matrix4_Isotropic.H, Matrix4_Diagonal.H) |
| Parallel-safe| no (blocks GPU suite run for the push)                       |

## Context budget

Read: src/Set/Matrix4_Isotropic.H:44-60, src/Set/Matrix4_Diagonal.H:30-45,
      src/Operator/Elastic.cpp:440-472 (call site)
Forbidden: docs/archive/*

## Objective

chamber-gpu HEAD (dd0f86c55) does not compile for CUDA: 5ae7dffa3's
conservative-face-flux Diagonal() branch calls
Matrix4<DIM,Isotropic>::operator() and Matrix4<DIM,Diagonal>::operator()
inside an AMREX_GPU_DEVICE lambda, but both element-access operators lack
AMREX_GPU_HOST_DEVICE (nvcc: "calling a __host__ function from a __device__
function", Elastic.cpp:459). Discovered 2026-07-18 building the 2D CUDA
binary for the GPU suite. NOTE: no CUDA binary anywhere has compiled the 14
unpushed commits (newest CUDA build = 3D, Jul 8, predates the Jul 17
recovery merge; the a100-sanitizer PASS gate evidence is stale).

## Oracle

Command(s): DIM=2 CUDA_FP=fast bash benchmark/build_alamo_local_gpu.sh
            (exit 0); then python3 tests/GPU/run_gpu_tests.py; then
            bash benchmark/ci_golden_compare.sh (CPU golden unchanged).
Covers: device compile, GPU runtime correctness (C*/F* tests), CPU parity.
Does NOT cover: A100 validation (local A1000 only; P* perf informational).

## Steps

### Step 1 - annotate
DO: add AMREX_GPU_HOST_DEVICE to
    - Matrix4_Isotropic.H operator () (line ~51, retrieval-only, pure
      arithmetic: mu/lambda sums; no abort path)
    - Matrix4_Diagonal.H operator () (line ~38, pure array read)
CHECK: 2D CUDA build exits 0. Iterate: build stopped at first error --
    further hidden device-compile breaks may surface; annotate same-pattern
    failures, STOP and reassess if any failure requires logic change.

### Step 2 - suite + gates
DO: rebuild 3D CUDA binary; run tests/GPU/run_gpu_tests.py; re-run
    ci_golden_compare.sh (host codegen unaffected -- macro empty on CPU
    builds -- but verify, don't assume).
CHECK: suite C*/F* pass (P* informational on 50W A1000); golden gate green.

## Checkpoints (tier 2)

- [ ] Before commit: diff summary (must be annotation-only) + oracle output.

## Closeout

- [ ] results/RESULT.md, touch results/DONE, SESSION_LOG.tsv line
- [ ] changelog note: CUDA-compile regression in unpushed set found+fixed
      before push (gate-evidence staleness noted)
