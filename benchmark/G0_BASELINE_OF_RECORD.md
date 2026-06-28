# G0 Baseline of Record

Date: 2026-06-20

Roadmap: `/home/jackplum/Desktop/GPU-OPT-ROADMAP.txt`

## Correctness/build status

- CUDA fast build present: `bin/alamo_gpu-2d-cuda86-g++`
- CUDA strict/no-fast build present: `bin/alamo_gpu-2d-nofast-cuda86-g++`
- CPU baseline build present: `bin/alamo-2d-clang++`
- Baseline references recorded under `benchmark/baseline_references/`
- Baseline check passed for all recorded references.

Commands:

```bash
CPU_NP=8 GPU_FAST_NP=1 GPU_STRICT_NP=1 python3 benchmark/baseline_suite.py record
CPU_NP=8 GPU_FAST_NP=1 GPU_STRICT_NP=1 python3 benchmark/baseline_suite.py check
CPU_NP=8 GPU_FAST_NP=1 GPU_STRICT_NP=1 python3 benchmark/baseline_suite.py report
```

## CPU/GPU timing baseline

CPU is the requested `-np 8` best-case local baseline. GPU profiles are single-rank.

| Case | Profile | MPI ranks | Elapsed s | vs CPU | Correctness |
| --- | --- | ---: | ---: | ---: | --- |
| canonical_step1 | cpu | 8 | 1.732 | 1.000x | reference |
| canonical_step1 | gpu_fast | 1 | 1.529 | 0.883x | ok |
| canonical_step1 | gpu_strict | 1 | 1.431 | 0.827x | ok |
| canonical_step2 | cpu | 8 | 1.954 | 1.000x | reference |
| canonical_step2 | gpu_fast | 1 | 1.865 | 0.954x | ok |
| canonical_step2 | gpu_strict | 1 | 1.849 | 0.946x | ok |
| eta_expression_step1 | cpu | 8 | 1.082 | 1.000x | reference |
| eta_expression_step1 | gpu_fast | 1 | 1.222 | 1.130x | ok |
| eta_expression_step1 | gpu_strict | 1 | 1.265 | 1.170x | ok |

## GPU-safe IC/BC matrix

See `docs/gpu_safe_ic_bc_matrix.md`.

Current canonical G0 surface:

- BMP IC: GPU-safe after the BMP max fix.
- Constant elastic BC: GPU-safe.
- Expression/PNG/PSRead/Trig/Laminate/vector-expression host-loop paths: guarded to fail fast under CUDA instead of silently writing device-arena data from host loops.

## Kernel resource snapshot

Static resource usage was captured with:

```bash
source benchmark/local_cuda_env.sh
cuobjdump --dump-resource-usage bin/alamo_gpu-2d-cuda86-g++ > benchmark/g0_cuda_resource_usage.txt
```

Representative static register counts from `benchmark/g0_cuda_resource_usage.txt`:

| Kernel family | Static registers/thread | Stack | Theoretical occupancy note, sm_86, 256 threads/block |
| --- | ---: | ---: | --- |
| `Flame::Advance` phase/thermal kernels | 55-56 | 0 B | about 66.7 percent register-limited occupancy |
| `Operator::Elastic::Fapply` | 87-101 | 0 B | about 33.3 percent register-limited occupancy |
| `Operator::Elastic::Diagonal` | 48-56 | 240 B | about 66.7 percent register-limited occupancy, plus stack spill pressure |
| `Newton::prepareForSolve` | 40-80 | 0-32 B | about 50-100 percent depending on kernel variant |
| AMReX `MLTensorOp::prepareForSolve` | 30 | 0 B | not register-limited |

These are static/theoretical numbers, not achieved occupancy.

## Nsight Compute status

Nsight Compute was unpacked locally without system installation:

```bash
apt-get download nsight-compute nsight-compute-target
dpkg-deb -x nsight-compute_2022.4.1.6~12.0.1-4build4_amd64.deb .local/nsight-compute
dpkg-deb -x nsight-compute-target_2022.4.1.6~12.0.1-4build4_amd64.deb .local/nsight-compute
```

The CLI runs:

```bash
.local/nsight-compute/usr/bin/ncu --version
```

But actual counter capture is blocked on this workstation by driver permissions:

```text
ERR_NVGPUCTRPERM - The user does not have permission to access NVIDIA GPU Performance Counters
```

The prepared capture command is:

```bash
benchmark/g0_ncu_capture.sh
```

Run it after enabling NVIDIA performance counter access or on a cluster node where Nsight Compute counters are permitted.

## G0 verdict

Local status: correctness/build/baseline timing requirements are satisfied; static register data is captured; achieved occupancy/register profiling is prepared but blocked by local `ERR_NVGPUCTRPERM`.

Strict roadmap status: G0 is not fully passable on this workstation until Nsight Compute performance counters are enabled or the same capture is run on NOVA/another permitted host.
