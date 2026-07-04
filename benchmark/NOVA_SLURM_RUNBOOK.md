# NOVA Slurm Runbook

This document records how to discover NOVA's Slurm shape, build Alamo on NOVA,
submit CPU/GPU jobs, and collect enough evidence to generate new `.slurm`
scripts without relying on stale cluster assumptions.

Status of the 2026-06-30 capture: `nova_doc_capture_20260630_143346.tar.gz`
was only 45 bytes and contained an empty tar stream. No NOVA inventory values
from that capture were available to analyze. The commands below are therefore
the canonical capture/run procedure; fill the inventory table from a non-empty
capture before changing resource defaults.

## Known Assumptions

- Account: `brunnels`
- Main partition: `nova`
- Optional partition seen in older examples: `scavenger` for preemptible GPU
  diagnostics
- Build scripts load `cuda`, `gcc`, and `openmpi`, with fallbacks for
  `cuda/12`, `gcc/12`, and `openmpi4`
- 3D scripts currently support `GPU_TYPE=v100|a100|h200`, mapped to CUDA
  architectures `70|80|90`
- NOVA Phase 3 runs default `ASYNC_OUT=0` because earlier NOVA `srun --mpi=pmix`
  launches reported `MPI_THREAD_SINGLE`; AMReX async output requires
  `MPI_THREAD_MULTIPLE`

Treat these as script defaults, not verified cluster inventory. Confirm them
with the capture commands before adding or regenerating scripts.

## Inventory To Capture Before Generating Scripts

Collect these facts any time NOVA changes, a new GPU class appears, or a script
needs new resource defaults:

- Slurm accounts, partitions, QOS values, wall limits, job limits, and preemption
  rules
- Accepted GPU GRES names, for example `gpu:v100`, `gpu:a100`, or `gpu:h200`
- GPU count per node, GPU memory, compute capability, driver version, and CUDA
  toolkit version
- CPU cores and memory per CPU node and GPU node
- Correct MPI launch mode from `srun --mpi=list`; current scripts use `pmix`
- Runtime MPI thread level reported by AMReX, especially whether async output is
  allowed
- Exact module names for CUDA, GCC, OpenMPI, Nsight Systems, Nsight Compute, and
  `compute-sanitizer`
- Whether `ncu` performance counters are permitted or fail with
  `ERR_NVGPUCTRPERM`
- Login-node policy: what may run on login nodes, and whether compile jobs must
  be submitted through Slurm
- Work/scratch paths and quotas for builds, plot files, and profiler traces

## Read-Only Slurm Capture

Run this on a NOVA login node from the directory where the capture should be
created. These commands do not submit jobs.

```bash
mkdir -p nova_doc_capture
cd nova_doc_capture

date > 00_date.txt
hostname -f > 01_hostname.txt
whoami > 02_user.txt
pwd > 03_pwd.txt

sacctmgr show user "$USER" withassoc \
  format=User,Account,Partition,QOS,DefaultAccount%25,MaxJobs,MaxSubmit,MaxWall,MaxTRESPU%60 -P \
  > 10_sacctmgr_user_assoc.txt 2>&1

sacctmgr show assoc where user="$USER" \
  format=Cluster,Account,User,Partition,QOS,GrpTRES%80,MaxTRES%80,MaxWall,Priority -P \
  > 11_sacctmgr_assoc.txt 2>&1

sinfo -Nel > 20_sinfo_nodes_long.txt 2>&1
sinfo -o "%P|%a|%l|%D|%t|%N|%G|%c|%m|%f" > 21_sinfo_partitions_nodes.txt 2>&1
scontrol show partition -o > 22_scontrol_partitions.txt 2>&1
scontrol show config > 23_scontrol_config.txt 2>&1

sinfo -h -N -o "%N" | sort -u | while read -r node; do
  scontrol show node -o "$node"
done > 24_scontrol_nodes_one_line.txt 2>&1

squeue -a -o "%.18i %.9P %.30j %.12u %.2t %.10M %.10l %.6D %.20R %.60b" \
  > 30_squeue_all.txt 2>&1
squeue -u "$USER" -o "%.18i %.9P %.30j %.2t %.10M %.10l %.6D %.20R %.60b" \
  > 31_squeue_user.txt 2>&1

module --version > 40_module_version.txt 2>&1
module avail > 41_module_avail.txt 2>&1
module spider cuda > 42_module_spider_cuda.txt 2>&1
module spider gcc > 43_module_spider_gcc.txt 2>&1
module spider openmpi > 44_module_spider_openmpi.txt 2>&1
module spider cmake > 45_module_spider_cmake.txt 2>&1
module spider nsight > 46_module_spider_nsight.txt 2>&1
module spider nvhpc > 47_module_spider_nvhpc.txt 2>&1
```

From the repo checkout, also capture the resource selector's interpretation:

```bash
cd /path/to/alamo
bash benchmark/select_nova_resources.sh --explain \
  > /path/to/nova_doc_capture/50_select_nova_resources.txt 2>&1
```

Package from the parent directory, not from inside `nova_doc_capture`, so the
archive contains the files:

```bash
tar -czf nova_doc_capture_$(date +%Y%m%d_%H%M%S).tar.gz nova_doc_capture
tar -tzf nova_doc_capture_*.tar.gz | head
```

If `tar -tzf` prints nothing, the archive is empty and must be recreated.

## Avoid Long Paste Failures

Long one-line `sbatch --wrap=...` commands are fragile in web terminals because
line wrapping can insert real newlines. Prefer batch files for probes and
benchmarks. Create the file with an editor, then submit the file with `sbatch`.

## CPU Probe Job

This verifies basic Slurm launch behavior, MPI modes, and exported environment
variables.

```bash
cat > cpu_probe.sbatch <<'EOF'
#!/usr/bin/env bash
#SBATCH --account=brunnels
#SBATCH --partition=nova
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --job-name=nova_cpu_probe
#SBATCH --output=60_cpu_probe.%j.out
#SBATCH --error=60_cpu_probe.%j.err

hostname
module purge
module load gcc 2>/dev/null || true
module load openmpi 2>/dev/null || true
module list
which mpirun srun || true
srun --mpi=list
srun -n 2 bash -lc 'hostname; env | sort | grep -E "SLURM|PMI|PMIX|OMPI"'
EOF

sbatch --parsable cpu_probe.sbatch > 60_cpu_probe_jobid.txt
squeue -j "$(cat 60_cpu_probe_jobid.txt)"
```

## GPU Probe Job

Run once per GPU GRES name reported by `sinfo`/`scontrol`. Change both
`--gres=gpu:a100:1` and the output filenames when probing another type.

```bash
cat > gpu_probe_a100.sbatch <<'EOF'
#!/usr/bin/env bash
#SBATCH --account=brunnels
#SBATCH --partition=nova
#SBATCH --gres=gpu:a100:1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=00:05:00
#SBATCH --job-name=nova_gpu_probe_a100
#SBATCH --output=61_gpu_probe_a100.%j.out
#SBATCH --error=61_gpu_probe_a100.%j.err

hostname
module purge
module load cuda 2>/dev/null || module load cuda/12 2>/dev/null || true
module load gcc 2>/dev/null || true
module load openmpi 2>/dev/null || true
module list
which nvcc nvidia-smi nsys ncu compute-sanitizer || true
nvcc --version || true
nvidia-smi -L
nvidia-smi --query-gpu=name,compute_cap,memory.total,driver_version --format=csv
srun --mpi=list
srun --gpus-per-task=1 -n 1 bash -lc 'hostname; nvidia-smi -L; env | sort | grep -E "SLURM|CUDA|GPU|PMI|PMIX|OMPI"'
EOF

sbatch --parsable gpu_probe_a100.sbatch > 61_gpu_probe_a100_jobid.txt
squeue -j "$(cat 61_gpu_probe_a100_jobid.txt)"
```

If Slurm rejects the GRES, identify accepted names with:

```bash
grep -i 'gpu' 21_sinfo_partitions_nodes.txt 24_scontrol_nodes_one_line.txt | head -120
```

## Build On NOVA

Use a work path with enough quota for AMReX, binaries, plot files, and profiler
traces.

```bash
mkdir -p /work/brunnels/$USER/alamo
cd /work/brunnels/$USER/alamo

# 2D build: produces bin/alamo_gpu-2d-profile-cuda<arch>-g++
sh /path/to/alamo/benchmark/build_alamo_nova.sh

# 3D build: produces bin/alamo_gpu-3d-profile-cuda<arch>-g++
sh /path/to/alamo/benchmark/build_alamo_nova_3d.sh
```

Useful overrides:

```bash
ACCOUNT=brunnels BUILD_PARTITION=nova ARCHES="70 80 90" BUILD_JOBS=64 BUILD_MEM=64G \
  sh benchmark/build_alamo_nova_3d.sh
```

Watch the build:

```bash
squeue -u "$USER"
tail -f alamo_build_3d.<jobid>.out
```

## Run 2D Jobs

The 2D GPU script expects the enclosing `sbatch` allocation to request GPU GRES.

```bash
sbatch --gres=gpu:a100:1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=8 --mem=32G \
  benchmark/nova_flame_gpu.slurm

GPU_TYPE=h200 sbatch --gres=gpu:h200:1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=8 --mem=32G \
  benchmark/nova_flame_gpu.slurm

MODE=fast sbatch --gres=gpu:a100:1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=8 --mem=32G \
  benchmark/nova_flame_gpu.slurm

sbatch benchmark/nova_flame_cpu.slurm
```

## Run 3D Jobs

Single GPU:

```bash
GPU_TYPE=v100 INPUT=input_3d_centre_bore_128 \
  sbatch --gres=gpu:v100:1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=8 --mem=32G \
  benchmark/nova_flame_gpu_3d.slurm
```

Multi GPU, one MPI rank per GPU:

```bash
GPU_TYPE=v100 NGPUS=2 INPUT=input_3d_centre_bore_128 \
  sbatch --gres=gpu:v100:2 --ntasks=2 --ntasks-per-node=2 --cpus-per-task=8 --mem=64G \
  benchmark/nova_flame_gpu_3d_multi.slurm
```

CPU baseline:

```bash
INPUT=input_3d_centre_bore_128 sbatch benchmark/nova_flame_cpu_3d.slurm
```

Keep `ASYNC_OUT=0` unless a probe or AMReX startup log confirms
`MPI_THREAD_MULTIPLE`.

## Queue And Failure Diagnostics

```bash
squeue -u "$USER" -o "%.18i %.9P %.30j %.2t %.10M %.10l %.6D %.20R %.60b"
scontrol show job -dd <jobid>
sacct -j <jobid> --format=JobID,JobName%30,Partition,Account,AllocTRES%80,State,ExitCode,Elapsed,MaxRSS,ReqMem
bash benchmark/slurm_pending_reason.sh <jobid>
```

Search logs for:

```bash
grep -E 'amrex::Abort|CUDA|Out of memory|Invalid generic resource|CANCELLED|TIME LIMIT|SIGABRT|SIGKILL|MPI_THREAD' \
  *.out *.err
```

Success markers:

```text
STEP <n> ends.
TinyProfiler total time
=== done ===
```

## Inventory Table

Fill this table from `nova_doc_capture` and probe logs. Do not regenerate Slurm
defaults until this is populated.

| Item | Observed value | Evidence file |
|------|----------------|---------------|
| Default account | pending | `10_sacctmgr_user_assoc.txt` |
| Valid partitions | pending | `21_sinfo_partitions_nodes.txt` |
| GPU GRES names | pending | `21_sinfo_partitions_nodes.txt`, `24_scontrol_nodes_one_line.txt` |
| GPU node CPU count | pending | `24_scontrol_nodes_one_line.txt` |
| GPU node memory | pending | `24_scontrol_nodes_one_line.txt` |
| CPU node CPU count | pending | `24_scontrol_nodes_one_line.txt` |
| CPU node memory | pending | `24_scontrol_nodes_one_line.txt` |
| CUDA module | pending | `42_module_spider_cuda.txt`, GPU probe output |
| GCC module | pending | `43_module_spider_gcc.txt`, probe output |
| OpenMPI module | pending | `44_module_spider_openmpi.txt`, probe output |
| `srun --mpi` mode | pending | CPU/GPU probe output |
| MPI thread level | pending | Alamo job output |
| Nsight Systems module/path | pending | `46_module_spider_nsight.txt`, GPU probe output |
| Nsight Compute counter access | pending | `ncu` diagnostic job output |
