# Phase C1 — NOVA A100 before/after for the `Fapply` source edits (TO RUN)

This is the **only remaining C1 gate item**. The CPU golden compare is done and
bit-identical (`PHASE_C1_cpu_golden_compare.md`); the source edits are described in
`PHASE_C1_fapply_occupancy.md`. What's owed is the roadmap rule "never land a kernel
change without an A100 before/after": measure `Fapply` **achieved occupancy +
registers/thread** and wall/step, modified vs baseline, on one A100.

All commands run from the worktree root on a NOVA login/build node unless noted. The
A/B is a **binary swap** (modified vs baseline binary, *identical* input) — not an
input-lever sweep like `phase_c_elastic_ab.sh`.

## 0. Sync this worktree to NOVA

The edits live only here (`chamber-gpu-elastic-opt`, uncommitted). Either commit to a
scratch ref and pull on NOVA, or rsync the tree. The single source delta is
`src/Operator/Elastic.cpp`; everything else matches `chamber-gpu` tip `7e972f1e8`.

## 1. Build BOTH arms (sm_80 A100), holding everything but `Elastic.cpp` fixed

```bash
# MODIFIED arm (edits present)
ARCHES=80 benchmark/build_alamo_nova_3d.sh          # -> bin/alamo_gpu-3d-profile-cuda80-g++
cp bin/alamo_gpu-3d-profile-cuda80-g++  bin/alamo_gpu-3d-cuda80-g++.MODIFIED

# BASELINE arm (revert only Elastic.cpp, rebuild)
git stash push -- src/Operator/Elastic.cpp         # or: git checkout -- src/Operator/Elastic.cpp
ARCHES=80 benchmark/build_alamo_nova_3d.sh
cp bin/alamo_gpu-3d-profile-cuda80-g++  bin/alamo_gpu-3d-cuda80-g++.BASELINE
git stash pop                                       # restore the edits
```

(If `build_alamo_nova_3d.sh` is wrapped in its own sbatch, build each arm in its own
job and copy the artifact out before reverting.)

## 2. Fapply occupancy / registers — the headline metric (ncu)

`benchmark/g0_ncu_capture.sh` already does the portable, version-proof capture
(auto-picks the ncu metric set, profiles a bounded launch window, prints a per-kernel
SoL table incl. `launch__registers_per_thread` and `sm__warps_active` ≈ achieved
occupancy). Point it at each saved binary and at the elastic-heavy input:

```bash
for arm in MODIFIED BASELINE; do
  GPU_BIN=bin/alamo_gpu-3d-cuda80-g++.$arm \
  INPUT=input_3d_centre_bore_256_a2 \
  OUT=benchmark/profiles_c1/${arm,,}_ncu \
  NCU_COUNT=80 \
  benchmark/g0_ncu_capture.sh 2>&1 | tee benchmark/profiles_c1/${arm,,}_ncu.log
done
```

Then read off, for the `Operator::Elastic<…>::Fapply` row in each run:
- `launch__registers_per_thread`  (baseline ≈ **255**, the cap — does it drop?)
- `sm__warps_active.avg.pct_of_peak_sustained_active`  (baseline ≈ **12.5%** occ)
- `gpu__time_duration.sum`  (per-launch Fapply time)

**Pass condition:** registers/thread falls below 255 (i.e. spills reduced) and/or
achieved occupancy rises from ~12.5%. Expected from `PHASE_A_FINDINGS.md` §8:
1.5–3× on Fapply if occupancy moves 12.5% → 25–50%. (Occupancy may *not* move if 255
is still the binding limit even after the live-double cut — in which case the honest
result is "register pressure reduced but occupancy unchanged," and the deferred
`__launch_bounds__` lever A1a becomes the next step. Record whatever ncu actually says.)

## 3. Wall/step A/B (TinyProfiler) — secondary

Short bench run of each arm on the same input; compare `MLMG::solve` inclusive time
and the `Fapply` TinyProfiler region. Reuse the A2 slurm by making each saved binary
the newest `cuda80` match it globs (it picks `ls -t bin/alamo_gpu-3d*cuda80*-g++`):

```bash
for arm in MODIFIED BASELINE; do
  cp bin/alamo_gpu-3d-cuda80-g++.$arm bin/alamo_gpu-3d-cuda80-g++   # make it newest
  INPUT=input_3d_centre_bore_256_a2 GPU_TYPE=a100 MODE=bench \
  EXTRA_ARGS="stop_time=0.06 elastic.solver.verbose=4 amr.plot_int=30 plot_file=out_c1_${arm,,}/plot" \
    sbatch --partition=nova --nodes=1 --gres=gpu:a100:1 --ntasks=1 \
           --cpus-per-task=72 --mem=0 --time=01:00:00 \
           -J c1_ab_${arm,,} benchmark/nova_flame_gpu_3d_a2.slurm
done
# after both finish:
grep -E 'MLMG::solve|Fapply|ELASTIC SOLVE' c1_ab_*.out
```

## 4. GPU stress-parity gate (reconfirm bit-identity on A100)

The CPU compare already proved bit-identity; this re-checks it on-device (HMM-free
A100, the real target). With the step-30 plotfiles written above:

```bash
python3 benchmark/compare_thermo.py out_c1_modified/ out_c1_baseline/   # stress parity
# or a hard byte compare of the node field, as the CPU compare did:
cmp out_c1_modified/plot/00030node/Level_0/Cell_D_00000 \
    out_c1_baseline/plot/00030node/Level_0/Cell_D_00000
```

GPU fast-math may make this not *bit*-exact run-to-run (atomics/reduction order); a
no-fast-math (`--cuda-fp strict`) build is the strict check if needed. Either way the
stress field must match the baseline to tolerance, else the edit is rejected.

## 5. Record the result

Write the ncu before/after table (regs + occupancy + Fapply µs) and the wall A/B into
`PHASE_C1_fapply_occupancy.md` and `PHASE_A_FINDINGS.md` §8, and update memory
`gpu_c1_fapply_source_opt.md`. Only after step 2 passes is the C1 source lever
"measured" per the roadmap; until then it remains correctness-verified but
perf-unproven.
