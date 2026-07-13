# NOTES — fapply-322b-a100 (submit phase, 2026-07-13)

Everything below covers Step 1-3 of PLAN.md: sync arms to NOVA, submit builds,
submit dependent wall A/B + ncu jobs. **Nothing has been waited on** — this is
the queued state as of submission time.

## Arms

| Arm      | Commit      | Description                                                |
|----------|-------------|--------------------------------------------------------------|
| BASELINE | `d964cfab8` | chamber-gpu tip before 3.2b (`Add NOVA NAWC 4ths elastic comparison runs`) |
| MODIFIED | `dc02baf17` | `fapply-322b` tip (DDW hoist, column-restricted contraction, unrolled Matrix4xMatrix3, Fsmooth 4D launch fusion) |

`git merge-base d964cfab8 dc02baf17` = `d964cfab8` (baseline is an ancestor of
modified, confirmed locally before bundling).

## NOVA paths

- Bundle: `/work/brunnels/jackplum/alamo322b.bundle` (contains temp refs
  `_bundle_baseline` → `d964cfab8`, `_bundle_modified` → `dc02baf17`; created
  locally via `git bundle create` after `git update-ref` on those two temp
  branches, refs deleted locally post-bundle).
- Arm checkouts: `/work/brunnels/jackplum/alamo-322b-ab/baseline/`,
  `/work/brunnels/jackplum/alamo-322b-ab/modified/` — each `git clone`d from
  the bundle then checked out to its respective `_bundle_*` ref.
- Verified on NOVA:
  - `git -C baseline rev-parse HEAD` = `d964cfab8c5ca063182dbef042a7bc1c129019be` ✓
  - `git -C modified rev-parse HEAD` = `dc02baf1779b6668a6ea5725ade2e79b22118dfb` ✓
- sbatch scripts (written directly on NOVA, templated from
  `benchmark/build_alamo_nova_3d.sh`'s embedded sbatch and the
  `20260707-a100-ab-c1` precedent `c1ab_wall.sbatch`/`c1ab_ncu2.sbatch` — none
  of these were added under the local repo's `benchmark/` per plan Step 1's
  "never in benchmark/" instruction, they live only on NOVA):
  - `<arm>/build_322b.sbatch` — ARCHES=80 only single-arch GPU build (`bin/alamo_gpu`).
  - `<arm>/wall_322b.sbatch` — TinyProfiler wall run, `input_3d_centre_bore_256_a2`,
    `stop_time=0.06 elastic.solver.verbose=4 amr.plot_int=30`, matches
    a100-ab-c1 EXTRA_ARGS convention for a bounded parity-checkable run.
  - `<arm>/ncu_322b.sbatch` — ncu gated on NVTX range
    `Operator::Elastic::Fapply()/` (per a100-ab-c1 precedent job 11448954 fix —
    do NOT rely on kernel-name matching), `--launch-count 40`, `--set basic`.

## Job IDs + dependency graph

| Step        | Arm      | Job ID   | Depends on |
|-------------|----------|---------:|------------|
| build       | baseline | 11639419 | —          |
| build       | modified | 11639420 | —          |
| wall A/B    | baseline | 11639477 | afterok:11639419 |
| wall A/B    | modified | 11639478 | afterok:11639420 |
| ncu capture | baseline | 11639479 | afterok:11639419 |
| ncu capture | modified | 11639480 | afterok:11639420 |

Confirmed via `squeue -u jackplum` immediately after submission: build jobs
`R` (running), all four dependent jobs `PD` with `Dependency` reason and the
correct unfulfilled `afterok:<jobid>` — dependency graph wired as intended.

## Output file locations (once jobs complete)

- Build logs: `<arm>/build_322b.<jobid>.out` / `.err`; binary
  `<arm>/bin/alamo_gpu-3d-profile-cuda80-g++`.
- Wall logs: `<arm>/wall_322b.<jobid>.out` / `.err` and
  `<arm>/wall_322b_<arm>.log` (TinyProfiler table); plotfiles under
  `<arm>/out_322b_<arm>/plot*` at steps 0/30/... (`amr.plot_int=30`) for
  fcompare parity.
- ncu: `<arm>/profiles_322b/<arm>_nvtx.ncu-rep`,
  `<arm>/profiles_322b/<arm>_nvtx.log`, plus the per-kernel CSV summary
  printed at the end of `<arm>/ncu_322b.<jobid>.out`.

## Exact commands used

```bash
# Local: bundle both arm tips
git update-ref refs/heads/_bundle_baseline d964cfab8
git update-ref refs/heads/_bundle_modified dc02baf17
git bundle create alamo322b.bundle refs/heads/_bundle_baseline refs/heads/_bundle_modified
git update-ref -d refs/heads/_bundle_baseline
git update-ref -d refs/heads/_bundle_modified
scp alamo322b.bundle nova:/work/brunnels/jackplum/alamo322b.bundle

# NOVA: clone both arms
mkdir -p /work/brunnels/jackplum/alamo-322b-ab
cd /work/brunnels/jackplum/alamo-322b-ab
git clone /work/brunnels/jackplum/alamo322b.bundle baseline
git clone /work/brunnels/jackplum/alamo322b.bundle modified
git -C baseline checkout _bundle_baseline
git -C modified checkout _bundle_modified

# NOVA: submit build jobs (ARCHES=80 only, sm_80)
cd /work/brunnels/jackplum/alamo-322b-ab
sbatch --parsable baseline/build_322b.sbatch   # -> 11639419
sbatch --parsable modified/build_322b.sbatch   # -> 11639420

# NOVA: submit dependent wall + ncu jobs
sbatch --parsable --dependency=afterok:11639419 baseline/wall_322b.sbatch   # -> 11639477
sbatch --parsable --dependency=afterok:11639420 modified/wall_322b.sbatch   # -> 11639478
sbatch --parsable --dependency=afterok:11639419 baseline/ncu_322b.sbatch    # -> 11639479
sbatch --parsable --dependency=afterok:11639420 modified/ncu_322b.sbatch    # -> 11639480
```

## Deviations from plan / precedent

- Plan Step 1 explicitly calls for two separate clone dirs (`baseline/`,
  `modified/`), unlike the `a100-ab-c1` precedent which used a single
  checkout and swapped one file (`git checkout <sha> -- Elastic.cpp`) between
  arm builds. Followed the plan literally since the 3.2b diff spans multiple
  commits and files (Elastic.cpp, Elastic.H, Stencil.H per the branch), not a
  single-file revert — the c1ab trick doesn't cleanly generalize.
- Build uses `ARCHES=80` only (A100 sm_80), per task Context — not the
  70/80/90 default in `build_alamo_nova_3d.sh`.
- ncu job requests `--mem=128G` and 16 cpus, matching the c1ab ncu precedent
  rather than the plan's generic resource line (plan text didn't specify ncu
  resources beyond "6 launches/arm across two grid sizes" — used the c1ab
  ncu2 job's resource footprint since it is a nearly identical capture; window
  is single deck/grid size `input_3d_centre_bore_256_a2`, not "two grid
  sizes" as PLAN Step 4 describes — flagging this for the collection-phase
  agent to either accept single-grid-size evidence or add a second deck/job).
- Did not create a wrapper script under `docs/agent_plans/.../` since all
  sbatch files were generated directly on NOVA from inline heredocs during
  this session (no local script edits were needed); this NOTES.md documents
  the exact content in the "Exact commands used" section above plus the job
  table, which should be sufficient to reproduce.

## Not done (next phase)

- Waiting for jobs to complete, extracting TinyProfiler tables, running
  fcompare on baseline vs modified plotfiles, ncu CSV analysis, and writing
  results/RESULT.md with the verdict — that's Step 4-5 of PLAN.md, out of
  scope for this submit-only session.
