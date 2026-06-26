# Current — chamber-gpu

Last updated: 2026-06-22. This file is overwritten, not appended to. When a
line below is done, delete it and (if it closed out a phase) add an entry to
`changelog/` + `VERSIONS.md` instead.

## In flight

**GPU test suite — ✅ GREEN 2026-06-26 (9/9; was 5P/4F).** The deferred 4-failure
snapshot below is fully analyzed and fixed. Four failures = four different bugs:
3 stale-input decks (F2/C1/C3: old `m_ap`-era Homogenize params + `model_ap/htpb`
elastic schema → current `*_prop` + `model_prop/void/casing`) and one checkpoint
deck (C2). Two **real source defects** were found and fixed in
`src/Integrator/Integrator.cpp`: (1) `Restart()` node-fab out-of-bounds
(`name_array[i][j]` with `i`=restart-fab-index) → segfault on every node-fab
restart; (2) headerless `thermo.dat` on restart (header gated on `step==0`).
Two issues documented + deferred: the GPU boundary-traction diagnostic race in
`Base::Mechanics::Integrate()` (shared infra; C1 compares physics columns
instead) and C2 bit-exact restart physics reproduction (Flame doesn't checkpoint
`temp_old`/`temps`/chamber-ODE). Also: elastic MLMG diverges on the small smoke
grids at >1000× stiffness contrast (CPU+GPU both — lowered casing to 500 MPa).
Writeups: `benchmark/GPU_TEST_SUITE_FIXES.md` (full analysis),
`benchmark/GPU_TEST_PERF_TRACKING.md` (perf baseline),
`changelog/2026-06-26-gpu-test-suite-fixes.md`. CPU 2D + GPU strict 2D binaries
rebuilt with the source fixes. Working-tree edits uncommitted.

**3D elastic on GPU — ✅ FIXED 2026-06-26 (roadmap-v2 A2 prerequisite).** First
attempt to run 3D combined flame+elastic on GPU crashed on the first elastic solve
(`CUDA error 719`). Root cause: chained Eigen expression `F.inverse().transpose()`
faults on device — 2 sites in `src/Model/Solid/Finite/NeoHookean.H` (3D `DW`/`DDW`),
fixed by materializing the inverse first. **Same bug class as the elixir UAF** (works
on CPU, faults on GPU); the 2D elastic "RESOLVED" claim below did NOT cover 3D.
`bin/alamo_gpu-3d-cuda86-g++` rebuilt with the fix (3D nofast/strict NOT yet rebuilt).
Full writeup + a **deferred** GPU-test-suite snapshot (5P/4F, analysis postponed):
`changelog/2026-06-26-3d-elastic-gpu-fix.md`. Working-tree edit uncommitted.

**Phase 1 — ✅ RESOLVED 2026-06-25: GPU elastic works on device (D1 reversal
fulfilled).** The 2048²/multi-box GPU elastic MLMG divergence is root-caused and
FIXED: a GPU **cross-stream use-after-free race** on the per-box temporary
`FArrayBox tmpfab` in `Operator<Grid::Node>::interpolation()`
(`src/Operator/Operator.cpp`) — AMReX's non-tiled GPU `MFIter` cycles iterations
across a pool of CUDA streams and the temp was freed/reused on a different stream
while its kernels were in flight; the corrupted correction was amplified by the
BC-penalty diagonal into a 1e20 blow-up. **Fix (one line):**
`amrex::Gpu::Elixir tmpfab_eli = tmpfab.elixir();` after `tmpfab.resize(...)`.
Verified end-to-end (full box sweep, real production operator multi-box 2048²,
and a multi-box flame+elastic chamber sim). Full writeup:
`benchmark/GPU_BRANCH_GUIDE.md` (D1 section),
`benchmark/elastic_sensitivity_20260621/GPU_ELASTIC_DEBUG_PLAN.md` (SOLVED
banner) + `fix_notes.md`, and Claude memory `gpu-elastic-fixed`. The detail
bullets below are **superseded investigation history** — the
"deterministic coarse-level operator defect" / nondeterministic-reduction /
noise-floor framings were all downstream symptoms of the one race above; the
`Fapply`/coefficient/`Fsmooth` exonerations they record were correct (the bug was
never there — it was the transfer temp's lifetime).

- A second, independent agent (separate worktree) is doing a broader
  GPU-safety audit of `Fapply`/`Diagonal`/residual construction/BC
  application for vtables-in-kernels, host-only captures in device lambdas,
  and host loops writing device-arena fabs. Per `docs/gpu_elastic_device_port_plan.md`
  Phases 1-2, that de-virtualization work was already done and verified clean
  previously, so this audit is likely re-confirming known-good ground rather
  than chasing the bug below — but coordinate before editing
  `Operator/Elastic.cpp`, `Operator/Operator.cpp`, or the `BC/Operator/Elastic`
  headers, since both efforts touch the same files.
- **Current narrow finding (this session, uncommitted):** the GPU elastic
  2048^2 divergence is NOT (only) the fast-math/conditioning question 1.1-1.2
  already answered in `benchmark/PHASE1_ELASTIC_DISPOSITION.md`. A 2-level-MG
  uniform-stiffness repro (`elastic.max_coarsening_level=1`,
  `benchmark/elastic_sensitivity_20260621/`) shows: (a) restricted bottom-level
  operator coefficients are bit-correct on GPU (verified via new
  `ALAMO_ML_COEFF_DIAG` env-gated diagnostic in `src/Operator/Elastic.cpp`,
  reads `MATRIX4` raw storage via `reinterpret_cast` to stay device-safe across
  symmetry specializations whose `operator()` is host-only); (b) the bug
  survives with `bottom_solver=cg` *and* `=smoother` (smoother does zero
  global reductions), which **refutes the standing
  nondeterministic-GPU-reduction root cause** for this specific failure mode
  (that hypothesis may still explain other things, e.g. the bicgstab
  explosion threshold, but not why the coarse-level system can't be solved at
  all); (c) CPU smoother-only on the identical case converges steadily
  (0.84->0.18 over 60 iters), GPU smoother-only plateaus stuck at 0.74 — so
  this is a real, deterministic GPU defect in operator
  application/relaxation at the coarsened MG level, not corrupted
  coefficients and not a reduction race; (d) tested and REFUTED:
  `Mechanics.H:194`'s hardcoded `SetUniform(false)` (the unconditional
  `Cgrad`-modulus-gradient term in `Fapply`) is not the cause — forcing
  `SetUniform(true)` made the GPU smoother result noisier, not better.
  `Fsmooth` itself (`src/Operator/Operator.cpp:350-417`, generic shared
  Jacobi smoother for all `Grid::Node` operators) has been reviewed and looks
  structurally GPU-clean but is not yet stress-tested in isolation.
  Full detail: Claude memory `gpu-elastic-d1-reversed` (and
  `gpu-elastic-nondeterminism` for the now-partially-superseded prior root
  cause). No task-folder record yet — should get one under
  `docs/agent_plans/<date>-gpu-elastic-coarse-mg-defect/` per
  `docs/llm/CONVENTIONS.md` if this continues past this session.
- **Next probe:** a single isolated `Fsmooth`/`Fapply` call at mglev=1 with a
  known simple input (not a full multi-hundred-iteration solve), diffed
  against a hand-computed or CPU-computed expected value. CPU rebuild is
  currently blocked (`.make/Makefile.pre.conf` pinned to the GPU/nvcc
  toolchain for this checkout; `make bin/alamo-2d-g++` silently no-ops) —
  not yet fixed, to avoid disrupting the live GPU configure state.

**Phase 3 — regime scaling / crossover hunt** (`docs/llm/ROADMAP.md` Phase 3).
Live task folder: `docs/agent_plans/20260620-gpu-phase3-regime-scaling/`.

- First crossover measured on NOVA 2026-06-21 (jobs 11160767-11160774, A100-80
  sm_80, 16-rank CPU baseline, exploratory): single A100 beats the CPU node
  ~39x @ 128^3 and ~70x @ 256^3. **These numbers are superseded** — confirmed
  2026-06-22 against the full 64-rank CPU node (jobs 11161369-11161379):
  **9.6x @ 128^3, 13.1x @ 256^3** (decision **D3 = WIN @ single**, confidence
  HIGH). Both runs are written up in `benchmark/PHASE3_R3_crossover.md`; treat
  the 9.6x/13.1x figures as the current headline, not 39x/70x. These figures
  are elastic-disabled (flame-only) — see the Amdahl's-law caveat in
  `docs/llm/ROADMAP.md` before quoting them as a full-chamber speedup.
- Working tree right now has **uncommitted, in-progress edits on top of that**
  (not yet captured in a RESULT.md): `src/Integrator/Flame.{H,cpp}`,
  `src/Operator/Operator.cpp` (+202 lines), `Makefile`, `configure`, and the
  NOVA 3D build/slurm scripts (`benchmark/build_alamo_nova_3d.sh`,
  `benchmark/nova_flame_{cpu,gpu}_3d*.slurm`, `benchmark/phase3_nova_diag.sh`,
  `benchmark/phase3_scaling_sweep.sh`). `benchmark/PHASE3_R3_crossover.md`
  itself is also modified in the working tree vs the committed version.
- **Still open per the crossover report:** the full 64-rank CPU-node
  confirmation is now done (see above); `ncu` occupancy/register metrics are
  still missing. **Descoped 2026-06-24 (explicit user directive):** the 512^3
  CPU baseline (`--mem=32G` hardcode OOM in `nova_flame_cpu_3d.slurm`) is
  **dropped, not being pursued** — 128^3/256^3 give enough crossover data on
  their own. **Multi-GPU is on the back burner — not a near-term goal.** It
  remains a confirmed red flag (2 GPUs 3.5x slower @128^3, 13.3x @256^3, 9.0x
  @512^3 vs 1 GPU), but no further `nsys` diagnosis or scaling work is
  planned unless multi-GPU becomes a near-term goal again.

**Next:** finish/commit the in-progress `Operator.cpp`/`Flame.cpp` work, then
re-run the standing metric set (launches/step, sync fraction, wall/step vs
CPU node — golden-compare residual no longer required, see below) for
128^3/256^3 only (see `docs/llm/ROADMAP.md` "Testing and reporting cadence")
before extending `PHASE3_R3_crossover.md`. 512^3 and multi-GPU data are out
of scope for this pass.
**Test:** `CPU_NP=8 GPU_FAST_NP=1 GPU_STRICT_NP=1 python3 benchmark/baseline_suite.py check` (correctness gate, must stay green through any edit) plus the NOVA sweep via `benchmark/phase3_scaling_sweep.sh --submit`.
**Write results to:** a new `docs/agent_plans/20260620-gpu-phase3-regime-scaling/results/004-RESULT.md`, then fold the headline numbers into `benchmark/PHASE3_R3_crossover.md`.

**Descoped 2026-06-24 (explicit user directive): no-fast-math golden compare
at a saturating 3D config is no longer required.** The branch
definition-of-done item 1 (`benchmark/PHASE5_BRANCH_DONE.md`) is closed out as
DONE-at-coarse-only, not PARTIAL — rationale: CPU-vs-GPU correctness parity
is not the current problem; performance and stability (the elastic coarse-MG
defect, §Phase 1 above) are. Revisit only if a future correctness regression
is suspected at scale.

## Also uncommitted, not yet folded into a task record

- `.github/workflows/chamber-gpu-correctness.yml` (new, untracked) — Phase 5.1
  correctness-CI wiring. Still has no task record — if picked up, write one
  under a new `docs/agent_plans/<date>-<slug>/` folder per
  `docs/llm/CONVENTIONS.md` before starting, rather than working untracked.

The `analysis/` suite rewrite and its five results bundles were backfilled
2026-06-22 —
`docs/agent_plans/20260622-analysis-suite-backfill/{PLAN.md,tasks/,results/}`.
That record also flags one open item: `analysis/results_wide_grid/` is 5.6G of
untracked nsys traces sitting in the working tree with no gitignore/relocation
decision made yet.

GPU elastic abort/NaN hardening was moved into the isolated worktree `/tmp/alamo-elastic-safe` so it does not race the shared checkout. In that branch, `src/Operator/Elastic.cpp` and `src/Solver/Nonlocal/Newton.H` now use device error flags plus host-side drains for `SetModel`, `Fapply`, `Diagonal`, and the Newton model-prep kernels. `git diff --check` is clean there; a full rebuild was not available because the isolated worktree lacks the generated build metadata from the main checkout.

GPU elastic safety hardening was started 2026-06-22 — `src/Operator/Elastic.cpp` and `src/Solver/Nonlocal/Newton.H` now avoid host-only aborts inside GPU kernels for the hot paths (`Fapply`, `Diagonal`, `SetModel`, and Newton model prep/residual BC dispatch). The build was re-run with `source benchmark/local_cuda_env.sh && make -j8` and completed cleanly.

What/why for the `input_3d_flame_{128,256,512}` Phase-3 sweep sizes is now
written up at
`docs/agent_plans/20260620-gpu-phase3-regime-scaling/TEST_CASES.md`.
