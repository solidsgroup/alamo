# RESULT — phase05-two-rank-probe

Campaign: `docs/agent_plans/20260727-gpu-memory-strategy/PLAN.md` §5.
Branch: `chamber-gpu-mem`. Started 2026-07-27.

Status: **IN PROGRESS** — Steps 1-2 done, Steps 3-6 pending.

---

## §1 Source confirmations (Step 1) — DONE 2026-07-27

### Q: is the `extensive` flag set for volume / area / mass_flux?

**Yes.** Risk-register row 1 ("`extensives` flag unset → reduction silently
rank-local after all") is **closed as not-a-risk**.

- `src/Integrator/Integrator.H:241` —
  `void RegisterIntegratedVariable(Set::Scalar*, std::string, bool extensive=true);`
  The flag **defaults to true**.
- `src/Integrator/Flame.cpp:195-197` — `volume`, `area`, `mass_flux` are
  registered with **two** arguments, i.e. they take the default `extensive=true`.
- Contrast `Flame.cpp:198-202` and `:228` — the five `thermo_*` max/min
  diagnostics and `chamber_pressure` pass an explicit `false`, which is correct:
  they are intensive and must not be summed across ranks.
- `Flame.cpp:233-235` — the `!thermal.on && variable_pressure` branch registers
  the same three variables the same way. Both registration paths are consistent.
- `src/Integrator/Integrator.cpp:1271-1275` — `ReduceRealSum` is applied to
  exactly the variables whose `extensives[i]` is true. Therefore volume, area and
  mass_flux **are** globally allreduced.

Ordering confirmed: `Integrator.cpp:1187-1190` calls `IntegrateVariables` (which
ends with the allreduce at `:1271-1275`) **before** `TimeStepComplete`, and the
chamber ODE advances at `Flame.cpp:706` inside `TimeStepComplete`. The scalar the
ODE consumes is post-allreduce.

**Verdict: the burning-surface reduction is structurally global. Campaign PLAN
§14 Q2 is fully answered; its "one residual check" is discharged.**

### Q (new, not in the campaign plan): can the chamber scalars go stale?

`IntegrateVariables` is gated at `Integrator.cpp:1229` on
`thermo.interval > 0 && step % thermo.interval == 0` (or a `thermo.plot_dt`
window). `chamber.model.Advance(...)` at `Flame.cpp:706` runs **every** step
regardless. So at `thermo.interval > 1` the ODE would advance on stale
`chamber.{mdot,volume}`.

Not currently a live defect:
- `Integrator.cpp:122-123` — `thermo.interval` **defaults to 1**, and
  `amr.thermo.int` defaults to 1 when queried.
- `input:17` — `amr.thermo.int = 1` explicitly.

Logged as a latent coupling, not a defect: the freshness of the chamber ODE input
is coupled to a *diagnostics-output* interval knob. Any future deck that sets
`amr.thermo.int > 1` with `variable_pressure = 1` silently integrates stale mass
flux. Recorded in the campaign `NOTES.md` (N6); no code change authorized here.

### Sync-site detail refined (feeds campaign PLAN §9)

`Flame.cpp:1127` `reduce_data.value(reduce_op)` is **not** once per box. `Integrate`
is invoked per `MFIter` box *and*, for every level below the finest, once per
element of `amrex::complementIn(box, cfba)` (`Integrator.cpp:1249-1255`). Each
invocation constructs its own `ReduceData` and syncs. The per-step sync count is
therefore `sum over levels of (boxes x complement-pieces)`, not `boxes`. Campaign
PLAN §9's "per box, per level" understates it.

---

## §2 GPU-aware MPI (Step 2) — local leg DONE 2026-07-27

**Local (kermit): GPU-aware MPI is NOT active.**

```
$ ompi_info --parsable --all | grep mpi_built_with_cuda_support
mca:mpi:base:param:mpi_built_with_cuda_support:value:false
mca:mpi:base:param:mpi_built_with_cuda_support:source:default
```
Open MPI 4.1.6, `/usr/bin/mpiexec` (distro build, no CUDA support compiled in).

Consequence: every device buffer handed to MPI is staged through host memory.
This confirms prior evidence (campaign PLAN §5 Q3, §11) at runtime rather than by
inference. It means:

1. A local two-rank GPU probe (Step 5) exercises the **host-staged** comm path.
   That is correct for a *correctness* probe and inadmissible as any Phase 2
   traffic measurement — consistent with §2's kermit-is-correctness-only policy.
2. Campaign PLAN §9's "device-buffer Allreduce" cannot be written as such
   unconditionally. Either NOVA provides a CUDA-aware MPI module, or the staged
   path is accepted and marked as debt in source.

**NOVA re-check (run on a compute node, inside the job, not the login node —
the MPI module differs):**

```bash
ompi_info --parsable --all 2>/dev/null | grep -E 'mpi_built_with_cuda_support|accelerator'
# and, module-agnostic:
srun -n1 --gpus=1 bash -lc 'echo "MPI=$(which mpiexec)"; ompi_info --version 2>/dev/null | head -1; \
  ompi_info --parsable --all 2>/dev/null | grep mpi_built_with_cuda_support'
# MPICH/Cray sites report differently:
echo "$MPICH_GPU_SUPPORT_ENABLED"
```
Fold this into the Phase 0 NOVA batch (campaign PLAN §2: batch by phase boundary,
never a job per question).

---

## §3 Probe harness (Step 3) — DONE 2026-07-27

`benchmark/two_rank_probe.sh` added. Runs the same deck at `-np 1` and `-np N`
into `benchmark/_two_rank_probe_*/` (gitignored, never `baseline_runs/`), then
compares with `benchmark/compare_thermo.py`. Deck overrides mirror
`baseline_suite.py:184-198` so the probe and the golden gate exercise the same
configuration. Knobs: `MGS`, `NP_HI`, `REL_TOL`, `ABS_TOL`, `EXTRA`, `DRYRUN`,
`ARENA_INIT_SIZE`, `DECK`, `BIN`.

**One design correction, caught mid-task.** The first CPU run used the deck's
default grid. The campaign deck sets no `amr.max_grid_size`, so AMReX's 2D
default (128) leaves level 0 as a **single 64×64 box owned entirely by rank 0**.
The probe would have passed partly vacuously at the coarse level. `MGS` now
defaults to 32, which puts ~4 boxes on level 0 and dozens on levels 2-3 across
both ranks. All results below are with `MGS=32`. Same reasoning as `MGS` in
`benchmark/local_a100_gate.sh`.

## §4 Probe results (Steps 4-5) — DONE 2026-07-27

Deck `input`, `amr.max_grid_size=32`, `amr.thermo.int=1`. All legs: **every
thermo.dat column bit-identical, `max_abs = max_rel = 0.0`**, exit 0 at
`rel_tol=1e-6`.

| Leg | Build | Ranks | Steps | Extra | Result | Artifacts |
|---|---|---|---|---|---|---|
| A | cpu `bin/alamo-2d-g++` | 1 vs 2 | 10 | — | PASS, bit-identical | `_two_rank_probe_20260727_181012_cpu` |
| B | cpu | 1 vs 4 | 10 | — | PASS, bit-identical | `_two_rank_probe_20260727_181047_cpu` |
| C | gpu `bin/alamo_gpu-2d-cuda86-g++` | 1 vs 2 | 10 | — | PASS, bit-identical | `_two_rank_probe_20260727_181057_gpu` |
| D | cpu | 1 vs 2 | 2 | `elastic.interval=1` | PASS, bit-identical | `_two_rank_probe_20260727_181125_cpu` |
| E | gpu | 1 vs 2 | 2 | `elastic.interval=1` | PASS, bit-identical | `_two_rank_probe_20260727_181147_gpu` |

Cross-check (not required by the exit gate): CPU np=1 vs GPU np=1 at 10 steps is
also bit-identical on every column.

**Legs D/E exist because the campaign deck sets `elastic.interval = 50`**
(`input:128`). A 10-step probe therefore never enters MLMG at all — legs A-C say
nothing about the elastic path under decomposition. With `elastic.interval=1` the
step-2 row carries non-zero boundary tractions (`trac_xhi_x = -7073.74`,
`trac_yhi_y = -6341.71`), and those match bit-exactly across rank counts. The
elastic/MLMG path is therefore covered, not assumed.

**On the strength of "bit-identical".** Do not over-read it. Most boxes
contribute exactly `0.0` to the burning-surface sums (`dvol`, `darea`, `dmdot`
are non-zero only near the interface), so the floating-point association across
ranks coincides more often than it would for a dense reduction. The load-bearing
evidence is not the last bit: **a rank-local reduction would halve `volume`,
`area` and `mass_flux` at 2 ranks and quarter them at 4** — an O(1) relative
error, impossible to miss. Nothing of the sort appears.

## §5 Verdict (Step 6)

**Phase 0.5 exit gate: PASS on kermit.** Two-rank (and four-rank) runs reproduce
single-rank chamber pressure history exactly, on both the CPU and CUDA builds,
with and without the elastic solve in the loop. The campaign's headline
multi-rank correctness risk is closed. No defect found; no tier-3 folder spawned.

Answers delivered to the campaign plan:

- §5 Q1 / §14 Q2 residual — **closed**, `extensive` defaults true (§1).
- §5 Q2 (domain-decomposition survival: ghost fill, BC at rank boundaries) —
  **closed** at the integrated-quantity level, CPU and GPU, elastic on and off.
- §5 Q3 (GPU-aware MPI) — **closed locally: inactive** (§2). NOVA leg outstanding,
  folded into the Phase 0 batch.
- Risk-register row 1 (`extensives` unset) — **closed as not-a-risk**.

### What this does NOT establish

1. **Per-cell field agreement.** thermo.dat is a reduction; a sign-cancelling
   field error at a rank boundary would survive it. If Phase 1's arena flip or
   Phase 2's reduction rewrite touches ghost handling, add a plotfile-level
   multi-rank compare rather than leaning on this result.
2. **Long-horizon behaviour.** 10 steps (2 with elastic). Nothing here speaks to
   divergence over a 6.5 s burn.
3. **Multi-*GPU*.** Legs C/E put both ranks on the *same* A1000. One-rank-per-GPU
   is a different code path (peer transfers, per-device arenas) and is a NOVA
   question.
4. **Multi-node**, and anything at all about performance (§2: kermit is
   correctness-only, shared and 50 W-capped).
5. **CUDA-aware MPI paths.** Local MPI stages through host (§2), so legs C/E
   never exercised a device-buffer transfer.

### Carried into Phase 0

- Run the §2 GPU-aware MPI probe on a NOVA **compute** node inside the Phase 0
  batch. It gates whether campaign PLAN §9's "device-buffer Allreduce" is
  literal or debt.
- Re-run leg C/E shape on NOVA as **one rank per GPU** before Phase 4 planning.
- Campaign PLAN §9's sync-site table understates `Flame.cpp:1127`: the sync count
  is `Σ_levels (boxes × complement-pieces)` per step, not `boxes` (§1).
- `NOTES.md` N6 (chamber ODE freshness coupled to `amr.thermo.int`) is open and
  unactioned.
