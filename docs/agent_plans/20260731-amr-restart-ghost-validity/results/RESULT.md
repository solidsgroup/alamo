# AMR restart ghost-validity result

Status: diagnosis complete; source fix and verification are waiting at the
Tier-2 plan checkpoint.

## Root cause

The checkpointed primary fields and their required halos are finite.  The
missing invariant is the persistent solid-temperature field `temps_mf`:

- `Flame::Parse` registers `temps_mf` with `writeout=false`.
- `Integrator::Restart` allocates every cell field, but copies only fields
  present in the checkpoint and does not initialize unmatched cell fields.
- The thermal kernel reads `temps_mf` before writing it:
  `Tsolid = dTdt + temps * (etanew - eta) / dt`.
- Consequently the first post-restart thermal update consumes uninitialized
  state.  Its output remains finite in the ordinary CPU allocation observed
  here, so that substep completes; the corrupted temperature then causes the
  Arrhenius mobility to overflow on the next fine substep.

This is not a regrid failure and not a bad checkpoint copy.  Disabling both the
generic runtime regrid and Flame's one-time refinement hook leaves the failure
unchanged.

## Evidence

### Exact failure location

A one-rank CPU restart with `amr.verbose=1` shows:

```text
[Level 0 step 11] Advanced 1024 cells
[Level 1 step 21] Advanced 2304 cells
[Level 1 step 22] ADVANCE with dt = 5e-06
non-finite value detected in Flame::Advance phase-field kernel at lev=1
```

Log:
`/tmp/alamo-c2-reproverbose.EgBdOw/run.log`.

Thus the first fine substep succeeds and the second fails; the failure is not
on the first level-1 phase-field evaluation after loading the checkpoint.

### Uninitialized-state tripwire

Running the same restart with `amrex.init_snan=1`, with both regrid paths
disabled, moves the failure to the earliest consumer:

```text
[Level 0 step 11] ADVANCE with dt = 1e-05
non-finite value detected in Flame::Advance thermal kernel at lev=0
```

Log:
`/tmp/alamo-c2-snan.1gnHUj/run.log`.

This is direct evidence that a cell field omitted from the checkpoint is read
before initialization.  Static tracing identifies that field as `temps_mf`;
the other non-checkpoint cell fields used in the phase/thermal pair are written
before being read or are not on this path.

### One-substep isolation

With `amr.nsubsteps=1`, the restart completes one level-1 thermal update and
writes a plot before a second phase update can consume the result:

```text
[Level 1 step 21] Advanced 2304 cells
STEP 11 ends. TIME = 0.0001099999975 DT = 1e-05
```

The checkpoint temperature is approximately 300 K:

| Level | Checkpoint minimum | Checkpoint maximum |
|-------|-------------------:|-------------------:|
| 0 | 300.000005 | 300.866286 |
| 1 | 300.000001 | 300.892429 |

After the isolated post-restart update:

| Level | Updated minimum | Updated maximum |
|-------|----------------:|----------------:|
| 0 | -0.0156235 | 299.731208 |
| 1 | -0.0377732 | 299.843516 |

Run and plot:
`/tmp/alamo-c2-onefine.1M064x`.

At the next fine substep, `Homogenize::get_L` evaluates
`exp(-E_prop / T)` with the corrupted near-zero negative temperature.  The
exponent becomes large and positive, producing the non-finite mobility caught
by the phase-field guard.

### Checkpoint and halo checks

- Cell valid regions round-trip exactly at both levels for all eight stored
  fields.
- Nodal `phi` and `chi` valid regions round-trip exactly.
- A fresh step-10 run and a restart-only step-10 dump show finite `eta` and
  `temp` through all three allocated ghost layers, and finite `phi`/`chi`
  through both allocated nodal ghost layers.
- Offline evaluation of every phase-kernel expression over all level-0 and
  level-1 valid cells, using the restored fields and the required `eta` halo,
  produces only finite values.

Fresh/restart ghost comparison:
`/tmp/alamo-c2-freshghost.yMCMQb` and
`/tmp/alamo-c2-ghosts.tyHs77`.

## Proposed repair

Treat `temps_mf` as restart state, not disposable scratch.  The preferred
minimal repair is to include it in the cell checkpoint and verify a newly
written two-level checkpoint through the CPU and strict two-rank GPU restart
oracles.  The implementation must also make the behavior for older
checkpoints that lack `temps` explicit; silently reading uninitialized memory
is not acceptable.

The same registration audit found one adjacent latent state gap:
`thermal.has_exceeded_Tcutoff` is persistent refinement history but is also
registered `writeout=false`.  It is not the C2 failure because C2's cutoff is
zero, so its restart-time read cannot take the reset branch while temperature
is positive.  A generic repair should either persist that field too or
deliberately reconstruct it; leaving it allocator-dependent would preserve the
same lifecycle defect on decks with a positive cutoff.

No application source was changed during this diagnosis.

## Checkpoint

- [x] Reproduced on CPU and strict GPU.
- [x] Identified the exact missing field and earliest invalid lifecycle point.
- [ ] Human approval to edit source.
- [ ] Focused CPU/GPU restart verification.
- [ ] `FULL=1 bash benchmark/status.sh`.
- [ ] Human pre-commit approval.
