# RESULT — a100-ab-c1 (2026-07-07)

## What changed

PLAN.md task 3.1 (A100 A/B for the committed C1 `Fapply` register-pressure
edits, `0bb893acc` on `chamber-gpu-elastic-opt`) was blocked last session
(2026-07-06) on non-interactive NOVA ssh from this machine. This session:
set up SSH ControlMaster multiplexing so the user's one interactive login
persists for tool use; confirmed prior jobs from 2026-07-06
(`/work/brunnels/jackplum/alamo-c1ab`, branch tip `123de00a2`) completed;
ran the still-outstanding parity check; collected the retargeted occupancy/
register ncu capture (job `11448954`, gated directly on the
`Operator::Elastic::Fapply()` NVTX range — the first attempt, job
`11433932`, mis-scoped onto init kernels); recorded everything in
`benchmark/PHASE_C1_fapply_occupancy.md` on `chamber-gpu-elastic-opt`
(commit `0ba677e37`).

## Evidence

**Registers / occupancy** (A100, `input_3d_centre_bore_256_a2`, 6 sampled
launches/arm across two AMR-level grid sizes):

| Arm      | registers/thread | achieved occupancy | Fapply/launch (large) | Fapply/launch (small) |
|----------|------------------:|--------------------:|------------------------:|------------------------:|
| BASELINE | 255               | ~12.10%             | 4.667 ms                | 2.286 ms                |
| MODIFIED | 244               | ~11.91%              | 4.256 ms                | 2.077 ms                |
| Δ        | −11 (−4.3%)       | flat (noise)         | **−8.8%**                | **−9.1%**                |

**Wall A/B** (TinyProfiler, job `11433933`, both arms did identical work —
105,920 `Fapply` calls, 13 MLMG solves):

| Region                       | BASELINE | MODIFIED | Δ              |
|-------------------------------|---------:|---------:|----------------|
| `Operator::Elastic::Fapply()` | 300.1 s  | 256.6 s  | −14.5% (1.17×) |
| `MLMG::solve()` (inclusive)   | 373.2 s  | 329.7 s  | −11.7%         |
| Everything else               | —        | —        | within ±0.5%   |

(Raw wall total was confounded by Lustre I/O noise on plotfile writes;
compute-minus-I/O matches the table above.)

**GPU stress parity** (`fcompare.gnu.ex`, steps 0/300/600): cell (flame)
fields bit-identical at every step. Node (elastic) fields: step 0 exact
zero; steps 300/600 relative error ~1e-7–1e-8 on stress/strain/disp,
consistent with GPU fast-math/atomic reduction-order noise — not a physics
regression.

**Verdict: PASS.** Pass condition (registers < 255 and/or occupancy >
~12.5%, stress matches to tolerance) is met via the register axis; occupancy
did not cross a new tier (244 reg/thread × 256 threads/block is still short
of the ≤128 reg/thread needed for a second block/SM on A100's 65536-register
file). The ~9% per-launch and ~14.5% wall Fapply win therefore comes from
reduced register-spill traffic, not from an occupancy increase — consistent
with the edit's by-construction removal of 90 live doubles (135→45 in the
grad(C) accumulation) plus the boundary-only `sig` sink.

## Deviations from plan

None — plan executed as written. The one judgment call: PHASE_A_FINDINGS.md
is referenced by the NOVA A/B procedure doc as a second place to record
results, but it now lives under `docs/archive/` on `chamber-gpu` (moved in
the 2026-06-30 almanac reorg) and CLAUDE.md marks `docs/archive/` forbidden
for edits; skipped it and recorded only in the live doc
(`PHASE_C1_fapply_occupancy.md`) and PLAN.md.

## Next step

Occupancy staying flat (~12%) after the register cut is exactly the trigger
PLAN.md task 3.2 (`__launch_bounds__` sweep) was designed to catch — that is
now the correct next move on this branch, not merely next-in-queue.
