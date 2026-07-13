# RESULT: fsmooth-launch-fusion

Status: **done, all gates PASS, committed on `fapply-322b`** (stacked on
20260709-fapply-kernel-surgery; same worktree, same merge path to chamber-gpu).

## What changed

`src/Operator/Operator.cpp`, `Operator<Grid::Node>::Fsmooth`: the Jacobi
update launched one kernel per component per box per sweep
(`for (n) { ParallelFor(bx, 3D) }`). Replaced with a single fused
`amrex::ParallelFor(bx, ncomp, ...(i,j,k,n))` launch — the AMReX-documented
GPU form ("component loop is moved to the innermost loop"). Loop-invariant
`auto m_omega = this->m_omega;` hoisted above the kernel (still a local
scalar copy — no implicit `this` capture). Update expression and branch
structure byte-identical. Bit-exact by construction: each (i,j,k,n) update
touches only its own component; launch partitioning cannot affect results.

Effect: 2x (2D) / 3x (3D) fewer kernel launches in the MLMG smoother — the
launch-latency win concentrates at coarse MG levels where boxes are tiny.

## Gates (all PASS, logs in this dir)

- device lint: `non-allowlisted violations: 0` (lint.log)
- strict CPU golden compare: all 4 decks `ok`, PASS (golden_compare.log)
- compute-sanitizer memcheck, 2D full-solve (max_step=55): ran to normal
  completion, `ERROR SUMMARY: 0 errors` (memcheck.log)

## Timing (local A1000, indicative)

Deck `input`, max_step=251, 3 runs, GPU idle. After-fusion median: total
63.83s, Fsmooth excl 1.239s, Fsmooth incl 8.908s, Fapply NCalls=63133 —
identical to pre-fusion baseline NCalls (solve path unperturbed). Total flat
vs baseline 63.33s (within noise; 50 W-capped shared A1000). Known artifact
gap: prior task saved only grep'd Fapply/MLMG/total rows, so no pre-fusion
Fsmooth row exists for a direct row diff (timing_summary.txt). Launch-count
reduction is structural and architecture-independent; wall judgment on A100
deferred with the rest of the branch.

## Deviations

None beyond the missing pre-fusion Fsmooth row noted above.
