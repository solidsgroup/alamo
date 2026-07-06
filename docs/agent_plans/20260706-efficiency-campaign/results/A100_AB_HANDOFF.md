# Step 3 handoff — A100 A/B for C1 (PLAN.md task 3.1)

**Status: BLOCKED from this machine** — `ssh nova.its.iastate.edu` rejects
non-interactive auth (gssapi/password only, no batch keys), and the runbook
(`benchmark/NOVA_SLURM_RUNBOOK.md`) requires a NOVA login node. Everything
else is staged; this is a copy-paste session for whoever holds NOVA creds.

## What to run

1. Log in to NOVA. Get the code there — the edits are **committed** (not
   uncommitted as the procedure doc's §0 says, which predates commit
   `0bb893acc`): push branch `chamber-gpu-elastic-opt` (local worktree
   `/home/jackplum/Projects/alamo-elastic-opt`, tip `123de00a2`) to a reachable
   remote and pull on NOVA, or rsync that worktree.
2. Follow `benchmark/PHASE_C1_nova_ab.md` (exists on the elastic-opt branch)
   §§1–5 verbatim: build both arms sm_80, ncu occupancy capture
   (`g0_ncu_capture.sh`, input `input_3d_centre_bore_256_a2`), TinyProfiler
   wall A/B via `nova_flame_gpu_3d_a2.slurm`, on-device stress parity, record.

## Pass condition (from the procedure)

`Fapply` registers/thread < 255 and/or achieved occupancy > ~12.5%; stress
field matches baseline to tolerance. If registers drop but occupancy doesn't
move, record honestly and proceed to the `__launch_bounds__` sweep (PLAN.md
task 3.2), which is the designed follow-up for exactly that outcome.
