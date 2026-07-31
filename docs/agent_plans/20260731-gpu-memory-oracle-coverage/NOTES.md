# Notes

## 2026-07-31 — Step 1 blocked by AMR restart state

The proposed two-rank, level-1 C2 continuous leg completed 20 steps and wrote
both `00010cell/Level_1` and `00010node/Level_1`.  The restart leg did not reach
step 11:

- an ordinary two-rank GPU launch reported CUDA error 700 while completing
  restart;
- `CUDA_LAUNCH_BLOCKING=1` changed the observable failure to the existing
  finite-value tripwire at `Flame::Advance`, level 1, step 11;
- a one-rank CPU restart from the same checkpoint failed at the same tripwire.

This rules out MPI rank count and GPU-only execution as the primary cause.  It
also means weakening the oracle to omit refined-state restart would conceal a
real defect.  Step 1 is paused pending the source-level prerequisite task
`20260731-amr-restart-ghost-validity`.
