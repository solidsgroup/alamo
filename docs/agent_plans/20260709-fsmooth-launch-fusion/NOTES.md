# NOTES — fsmooth-launch-fusion (2026-07-09)

Backlog findings from the same Fsmooth read, NOT in this task's scope:

1. **Per-call MultiFab temporaries.** `Operator<Grid::Node>::Fsmooth`
   allocates Ax, Dx, Rx (3 full MultiFabs) on EVERY smoother call
   (Operator.cpp:361-363), inside the hottest MLMG loop. AMReX guidance is
   stream-ordered temporaries (The_Async_Arena) or operator-owned scratch.
   Caching them per (amrlev, mglev) like m_diag would kill the alloc/free
   churn. Not bit-exactness-sensitive, but structural — needs its own task.

2. **Fusible elementwise chain.** Dx=Copy(x); Dx*=diag; Rx=Copy(Ax); Rx-=Dx;
   then the update kernel reads Rx. Four full-fab passes + the update could
   fold to one kernel computing (b - (Ax - x*diag)) inline. NOT bit-exact-
   trivial: nvcc/gcc FP contraction (FMA) can change rounding vs the stored-
   intermediate version — needs the golden gate and possibly
   -ffp-contract pinning. Bigger win than launch fusion; tier 3.

3. **Dead branch in update kernel**: `if (!domain.contains(i,j,k)) {}` empty
   body (Operator.cpp:396-399) — harmless, left untouched to keep this task
   purely mechanical.

4. **Linear.H:83-93 per-component ParallelFor** (solve setup, negation/BC
   zeroing) — cold path (once per solve), plus 4x `Util::Message(norm0())`
   diagnostic reductions around it. Skipped: no measurable win.

5. **IC/Constant.H:62 per-component ParallelFor** — init only. Skipped.
