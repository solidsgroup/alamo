# RESULT — elastic-psi-stencil-oob (N11)

## §1 Provenance — PRE-EXISTING. Step 1 done 2026-07-28.

Question: does the `Elastic::Diagonal` OOB reproduce without the uncommitted
`Flame.{cpp,H}` changes?

It matters because that diff is not inert with respect to this defect — it
touches `elastic.chi_refinement_criterion` / `casing_support_refinement_criterion`,
i.e. AMR tagging, which changes the BoxArray, which changes the box geometry a
low-corner OOB depends on.

Method: `git stash push src/Integrator/Flame.{cpp,H}` → `./configure --dim=3
--comp=g++ --cuda 86` → `make -j` → `TIERS=2 bash benchmark/local_a100_gate.sh`
→ `git stash pop`. Stash restored under an EXIT trap; `git status --porcelain
src/` confirms both files back.

| tree | invalid reads | ERROR SUMMARY | fault site |
|---|---:|---:|---|
| dirty (Flame changes present) | 2,497 | — | `Elastic.cpp:461` |
| **clean (HEAD)** | **2,337** | **2,338 errors** | **`Elastic.cpp:461`** |

Logs: `benchmark/_a100_gate_20260728_125617/` (dirty),
`benchmark/_a100_gate_20260728_131743/` (clean).

**Verdict: pre-existing on `chamber-gpu-mem`.** The uncommitted work is not the
cause. The ~6% difference in error count is consistent with the tagging diff
shifting the box layout slightly; it does not change the conclusion.

No checkpoint triggered — Step 1's stop condition was a GREEN clean tree, and it
is red.

**Note for later steps:** `bin/alamo_gpu-3d-cuda86-g++` is now built from the
CLEAN tree and is therefore stale with respect to the restored working
directory. Rebuild before any measurement that depends on the Flame changes.

## §2-§5

Not started.
