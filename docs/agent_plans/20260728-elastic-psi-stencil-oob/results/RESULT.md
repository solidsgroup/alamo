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

## §2 Read geometry — characterized. Checkpoint pending 2026-07-28.

`psi` and the operator diagonal start from the same cell-centered
`m_grids[amrlev][mglev]` and the same distribution map, so an `MFIter` index
selects corresponding fabs. Their index types and halo requirements differ:

- `Elastic.cpp:63,76-77` allocates cell-centered `m_psi_mf` with two ghost
  cells.
- `Operator.cpp:467,476-477` allocates node-centered `m_diag` with two ghost
  nodes.
- `Elastic.cpp:406-425` grows each node-centered diagonal valid box by all two
  diagonal ghosts before launching the kernel.
- `Stencil.H:1584-1592` defaults `CellToNodeAverage` to a central stencil and
  reads both the current cell index and one index lower in every dimension.

For a cell-centered fab with valid interval `[L,H]`, the allocated `psi`
interval is `[L-2,H+2]`. The corresponding nodal valid interval is
`[L,H+1]`; `Diagonal` grows it by two to `[L-2,H+3]`. Therefore the default
cell-to-node average can request `L-3` on the low halo and `H+3` on the high
halo. Both are one cell outside the `psi` allocation. The global `stencilbox`
clips physical boundaries, but it does not clip an interior fab boundary; an
interior low halo remains `Central`. This explains why the defect appears on
the multi-box 3D layout. At a non-periodic physical low boundary, the clipped
node is legal in storage because `i-1` is only the first `psi` ghost.

The other two Elastic call sites do not have the same launch geometry:

| Site | Node launch extent relative to cell valid `[L,H]` | `psi` reads | Verdict |
|---|---|---|---|
| `Fapply`, `Elastic.cpp:207-210,255` | At most `[L-1,H+2]` (`validbox().grow(1)`) | `[L-2,H+2]` | In bounds with the existing two ghosts |
| `Diagonal`, `Elastic.cpp:406-425,461` | `[L-2,H+3]` (two diagonal ghosts) | `[L-3,H+3]` | **Out of bounds by one on both fab sides** |
| `Stress`, `Elastic.cpp:645-650,669` | Ungrown nodal tile `[L,H+1]` | `[L-1,H+1]` | In bounds with the existing two ghosts |

This narrows N11 to a halo-width mismatch in `Diagonal`; it is not evidence
that every stencil-less `CellToNodeAverage` call is unsafe. No fix has been
selected or applied. Per the tier-3 plan, Step 2 stops here for the required
user checkpoint before Step 3 evaluates the alternatives.

## §3 Fix selection — design recorded. Checkpoint pending 2026-07-28.

**Selected: grow only the cell-centered `psi` coefficient from two to three
ghost cells.** Keep `m_ddw_mf` at two ghosts by splitting the current shared
`model_nghost` allocation constant into separate model and `psi` widths. Do
not change `CellToNodeAverage` or any of its three Elastic call sites.

The required width follows directly from the access contract: `Diagonal`
computes two node ghosts, and a cell-to-node average needs one additional
backing cell. Three `psi` ghosts cover `[L-3,H+3]` while preserving the exact
existing interpolation at every node. `Elastic.cpp:1094-1097` calls
`FillBoundaryAndSync` after coefficient averaging; AMReX's no-width overload
uses the MultiFab's full `nGrowVect()` (`AMReX_FabArray.H:3505-3510`), so the
third interior/periodic ghost layer receives real neighboring data. The kernel
is clipped to the physical domain in non-periodic directions, so unfilled
physical exterior ghosts are not consumed.

Only `psi` grows. Growing `m_ddw_mf` would add a third halo to its
`AMREX_SPACEDIM + 1` `Matrix4` components without satisfying any demonstrated
need, creating avoidable memory traffic and footprint in a memory-strategy
campaign.

Rejected alternatives:

1. **Pass `sten` to `CellToNodeAverage`.** `sten` is derived from the global
   `stencilbox`, not the local fab bounds. It remains `Central` on interior fab
   halos and therefore does not prevent the demonstrated `L-3`/`H+3` accesses.
   At physical boundaries it would also change interpolation values, making a
   numerics change without fixing the multi-box storage defect.
2. **Clamp the read.** Clamping to the local fab allocation would duplicate a
   box-edge value instead of reading the neighboring box's valid coefficient.
   The diagonal would become decomposition-dependent. Clamping to the global
   domain would leave the interior-fab OOB unchanged.

Per-site disposition:

- `Elastic.cpp:461` (`Diagonal`): no call-site edit; its required backing
  extent is supplied by the third `psi` ghost.
- `Elastic.cpp:255` (`Fapply`): explicitly unchanged because its one-node
  grown launch is already covered by two `psi` ghosts.
- `Elastic.cpp:669` (`Stress`): explicitly unchanged because its ungrown
  nodal launch is already covered by two `psi` ghosts.

This design preserves boundary-stencil numerics and changes only coefficient
storage. Per the tier-3 plan, no source edit occurs until the user confirms
this Step 3 choice.

## §4 Implementation — source applied; rebuilt candidate green.

The selected allocation change is present in `src/Operator/Elastic.cpp`: model
coefficients retain two ghosts and only `m_psi_mf` receives three.

First `FULL=1 bash benchmark/status.sh` result:

| Leg | Result |
|---|---|
| device lint | PASS |
| `golden-gpu-strict` | PASS |
| Tier-2 memcheck | FAIL, but **not a candidate-binary result** |

Tier-2 used `bin/alamo_gpu-3d-cuda86-g++` timestamped 13:17, while the source
edit is timestamped 15:35. The sanitizer frame still identifies
`Elastic.cpp:461`; in the edited source the interpolation moved to line 464.
This is the stale clean-tree binary already warned about in §1, not evidence
that the three-ghost candidate was compiled and failed.

Logs:
`benchmark/_a100_gate_20260728_153809/tier2_memcheck.log` and
`benchmark/_gate_logs/local_a100_gate.log`.

No correctness verdict is claimed. Per the task's stop-on-failed-VERIFY rule,
execution stops before rebuilding and rerunning.

The user then authorized the required rebuild. `./configure --dim=3 --comp=g++
--cuda 86 && make -j` rebuilt `bin/alamo_gpu-3d-cuda86-g++` at 15:44, newer
than the 15:35 source edit. A solo `FULL=1 bash benchmark/status.sh` on that
candidate passed:

| Leg | Result |
|---|---|
| device lint | **PASS** |
| `golden-gpu-strict` | **PASS** |
| Tier-1 runtime-strict smoke | **PASS** |
| Tier-2 compute-sanitizer memcheck | **PASS — `ERROR SUMMARY: 0 errors`** |

Valid candidate logs:
`benchmark/_a100_gate_20260728_154505/tier1_runtime_strict.log`,
`benchmark/_a100_gate_20260728_154505/tier2_memcheck.log`, and
`benchmark/_gate_logs/local_a100_gate.log`.

**Verdict:** the three-ghost `psi` allocation removes the demonstrated
`Elastic::Diagonal` invalid reads without moving the strict GPU golden values.

## §5

Not started.
