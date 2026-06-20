# Task A: Lever 2.5 — field packing (cut FillBoundary launches)

## Goal
Reduce the phase-field GPU launch count by packing co-evolving, ghost-bearing,
cell-centered fields into multi-component MultiFabs so a single FillBoundary /
FillPatch covers them instead of one per field. FillBoundary is **56% of all
cudaLaunchKernel** in the parity nsys capture. Target: a measurable launches/step
drop AND a wall/step improvement, with the golden compare preserved bit-for-bit.

This lever is **low-ROI and high refactor risk** (best case ~25% fewer launches;
the 30% AMR-interp launches and everything else stay subcycle-multiplied). The
user wants it attempted for diligence. **A well-documented negative result
(no ≥5% wall/step gain, reverted to keep numerics safe) is a valid, expected
outcome.** Do not force a fragile change to manufacture a win.

## Context (only what this task needs)
- Repo `/home/jackplum/Projects/alamo`, branch `chamber-gpu` @ `fe3844f31`, main tree.
- Read `../PLAN.md` for ownership rules. You OWN `src/`, all builds, `bin/`,
  `tmp_build_dir/`, `.alamo`, `Makefile`, `benchmark/phase2_box_sweep/`.
- **Evolving ghost-bearing fabs** (declared in `Flame.H`, set up in `Flame.cpp`
  ~122–191): `eta_mf` (BC `bc_eta`, nghost 3), `psi_mf` (BC `bc_eta`, nghost 2),
  `temp_mf` (BC `bc_temp`, nghost 3). Swap partners `eta_old_mf` (g2) /
  `temp_old_mf` (g3) are `evolving=false`, swapped via `std::swap` at
  `Flame.cpp:611`. nghost=0 fabs (mdot, heatflux, L, alpha, temps) generate **no**
  FillBoundary — packing them saves nothing.
- **Key opportunity / blocker:** `eta` and `psi` BOTH use `bc_eta` → packing
  eta+psi into one multi-component fab (nghost=max(3,2)=3, BC `bc_eta`) is the
  **low-risk pack** (same BC). `temp` uses a different BC (`bc_temp`) → packing
  temp with eta needs a **combined multi-component BC** (eta-cond on comp0,
  temp-cond on comp1); `BC::Constant` supports component ranges. Do the same-BC
  pack first; attempt the cross-BC pack only if the first stays green and time allows.
- Framework: `RegisterNewFab(Set::Field&, BC*, ncomp, nghost, name, writeout,
  evolving)` at `Integrator.cpp:340`; per-fab FillPatch via `physbc_array`
  (`Integrator.cpp:394`, `1218`). Packing touches: field decls (`Flame.H`),
  registration, **every access** (`.array(mfi)` / `.Patch(lev,mfi)` / `[lev]`),
  BC construction, and the swap logic.

## Files to read first
- `src/Integrator/Flame.H`, `src/Integrator/Flame.cpp` (field decls/registration/
  accesses/swap), `src/Integrator/Integrator.cpp` (RegisterNewFab, physbc_array).
- `../PLAN.md`. `benchmark/baseline_suite.py` and `benchmark/phase2_box_sweep.py`
  (to know exactly what the gates run).

## Files allowed to modify
- `src/Integrator/Flame.cpp`, `src/Integrator/Flame.H`
- BC construction code these reference, ONLY if a combined BC is needed and ONLY
  within Flame's setup (do not change shared `BC::Constant` semantics for others).

## Files NOT allowed to modify
- Anything Track B owns (`input`, `input_phase1_*`, `benchmark/phase1_elastic_2048/`,
  `analysis/results_phase1*/`).
- Shared framework headers beyond what the pack strictly requires. Do NOT broaden
  scope or do unrelated cleanup.

## Implementation steps
1. **Baseline measurement first.** `source benchmark/local_cuda_env.sh`. Capture
   the parity nsys run (command in PLAN / handoff §5) and record per-step
   `cudaLaunchKernel`, `cudaStreamSynchronize`, and the FillBoundary share. Confirm
   which of eta/psi/temp actually generate **per-step** FillBoundary (check the
   `evolving` flag and FillPatch call sites) — if `psi` is effectively static,
   packing it saves nothing; adjust the plan to whatever truly fires per step.
2. **Same-BC pack (eta+psi) first** if both fire per step: combine into one
   2-component `Set::Field` with `bc_eta`, nghost 3. Update declaration,
   `RegisterNewFab`, every access (component-index the formerly-separate fields),
   and the swap logic. Build CPU + gpu_fast + gpu_strict (all three — the change
   affects every build). Re-run **both gates**: `baseline_suite.py check` (expect
   9/9 ok, strict must be bit-identical) and `phase2_box_sweep.py` (5/5). If the
   golden compare breaks → revert this step, record why, stop.
3. **Re-capture nsys** identically to step 1. Record launches/step before→after
   and the FillBoundary-share delta. Measure wall/step vs the CPU node on the best
   box-sweep config (`wide_512_bf32_mgs128`).
4. **Cross-BC pack (add temp)** ONLY if step 2 is green and time allows: build the
   combined multi-component BC, pack temp in, repeat the full build + both gates +
   nsys + wall/step. Revert on any golden break.
5. **Decision (D2):** if marginal wall/step improvement < 5%, stop pursuing further
   packing. Record the attribution table.

## Invariants (device, concurrency, architecture)
- Golden compare (strict/no-fast-math) must stay **bit-for-bit** vs
  `benchmark/baseline_references/`. Any divergence = revert.
- No new host-loop writes into device-arena fabs.
- Do not change observable thermo.dat output ordering/semantics.
- You may rebuild bin/ freely; Track B uses a private binary copy, so your
  rebuilds will not disturb it. Never touch the shared `input` file's contents.

## Build and test commands
```bash
source benchmark/local_cuda_env.sh
./configure --comp=clang++ --dim 2 && make -j"$(nproc)"
COMP=g++ CUDA_FP=fast   SMOKE=0 ./benchmark/build_alamo_local_gpu.sh
COMP=g++ CUDA_FP=strict SMOKE=0 ./benchmark/build_alamo_local_gpu.sh
CPU_NP=8 GPU_FAST_NP=1 GPU_STRICT_NP=1 python3 benchmark/baseline_suite.py check
python3 benchmark/phase2_box_sweep.py
```

## Expected result
Either a green pack with launches/step + wall/step deltas quantified (and a
recommendation to commit), OR a documented negative result reverted to keep
numerics safe. Leave your final src state in the working tree for lead review.
**Do NOT git commit.**

## Non-goals
- 3D, multi-GPU, CUDA Graphs (2.6), elastic, or any roadmap phase other than 2.5.
- Performance work beyond field packing.

## Stop conditions
- Golden compare breaks and can't be made transparent within the pack → revert, stop.
- Marginal wall/step gain < 5% after the achievable packs → stop, document.
- Build can't be made to link → record the exact error, stop.

## Final report: write `../results/A-RESULT.md`
```
# Result: Task A — Lever 2.5 field packing
## Summary (green-with-win | green-no-win-reverted | blocked)
## Files changed (left in working tree, uncommitted)
## Measurements: launches/step before→after, FillBoundary share, sync fraction,
   wall/step vs CPU node, golden-compare result (9/9? bit-identical?)
## What was packed and why (and what was left unpacked)
## Tests run and results
## Issues found
## Deviations from task
## Recommendation: commit or revert, with the numbers behind it
```
