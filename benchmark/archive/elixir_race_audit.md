# Elixir-Race Audit (A4)

**Date:** 2026-06-25
**Auditor:** automated grep + manual review
**Scope:** all ALAMO custom GPU `MFIter` loops in `src/Operator/`, `src/Solver/`, `src/Integrator/`
**Trigger:** `GPU_ROADMAP_V2.md` task A4 — confirm the `interpolation()` race was the only instance

---

## The pattern being audited

A cross-stream use-after-free occurs when:
1. A `FArrayBox` (or `BaseFab`) is created as a local variable inside an `MFIter` loop.
2. Its device memory is passed to an async GPU kernel (`amrex::ParallelFor` / `AMREX_GPU_DEVICE`).
3. No `elixir()` is called — so the device memory is released at the end of the loop body
   while the kernel on a different CUDA stream may still be reading/writing it.

Safe-by-construction: loops that obtain a fab via `multifab[mfi]` or `.array(mfi)` on a
pre-existing `MultiFab` do NOT create new allocations inside the loop — the device memory
is owned by the `MultiFab` and lives beyond the loop body. No `elixir()` needed.

---

## Search methodology

```bash
# Find all local FArrayBox value declarations (not references, not via [mfi])
grep -rn "FArrayBox [a-zA-Z_][a-zA-Z0-9_]*;" src/ --include="*.cpp" --include="*.H"
grep -rn "amrex::FArrayBox [a-zA-Z_][a-zA-Z0-9_]*;" src/ --include="*.cpp" --include="*.H"
grep -rn "BaseFab<[^>]*> [a-zA-Z_][a-zA-Z0-9_]*;" src/ --include="*.cpp" --include="*.H"

# Find .resize() calls inside MFIter-containing files (could indicate FArrayBox resize)
grep -rn "\.resize(" src/Operator/ src/Solver/ src/Integrator/ --include="*.cpp" --include="*.H"
```

---

## Audit table

| File | Line | Variable | Inside MFIter? | Feeds device kernel? | Verdict |
|------|------|----------|----------------|----------------------|---------|
| `src/Operator/Operator.cpp` | 715 | `FArrayBox tmpfab` | **YES** (line 709 loop) | **YES** (`ALAMO_OPERATOR_FOR` at 739, `fine[mfi].plus` at 780) | **FIXED** — `elixir()` at line 728 |

No other local `FArrayBox` or `BaseFab` value declarations exist anywhere in the codebase.

---

## Per-function analysis (all MFIter+device-kernel loops)

### `Operator<Grid::Node>::interpolation()` — `Operator.cpp:709`
**Pattern:** creates `FArrayBox tmpfab` inside the loop, resizes it, passes to `ParallelFor`.
**Status: FIXED** — `amrex::Gpu::Elixir tmpfab_eli = tmpfab.elixir()` at line 728.
This was the root cause of the multi-box MLMG divergence (commit `c00f69086`).

### `Operator<Grid::Node>::restriction()` — `Operator.cpp:920`
**Pattern:** operates directly on `fine_res_for_coarse.array(mfi)` and `fine_res.array(mfi)`.
No temporaries created inside the loop. **SAFE by construction.**

### `Elastic<SYM>::averageDownCoeffsSameLevelAMRLevel()` — `Elastic.cpp:1306`
**Pattern:** `fine_on_crseba` is a `MultiTab` defined *outside* the `MFIter` loop
(line 1301–1303). Inside the loop, access is via `fine_on_crseba.array(mfi)` and
`crse.array(mfi)` — no new allocations. **SAFE by construction.**

### `Elastic<SYM>::averageDownCoeffsSameLevelAMRLevel()` psi loop — `Elastic.cpp:1385`
**Pattern:** `fine_psi_on_crseba` defined outside the loop (line 1381–1383). Inside
the loop, access via `.array(mfi)` only. **SAFE by construction.**

### `Elastic<SYM>::restrictCoeff()` — `Elastic.cpp:1193`
**Pattern:** operates on `fine_ddw_for_coarse.array(mfi)` and `fine_ddw.array(mfi)`.
No temporaries inside loop. **SAFE by construction.**

### `Newton::prepareForSolve()` — `Newton.H:185`
**Pattern:** `model`, `u`, `dw`, `ddw` all accessed via `.array(mfi)` on pre-existing
`Set::Field` members. No new allocations inside loop. **SAFE by construction.**

### `Newton::linearUpdate()` — `Newton.H:80`
**Pattern:** same — all `.array(mfi)` accesses on pre-existing fields. **SAFE by construction.**

### `Newton::computeRHS()` — `Newton.H:134`, `Newton.H:243`
**Pattern:** same — all `.array(mfi)` accesses. **SAFE by construction.**

### `Newton::Solve()` residual/update loops — `Newton.H:580`, `Newton.H:623`
**Pattern:** same — all `.array(mfi)` accesses. **SAFE by construction.**

### `Elastic<SYM>::Fapply()` / `Diagonal()` — `Elastic.cpp:391–475`, `808`
**Pattern:** all accesses via `a_f.array(mfi)`, `a_u.array(mfi)`, `a_diag.array(mfi)`.
No temporaries. **SAFE by construction.**

---

## Conclusion

**The `interpolation()` `tmpfab` (now fixed) was the only instance of the race pattern
in the entire ALAMO custom GPU codebase.** No further elixir fixes are needed.

All other MFIter+device-kernel loops in `src/Operator/`, `src/Solver/Nonlocal/Newton.H`,
and `src/Integrator/` access pre-existing MultiFab data via `.array(mfi)` references —
these hold device pointers into memory owned by the MultiFab, which outlives the loop
body regardless of stream scheduling.

---

## Re-verification

To re-run this audit after future changes:
```bash
grep -rn "FArrayBox [a-zA-Z_][a-zA-Z0-9_]*;" src/ --include="*.cpp" --include="*.H"
```
Any new hit inside an MFIter loop that feeds a ParallelFor requires an `elixir()` before
the first kernel that reads the temporary.
