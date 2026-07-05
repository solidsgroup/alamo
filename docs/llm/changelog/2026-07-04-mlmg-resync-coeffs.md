# MLMG stale-hierarchy fix — resync coarse operators after Newton relinearization — 2026-07-04

Branch `chamber-gpu`, committed as `1cc3c32db` ("Improve high-contrast Newton
robustness"). Fixes the high-modulus-contrast elastic divergence that was the
real stability cliff in the 2026-07-02..07-04 MLMG campaign. Full taxonomy:
`benchmark/MLMG_HIGH_CONTRAST_FINDINGS.md`; source-level writeup:
`benchmark/mlmg_high_contrast_20260702/RESYNC_COEFFS.md`.

## TL;DR

| Finding | Class | Status |
| --- | --- | --- |
| Persistent `MLMG` reused across Newton iters keeps stale coarse operators + stale smoother/normalize diagonal after relinearization | correctness defect (explosive divergence at high contrast) | **fixed** — gated `Operator::SyncCoefficients()` call in `Newton` after `prepareForSolve()` on iters > 0 |
| Default input decks | unchanged | gate `elastic.solver.resync_coeffs` defaults off; default behavior bit-identical (golden PASS) |

## Root cause

`Solver::Nonlocal::Newton` reuses one persistent `amrex::MLMG` across Newton
iterations. AMReX runs `linop.prepareForSolve()` — which derives every coarse
MG-level coefficient field and the smoother/normalize diagonal from `mglev=0` —
**only on the MLMG object's first solve**. Newton's own `prepareForSolve()`
calls `SetModel()`, which updates the `mglev=0` coefficients only. So on Newton
iterations >= 2 the V-cycle smooths the freshly relinearized fine-level operator
against **stale coarse operators and a stale diagonal**. At high modulus
contrast this produces explosive divergence within 1-2 V-cycles.

## Fix

`Operator<Grid::Node>` exposes (`src/Operator/Operator.H:66`):

```cpp
void SyncCoefficients() { averageDownCoeffs(); Diagonal(true); }
```

`Solver::Nonlocal::Newton` calls it after `prepareForSolve()` on iterations
after the first, gated on `elastic.solver.resync_coeffs`
(`src/Solver/Nonlocal/Newton.H:398` and `:561`; member `:791`; parse `:946`).
The gate defaults to `false`, so historical decks keep their old trajectory
unless they opt in. This is the "force operator/MLMG regeneration" candidate
(correct, possibly slow); no incremental-update scheme was attempted.

## Evidence (instrumented A/B, not narrative)

Anchor deck, `psi_floor=0`, same binary/config except the gate
(`benchmark/mlmg_high_contrast_20260702/validation_20260704/`):

| run | result |
|---|---|
| `resync_coeffs=1` | PASS — Newton linear solves converged 240/248/228 iters |
| `resync_coeffs=0` (control) | FAIL — Newton iter 2 reproduced stale-hierarchy divergence after 9 MLMG iters, `resid/resid0 ~= 2.56e22` |

Supporting: diag-probe confirmed no `Diagonal()` rebuild at iter 2 pre-fix; the
`nriters=1` two-timestep control converges (isolating relinearization, not the
timestep, as the trigger); the J-collapse hypothesis was instrumented and
refuted for this within-solve mechanism. Strict-3D-CUDA parse and step smoke
tests finalized (`codex_gpu3d_resync_*`).

## Companion input knobs (same commit / campaign)

- `elastic.solver.nr_convergence=psi_update` — gates Newton convergence on
  psi-weighted accepted updates instead of raw max update in near-void
  displacement (the performance win: 13 solves / ~27.5 s vs 178 / ~215 s on the
  anchor first-elastic window).
- `elastic.zero_out_displacement=1` — cold-starts each solve to avoid
  warm-started void displacement / J-collapse accumulation on anchored runs.
- `elastic.solver.line_search=1` — kept enabled (undamped-Newton overshoot fix).

## Phase 5.4 enshrinement + validation (done 2026-07-05, local A1000 sm_86)

- **Permanent regression case** `rod_and_tube_step2` added to
  `benchmark/baseline_suite.py` (driven by `benchmark/ci_golden_compare.sh`):
  full 2D rod-and-tube (`input_rod_and_tube_2d`, expression IC, no quarter
  symmetry) with a near-floating stiff rod coupled to a stiff tube only through
  a thin soft void seam. `elastic.interval=1` fires the elastic solve at step 2;
  Newton reaches iteration 2, so the resync `SyncCoefficients` path is exercised;
  thermo.dat's `disp_*`/`trac_*` columns capture the solve. A resync/MLMG-recipe
  regression trips it via divergence->abort or changed boundary tractions.
  Note this is the hardest-geometry full-recipe stress case, not a pure Mode-C
  A/B isolator (below psi_floor=0.01 the 2D floating rod is near-singular and
  detonates even with resync; at 0.01 it is Mode-D-safe and Mode-C-benign in the
  first window). The 3D extrusion is `input_rod_and_tube_3d` (periodic z).
- **Full golden compare passes** (rebuilt CPU + sm_86 CUDA binaries, all current
  with the fix): `ci_golden_compare.sh` cpu PASS; `baseline_suite.py check`
  cpu / gpu_fast / gpu_strict all 4 cases green (12/12). The 3 pre-existing cases
  reproduce their committed references with the freshly-rebuilt binaries,
  re-confirming the fix is bit-identical when off. Boundary traction on the new
  case is bit-identical across cpu/gpu_fast/gpu_strict.
- **compute-sanitizer clean**: memcheck on the rod-and-tube repro (strict sm_86
  build, full converged solve incl. the resync path) => `ERROR SUMMARY: 0 errors`.
- VERSIONS.md bumped (`gpu-v0.1.1`, patch) and this changelog entry are the
  Phase 5.4 doc deliverables.

Refuted / not committed (recorded so they are not revisited): diagonal
inflation, residual-only / stagnation-gated Newton exits, and the quick AMReX
`GMRES_MLMG` wrapper — see `benchmark/mlmg_high_contrast_20260702/NEGATIVE_RESULTS.md`.
