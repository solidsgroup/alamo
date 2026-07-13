# RESULT: fapply-kernel-surgery (task 3.2b)

Status: **edits complete, ALL ORACLE GATES PASS, verifier CONFIRMED, committed.**
Fresh-context adversarial verifier independently regenerated every hand-unrolled
term (script-generated from the accessor uid tables) — character-identical to the
diff, 2D and 3D; MulColMajor rows match the original operator* rows term-for-term;
gate logs real; diff scope clean. One verifier finding (Tier-2 deck abort) disclosed
in the gate table below and closed with a supplemental full-solve memcheck.
Work done in isolated worktree
`/home/jackplum/Projects/alamo-fapply-322b` (branch `fapply-322b`, pinned at
HEAD `d964cfab8`) because a concurrent session is live-editing the elastic
solver in the main tree.

## What changed

Two files, four bit-exact-intent edits (a,b,c,d):

### `src/Set/Matrix4_Major.H` (edits b, c)
- **(c) `operator*(Matrix4<D,Sym::Major>, Set::Matrix3)`** — replaced the branchy
  `(i,J,k,L)` accessor loop with a hand-unrolled `data[]`-indexed expression for
  each `ret(i)`, 2D and 3D (`#if AMREX_SPACEDIM`). The `(i,J,k,L)->data[]` map was
  derived mechanically from the class's own `operator()` and verified to reproduce
  the existing `operator*(Matrix4,Set::Matrix)` byte-for-byte. Accumulation order
  (J outer, k, L inner; `a*b` term order) preserved exactly.
- **(b) `MulCol` / `MulColMajor`** — new column-restricted `Matrix4*Set::Matrix`.
  Public entry is a single template `MulCol<D,S>` that dispatches (via `if constexpr
  S==Major`) to friend helpers `MulColMajor` (2D/3D, `data[]`-unrolled, computing
  only the requested column with the exact per-entry summation of the matching
  `operator*` row) and, for every other Sym, to `(a*b).col(col)` unchanged. Single
  template deliberately keeps the Major fast path out of the overload set so a
  `Matrix4<D,OtherSym>` argument never becomes a conversion candidate.
- Added the required `friend` declarations for the Matrix3 `operator*` and
  `MulColMajor` (the hand-unrolled forms now touch private `data[]`; the original
  accessor-based versions did not need friendship).

### `src/Operator/Elastic.cpp` (edits a, d, and b call-site)
- **(a) Fapply DDW hoist** — `MATRIX4 const ddw = DDW(i,j,k);` loaded once; replaced
  the four center reads (sig ~532, C·gradgradu ~613, sine-probe ~607, grad(psi)
  correction ~629). Debug-only neighbor reads left as `DDW(...)`.
- **(b) Fapply column contraction** — `(Cgrad1*gradu).col(0)+...` replaced with
  `Set::MulCol(Cgrad1,gradu,0)+...`, same vector-sum and `*psi_avg` order.
- **(d) Diagonal DDW hoist** — `MATRIX4 const ddw = DDW(i,j,k);` hoisted above the
  per-component `for(p)` loop; replaced the two reads at ~860 (boundary sig) and
  ~868 (interior).

Only these two files are modified (`git status --porcelain`: `M src/Operator/Elastic.cpp`,
`M src/Set/Matrix4_Major.H`). No unrelated files touched.

## Oracle gates — all PASS

| Gate | Command | Result |
|------|---------|--------|
| Device lint | `benchmark/lint_device_patterns.sh` | **PASS** (exit 0, non-allowlisted violations: 0) |
| CPU golden compare | `GOLDEN_MODE=cpu benchmark/ci_golden_compare.sh` | **PASS** — `canonical_step1/2`, `eta_expression_step1`, `rod_and_tube_step2` all `fresh vs reference: ok`; NaN smoke OK; `=== ci_golden_compare.sh (cpu): PASS ===` |
| compute-sanitizer memcheck | `TIERS=2 benchmark/local_a100_gate.sh` (BIN=3d profile, SAN/CUDA_HOME→main-tree .local) | **PASS** — `TIER 2 ... PASS wall=158s`; `ALL TIERS PASSED -- cleared for NOVA`. **Disclosure (verifier finding):** the Tier-2 deck (`input_3d_centre_bore_128_a2`) aborts at step 6 with `MLMG failed` — a pre-existing divergence at the pinned clean HEAD (deck is tuned for a concurrent session's in-flight solver fixes), not caused by this diff. The gate's Tier-2 criterion checks only `ERROR SUMMARY: 0 errors` and legitimately ignores the abort, so this PASS covers device-safety only for the ~6 steps that ran. |
| compute-sanitizer memcheck, supplemental full-solve | `compute-sanitizer --tool memcheck alamo_gpu-2d-profile input max_step=55` (2D converging deck, ≥1 complete elastic solve, run to normal completion) | **PASS** — run finished all 55 steps, AMReX finalized normally, `ERROR SUMMARY: 0 errors` (`gate_memcheck_2d_fullsolve.log`). Added by orchestrator to close the coverage gap above. |
| Standalone bit-exactness | nvcc host build, strict IEEE (no fast-math), 20000 randomized trials + `Increment()` signed-zero corner | **PASS** — `[DIM=2] 0 mismatches`, `[DIM=3] 0 mismatches` |

The CPU golden compare is the authoritative FP-accumulation-order oracle (plain g++,
strict IEEE); its bit-match confirms the rewritten expressions are bit-identical to the
originals on real physics. Logs in this dir: `gate_lint.log`, `gate_golden_cpu.log`,
`gate_sanitizer_tier2.log`, `gate_standalone_bitexact.txt` (+ test source
`matrix4_bitexact_test.cpp`).

## Static register / stack (cuobjdump, sm_86, 3D binary)

Chamber path is `Elastic<1>` (Sym::Major, NeoHookeanPredeformed). Full tables in
`regs_baseline.txt` / `regs_after.txt`.

| Kernel (Elastic<1>) | baseline REG / STACK | after REG / STACK |
|---------------------|----------------------|-------------------|
| Fapply (main)       | 255 / **192**        | 254 / **48**      |
| Diagonal            | 189 / 936            | 184 / 936         |

Headline: the hot `Elastic<1>::Fapply` static **stack spill drops 192→48 bytes** (reg
count was already pinned at the 255 cap → 254). Other Sym rows for reference: Fapply
`<3>` 200→176, `<4>` 122→149, `<6>` 122→129 (the non-Major rows use the unchanged
`(a*b).col()` fallback; small ±deltas are inlining side-effects of the ddw hoist, not
the chamber path). A100 occupancy/reg judgment deferred (NOVA off-limits this session).

## Wall time (local A1000, indicative only)

Deck: `input` (2D, 64², star geometry, elastic on), `max_step=251` (5 elastic solves),
`plot/thermo` off. 3 runs each. GPU idle before every run (compute-apps snapshots empty
— no contention from the concurrent session). Binary: `alamo_gpu-2d-profile-cuda86-g++`.

| Metric | baseline runs | median | after runs | median | Δ |
|--------|---------------|--------|------------|--------|---|
| Fapply excl (s) | 5.530, 5.572, 5.535 | **5.535** | 5.294, 5.320, 5.278 | **5.294** | **−4.4%** |
| TinyProfiler total (s) | 62.98, 63.96, 62.86 | 62.98 | 62.17, 63.33, 63.64 | 63.33 | ~flat (noise) |

`Fapply() NCalls = 63133` in every run, baseline and after — identical solve path, i.e.
the edits did not perturb convergence (consistent with bit-exactness).

**Clock / power caveat:** RTX A1000, 8 GB, **50 W power cap**, shared desktop. SM clock
idles at 210 MHz and boosts to ~1882–1897 MHz under load in *both* baseline and after
(matched), but the 50 W cap makes absolute wall numbers indicative only. Total wall is
dominated by non-elastic flame stepping (Fapply ~8.4% of total), so the ~4.4% Fapply
improvement does not move total wall out of run-to-run noise. **A100 wall judgment is
deferred** — NOVA was off-limits this session.

`Diagonal()` is not invoked by this deck's smoother configuration, so it produced no
TinyProfiler wall row; edit (d) is validated by the register table + golden + sanitizer.

## Deviations from PLAN.md

1. **Worktree instead of main tree.** Step 0 VERIFY (`git status --porcelain -- src/Operator
   src/Set` clean) could not be met in the main tree — a concurrent session had an
   uncommitted `src/Operator/Elastic.H` change (m_psi_small default 1E-8→0.0). Orchestrator
   provided isolated worktree `alamo-fapply-322b` at pinned HEAD; all work done there.
2. **Timing deck.** PLAN preferred `input_nova_centre_bore`; at the pinned clean HEAD its
   first elastic solve **diverges** (MLMG resid→1e20, aborts — deck is tuned for the
   concurrent session's in-progress solver fixes) and the 3D `input_3d_centre_bore_128_a2`
   **OOMs** the 8 GB A1000 during the full elastic solve. Fell back to the canonical
   golden-suite deck `input` (2D, elastic-on, converges cleanly at HEAD), which is the
   correct "smallest converging elastic deck" per the PLAN's fallback intent.
3. **compute-sanitizer paths.** The worktree has no `.local/`; pointed `SAN`/`CUDA_HOME`/`BIN`
   at the main-tree `.local` and the 3D profile binary. Gate ran clean.
4. **Standalone test build.** Compiled with nvcc (AMReX headers are GPU-configured;
   host g++ can't parse `__host__ __device__`), strict IEEE (no `--use_fast_math`) to
   match the CPU-golden semantics, linked against libamrex with three minimal `Util::`
   symbol stubs (Random/Abort/globalprefix) to avoid dragging in the IO/ParmParse graph.
   This test **caught a real bug** mid-implementation: the hand-unrolled Matrix3 operator
   accessed private `data[]` without a friend declaration (the original used the public
   accessor) — fixed before any gate ran.

## Landing status (orchestrator, 2026-07-09 evening)

Committed as **9470889b1** on branch **`fapply-322b`** (worktree
`/home/jackplum/Projects/alamo-fapply-322b`, based on d964cfab8).
**NOT cherry-picked onto chamber-gpu**: by commit time the concurrent session
had a wholesale uncommitted rewrite of `src/Operator/Elastic.cpp`
(−461/+218 over the exact Fapply/Diagonal regions this task edits), so
landing now is impossible without destroying its in-flight work.

Merge instructions once the concurrent session's Elastic work lands:
1. `git merge fapply-322b` (or `git cherry-pick 9470889b1`) on chamber-gpu.
2. `src/Set/Matrix4_Major.H` will merge clean (untouched by the other
   session) — this carries edit (c) unrolled `operator*(Matrix4,Matrix3)`
   plus `MulCol`/`MulColMajor`.
3. `src/Operator/Elastic.cpp` WILL conflict. If the conflict is gnarly,
   re-apply the three call-site intents by hand onto the new code — each is
   a 1-5 line, bit-exact change:
   (a) hoist `MATRIX4 const ddw = DDW(i,j,k);` once in the Fapply lambda,
       use for all center-node reads;
   (b) replace `(Cgrad_c * gradu).col(c)` terms with `Set::MulCol(Cgrad_c, gradu, c)`;
   (d) hoist the same `ddw` above Diagonal's per-component `p` loop.
4. Re-run the full gate set (lint, CPU golden, sanitizer) after merge —
   the other session changes solver semantics (psi regularization), so the
   combined state needs its own pass.
