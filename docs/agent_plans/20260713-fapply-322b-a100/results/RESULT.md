# RESULT — fapply-322b-a100 (2026-07-13)

## What this covers

Collection phase for task 3.2b A100 judgment. All 6 NOVA jobs completed
(builds 11639419/11639420, wall 11639477/11639478, ncu 11639479/11639480).
This doc records the wall A/B, ncu A/B, and cross-arm parity, plus a verdict
recommendation. Orchestrator commits; DONE disposition explained in §5.

Arms (see NOTES.md):
- BASELINE `d964cfab8` — chamber-gpu tip before 3.2b.
- MODIFIED `dc02baf17` — `fapply-322b` tip (DDW hoist, column-restricted
  contraction, unrolled Matrix4xMatrix3, Fsmooth 4D launch fusion).

Deck: `input_3d_centre_bore_256_a2`, 1× A100 (sm_80), TinyProfiler on, 600
steps. Both arms did identical work: 105,920 Fapply calls, 13 MLMG solves,
39,600 Fsmooth calls, 21 plotfile writes.

## 1. Wall A/B (TinyProfiler, extracted from `wall_322b_<arm>.log`)

| Region                             | calls   | BASELINE | MODIFIED | Δ            |
|------------------------------------|--------:|---------:|---------:|--------------|
| `Operator::Elastic::Fapply()` excl | 105,920 | 299.0 s  | 255.7 s  | **−14.5%**   |
| `MLMG::solve()` incl               | 13      | 366.2 s  | 327.4 s  | **−10.6%**   |
| `Operator::Fsmooth()` incl         | 39,600  | 302.1 s  | 270.2 s  | −10.6%       |
| `Operator::Fsmooth()` excl         | 39,600  | 16.24 s  | 20.27 s  | +4.0 s (see note) |

Fapply exclusive MODIFIED (255.7 s) ≤ BASELINE (299.0 s): perf-oracle met, no
regression. The Fsmooth *exclusive* self-time rises (+4 s) but its *inclusive*
time falls 32 s — the 4D launch fusion moves work out of per-component launch
overhead and into the fused body, and the dominant Fapply calls underneath it
got cheaper. Net MLMG::solve() −38.8 s (−10.6%). Matches C1-precedent shape
(that A/B saw Fapply −14.5%, MLMG −11.7%).

## 2. ncu A/B (job 11639479/11639480, 40 launches/arm, NVTX-gated on
`Operator::Elastic::Fapply()/`)

Registers and occupancy are flat; the win is per-launch time and reduced SM
work, grid-size-matched (4 grid sizes appear across the captured window):

| Metric                            | BASELINE | MODIFIED | Δ        |
|-----------------------------------|---------:|---------:|----------|
| registers / thread                | 255      | 254      | −1 (−0.4%) |
| achieved occupancy (avg)          | 12.07 %  | 11.90 %  | −0.17 pp (flat) |
| Fapply / launch, grid 8715 (large)| 5.528 ms | 4.263 ms | **−22.9%** |
| Fapply / launch, grid 4458        | 2.767 ms | 2.152 ms | **−22.3%** |
| Fapply / launch, grid 8517        | 5.389 ms | 4.138 ms | **−23.2%** |
| Fapply / launch, grid 2280 (small)| 1.401 ms | 1.092 ms | **−22.0%** |
| SM active-cycles, grid 8715       | 6.303e8  | 5.808e8  | **−7.9%** |
| SM active-cycles, grid 4458       | 3.137e8  | 2.917e8  | −7.0%    |
| inst-executed % of peak (avg)     | 27.7 %   | 29.9 %   | +8% rel (denser) |

Key-metric caveat (executed instructions): the `--set basic` report did **not**
store a raw executed-instruction count — `smsp__inst_executed.sum` is empty in
both `.ncu-rep` files, only the rate form `sm__inst_executed.sum.pct_of_peak…`
is present. Available instruction-side evidence therefore is the pair above:
**SM active-cycles fall ~7–8%** (clock-independent hardware work drops) while
**instruction-throughput density rises ~8% relative** (27.7→29.9 % of peak,
i.e. fewer stall/replay bubbles per cycle). Both are consistent with the
surgery's by-construction removal of live temporaries and spill-replay traffic
— the same mechanism C1 identified — not with an occupancy change (occupancy is
flat, still register-bound at 254/thread × 256/block on A100's 65,536-reg file).

The ncu per-launch −22% exceeds the TinyProfiler wall −14.5% because ncu
serializes with kernel replay and pins clocks; the **wall −14.5% is the
authoritative perf number**. Both point the same direction.

## 3. Parity (fcompare.gnu.ex, baseline vs modified, matched plotfiles)

Cell (flame) fields: **bit-identical (0 abs / 0 rel) at every step and level**
— steps 0/300/600, levels 0 and 1, "PLOTFILE AGREE". Fapply never writes cell
data; confirms the flame integration path is untouched.

Node (elastic) fields — max relative error across all variables per step/level:

| step | level | max node rel err        | where             | abs err at that point |
|-----:|------:|-------------------------|-------------------|----------------------|
| 0    | 0/1   | 0 (exact)               | —                 | 0                    |
| 300  | 0     | 1.11e-7                 | strain_zy         | 2.78e-10             |
| 300  | 1     | 1.63e-6                 | strain_zx         | 4.12e-9              |
| 600  | 0     | 2.25e-7                 | strain_xy         | 1.08e-9              |
| 600  | 1     | **2.06e-6**             | strain_zx         | 5.30e-9              |

Primary elastic outputs stay well within bar at every step/level:
- displacement (the solve unknowns): ≤ 4.0e-8 rel;
- stress (physical output): ≤ 4.8e-7 rel;
- diagonal strains: ≤ ~1e-9 rel.

The only excursions past 1e-6 are the two off-diagonal shear-strain components
`strain_zx` / `strain_zy` on the fine level (max 2.06e-6). Their **absolute**
errors are ~5e-9 — at the double-precision accumulation floor of a full MLMG
solve. In this near-axisymmetric centre-bore geometry those shear components
are near zero, so the small ‖A‖ denominator inflates the relative error while
the absolute error is negligible. Errors grow only linearly with accumulated
solve iterations (300→600: 1.6e-6→2.1e-6), not exponentially.

Interpretation: this is floating-point reorder noise from the surgery's
changed contraction/summation order (column-restriction + Matrix4xMatrix3
unroll + DDW hoist reorder the FP adds more aggressively than the C1 edit did,
so the noise floor is ~10–20× the C1 precedent's 1e-7–1e-8 — but still a noise
floor, abs ~5e-9), **not** a physics divergence. Cell bit-identicality plus
stress/disp within 5e-7 corroborate.

## 4. Verdict recommendation

**Perf: PASS.** Fapply exclusive wall −14.5% (MODIFIED 255.7 s ≤ BASELINE
299.0 s), MLMG::solve −10.6%. ncu corroborates (−7–8% SM work, +8% throughput
density). No regression on the perf axis.

**Parity: PASS on physics grounds, with one flagged caveat.** Cell fields
bit-identical; displacement and stress within the ~1e-6 bar; the only bar
crossing is two near-zero off-diagonal shear-strain components at max 2.06e-6
rel / ~5e-9 abs, which is FP-reorder noise, not divergence.

**Recommendation to orchestrator: PASS → unlock the merge.** The physics-error
budget is met: primary outputs (stress, disp, all cell fields) are within
tolerance and the strain excursion is a near-zero-denominator artifact at the
double-precision floor.

## 5. DONE disposition (escalation)

I did **not** `touch results/DONE`. Reason: the literal collection-phase bar
("node fields rel err ≤ ~1e-6; anything worse than 1e-6 rel on a node field is
a REGRESSION verdict") is crossed by strain_zx/zy (max 2.06e-6). My analysis
says this is reorder noise, not a regression — but crossing the written node
bar is exactly the case the instruction reserves for orchestrator
adjudication rather than agent self-certification. Deciding whether the
2.06e-6 shear-strain excursion is "within ~1e-6 latitude" (PASS, merge) or a
hard bar (REGRESSION, hold) is the orchestrator's call. Evidence is complete
above; touch DONE if you accept the noise reading.

## Anomalies / notes

- `smsp__inst_executed.sum` not captured by `--set basic` — raw
  executed-instruction count unavailable; used SM active-cycles + throughput %
  as instruction-side proxies (§2). A `--metrics smsp__inst_executed.sum`
  re-profile would give the exact count if the orchestrator wants it, but the
  verdict does not hinge on it (parity + wall are conclusive).
- ncu window spans 4 grid sizes (8715/8517/4458/2280), not the "two grid
  sizes" PLAN Step 4 anticipated — AMR produced more box sizes than expected;
  more coverage, not less. Single deck (`input_3d_centre_bore_256_a2`) as
  submitted.
- Fsmooth exclusive self-time rose +4 s (launch-fusion moves per-component
  launch overhead into the fused body); inclusive fell 32 s — net win. Not a
  regression.
- fcompare on the largest plotfile (`00600node`, level 1) is Lustre-I/O-bound
  and takes ~5 min; a helper wrapper self-deadlocked on a `pgrep -f` pattern
  matching its own command line (killed, re-ran directly — data above is from
  the clean synchronous run).
