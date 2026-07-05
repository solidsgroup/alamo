# Physics-Quantity Registry & Error Budget — Roadmap Task 1.A

This is the contract that defines "correct enough" for chamber-gpu (originally
GPU Roadmap v3, `docs/archive/GPU_ROADMAP_V3.md`, historical; current status
`docs/llm/PLAN.md`). Every observable the validation suite tracks is
listed here with its extraction source and the tolerance class it must satisfy.
The machine-readable mirror is `physics_budget.yaml`; `compare_validation.py`
(task 1.D) reads the YAML, this file is the human-reviewable rationale.

No optimization (input lever or kernel rewrite) ships unless it passes every
CORRECTNESS observable and stays inside every ENGINEERING-TRAJECTORY tolerance
below. SOLVER-HEALTH observables are tracked for regression visibility but do
not gate a verdict on their own.

## How sources were determined

Every observable below is tied to a real, currently-registered field or log
line in this checkout, confirmed by grep against `src/`:

- `thermo.dat` columns come from `RegisterIntegratedVariable` calls in
  `src/Integrator/Flame.cpp` (chamber quantities) and
  `src/Integrator/Base/Mechanics.H` (corner displacement).
- Plotfile fields come from `RegisterNodalFab`/`RegisterGeneralFab` calls in
  `src/Integrator/Mechanics.H` (`eta`), `src/Integrator/Flame.cpp` (`phi`), and
  `src/Integrator/Base/Mechanics.H` (`disp`, `stress`, `strain`, `model`).
- Solver-health lines come from AMReX's own verbose output: `MLMG: Iteration N
  Fine resid/<bnorm|resid0> = ...` (internal MLMG `verbose>=2`) / `MLMG: Final
  Iter. N resid, resid/<bnorm|resid0> = <abs>, <rel>` (internal `verbose>=1`) in
  `ext/amrex/Src/LinearSolvers/MLMG/AMReX_MLMG.H:508,532`, bottom-solver
  `MLCGSolver_BiCGStab: Final: Iteration N rel. err. X` in
  `AMReX_MLCGSolver.H:260` (needs `mlmg.setBottomVerbose()>0`), and `NR
  iteration N ... relnorm = ...` / `Newton Iteration N of M` in
  `src/Solver/Nonlocal/Newton.H:382,455`. **Important indirection:**
  `Linear.H:268-272` (`PrepareMLMG`) maps the single input knob
  `elastic.solver.verbose` to these AMReX internals as `mlmg.setVerbose(verbose
  - 1)`, and only calls `mlmg.setBottomVerbose(verbose)` when `verbose > 4` —
  there is **no separate `bottom_verbose` key**. So `elastic.solver.verbose=5`
  is the single setting that satisfies all three thresholds (per-V-cycle lines
  need raw `verbose>=3`, the `Final Iter.` line needs raw `verbose>=2`, the
  bottom-solver line needs raw `verbose>4`), captured in `run.log`.

No observable below is aspirational — every one already exists as data the
binary produces today.

---

## Tolerance class definitions

### CORRECTNESS
Must match the CPU-strict (no-fast-math) golden to **~1e-6 relative**. Evaluated
once, on the **final** plotfile/checkpoint of a short, deterministic run (small
`max_step`, no stochastic regrid timing). Any CORRECTNESS fail ⇒ overall FAIL,
unconditionally — these are the observables a register-shaved kernel or a
loosened solve tolerance could silently corrupt while still "running fine."

### ENGINEERING-TRAJECTORY
Must match within a stated **physical %** over the **whole run** (not just the
endpoint), via two metrics: (a) max relative deviation of the time series, and
(b) a phase-lag / cross-correlation check, so a solve that converges to the
right final state but drifts in time (a stale-but-converged field) cannot pass
on its endpoint alone. A fail here is surfaced with magnitude, not an automatic
hard-fail of the whole budget — magnitude and which observable failed decide
whether the change is acceptable (see `compare_validation.py` roll-up logic).

### SOLVER-HEALTH
Regression-tracked, **non-gating**. A SOLVER-HEALTH delta does not fail the
budget by itself, but a large jump (e.g. Newton iterations doubling, MLMG
V-cycles departing from the historical count, bottom-solver iterations
exploding) is the first symptom of a stale or under-resolved solve — Phase 4.F
ties these into the gate as a guard condition once that task lands.

---

## CORRECTNESS observables

| Name | Source | Extraction | Threshold |
| --- | --- | --- | --- |
| `eta_field_l2` | plotfile `eta` (node, `Mechanics.H:136`) | per-component L2 norm, final plotfile, GPU-strict vs CPU-strict | rel ≤ 1e-6 |
| `eta_field_linf` | plotfile `eta` | per-component L∞ norm, final plotfile | rel ≤ 1e-6 |
| `phi_field_l2` | plotfile `phi` (node, `Flame.cpp:250`) | per-component L2 norm, final plotfile | rel ≤ 1e-6 |
| `temp_field_l2` | plotfile `temp` (`RegisterNewFab`, `Flame.cpp:180`) | per-component L2 norm, final plotfile | rel ≤ 1e-6 |
| `disp_field_l2` | plotfile `disp` (vector, `Base/Mechanics.H:71`) | per-component (x/y[/z]) L2 norm, final plotfile | rel ≤ 1e-6 |
| `disp_field_linf` | plotfile `disp` | per-component L∞ norm, final plotfile | rel ≤ 1e-6 |
| `stress_field_l2` | plotfile `stress` (symmetric matrix, `Base/Mechanics.H:73`) | per-component (xx/yy[/zz]/xy[/xz/yz]) L2 norm, final plotfile | rel ≤ 1e-6 |
| `stress_field_linf` | plotfile `stress` | per-component L∞ norm, final plotfile | rel ≤ 1e-6 |
| `strain_field_l2` | plotfile `strain` (`Base/Mechanics.H:74`) | per-component L2 norm, final plotfile | rel ≤ 1e-6 |
| `elastic_residual_at_convergence` | `run.log`, `MLMG: Final Iter. N resid, resid/<bnorm\|resid0> = <abs>, <rel>` | parse the **relative** residual (`<rel>`, second value after `=`) — not the absolute one: it scales with the problem's physical units (Pa) and stays large (observed 0.125 on a real converged local run, vs `<rel>`≈5e-9) even when AMReX itself reports convergence, since the solve's actual gate is `resid ≤ max(tol_rel·bnorm, tol_abs)` | ≤ 1e-8 (matches `elastic.tol_rel` default, `Base/Mechanics.H:142`) |

CORRECTNESS is evaluated **only** on strict (no-fast-math) builds, per v3
guiding principle 3 — fast-math binaries are never used to make a correctness
claim, only a perf-labelled smoke row.

## ENGINEERING-TRAJECTORY observables

| Name | Source | Extraction | Threshold |
| --- | --- | --- | --- |
| `chamber_pressure` | `thermo.dat` column `chamber_pressure` (`Flame.cpp:224`) | full time series | max rel dev ≤ 2%; phase-lag ≤ 1 output interval |
| `chamber_volume` | `thermo.dat` column `volume` (`Flame.cpp:189`) | full time series | max rel dev ≤ 2% |
| `burn_area` | `thermo.dat` column `area` (`Flame.cpp:190`) | full time series | max rel dev ≤ 2% |
| `total_mdot` | `thermo.dat` column `mass_flux` (`Flame.cpp:191`) | full time series | max rel dev ≤ 2%; phase-lag ≤ 1 output interval |
| `mdot_max` | `thermo.dat` column `mdot_max` (`Flame.cpp:193`) | full time series | max rel dev ≤ 3% |
| `max_temp` | `thermo.dat` column `max_temp` (`Flame.cpp:192`) | full time series | max rel dev ≤ 2% |
| `interface_extent` | `thermo.dat` column `eta_min` (`Flame.cpp:196`) proxy for burn-front progress | full time series | max abs dev ≤ 0.02 (eta is in [0,1]) |
| `corner_displacement` | `thermo.dat` columns `disp_xhi_x/y`, `disp_yhi_x/y` (`Base/Mechanics.H:114-117`) | full time series | max rel dev ≤ 3% |
| `max_von_mises_stress` | derived from plotfile `stress` components per case dimensionality | computed per snapshot: `sqrt(0.5*((sxx-syy)^2+(syy-szz)^2+(szz-sxx)^2+6*(sxy^2+syz^2+szx^2)))` (2D: `szz=syz=szx=0`); max over domain, tracked per output | max rel dev ≤ 5% |
| `mean_von_mises_stress` | same derivation as above | domain mean per output | max rel dev ≤ 5% |
| `max_principal_stress` | derived from plotfile `stress` | largest eigenvalue of the stress tensor per node, max over domain per output | max rel dev ≤ 5% |
| `elastic_energy` | derived from plotfile `stress`/`strain` | `0.5 * sum(stress : strain) * cell_volume` per output | max rel dev ≤ 3% |

Von Mises / principal stress / elastic energy are **derived**, not raw
registered fields — `field_norms.json` (bundle schema, task 1.C) stores the raw
per-component stress/strain norms; the derived scalars above are computed by
the comparator from those norms at compare time, not by the simulation binary.

`max rel dev` = `max_t |candidate(t) - reference(t)| / max(|reference(t)|, eps)`
over the matched/interpolated time series. Phase-lag is reported wherever
chamber pressure or mdot exhibits a recognizable transient (the ignition
spike); it is a cross-correlation-offset metric, not used where the trajectory
is monotone (no extremum to lag).

## SOLVER-HEALTH observables (regression-tracked, non-gating)

| Name | Source | Extraction |
| --- | --- | --- |
| `newton_iters_per_solve` | `run.log`, count of `NR iteration N ...` lines between successive `Newton Iteration` blocks (`Newton.H:382,455`) | per elastic solve, requires `verbose>0`; `Newton: public Linear` (`Newton.H:27`) so this is the *same* `elastic.solver.verbose` field, already satisfied by `=5` |
| `newton_final_relnorm` | `run.log`, trailing `relative norm(ddisp)` value on the last `NR iteration` line per solve | per elastic solve |
| `mlmg_vcycles_per_solve` | `run.log`, count/max-index of `MLMG: Iteration N ...` lines, terminal `MLMG: Final Iter. N` (`AMReX_MLMG.H:508,532`) | per elastic solve |
| `bottom_solver_iters` | `run.log`, `MLCGSolver_BiCGStab: Final: Iteration N` (`AMReX_MLCGSolver.H:260`), requires `elastic.solver.verbose=5` (see indirection note above) | per V-cycle |
| `residual_monotonic` | `run.log`, sequence of `MLMG: Iteration N Fine resid/...` residual values | boolean: non-increasing across the V-cycle sequence, per solve |

These require `elastic.solver.verbose=5` (`Linear.H:315`, see indirection note
above — this single knob is the only way to get the bottom-solver line at all)
in the validation case overrides — the case manifest (task 1.B) sets this for
every case so `run.log` always carries the lines above.

---

## Review checklist (the "Done-when" for 1.A)

Before this budget is treated as final, it must be reviewed against the
chamber model's engineering intent — does a PASS on every observable above
actually mean the burn/pressure/stress prediction the chamber model exists to
produce is preserved? Specifically:

- [ ] Are the 2% / 3% / 5% ENGINEERING thresholds tight enough to catch a
      change that would alter a chamber-design decision, and loose enough that
      legitimate floating-point/reduction-order noise (GPU non-associative
      sums) doesn't generate false FAILs? *(needs a domain-expert sanity pass —
      not a number I can derive from source alone.)*
- [ ] Is von Mises / principal stress the right structural-margin proxy for
      this propellant grain's failure mode, or should burn-front curvature /
      local stress concentration at the bore also be tracked?
- [ ] Is 1-output-interval phase-lag tolerance appropriate, or too coarse to
      catch a meaningfully shifted ignition transient?

These three are flagged for human review rather than decided unilaterally —
the thresholds above are a defensible starting draft grounded in existing
`abs_tol`/`rel_tol` conventions already used by `baseline_suite.py`
(1e-8 abs / 1e-6 rel for CORRECTNESS-class thermo comparisons), not yet
validated against real GPU-vs-CPU noise floors. Task 1.E (golden references)
will surface the *actual* current GPU-vs-CPU deltas, which should be used to
sanity-check these thresholds before 1.F makes them gating.
