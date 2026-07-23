# Checkpoint before the frozen-mechanics test hook

## Outcome so far

The existing synchronized final output strongly supports H0, but Step 1 cannot
yet pass because the normal restart/evolve path does not preserve the exact
constitutive and load state used by the saved mechanics solve.

No production source was changed by this investigation.

## Pressure-flux and corrected-traction evidence

The analysis oracle recovered the lagged mechanics pressure as
`4,316,386.252 Pa`, consistent with the original log value of approximately
`4,316,390 Pa`.  The 2-D `CellGradientOnNode` stencil was represented exactly
as the divergence of transverse eta-face averages; the maximum identity error
was `5.12e-13 1/m`.

Across both circular interfaces on both symmetry axes:

| Metric | Raw elastic traction | Pressure-corrected traction |
|---|---:|---:|
| Nyquist projection | 18.53-19.08 kPa | 0.548-0.902 kPa |
| Detrended RMS | 223.7-247.4 kPa | 10.5-29.2 kPa |
| Historical max second difference | 146.5-146.9 kPa | 13.2-24.2 kPa |

Thus adding the exact pressure flux removes 95.3-97.1% of the alternating
component measured in raw elastic traction.  This is the predicted signature
when raw stress contains the diffuse transition needed to balance
`-p grad(eta)`.  It is provisional because the frozen repeat, planar control,
and refinement checks have not yet run.

The RHS fit to the final cell eta has relative Linf error `4.58e-4`.  This is
consistent with the known time ordering: mechanics used the state at
`t=1.0002`, while the final cell plot was written after the phase/thermal
advance at `t=1.0004`.

## Why the input-only frozen restart failed

The initial dual-centering restart output proves:

- `disp_x`, `disp_y`, `rhs`, stress, strain, phi, casing support, and the
  cell-centered eta are restored bitwise on all three levels;
- the `model_mu`, `model_kappa`, and `model_F0*` fields are all zero after
  restart, despite their names being matched; and
- `Set::Field<NeoHookeanPredeformed>` implements plot serialization but no
  `CopyFrom`, so restart deserialization is a no-op.

Normal `TimeStepBegin` then zeros RHS and calls `UpdateModel`, rebuilding RHS
and F0 from the later plotted eta/temperature rather than retaining the state
that produced the saved displacement.  In the serial diagnostic run the
nonlinear residual remained `1.46e12 Pa/m` after three damped iterations,
instead of the saved converged value `2.32284e6 Pa/m`; that invalid run was
stopped.

The eight-rank dual restart also segfaulted while rebuilding nodal fields.  The
serial path avoids that independent restart/distribution issue, so the first
valid frozen check should remain serial.

## Proposed experimental hook

Add one narrowly scoped, default-off diagnostic option:

`elastic.frozen_restart_fields=1`

Its behavior would be limited to a dual-centering restart:

1. deserialize all six `NeoHookeanPredeformed` components (`mu`, `kappa`, and
   `F0`) from the nodal plot;
2. retain the restored RHS and model instead of zeroing/recomputing them in
   `TimeStepBegin`;
3. set the boundary-condition time, solve once, compute the composite residual,
   and use the already configured zero-mobility/zero-energy advance solely to
   reach the normal plot write; and
4. assert that eta, temperature, phi, casing support, RHS, model, hierarchy,
   and pressure are unchanged before and after the solve/write cycle.

The discriminating result is precise: starting from the restored displacement,
Newton should reproduce the saved `2.32284e6 Pa/m` residual scale with a zero or
near-zero update, two identical serial runs should agree bitwise, and the
pressure-corrected Nyquist feature should remain below 1 kPa.  A materially
different result rejects the current H0 evidence and returns the investigation
to the source/model coupling hypotheses.

The hook will be kept as a separate experimental diff.  After the frozen,
planar, and refinement experiments, it will be reverted unless its model
deserialization portion is independently selected and tested as a production
restart fix.  No operator, stencil, material law, or production default will
change in this experiment.

## Reproduction

```bash
python3 docs/agent_plans/20260722-eta-ripple-rootcause/analyze_eta_ripple.py --self-test
python3 docs/agent_plans/20260722-eta-ripple-rootcause/analyze_eta_ripple.py \
  --case frozen-baseline --check-pressure-flux-identity
NPROCS=1 docs/agent_plans/20260722-eta-ripple-rootcause/run_frozen_baseline.sh frozen-single-c
```

The third command is the deliberately stopped, invalid pre-hook reproduction;
it must not be used as a physical result.
