# Feature and numerical evidence

Classification: chamber-gpu corpus evidence. This material illustrates why a
decision must be surfaced; the normative workflow is
`../ARCHITECTURE_POLICIES.md#feature-and-num-surfacing`.

- The high-contrast elasticity, AMR recovery, traction, and Newton records describe device-independent solver behavior. They support `[NUM]` entries in `FEATURES.md`, never GPU conversion Transforms.
- Relevant sources include `benchmark/MLMG_HIGH_CONTRAST_FINDINGS.md`, `benchmark/mlmg_high_contrast_20260702/`, `docs/agent_plans/20260709-elastic-void-robustness/`, and `docs/agent_plans/20260713-elastic-regression-recovery/`.
- Retired NUM-001 through NUM-004 cover conservative face flux, coarse-fine policies, component layouts required by that formulation, and damped Newton updates. Each requires explicit task-level opt-in.
- Golden tolerances, convergence criteria, residual definitions, and AMR policies are configuration/physics evidence. A port worker must surface them to the user rather than importing them from chamber-gpu.
