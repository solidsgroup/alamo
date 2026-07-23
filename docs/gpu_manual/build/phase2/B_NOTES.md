# Partition B refinement notes

- `B_HUNKS.csv` has 218 rows: 216 distinct base contextual parents plus four explicit subentries (`src/Operator/Elastic.cpp:173-264#1/#2` and `src/Operator/Operator.cpp:368-380#1/#2`). Parent denominator remains 216; cluster counts deduplicate parent IDs.
- Broad rejected clusters were split into recognizable transforms: host/device annotations, include-only oneoffs, host-loop guards/fallbacks, static BC dispatch, device-safe vector capture, device error flag/API, matrix unrolling, component layout/interpolation, elixir lifetime, and narrow operator/Newton policies.
- `ELASTIC_COMPONENT_LAYOUT` and `ELASTIC_COMPONENT_INTERPOLATION` are repeated component-shape transforms; `OPERATOR_GHOST_POLICY` is repeated ghost-row handling. Composite norm, coefficient resync, quadratic interpolation, diagnostics, reflux, line-search configuration, and other unique behavior are ONEOFF.
- Parent 62 contains capture plumbing (#1) and elastic dispatch behavior (#2). Parent 100 contains diagnostics (#1) and coarse-ghost sanitization (#2).
- Distinct-parent classifications: CANDIDATE 158, ONEOFF 40, NOISE 18. No base parent spans multiple classifications; two candidate parents span multiple candidate clusters via their explicit subentries.
- Include-only compile fixes are ONEOFF; comments/blank/formatting remain NOISE. No markdown or Tier 1 files were touched.
