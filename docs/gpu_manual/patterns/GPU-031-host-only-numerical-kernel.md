# GPU-031: Port host-only numerical kernels as explicit device algorithms
Transform status: draft
Class: correctness
Detection: Advisory compiler diagnostic plus call-graph and inspection worklist; regex cannot prove transitive callability, container ownership, iteration bounds, or device reachability.
Invariant: A device launch may reach only device-callable code and state; numerical termination, error, and ownership semantics must remain explicit.
Port contract: Inventory mathematical inputs/outputs, conserved or bounded quantities, iteration/convergence rules, branch map, temporaries, dispatch alternatives, container storage, transitive callees, and CPU behavior. Complete `GPU_NATIVE_SHAPE.md` for kernel boundary, divergence, and resource evidence.
Transform:
  Before: A device lambda calls a host-oriented numerical solver/model with virtual dispatch, host containers, exceptions/diagnostics, or unbounded/hidden iteration.
  After: Treat the algorithm as a sub-port: stage device-safe values/views, use static selection, make the full call chain callable, bound or justify iteration, propagate errors by policy, and choose launch boundaries from recorded resource evidence.
Corpus example: The non-chamber current-target Hydro Riemann and Fracture crack paths in `evidence/host-only-numerical-kernels.md` are unported corpus candidates, not successful recipes.
Constraints: Annotation alone is insufficient. Preserve the numerical algorithm unless FEATURE/[NUM] is approved; do not flatten branches, cap iterations, or replace containers without matching validation evidence.
Verify: Instantiate `VALIDATION.md`; require `strict-build`, `analytic-exact`, `conservation`, `golden-regression`, `multi-box`, and `sanitizer` rows, plus `GPU_NATIVE_SHAPE=pass` and a written tolerance rationale.
Failure modes: Partial annotation leaves deeper host calls; copied host containers retain invalid addresses; virtual dispatch, exception paths, data-dependent loops, divergence, or register spills can fail compilation, corrupt results, or produce a nominally GPU kernel that is slower than the CPU path.
Evidence: Primary: `evidence/primary-sources.md` (device callability and workload-shape anchors). Corpus: current-target inspection at commit 13342c7cdf3141d86011fa672ee26c658fb23a7f; `evidence/host-only-numerical-kernels.md`.
