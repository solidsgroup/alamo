# GPU manual hostile-review result

Outcome: the five live blocking defects are repaired in documentation and
recognizer tooling. No Hydro, Fracture, or other source file was changed, and
no claim is made that either integrator is ported, correct, or performant.

## Delivered

1. Recognizer schema v3 removes the named Flame/Elastic candidate literals,
   uses identifier-independent high-recall shapes, and leaves semantic facts to
   compiler/call-graph/ledger confirmation. Production-table fixtures exercise
   GPU-002/003/004/013/018/019/022/030; GPU-031 remains deliberately manual.
2. Static root `docs/gpu_manual/COVERAGE.csv` is retired. Each scan requires a
   `port_id` and immutable `source_revision`; reviewed dispositions carry only
   when both match. A revision change reopens candidates.
3. GPU-002 now covers device-reached host polymorphism and container-backed
   state generally, with a self-contained retained-buffer/value-view/static
   tuple/switch example.
4. `GPU_NATIVE_SHAPE.md` is a mandatory pre-baseline contract for field layout,
   component order, kernel/MFIter graph, transfers, numerical call chains,
   compiler resources, divergence, shared memory, atomics/reductions, and
   arena/allocation lifetime. Three templates make the evidence instantiable.
5. GPU-031 treats iterative, branch-heavy, host-oriented numerical algorithms
   as explicit sub-ports. Current Hydro/Riemann and Fracture sites are recorded
   as unported candidate evidence, not recipes.
6. GPU-007 and GPU-016 are now `file-verified`, not transfer-verified. Promotion
   requires a frozen unseen target; repair/retest makes the target authoring
   evidence and cannot establish transfer.

The old review's statements about INDEX-driven static coverage and
physics-specific Verify commands were already stale after generalization and
were not reintroduced. The historical 103-FEATURE denominator issue is retained
as open blind spot BS-008 pending a bounded independent re-cluster.

## Evidence

- Permanent manual validator: 26 patterns; 2 file-verified, 0
  transfer-verified, 0 cross-family.
- Scanner tests: 10/10 pass, including revision-change reopening and generic
  production-rule fixtures.
- Fresh full-source advisory scan: 1,593 rows; exact GPU-013 sites include
  Hydro.cpp 466–469 and Fracture.H 696; regenerated output byte-matches the
  recorded revision-bound artifact.
- `git diff --check`: pass.
- Source-write guard: no path under `src/` changed.
- Independent hostile reviewer: initial evidence-scope/carryover findings fixed;
  recheck found no remaining high/medium defect.
- Project status at closeout: device-lint PASS. No build, sanitizer, numerical,
  or GPU performance result is claimed by this docs/tooling task.

## Remaining gates

The 1,593 audit rows are intentionally undispositioned task evidence, not a
completed-port worklist. A real port must scan its declared closure, disposition
every row, execute validation and GPU-native-shape/baseline contracts on named
hardware, run a frozen unseen-target gate, and harvest changes back into the
manual. Those pilot-dependent statuses remain open.
