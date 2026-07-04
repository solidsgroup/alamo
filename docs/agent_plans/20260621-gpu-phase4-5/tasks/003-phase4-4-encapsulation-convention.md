# Task 003: Phase 4.4 — encapsulation / device-capture convention

## Goal
Document the standard pattern that replaced ad-hoc shadowing in `Flame.cpp`:
"copy POD params into a local, capture by value (never `this`) in device lambdas."
Make it a referenceable convention so future GPU-ported integrators follow it
instead of making integrator internals `public` case by case.

## Context
nvcc rejects `this`-capture in `__device__` lambdas. The branch's idiom: name the
integrator's anonymous param structs (`pf_struct pf;` etc.), make them capturable,
and inside each `ParallelFor` do `auto pf = this->pf;` (shadow into locals) so the
lambda captures by value with no `this`. See `src/Integrator/Flame.cpp` (the
`auto X = this->X;` shadowing before `amrex::ParallelFor`) and `Flame.H` (named,
public param structs). An existing doc `docs/gpu_device_capture_conventions.md`
starts this — EXTEND it into the full convention; do not create a duplicate.

## Files to read first
- `docs/gpu_device_capture_conventions.md` (current content)
- `src/Integrator/Flame.cpp` (grep `auto .* = this->`; the ReduceOps lambda)
- `src/Integrator/Flame.H` (the named `*_struct` param structs)

## Files allowed to modify
- `docs/gpu_device_capture_conventions.md` ONLY (extend it).

## Files NOT allowed to modify
- Everything else (esp. do NOT touch `Flame.H`/`Flame.cpp` — doc only).

## Implementation steps
1. Read the current doc + the Flame.cpp/Flame.H patterns.
2. Extend the doc to cover, with short real code excerpts from Flame:
   - The rule: device lambdas capture by value, never `this`.
   - The mechanism: named param structs (POD) + `auto x = this->x;` shadow.
   - When you still need `public` (only to let a captured-by-value struct be read),
     and the preferred alternative (shadow local) so internals stay encapsulated.
   - A short checklist for porting a new integrator's hot loop.
3. Keep it concise and practical (a convention, not an essay).

## Invariants
Documentation only — describe existing code, do not change behavior.

## Expected result
`docs/gpu_device_capture_conventions.md` is a complete, reusable convention.

## Non-goals
Editing source; new files.

## Stop conditions
None expected.

## Final report: write results/003-RESULT.md
