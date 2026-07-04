# Result: Task 003 — Phase 4.4 encapsulation/device-capture convention doc

## Status
DONE

## What changed
Extended `docs/gpu_device_capture_conventions.md` (sole file modified) with a
new section, "The Encapsulation Convention: Named Param Structs + Local
Shadowing", appended after the existing "Required Pattern" / "Rules" content.
No other files were read for editing purposes beyond reference, and none were
modified: `src/Integrator/Flame.H` and `src/Integrator/Flame.cpp` were only
grepped/read to source real excerpts.

The new section covers, with short verbatim excerpts pulled from the current
working tree:

1. **Why `public`**: quotes the existing `Flame.H` comment block (lines
   ~50-55) explaining nvcc's restriction (extended `__device__` lambdas
   cannot live in private/protected member functions, and cannot capture a
   variable whose type is private/protected/unnamed).
2. **Named param structs**: excerpts of `pf_struct` and `elastic_struct` from
   `Flame.H` (lines 107-153) showing the POD-only, publicly-named shape.
3. **The shadow line**: three real excerpts of `auto X = this->X;` from
   `Flame.cpp` (lines 329-330 in `Initialize`'s pre-relax kernel, lines
   432-433 in the elastic/model-build loop, and lines 629-635 in `Advance`,
   which mixes whole-struct shadowing (`pf`, `propellant`) with individual
   scalar shadowing (`thermal_on`, `thermal_hc`, `thermal_Tcutoff`,
   `thermal_Tfluid`) — used to explain when to shadow a whole struct vs. just
   the fields a kernel needs.
4. **When you still need `public`, and the preferred alternative**: clarifies
   `public` is required only for struct *types* captured by value into a
   device lambda, not a blanket rule; non-captured members should stay
   `protected`; explicitly calls out that capturing `this` directly inside an
   extended device lambda (`[this](...){...}`) is the anti-pattern this
   convention replaces, and that complex non-POD members should get an
   explicit device-copyable view/accessor instead of being made public.
5. **Checklist for porting a new integrator's hot loop**: 7 concrete steps
   (identify member reads -> group into POD structs -> mark public only if
   captured -> shadow into locals before the kernel -> convert host-owned
   pointers to value types -> verify `[=]` capture and no remaining `this->`
   -> build under the CUDA POSTFIX and check for extended-lambda visibility
   errors).

## Files touched
- `docs/gpu_device_capture_conventions.md` (extended only; was already
  untracked/new in git, no prior history to preserve).

## Files explicitly NOT touched
- `src/Integrator/Flame.H`, `src/Integrator/Flame.cpp`, `Makefile`,
  `configure` — read-only for sourcing excerpts. Confirmed via
  `git diff --stat` before/after that none of these changed as a result of
  this task (their pre-existing uncommitted diffs from before this task
  started are untouched).

## Verification
- `git status --porcelain docs/gpu_device_capture_conventions.md` shows the
  file as the only artifact from this task.
- `git diff --stat -- src/Integrator/Flame.H src/Integrator/Flame.cpp
  Makefile configure` shows only the pre-existing (pre-task) diffs, unchanged
  by this task.
- No commit was made (per instructions).

## Open follow-ups / risks
- None. The doc describes existing, already-working code; no behavior was
  changed. Future GPU-ported integrators should follow the checklist before
  being added to `src/GPU/IntegratorPolicy.mk` (per the doc's closing line,
  unchanged from the original).
