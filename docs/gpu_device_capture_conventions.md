# GPU Device Capture Conventions

Device kernels must capture a self-contained, device-valid snapshot of the data
they need. Do not rely on implicit `this` capture or host-owned pointers inside
`AMREX_GPU_DEVICE` lambdas.

## Required Pattern

Before `amrex::ParallelFor`, copy every member, parameter bundle, geometry value,
and model object used by the lambda into a named local value:

```cpp
Set::Vector DX(geom[lev].CellSize());
const auto dx = DX.data();
const auto phase = phasefield;
const auto propellant_model = propellant.model;

amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k) {
    // Use dx, phase, and propellant_model here. Do not use this->... here.
});
```

This makes captures explicit, keeps the CUDA object closure small, and avoids
capturing host-only storage such as geometry-owned arrays.

## Rules

1. Capture only named local values in device lambdas. Avoid member access,
   implicit `this`, references to parser-owned state, and virtual dispatch.
2. Convert host-owned array pointers into owning value types before capture.
   Examples: `Set::Vector(geom[lev].CellSize())` for cell size and
   `amrex::GpuArray` for small fixed arrays.
3. Keep captured structs POD-like and publicly named. CUDA rejects some extended
   lambdas that capture private, protected, anonymous, or host-only types.
4. For complex objects, add an explicit device-copyable view/accessor rather than
   exposing internals ad hoc.
5. Keep comments local and short. Point unusual cases back to this document
   instead of restating the convention at every kernel.

`Integrator::Flame` currently follows this convention for chamber CUDA kernels.
New GPU-clean integrators should match it before being added to
`src/GPU/IntegratorPolicy.mk`.

## The Encapsulation Convention: Named Param Structs + Local Shadowing

The pattern above ("copy every member ... into a named local value") has one
recurring wrinkle: `nvcc` rejects an extended `__device__` lambda defined inside
a private/protected member function, and it also rejects a captured variable
whose *type* is private, protected, or unnamed. A lambda capturing
`this->elastic` by value is only legal if `elastic`'s type and the enclosing
method are both visible to the lambda's instantiation context. The fix used
throughout `Flame` is two-part: (1) group each subsystem's POD parameters into
a named, `public` struct, and (2) immediately before the kernel, shadow each
struct member you need into a local of the same name via `auto x = this->x;`.
The lambda then captures the *local*, never `this`, and never needs to know
the struct came from a class member at all.

`src/Integrator/Flame.H` documents why these members are `public` instead of
`protected` (the conventional default for an `Integrator` subclass):

```cpp
// NOTE: members below are public (rather than protected) because their
// override methods define extended __device__ lambdas, and the parameter
// structs are captured by value into device lambdas. nvcc forbids both an
// extended device lambda inside a private/protected member function and a
// captured variable whose type is a private/protected (or unnamed) member.
// See docs/gpu_device_capture_conventions.md for the capture convention.
public:
```

Below that comment, each subsystem gets one named struct instance, e.g.:

```cpp
// Phase field model variables
struct pf_struct {
    Set::Scalar eps = NAN;
    Set::Scalar lambda = NAN;
    Set::Scalar kappa = NAN;
    Set::Scalar w1 = NAN, w12 = NAN, w0 = NAN;
    Set::Scalar min_eta = 0.001;
    int relax_steps = 0;
} pf;

// Thermoelastic variables
struct elastic_struct {
    int on = 0;
    Set::Scalar Telastic = NAN;
    model_type model_prop, model_void, model_casing;
    Set::Scalar traction = NAN;
    int traction_from_chamber = 0;
    ...
} elastic;
```

`thermal_struct` and `chamber_struct` follow the same shape. Each struct is
POD (or POD-of-PODs, e.g. `model_type`): plain scalars, small fixed structs,
and no pointers to host-only storage, no virtual methods, no `Set::Field`
(field handles stay outside the struct and are turned into `Set::Patch`
views per-`MFIter` as usual).

### The shadow line

Inside `Flame.cpp`, right before the loop that builds device lambdas, the
struct is shadowed into a local of the same name:

```cpp
// GPU: copy members into locals so the model-build kernels below
// capture by value instead of the host `this` pointer.
auto       elastic     = this->elastic;
const bool homogenized = this->homogenized;
```

and, in `Advance`, several subsystems at once:

```cpp
// GPU: device lambdas may not capture the host `this` pointer. Copy the POD
// phase-field params, the propellant model, and the `small` floor into locals
// (shadowing the members), and pull the thermal scalars used inside kernels
// into named locals so the kernels below capture by value, not via `this`.
auto              propellant      = this->propellant;
auto              pf              = this->pf;
const Set::Scalar small           = this->small;
const bool        thermal_on      = thermal.on;
const Set::Scalar thermal_hc      = thermal.hc;
const Set::Scalar thermal_Tcutoff = thermal.Tcutoff;
const Set::Scalar thermal_Tfluid  = thermal.Tfluid;
```

and, for the pre-relax kernel in `Initialize`:

```cpp
// GPU: copy POD members into locals so the relax kernel captures by value.
auto pf = this->pf;
const Set::Scalar small = this->small;
```

Note the pattern in the `Advance` example: rather than shadowing the whole
`thermal_struct`, only the scalars actually used inside the kernel
(`thermal_on`, `thermal_hc`, `thermal_Tcutoff`, `thermal_Tfluid`) are pulled
out as individually named locals. Shadow the whole struct when the kernel
uses several of its fields (cheap, POD, one line); pull out individual
scalars when only one or two fields are needed, or when the struct mixes
device-safe scalars with host-only members (e.g. `thermal_struct` also holds
`Set::Field` and `IC::IC<Set::Scalar>*`, which must never be captured).

The variable name does not change across the shadow — `auto pf = this->pf;`
keeps every later use of `pf.eps`, `pf.kappa`, etc. inside the kernel
unchanged from a textual standpoint. Only `this->` resolution is removed at
the `pf` name's declaration point; everywhere else in the lambda body refers
to the new local. This makes the diff between "broken on device" and
"works on device" a one-line addition, not a rewrite of the kernel body.

### When you still need `public`, and the preferred alternative

`public` is required for a member only when a lambda will capture that
member's *struct type* by value (the nvcc restriction above bites on the
type, not just on `this`-access). It is not a blanket license to make
everything public:

- If a member is a POD struct that some kernel captures via shadowing
  (`pf`, `thermal`, `elastic`, `chamber` above), it must stay `public`, with
  a comment pointing back to this doc explaining why.
- If a member is **not** captured into any device lambda (host-only logic,
  bookkeeping, anything not read inside `[=] AMREX_GPU_DEVICE(...)`), prefer
  keeping it `protected`/private and reading it through a normal accessor or
  local copy on the host side. Do not widen visibility preemptively.
- Prefer the shadow-local idiom over capturing `this` and reaching through it
  inside the lambda (`[this](...){ ... this->elastic.on ... }`) even when it
  would compile on the host — it will not compile under `nvcc` for an
  extended device lambda, and it reintroduces the exact problem this
  convention exists to avoid.
- For something more complex than a POD struct (e.g. a model with internal
  invariants you don't want copied/mutated freely), prefer adding a small,
  explicit device-copyable view/accessor type over making the original type
  or its members public (see Rule 4 in "Rules" above).

### Checklist for porting a new integrator's hot loop

1. Identify every member the kernel body reads: scalars, POD structs,
   model objects, geometry-derived values (`DX`, domain boxes), and any
   `Set::Patch` views.
2. Group the kernel's POD inputs into one or a few named structs
   (`foo_struct foo;`) if they aren't already grouped; keep `Set::Field`
   handles and pointer/IC members out of those structs.
3. Mark the struct(s) `public` only if a lambda will capture the struct
   itself by value; add the same `// NOTE: ... public ...` comment pattern
   pointing back to this doc.
4. Immediately before the `MFIter`/`ParallelFor` block, shadow each needed
   struct (or, for mixed structs, each needed scalar field) into a local of
   the same name: `auto foo = this->foo;`.
5. Convert any host-owned array pointers (cell size, domain box) into owning
   value types before the shadow block, per Rule 2 above.
6. Confirm the lambda's capture clause is `[=]` (by value) and that no
   `this->` or bare member name resolving through `this` remains inside the
   lambda body — only the shadowed locals.
7. Build with the CUDA POSTFIX and confirm there is no "calling a
   `__host__` function from a `__device__` function" or extended-lambda
   visibility error before moving on.
