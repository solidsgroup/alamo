# Plan: Full Device-Native Elastic Solve (chamber-gpu, Option B)

**Status: controlled experimental branch, NOT the mainline architecture.**
`benchmark/PHASE1_ELASTIC_DISPOSITION.md` carries the standing, final D1
verdict: **CPU-resident elastic** — a fair `np8` CPU comparison beat GPU
elastic at every resolution that converges (2.27x-3.40x), and GPU elastic
diverges entirely above 1024^2 on a non-fast-math build. The mainline
`chamber-gpu` architecture is therefore **GPU phase-field + CPU-resident (or
disabled) elastic** until this plan produces a fair benchmark that overturns
D1. Do not let this plan's progress (e.g. "reaches `STEP 60 ends`") be read as
superseding D1 — running is necessary but not sufficient.

Before this plan's elastic path can be promoted to mainline, it must clear
three hard gates, on the actual target GPU (A100, not the sm_86 dev box this
plan was drafted against):
1. **No-fast-math convergence** at high resolution (the case that previously
   diverged at 2048^2 in `PHASE1_ELASTIC_DISPOSITION.md` must converge here).
2. **A fair N-rank, fully-subscribed CPU-node solve comparison** at matched
   resolution and matched iteration count — not a single-core baseline.
3. **`ncu` occupancy/register-pressure data on A100.** The prior
   register-pressure finding (up to ~255 registers/thread on sm_86, see
   `benchmark/PHASE3_3D_READINESS.md`) was never remeasured on A100 and may
   not transfer.

> **Gate 1 is the active blocker and is being worked separately.** The 2048²
> no-fast-math divergence is under root-cause in
> `benchmark/elastic_sensitivity_20260621/GPU_ELASTIC_DEBUG_PLAN.md` ("CURRENT
> PRIORITIES"). Standing cause (2026-06-22): a deterministic GPU coarse-level
> (mglev≥1) operator-application defect — NOT fast-math, NOT bottom-solver
> acceptance policy, NOT a nondeterministic reduction (all refuted). Interim:
> `elastic.max_coarsening_level=0` already converges at 2048² on GPU (single MG
> level) and may satisfy "elastic runs on GPU" pending a wall-clock-vs-CPU number.
> This port plan ("does it RUN", de-virtualization) and the divergence hunt
> ("does it CONVERGE at scale") are distinct; this plan's Gate 1 == that hunt's
> success.

**Goal:** Make the elastic (Mechanics/Newton/Operator) solve run *on the GPU*, so a
full chamber run completes on `bin/alamo_gpu-2d-cuda86-g++` and the elastic solve
(currently ~27% of CPU wall) is accelerated — not just the Flame phase-field path.

**Current status:** Phase 2 source builds for CPU clang++ and CUDA g++.
`bin/alamo_gpu-2d-cuda86-g++` now reaches `STEP 60 ends` on `input_copy` with
`model_void.kappa=20_MPa model_void.mu=20_MPa`, passing the first elastic solve.
Correctness validation is deferred until a chamber-branch 20 MPa CPU reference is
available; `cpu_devirt_s60` is not a valid baseline because it came from
`chamber-gpu` and used 4 MPa void stiffness metadata.

---

## Root cause (confirmed, two stacked problems)

1. **Polymorphic Solid model on device.** `Model::Solid::Solid<SYM>`
   (`src/Model/Solid/Solid.H:37`) declares `virtual ~Solid`, `virtual W/DW/DDW/
   Advance/ContainsNan/Print`. The Newton device kernels (`Newton.H:167,544,587`)
   call `model(i,j,k).DW(...)`. The model object lives in device memory but its
   **vtable pointer points to host memory**; the virtual dispatch derefs a host
   pointer on the device → illegal access → `CUDA error 719`. (`origin/gpu` ships
   the same virtual models and relies on the compiler *devirtualizing* — that is
   not happening in this build, and is not guaranteed by the standard anyway.)
   Note `Solid.H:91` already warns via `std::has_virtual_destructor<T>` — the
   codebase knows virtual dtors are a GPU hazard.

2. **Polymorphic Elastic BC operator on device.** `BC::Operator::Elastic::Elastic`
   (`src/BC/Operator/Elastic/Elastic.H:18`) has `virtual ~`, `virtual operator()`,
   `virtual getType`. It is invoked as `(*m_bc)(...)` *inside a device ParallelFor*
   at `src/Operator/Elastic.cpp:364` (and host loops at 168/205). Same vtable-on-
   device hazard. This is why `gpu` demoted some elastic kernels to host loops.

A naive "demote to host loops" fix fails because the solver's **scratch fabs
(`dw_mf/ddw_mf/rhs_mf/res_mf`) are in device-only memory**, so host loops segfault
writing them, and `amrex.the_arena_is_managed=1` did not change that. Option B
sidesteps this entirely by keeping everything on the device.

---

## Strategy

De-virtualize both hierarchies so every elastic kernel is a true device kernel
operating on concrete, statically-typed objects (all elastic code is already
templated on the concrete model type `T` via `Newton<T>` / `Mechanics<T>` /
`Operator::Elastic<SYM>`, so **runtime polymorphism is not actually required** —
this must be verified in Phase 1a). No managed memory needed: data already lives
on-device (proven — the model fab read succeeds on device; only the vtable deref
failed).

## Solver safety: what has to be true for MLMG/BiCGStab/CG

The Krylov method itself is not the GPU risk. The risk is the **entire solver
call chain** beneath it. A GPU MLMG/CG/BiCGStab path is only safe if every piece
it calls is device-safe:

1. `Fapply`, `Diagonal`, residual construction, smoother steps, and BC
   application must run as device kernels or device-callable functions.
2. No virtual dispatch may occur inside a device kernel unless the object is
   proven device-safe end to end; treat `Solid::DW/DDW` and elastic BC
   `operator()` calls as hazards until de-virtualized.
3. No host-only objects may be captured by device lambdas. This includes
   parser-owned state, `std::function`, `std::string`, `std::ostream`, and any
   host-owning pointer members.
4. No host loops may write into device-arena scratch fabs. If a host fallback is
   ever used, it must be backed by explicit managed-memory policy, not accidental
   writes into device-only allocations.
5. Device abort / NaN handling must be redesigned as a device-side flag plus a
   host-side reduction and abort. Do not call host-style `Util::Abort` or stream
   printing from inside kernels.
6. Numerical safety is separate from memory safety. A device-safe solver can
   still be ill-conditioned or convergence-fragile; the no-fast-math 2048^2
   convergence test remains mandatory before any promotion.

Practical reading:

- **CG** is only acceptable if the effective operator is actually SPD under the
  current masking, BC, and material state.
- **BiCGStab** is more tolerant of nonsymmetry, but it still depends on the same
  device-safe operator and smoother chain.
- **Smoothers** must be GPU-friendly variants (Jacobi / Chebyshev / parallel
  red-black / AMReX-supported device smoothers). Sequential lexicographic
  smoothers are not a good GPU target.

The decision criterion is therefore not “can CG/BiCGStab run on GPU?” but
“does the full elastic operator chain satisfy device-safety and convergence
requirements when driven by that solver?”


Work in phases; each phase ends with a build + correctness gate.

---

## Phase 0 — Harness & guardrails (do first, ~30 min)

- [ ] Confirm clean baseline: `git status` (Newton.H must be unmodified).
- [ ] Record the CPU golden reference for correctness diffing. Full run already
      exists: `out_cpu_star_20mpa.log` (15000 steps, 1.5 s, wall 1401.7 s). Also
      keep a SHORT golden: run CPU `max_step=60` (covers 1 elastic solve) with a
      plotfile, save the plt dir + `thermo.dat` as the bit-compare target.
      ```
      Current directive: do not use `cpu_devirt_s60` for this gate. If a short
      correctness baseline is needed, generate it from a chamber branch run with
      `model_void.kappa=20_MPa model_void.mu=20_MPa`.
      ./bin/alamo-2d-clang++ input_copy max_step=60 \
        model_void.kappa=20_MPa model_void.mu=20_MPa \
        plot_file=golden_cpu_s60 amr.plot_int=60 amr.thermo.plot_int=10
      ```
- [ ] Build helper: `source benchmark/local_cuda_env.sh` then
      `./configure --comp=g++ --dim 2 --cuda 86 && make -j8`
      (incremental `make` after header edits recompiles only dependents).
- [x] Fast GPU crash-repro / pass test (step past the first elastic solve):
      ```
      ./bin/alamo_gpu-2d-cuda86-g++ input_copy max_step=60 \
        model_void.kappa=20_MPa model_void.mu=20_MPa \
        amr.plot_int=-1 amr.thermo.plot_int=-1 \
        elastic.solver.verbose=0 elastic.print_model=0 plot_file=/tmp/gpu_s60
      ```
      Success = reaches `STEP 60 ends` with exit 0.
- [ ] Analysis suite is ready at `analysis/` (`analysis/run_all.sh`) for the final
      CPU-vs-GPU comparison once the GPU run completes.

## Phase 1 — De-virtualize the Solid model hierarchy (core enabler)

**1a. Audit for real runtime polymorphism (BLOCKING — do before editing):**
- [ ] `grep -rn "Solid<" src | grep -vE "class|: public|template"` and look for any
      `Solid<...>*` / `Solid<...>&` base-pointer/reference call sites, `dynamic_cast`,
      or containers of base pointers. Expectation: none in the solver hot path
      (everything is `Set::Field<T>` + templated operators). If found, those sites
      must be templated/specialized instead of relying on virtual dispatch.
- [ ] Check `src/Test/**` and `src/Model/Solid/Solid.H:91` (the
      `has_virtual_destructor` self-test) — update/keep as appropriate.

**1b. Remove virtual from the base (`src/Model/Solid/Solid.H`):**
- [ ] Drop `virtual` from `~Solid`, `W`, `DW`, `DDW`, `Advance`, `ContainsNan`,
      `Print`. Keep them `AMREX_GPU_HOST_DEVICE`. Make `~Solid` a defaulted
      non-virtual dtor (or rely on trivial). Keep the non-pure default bodies.

**1c. Remove `override`, drop `virtual ~`, optionally mark `final`, in EVERY model:**
- [ ] Files: `src/Model/Solid/Finite/NeoHookean.H`,
      `NeoHookeanPredeformed.H`, and ALL other `src/Model/Solid/**` models
      (Linear/*, Affine/*, etc. — enumerate with
      `grep -rln "virtual\|override" src/Model/Solid`).
- [ ] For each: remove `virtual` and `override` keywords on W/DW/DDW/Print/dtor;
      keep `AMREX_GPU_HOST_DEVICE`. Consider `class X final : public ...`.
- [ ] Watch `Print(std::ostream&)`: keep it guarded with `AMREX_IF_ON_HOST(...)`
      (already is in NeoHookeanPredeformed) since `std::ostream` is host-only.

**1d. Device-safe the in-kernel aborts (gpu idiom):**
- [ ] Guard host-only `Util::Abort`/`contains_nan()`-abort calls that sit *inside*
      device `ParallelFor` with `#ifndef ALAMO_GPU` (ALAMO_GPU is `-D`-defined by
      `configure` for `--cuda`). Sites incl. `Operator/Elastic.cpp:77,734,843` and
      any model `ContainsNan`. Mirror how `Flame.cpp` already does this.

**1e. Gate:**
- [ ] CPU build (`./configure --comp=clang++ --dim 2 && make`) compiles; a
      `max_step=60` CPU run is **bit-identical** to the golden (de-virtualization
      must not change numerics — it's pure dispatch).
- [ ] GPU build compiles; the Phase-0 GPU step-60 test now **passes the elastic
      solve** (the Newton model kernels at 167/544/587 run on device). The BC
      kernels are still host (Phase 2), so this should already complete.

## Phase 2 — Device-native the Elastic BC operator

The BC is only evaluated on domain-boundary nodes, but `gpu` left it on the host.
For full device-native, pick ONE approach (recommend 2a):

**2a. (Recommended) De-virtualize + capture-by-value the concrete BC.**
- [ ] Remove `virtual`/`override` from `BC::Operator::Elastic::Elastic` base and
      `Constant.H`/`Expression.H` (`operator()`, `getType`, dtor); keep
      `AMREX_GPU_HOST_DEVICE`. Make the BC a small POD-ish object (no heap-owning
      members; check `Constant`'s members are trivially device-copyable —
      `m_bc_*` arrays of scalars/Type enums are fine; flag any std::string/std::function).
- [ ] Where the operator holds `m_bc` as a base pointer, template
      `Operator::Elastic<SYM>` (or just the kernels) on the concrete BC type, OR
      store the concrete BC by value and capture it into the kernel `[=]`.
      Convert the host loops at `Operator/Elastic.cpp:168` (Stress) and `:516`,
      and any Newton BC loop, to device `ParallelFor` capturing the BC by value.
- [ ] `Expression.H` likely holds a parsed-expression object — verify it is
      device-constructible; if not, keep Expression-BC on host and only run
      Constant-BC on device (document the limitation).

**2b. (Fallback) Keep BC on host, models on device (hybrid).**
- Accept BC boundary loops on host *with managed scratch fabs* (allocate the BC
  rhs/stress scratch in managed memory). Smaller, but not "full device-native".
  Only use if 2a's BC object proves not device-copyable.

**Gate:**
- [ ] GPU full 1.5 s run completes, exit 0.
- [ ] GPU vs CPU correctness: `thermo.dat` and field tree within tolerance
      (de-virtualization is exact; expect bitwise or ~1e-12 from FP reassociation
      under `--use_fast_math`). If `--use_fast_math` causes drift, build a
      `--profile` (no fast-math) variant for the correctness check.

## Phase 3 — Benchmark & report

- [ ] Run `analysis/run_all.sh` (CPU log `out_cpu_star_20mpa.log` already exists;
      produce the GPU log `out_gpu_star_20mpa.log` from a full 1.5 s GPU run).
- [ ] Headline: GPU vs CPU wall clock, and specifically the elastic-solve cost now
      on device (parse MLMG timers; `analysis/lib/parse_log.py`).
- [ ] Note: 64^3 base / 3 levels is small — if GPU still trails, scale
      `amr.n_cell` / `max_level` to show the crossover (see port memory notes).

---

## Risk register / gotchas

- **Hidden runtime polymorphism** (Phase 1a) is the only thing that can make
  de-virtualization wrong. If a base `Solid*`/`Elastic*` is dispatched at runtime
  anywhere, that path must be templated instead. Audit hard before editing.
- **Solver safety is end-to-end, not algorithm-only.** BiCGStab / CG / MLMG are
  safe on GPU only if their operator, smoother, BC, scratch-fab, and abort/NaN
  paths are all device-safe. The solver name alone does not make the path safe.
- **`Expression` BC / models** may carry host-only members (`std::function`,
  parsed AST, `std::string`) — not device-copyable. Constant BC (used by
  `input_copy`) is the priority; gate Expression behind host or make it POD.
- **In-kernel `Util::Abort`/`Util::Message`/`contains_nan`** must be
  `#ifndef ALAMO_GPU`-guarded or they break the device compile/run.
- **`--use_fast_math`** can perturb the last bits → use a `--profile` build for
  the bit-compare gate; use fast-math build for timing.
- **Don't add managed memory** for Option B — it's unnecessary and slow; the whole
  point is keeping everything device-resident.
- Keep CPU build bit-identical after every phase (it's the regression oracle):
  `./configure --comp=clang++ --dim 2 && make` then diff a `max_step=60` run vs
  `golden_cpu_s60`.

## Key file map

| Concern | File(s) |
|---|---|
| Solid base virtuals | `src/Model/Solid/Solid.H` (37–93) |
| Models to de-virtualize | `src/Model/Solid/**` (Finite/NeoHookean*, Linear/*, Affine/*, …) |
| Model operator macros | `src/Model/Solid/{In,Ext}ClassOperators.H` |
| Newton device kernels (model calls) | `src/Solver/Nonlocal/Newton.H` (167, 544, 587) |
| Elastic operator kernels + BC call | `src/Operator/Elastic.cpp` (110,133,168,205,325-364,451,516,580,685,796,873) |
| BC operator virtuals | `src/BC/Operator/Elastic/{Elastic,Constant,Expression}.H` |
| Mechanics integration / fab alloc | `src/Integrator/Base/Mechanics.H` (71–102, 207) |
| Build | `configure`, `benchmark/local_cuda_env.sh`, `benchmark/build_alamo_local_gpu.sh` |
| Comparison harness | `analysis/run_all.sh` |

## Definition of done

- [ ] `bin/alamo_gpu-2d-cuda86-g++` completes a full 1.5 s `input_copy` (star geom,
      void=20 MPa) run, exit 0, no host loops in the elastic solve.
- [ ] GPU result matches CPU golden within tolerance.
- [ ] `analysis/run_all.sh` produces the CPU-vs-GPU report with the elastic solve
      timed on-device.
