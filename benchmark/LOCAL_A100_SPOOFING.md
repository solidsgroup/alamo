# Spoofing the NOVA A100 environment on the local RTX A1000

Goal: catch **A100/NOVA-class GPU defects locally** before spending NOVA time, on the
chamber-gpu high-performance test branch.

## The host-pointer-on-device bug class (instances found + fixed)

Every NOVA error-700 on this branch so far is the same class: a **host pointer
dereferenced inside a device kernel**. The local A1000 hides it (HMM migrates the page,
see below); the A100 faults. The compute-sanitizer gate is HMM-immune and catches all of
them. Known instances, all fixed on `chamber-gpu`:

| Site | What was wrong | Fix | Commit |
|---|---|---|---|
| `Base::Mechanics::Integrate` | host members `trac_hi`/`disp_hi` (array decays to host ptr) written in a device `ParallelFor` | `amrex::ReduceOps` sum-reduction, add to host members afterward | `76ae550e0` |
| `Solver::Nonlocal::Newton` (`prepareForSolve`, `rhs`, `DW`) | `const Real* dx = Geom().CellSize()` host pointer dereferenced in the elastic device kernels | copy into `Set::Vector DX(...CellSize())` (by-value, device-resident) and pass `DX.data()` | `54a941433` |
| `Operator::Elastic` (`Fapply`, `Diagonal`, `Strain`, `Stress`, `Energy`) | same `CellSize()` host-pointer deref across the MLMG operator kernels | same `Set::Vector DX` + `DX.data()` idiom (precedent: `Mechanics.H`, `TopOp.H`) | `54a941433` |

The `Set::Vector DX(ptr)` idiom works because the small Eigen vector is captured **by
value** into the `[=]` device closure, so `DX.data()` points at the device-resident copy
rather than the host `Geometry`'s cell-size array. These conversions are value-preserving
(CPU output unchanged) and were verified by the gate: memcheck = 0 errors through the
multi-box elastic solve, which exercises exactly these kernels.

Two companion changes shipped alongside: `src/alamo_gpu.cc` drops an obsolete diagnostic
probe (an env-gated CUDA device-stack bump that was chasing the now-explained
`Newton.H:196` "unspecified launch failure"); and the CUDA build's integrator closure
was refactored out of the `Makefile` into `src/GPU/IntegratorPolicy.mk` + a
`configure --gpu-integrator` flag (currently `flame` only) -- behavior-preserving (the
`flame` source list is identical to the old hardcoded one), just explicit and extensible.

## The problem we are solving

The `Base::Mechanics::Integrate` traction-diagnostic bug (a device kernel writing to
host members `trac_hi`/`disp_hi`) **hard-faulted on NOVA's A100** with
`CUDA error 700 (illegal memory access)` but **ran clean on the local A1000**. The
divergence is not random and not power/clock related. Root cause, confirmed locally:

```
$ nvidia-smi -q | grep 'Addressing Mode'                      ->  HMM
$ cat /sys/module/nvidia_uvm/parameters/uvm_disable_hmm       ->  N   (HMM enabled)
```

The A1000's driver has **HMM (Heterogeneous Memory Management)** enabled. A device
dereference of a *system-allocated* host pointer is transparently serviced by page
migration instead of faulting. The NOVA A100 stack does not migrate it, so the same
access is an illegal access. **Anything that relies on the hardware to fault will be
fooled by HMM locally** -- including ordinary runs and even AMReX `--debug` post-kernel
error checks (no CUDA error is ever raised, because from the driver's view the access
"succeeded").

## What can and cannot be spoofed

| Lever | Spoofs the A100? | Notes |
|---|---|---|
| **compute-sanitizer** (memcheck/racecheck/initcheck) | **Yes (best)** | HMM-immune: validates each access against the *intended* allocation, independent of migration. Catches the error-700 class today's bug belongs to. No root needed. Cost: 10-50x slowdown. |
| **Disable HMM** (`modprobe nvidia_uvm uvm_disable_hmm=1`) | **Yes** | Makes the A1000 hard-fault like the A100, full speed. Needs **root** + module reload, and **disrupts any GPU/display session** on the device. Gated behind `ALLOW_HMM_DISABLE=1` in the gate script. |
| **Runtime flags** (`CUDA_LAUNCH_BLOCKING=1`, `amrex.abort_on_out_of_gpu_memory=1`) | Partial | Launch-blocking reports async faults at the exact kernel; abort-on-OOM is a clean guard. Does **not** defeat HMM for plain host members. Cheap; always worth running. NB: pure device arena (`the_arena_is_managed=0`) is impractical here -- it OOMs reserving its init block on the shared 8 GB desktop, and the chamber **elastic solve only fits at all via managed/HMM page migration** -- so the gate keeps managed memory. |
| **AMReX `--debug` build** (`./configure --debug --cuda 86`) | Partial | Adds `Array4` bounds asserts + `AMREX_GPU_ERROR_CHECK` after every launch. Catches device-arena OOB and post-kernel CUDA errors immediately, faster than sanitizer -- but **not** the HMM-masked host-pointer class. Separate postfix (`3d-debug-cuda86-g++`), needs a full AMReX-debug rebuild. |
| **GPU architecture** (sm_80 vs sm_86) | No | The A1000 cannot execute sm_80-only code; this is a minor correctness factor. Build local sm_86; treat arch as out of scope. |
| **CUDA toolkit version** | Mostly N/A | Local build toolkit is CUDA 12.6 (`.local/cuda`); sanitizer is 12.0 (`.local/cuda-12.0-ubuntu`). Matches modern NOVA closely enough for correctness work. |
| **Device memory size** (8 GB vs 40/80 GB) | No -- shrink instead | Cannot enlarge local VRAM. Reproduce NOVA's **box decomposition** at smaller total size: use `input_3d_centre_bore_128_a2` + `amr.max_grid_size=64` to get a multi-box (2x2x1) layout like NOVA's 256^3 deck, which is what exercises the cross-box GPU paths. |
| **ECC / power cap / clocks** | Irrelevant | Affect perf and bit-error rates, not logic-bug reproduction. |

## The gate: `benchmark/local_a100_gate.sh`

A tiered, HMM-immune pre-NOVA correctness gate built around compute-sanitizer.

```
bash benchmark/local_a100_gate.sh                 # default: tiers 1 2 (smoke + memcheck)
TIERS="1 2 3" bash benchmark/local_a100_gate.sh   # thorough: + racecheck/initcheck (slow)
ALLOW_HMM_DISABLE=1 TIERS=4 bash benchmark/local_a100_gate.sh   # root: real hard-fault spoof
```

- **Tier 1** runtime-strict smoke (launch-blocking + abort-on-OOM, managed arena, stops
  before the elastic solve). Fast (~2 s). Surfaces any async fault at the exact kernel.
- **Tier 2** compute-sanitizer **memcheck** on a multi-box layout -- the load-bearing
  A100 spoof; this is the tier that would have caught the traction-diagnostic bug.
- **Tier 3** **racecheck** (non-atomic global `+=` race class) + **initcheck**
  (uninitialized device reads).
- **Tier 4** optional root-gated HMM-disabled hard run (no sanitizer overhead).

The gate runs ~6 steps of the 128^3 chamber deck with tiny elastic-solver iteration
counts (`nriters=1, max_iter=2`): the goal is to exercise every kernel (the traction
diagnostic runs every step; one elastic solve at step 5), not to converge the solve.

### Measured costs on this A1000 (128^3, mgs=64), all PASSING on the fixed binary
- Tier 1 runtime-strict smoke (5 steps, pre-elastic): **~2 s**.
- Tier 2 memcheck (6 steps incl. one multi-box elastic solve): **~150 s**, `ERROR SUMMARY: 0 errors`.
- Tier 3 initcheck (3 steps): **~33 s**. racecheck: the slowest tool by far (many minutes for
  even 3 steps -- Eigen-dense kernels are expensive to instrument); reserve for pre-NOVA only.

## Recommended workflow

1. Every source change touching device kernels: **Tier 1 + Tier 2** (smoke + memcheck).
2. Before a NOVA run: full `local_a100_gate.sh` (all sanitizer tiers).
3. The standing rule from this episode: **memory model is a performance/convenience
   lever, never a correctness substitute.** HMM/managed memory hide host/device pointer
   bugs; keep scalars in `amrex::ReduceOps` reductions and validate with sanitizer.
