# GPU Elastic Solver: Break-vs-Sensitivity Study

Date: 2026-06-21. Branch `chamber-gpu`. Local RTX A1000 (sm_86).

## Question
Does the GPU build (with/without fast-math) **break** the elastic MLMG solve, or
just make it **more sensitive**? Do `use_psi=1`/`psi_floor`, void stiffness, or
timestep fix the high-resolution divergence?

## Method
- Harness: `run_case.sh` — one elastic solve (`max_step=1`, `elastic.interval=1`,
  `elastic.tstart=0`, `solver.verbose=4`), star geometry, base = `input_base`
  (copy of `input_phase1_elastic`: `use_psi=1`, `psi_floor=0.05`, void κ=μ=4 MPa,
  prop κ=162/μ=113.6 MPa, stiff casing 5GPa/2GPa, all-clamped BCs).
- Binaries (Jun 20, NOT rebuilt — isolated from concurrent agent builds):
  - fast:    `bin/alamo_gpu-2d-cuda86-g++`        md5 63e65383...
  - no-fast: `bin/alamo_gpu-2d-nofast-cuda86-g++` md5 71e44c51... (`--fmad=false`)
  - CPU:     `bin/alamo-2d-g++` (np8)             md5 1ed0cec3...
- No source edits, no rebuild. All I/O confined to this directory.

## Results

### Resolution parity (baseline conditioning)
| Finest | CPU np8 | GPU no-fast | GPU fast |
|---|---|---|---|
| 512²  | conv 22+22 (Phase-1) | conv 22+22 → 8.381e-9 | conv 22+22 → 8.381e-9 |
| 1024² | conv (Phase-1)       | conv 16+16 → 4.655e-9 | conv 16+16 → 4.655e-9 |
| 2048² | **conv 13+13 → 4.04e-9** | **DIVERGE → 1.07e20** | **DIVERGE → 1.54e20** |

- fast vs no-fast agree to ~7 sig figs where they converge → **fast-math is not the cause**.
- fast-math buys only ~6–10% on the MLMG Solve timer (512²: 3.105→2.825 s; 1024²:
  3.497→3.293 s) — solve is launch/latency-bound, little arithmetic for it to speed up.

### Conditioning levers at 2048² (GPU no-fast) — NONE fix it
| Lever | Result (abort iters → resid/bnorm) |
|---|---|
| baseline psi=0.05, void=4 MPa | DIVERGE 270 → 1.07e20 |
| psi_floor=0.1  | DIVERGE 673 → 2.78e20 |
| psi_floor=0.2  | DIVERGE 312 → 1.27e20 |
| psi_floor=0.3  | DIVERGE 139 → 1.37e20 |
| void=40 MPa (10×, contrast ~4) | DIVERGE 134 → 1.03e20 |
| void=162 MPa (=prop, ~uniform stiffness, contrast ~1) | DIVERGE 540 → 1.32e20 |
| psi=0.2 + void=162 (both maxed) | DIVERGE 246 → 1.21e20 |

psi_floor/void only change *how long it limps before exploding*, never *whether*.
Even a near-uniform (best-conditioned) operator diverges → **not a conditioning problem.**

### Timestep at 2048² (GPU no-fast) — no effect
| dt | Result |
|---|---|
| 1e-4 | DIVERGE 270 → 1.07e20 |
| 1e-5 | DIVERGE 279 → 1.02e20 |
| 1e-6 | DIVERGE 200 → 1.38e20 |

All three have **bit-identical `resid0 = 33682118.18`** — the t=0 IC solve does not
depend on dt (no time has elapsed), so dt cannot be a lever for this failure. The
run never survives the first (cold) solve, so warm-start/dt robustness never engages.

### Determinism test at 2048² (identical command ×2) — ROOT CAUSE
| run | resid0 | iter 1..5 Fine resid/bnorm | abort |
|---|---|---|---|
| rep1 | 33682118.18 | 1.030, 1.288, 2.341, 1.267, 0.125 | 131 → 2.02e20 |
| rep2 | 33682118.18 | 0.254, 0.172, 0.108, 0.196, 0.146 | 187 → 1.02e20 |

Byte-identical input, identical `resid0`, **divergent V-cycle trajectories from
iteration 1** and different abort points → the GPU MLMG V-cycle is
**nondeterministic** (order-dependent reduction/atomic accumulation; intra-GPU,
np=1, so not MPI). On CPU the same problem is deterministic and converges 13+13.

## Conclusion
- **Does the GPU break the solver? No.** Both fast and no-fast converge *identically
  to CPU* (matching iteration counts, ~7 sig figs) at ≤1024². The operator is
  correct and even makes initial progress at 2048².
- **Is it "more sensitive"? Yes, but the mechanism is GPU reduction nondeterminism,
  not input conditioning.** The nondeterministic V-cycle is damped at ≤1024² and
  amplified by the marginally-stable V-cycle at 2048² → divergence. That is why no
  *input* lever (psi_floor, void stiffness, dt) can fix it.
- **fast-math is irrelevant** to correctness here (only ~6–10% perf at low res).
- The CPU "pickiness" the user has seen separately is a genuine conditioning/geometry
  issue (e.g. centre-bore, high void contrast) and is *distinct* from this GPU-only
  high-res nondeterministic divergence.

## UPDATE — conditioning RULED OUT; trigger is resolution/concurrency
The operator stacks two contrast sources: 3-material blend (casing 5 GPa = 31× prop,
1250× void; Flame.cpp:501-506) AND psi weight (psi_floor→1). So earlier "void=162"
runs were never uniform. Decisive test — a **genuinely uniform operator**
(`use_psi=0` + void=casing=prop, zero contrast):

| Uniform operator | GPU no-fast | CPU |
|---|---|---|
| 512²  | conv 9 it | conv 9 it |
| 1024² | conv 9 it | conv 9 it |
| 2048² | **DIVERGE → 1.53e20** | **conv 9 it** |

GPU @2048² divergence shape (uniform): contracts cleanly 0.066→0.028→0.013 for 3
iters (just like CPU), then abruptly reverses → 1.5e20. A problem the CPU solves in 9
iterations. Conditioning cannot be the cause; the threshold is set by resolution / MG
depth (8 levels @2048 vs 7 @1024), i.e. amount of GPU concurrency, not difficulty.
Pure FP reduction-reordering (~1e-9) is too small for an O(1) blowup → genuine
concurrency defect (race/sync/atomic) in the GPU nodal MLMG kernels.

**Two distinct stability boundaries:** (1) a conditioning boundary hitting BOTH
CPU+GPU at extreme contrast (user: CPU diverges for psi_floor ≪ 0.05 or void ≪ ~0.5
MPa — psi_floor/void are the right levers here); (2) the GPU-only resolution/
concurrency boundary that kills 2048² regardless of conditioning. Different problems.

## Implications / next steps
- Confirms the D1 decision (elastic CPU-resident) for the right reason.
- A real GPU fix would target deterministic reductions in the nodal MLMG smoother/
  residual (AMReX deterministic-reduce paths), not input tuning. Worth checking
  whether the nodal Gauss-Seidel/residual uses atomic accumulation across blocks.
- Practical mitigation if device-elastic is ever needed: cap finest resolution at
  1024² (last resolution that converges), or run elastic on CPU.
