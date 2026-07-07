# Phase C1 — Elastic `Fapply` register pressure / occupancy (DONE 2026-07-07)

**Branch / workspace:** `chamber-gpu-elastic-opt` (git worktree at
`/home/jackplum/Projects/alamo-elastic-opt`, forked from `chamber-gpu` tip
`7e972f1e8`, carrying the same pending working-tree GPU fixes the main checkout
had). Kept isolated so other agents can keep working the main `chamber-gpu`
checkout.

**Goal (roadmap C1 / `PHASE_A_FINDINGS.md` §8 lever A1):** raise the ~12.5%
achieved occupancy of `Operator::Elastic<SYM>::Fapply`, which Phase A measured at
**74.6% of all GPU kernel time** and **255 registers/thread** (the CUDA hard cap —
the compiler is spilling; the build already passes `-maxrregcount=255` globally).

---

## What changed in source so far (`src/Operator/Elastic.cpp`, `Fapply`)

Two **bit-identical, monotonic** edits to the hot per-node lambda. Neither reorders
any floating-point operation, so the operator output is unchanged to the last bit;
they only remove work and shrink live register state.

### Change 1 — accumulate `grad(C):grad(u)` one direction at a time (the A1b lever)

The `!m_uniform` branch previously named three live `Matrix4` derivative temps
(`Cgrad1/2/3`), whose lexical scope spanned the whole summation expression:

```cpp
MATRIX4 Cgrad1 = Stencil<…1,0,0>::D(DDW,…),
        Cgrad2 = Stencil<…0,1,0>::D(DDW,…),
        Cgrad3 = Stencil<…0,0,1>::D(DDW,…);
f += ((Cgrad1*gradu).col(0) + (Cgrad2*gradu).col(1) + (Cgrad3*gradu).col(2)) * psi_avg;
```

In 3D a `Matrix4<3,Major>` is **45 doubles**, so three live temps = **135 live
doubles (~1.1 KB)** — `PHASE_A_FINDINGS.md` §4 names this the dominant Fapply
register-spill source. The new form holds **one** `Matrix4` temp at a time (each
stencil result is consumed and dies at its statement boundary):

```cpp
Set::Vector graddc = Set::Vector::Zero();
graddc += (Stencil<…1,0,0>::D(DDW,…) * gradu).col(0);
#if AMREX_SPACEDIM > 1
graddc += (Stencil<…0,1,0>::D(DDW,…) * gradu).col(1);
#endif
#if AMREX_SPACEDIM > 2
graddc += (Stencil<…0,0,1>::D(DDW,…) * gradu).col(2);
#endif
f += graddc * psi_avg;
```

Peak live derivative state drops from 135 → ~45 doubles. Summation order is
unchanged (`0 + t0 + t1 + t2`, then `× psi_avg`), so `f` is bit-identical.

### Change 2 — sink the boundary-only stress tensor into the boundary branch

`sig = (DDW(i,j,k) * gradu) * psi_avg` was computed for **every** node but is only
read by the domain-boundary BC evaluation. Interior nodes — the vast majority, and
the 1–5 ms fine-level applies that are ~90% of Fapply time — paid a full
`Matrix4·Matrix` product (9 outputs × up to 9 terms) and a `Set::Matrix` (9 live
doubles) for a value they discard. Moved the computation inside `if (boundary)`.
Bit-identical for the boundary nodes that consume it; pure removal everywhere else.

---

## Verification done locally

- **Compiles for the real target.** Single-TU build of the modified `Elastic.cpp`,
  GPU 3D, sm_86, fast-math, `-maxrregcount=255`, exact production nvcc flags →
  **exit 0**, device object emitted, only the pre-existing host/device
  redeclaration warnings (`#20040-D`). Command archived in this branch's session
  notes; reproduce by compiling against
  `ext/AMReX-Codes/amrex/3d-cuda86-g++-26.06-dirty/include`.
- **Correctness argument is by construction** (no FP reordering), so CPU and GPU
  outputs are unchanged. **Now confirmed end-to-end on CPU — see next section.**

## ✅ CPU golden compare — DONE (2026-06-28, bit-identical)

The by-construction claim is now **verified end-to-end** with a held-everything-else-
identical A/B in this worktree:

- **Method.** Built the worktree CPU binary `bin/alamo-2d-g++` (2D, g++, `-O3 -flto`,
  IEEE-strict — *not* fast-math; AMReX `2d-g++-26.06`) **with** the two `Fapply`
  edits → "MODIFIED". Then `git checkout -- src/Operator/Elastic.cpp` (revert to the
  `chamber-gpu` baseline `Fapply`, **only** that file changed), rebuilt → "BASELINE".
  Worktree restored to the modified source + binary afterward.
- **Case.** `tests/GPU/C1_correctness_elastic/input`, `max_step=30` → **6 elastic
  solves** at `interval=5`, `mpiexec -np 1` (deterministic). This deck exercises
  **both** edited paths: `model_prop` vs `model_void` ⇒ non-uniform `C` ⇒ the
  `!m_uniform` grad(C) branch; `elastic.traction = 1.0_MPa` ⇒ the domain-boundary
  `sig` path (`trac_xhi_x = -280.721` etc., nonzero — so `Fapply` ran on real data,
  not a trivial zero field).
- **Result — byte-for-byte identical:**
  - `thermo.dat`: identical across all 31 rows, **including** the elastic
    `trac_*` / `disp_*` boundary diagnostics through all 6 solves.
  - Field plotfiles (`amr.plot_int=30`): **every** non-metadata data file identical
    via `cmp`, including the **node-centered elastic field** `00030node/Level_{0,1,2}/
    Cell_D_00000` on all three AMR levels (the direct `Fapply` solution output) —
    this is what proves the **interior**-node edits (sig-sink + one-Matrix4 grad(C))
    changed nothing. Only `metadata` and the embedded `diff.patch` differ (build
    provenance + the source diff itself — expected).
- **Record:** `benchmark/PHASE_C1_cpu_golden_compare.md` (full commands + output).

## ✅ A100 before/after — DONE (2026-07-06/07, NOVA `nova21-gpu-9`)

Ran `benchmark/PHASE_C1_nova_ab.md` §§1–5 end to end. Worktree `alamo-c1ab` on NOVA
(`/work/brunnels/jackplum/alamo-c1ab`), branch tip `123de00a2`. Jobs: build
`11433898`, wall A/B `11433933`, first ncu attempt `11433932` (mis-scoped — profiled
init kernels, not the solve), retargeted ncu `11448954` (gated directly on the
`Operator::Elastic::Fapply()` NVTX range — same fix as the June 27 ncu-empty
postmortem).

**Registers / occupancy** (`--nvtx-include "…/Operator::Elastic::Fapply()/"`, 6
launches sampled per arm, two grid sizes from AMR level split):

| Arm      | registers/thread | achieved occupancy | Fapply/launch (large, 8715 blk) | Fapply/launch (small, 4458 blk) |
|----------|------------------:|--------------------:|---------------------------------:|---------------------------------:|
| BASELINE | 255               | ~12.10%             | 4.667 ms                         | 2.286 ms                         |
| MODIFIED | 244               | ~11.91%              | 4.256 ms                         | 2.077 ms                         |
| Δ        | −11 (−4.3%)       | ~flat (−0.2pp, noise)| **−8.8%**                        | **−9.1%**                        |

Registers dropped off the 255 hard cap as predicted; occupancy did **not** cross a
new tier (still ~12%, same regime as Phase A) — 244 reg/thread × 256 threads/block
is still short of the ≤128 reg/thread needed to fit a second block/SM on A100's
65536-register file, so this is the "register pressure reduced but occupancy
unchanged" outcome the procedure flagged as possible. Per-launch time still dropped
~9%, consistent with reduced register-spill traffic rather than an occupancy-driven
win — matches the by-construction removal of 90 live doubles (135→45) and the
boundary-only `sig` sink.

**Wall A/B** (TinyProfiler, device-synced, job `11433933`, identical work both arms —
105,920 `Fapply` calls, 13 MLMG solves):

| Region                        | BASELINE | MODIFIED | Δ               |
|--------------------------------|---------:|---------:|-----------------|
| `Operator::Elastic::Fapply()`  | 300.1 s  | 256.6 s  | −14.5% (1.17×)  |
| `MLMG::solve()` (inclusive)    | 373.2 s  | 329.7 s  | −11.7%          |
| Everything else                | —        | —        | within ±0.5%    |

(Raw total-wall comparison is confounded by `VisMF::Write` plotfile I/O on Lustre,
904s vs 523s for the same 84 writes — filesystem noise, not compute; compute-minus-
I/O matches the Fapply delta above.)

**GPU stress parity** (`fcompare.gnu.ex`, steps 0/300/600 of
`input_3d_centre_bore_256_a2`, `stop_time=0.06`): cell (flame) fields bit-identical
at every step (`PLOTFILE AGREE`, exact 0). Node (elastic) fields: step 0 exact 0;
steps 300/600 relative error ~1e-7–1e-8 on stress/strain/disp — consistent with
GPU fast-math/atomic reduction-order noise (same class as the CPU-vs-GPU note in
`gpu_elastic_d1_reversed`), not a physics divergence.

**Verdict — PASS.** Pass condition was "registers/thread < 255 and/or achieved
occupancy > ~12.5%; stress field matches baseline to tolerance." Registers cleared
(255→244); occupancy did not move; parity holds. C1 is now perf-proven, not just
correctness-verified: **~9% per-launch Fapply speedup, ~14.5% Fapply wall, ~11.7%
MLMG wall**, zero physics regression. Occupancy being flat is the designed trigger
for the next lever — PLAN.md task 3.2 (`__launch_bounds__` sweep, deferred lever A1a
above) is now the correct next step, not optional follow-up.

## Deferred levers (intentionally not touched here)

- **A1a — per-kernel `__launch_bounds__` register cap.** The build already caps at
  255 globally; a kernel-local `__launch_bounds__(256, N)` could raise occupancy
  but **can regress** by forcing more spills. Stage it as a measured toggle once
  the A1 ncu SoL exists — do not land blind.
- **A2 / C2 — SetModel per-node "material-interface" mask + AoS→SoA Matrix4.**
  Biggest single kernel win but larger refactor; gate on the ncu memory-vs-compute
  SoL.
- **A4 / C4 — reuse the elastic operator across solves** (stop rebuilding
  `Elastic::define()` + `prepareForSolve` every solve). Independent setup-cost
  lever.

These two source edits are the safe, counter-justified first step; everything
above is sequenced behind a NOVA measurement.
