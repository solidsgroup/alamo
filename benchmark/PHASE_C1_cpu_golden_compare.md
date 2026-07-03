# Phase C1 — CPU golden compare for the `Fapply` source edits (DONE, bit-identical)

**Date:** 2026-06-28 · **Branch/worktree:** `chamber-gpu-elastic-opt` at
`/home/jackplum/Projects/alamo-elastic-opt`.

This is the correctness half of the C1 gate (`PHASE_C1_fapply_occupancy.md`): prove
the two `Operator::Elastic<SYM>::Fapply` edits (one-`Matrix4`-at-a-time grad(C);
boundary-only `sig`) are bit-for-bit equivalent to the `chamber-gpu` baseline. The
edits introduce no FP reordering, so the expectation was bit-identity; this confirms
it end to end. The A100 perf before/after remains owed — see `PHASE_C1_nova_ab.md`.

## Build (both arms identical except `src/Operator/Elastic.cpp`)

CPU 2D, g++, production `-O3 -flto`, **IEEE-strict (not fast-math)**, AMReX
`2d-g++-26.06` (reused from the main checkout). Configured with:

```bash
cd /home/jackplum/Projects/alamo-elastic-opt
./configure --dim 2 --comp g++ --no-debug \
    --amrex /home/jackplum/Projects/alamo/ext/AMReX-Codes/amrex/2d-g++-26.06 --eigen /usr
# configure left AMREX_TARGET empty and set EIGEN=/usr/include (the latter breaks the
# libc #include_next chain). Fix applied to .make/Makefile.pre.conf:
#   - set AMREX_TARGET = <the 2d-g++-26.06 path>
#   - drop the EIGEN line (system /usr/include/eigen3 is found by default)
make -j24      # -> bin/alamo-2d-g++   == MODIFIED
```

Baseline arm: `git checkout -- src/Operator/Elastic.cpp` (reverts **only** that file
to the `chamber-gpu` tip `Fapply`), `make -j24` → BASELINE. The worktree source +
binary were restored to MODIFIED afterward (`git diff --stat` again shows
`Elastic.cpp | 32 +++…`).

## Case — exercises BOTH edited paths

`tests/GPU/C1_correctness_elastic/input`, run from the worktree root (so the relative
`simple_circle.bmp` / `base_circle0.bmp` ICs resolve):

```bash
mpiexec -np 1 ./bin/alamo-2d-g++ tests/GPU/C1_correctness_elastic/input \
    max_step=30 amr.plot_int=30 plot_file=<out>/plot \
    elastic.solver.verbose=0 elastic.print_model=0
```

- `max_step=30`, `elastic.interval=5` ⇒ **6 elastic solves**; `np 1` ⇒ deterministic.
- `model_prop` (kappa 162 MPa) vs `model_void` (4 MPa) ⇒ **non-uniform `C`** ⇒ the
  `!m_uniform` grad(C) branch (edit 1) runs at every interior node.
- `elastic.traction = 1.0_MPa` ⇒ the **domain-boundary `sig`** branch (edit 2) runs;
  observed `trac_xhi_x = -280.721`, `trac_yhi_y = -278.203` (nonzero — `Fapply` ran on
  a real loaded field, not a trivial zero state).

## Result — byte-for-byte identical

- `thermo.dat`: `cmp` identical across all **31 rows** (header + 30 steps), including
  every elastic boundary diagnostic (`trac_*`, `disp_*`) through all 6 solves.
- Field plotfiles: `cmp` over **every** non-metadata file under the step-0 and step-30
  plotfiles passed — cell data and, critically, the **node-centered elastic field**
  `00030node/Level_{0,1,2}/Cell_D_00000` on all three AMR levels (the direct `Fapply`
  solution output). This is the part that proves the **interior**-node edits changed
  nothing — the boundary-only `thermo` diagnostics alone could not.
- Only `metadata` (build/timestamp provenance) and `diff.patch` (the embedded git
  diff, which *is* the source change) differ. Both expected; neither is physics data.

## Verdict

The `Fapply` source edits are **bit-identical** to the `chamber-gpu` baseline on CPU,
end to end (thermo + full AMR field, 6 elastic solves). The "correctness by
construction" argument is now empirically confirmed. The only remaining C1 gate item
is the **A100 occupancy/perf before/after** (`PHASE_C1_nova_ab.md`).
