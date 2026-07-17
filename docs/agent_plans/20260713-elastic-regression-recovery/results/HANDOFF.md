# Elastic regression recovery continuation handoff

## Start here

Status is **BLOCKED**, not complete.  Work only in:

~~~
/tmp/alamo-elastic-void-continue
~~~

The branch is `elastic-void-wip-20260713` at HEAD `a06c6f15d`; all recovery
changes are uncommitted.  Preserve the dirty worktree and retained output
directories.  Do not use a psi floor, artificial stiffness, tolerance
relaxation, a larger nonlinear cap, forced acceptance, or reference refresh.
Do not create `results/DONE` or commit unless every plan gate passes and the
user explicitly requests a commit.

Read, in order:

1. `RESULT.md` — concise outcome and verification matrix.
2. `REVIEW.md` — blocking findings and risks.
3. `../NOTES.md` — chronological experiment almanac and rejected hypotheses.
4. `../PLAN.md` — acceptance gates.
5. `../../20260709-elastic-void-robustness/results/REGRESSION_RECOVERY_HANDOFF.md`
   only for the original failure inventory.

Run `benchmark/status.sh` before doing anything else.  Hooks are already set
to `.githooks`.

## Exact current state

- Current configuration is a 2-D production GCC build.  `make -j4` and
  `bin/test-2d-g++` pass on the exact closeout source.
- The paired elastic rewrite has been replaced by trusted parent `d964cfab8`
  behavior.  Relative to that parent, the intentional source delta is limited
  to Newton hardening, Constant BC spacing, Flame periodic/full-halo model
  construction, two-ring diagonal construction, and valid-only normalization.
- Fresh final-source 2-D serial and two-rank runs complete and are numerically
  identical, but their checker fails force balance at `0.270767`.
- The retained final 3-D run completes with finite real MLMG work and valid
  periodic symmetry, but fails force balance at `0.270142`.
- Tightening the positive update tolerance to `1e-7` does not solve the issue:
  five 3-D iterations end at nonlinear residual ratio `0.139951` and abort.
- No final complete suites, CUDA correctness, sanitizer, comparison, or timing
  campaign has been run.  The GPU is visible, but no `nvcc` is currently
  available.

Canonical artifacts:

~~~
tests/ElasticSoftVoid/output_2026-07-16_15.26.33_kermit_2d-serial
tests/ElasticSoftVoid/output_2026-07-16_15.26.33_kermit_2d-parallel
tests/ElasticSoftVoid/output_2026-07-13_19.02.34_kermit_3d-serial
report/output_2026-07-16_15.26.33_kermit.json
report/output_2026-07-16_15.26.33_kermit.html
~~~

Ignore `report/output_2026-07-16_15.24.46_kermit.*`: an argparse ordering
mistake caused unrelated Dendrite/ThermoElastic sections to be interrupted.

## First investigation — do not tune parameters

The first final-source 2-D solve is the cheapest deterministic reproducer and
uses a linear constitutive model.  It starts from nonlinear RHS
`2.616043465e8`; after the first linear correction (11 MLMG iterations,
relative residual `2.962146737e-6`) the assembled nonlinear residual is still
`6.9295e7`.  After the second correction it remains `6.19397e7`, even though
the update criterion declares convergence.  A consistent linear residual and
tangent should not behave this way.

Before any source edit:

1. Locate the maximum final plotted `res_*` by AMR level and classify it as
   interior, physical-boundary, or outer coarse/fine-union row.
2. On Newton iteration 1, compare valid rows of the residual assembled by
   `prepareForSolve` with the corresponding `Fapply(u)`/boundary equation.
3. Verify `Fapply(dsol)` equals the current linear RHS to the requested MLMG
   tolerance, separately for interior, physical-boundary, and C/F rows.
4. Repeat minimal toggles only to localize the contract: AMR level 0 versus
   explicit AMR, uniform 140/150 MPa versus the 0.2 MPa mixture, and periodic
   extrusion versus 2-D.  These are diagnostics, not candidate acceptance
   settings.
5. Record the invariant and expected repair in `../NOTES.md` before editing.

Likely code surfaces are:

~~~
src/Solver/Nonlocal/Newton.H        prepareForSolve/residual and coefficient resync
src/Operator/Elastic.cpp            Fapply, Diagonal, coefficient hierarchy
src/Operator/Operator.cpp           smoother/normalize and C/F synchronization
src/BC/Operator/Elastic/Constant.H  physical boundary equation
src/Integrator/Flame.cpp            model/RHS construction and plot residual
~~~

Do not begin by adding `nr_convergence=residual`.  That can be a useful final
contract only after residual/tangent consistency is proved; with the current
five-iteration cap it would correctly abort rather than pass.

## Reproducer commands

Keep the positional test directory before `--sections`, because that option
uses `nargs=*`.

~~~
./configure --dim=2 --comp=g++ --no-debug --build-amrex --offline
make -j4
bin/test-2d-g++
scripts/runtests.py tests/ElasticSoftVoid --comp=g++ --no-backspace \
  --no-clean --sections 2d-serial 2d-parallel
~~~

The expected current result is two run passes and two check failures at
equilibrium ratio `0.270767`.  Once 2-D passes, rebuild 3-D:

~~~
./configure --dim=3 --comp=g++ --no-debug --build-amrex --offline
make -j4
bin/test-3d-g++
scripts/runtests.py tests/ElasticSoftVoid --comp=g++ --no-backspace \
  --no-clean --sections 3d-serial
~~~

Then rerun focused trusted mechanics before the full suites:

~~~
scripts/runtests.py tests/PlateHole --comp=g++ --no-backspace --no-clean \
  --sections 2D-serial
scripts/runtests.py tests/FracturePFCZM --comp=g++ --no-backspace --no-clean \
  --sections 2d-serial-notch
scripts/runtests.py tests/SCPSpheresElastic --comp=g++ --no-backspace \
  --no-clean --sections 2d-parallel-long
~~~

## Acceptance sequence after the repair

1. `ElasticSoftVoid` 2-D serial/MPI and 3-D all pass their unchanged
   equilibrium, sub-MPa, AMR, mixture, and symmetry gates.
2. Unit, PlateHole, FracturePFCZM, EshelbyFiniteKinematics,
   VoronoiSimplePeriodic, and SCPSpheresElastic pass on the exact source.
3. Complete 2-D and 3-D GCC regression suites and `make test` pass.
4. A fresh Tier-3 reviewer inspects the final diff for residual/tangent/
   diagonal agreement, boundary closure, halo lifecycles, GPU safety, and
   oracle weakening.
5. Make a CUDA toolkit available and build strict and fast Flame binaries for
   dimensions 2 and 3.  Run plotted correctness and the repository comparison
   gate before performance.
6. For performance, use one excluded warm-up plus three rank-one repetitions
   per backend/dimension, exactly five elastic solves, no plot I/O, and retain
   wall time, MLMG/Newton histories, hashes, peak GPU memory, medians, and
   speedups.  Report the reduced 3-D and local rank-one scope explicitly.

CUDA was previously planned with:

~~~
./configure --dim=2 --comp=g++ --no-debug --cuda=86 --cuda-fp=strict \
  --gpu-integrator=flame --build-amrex --offline
./configure --dim=3 --comp=g++ --no-debug --cuda=86 --cuda-fp=fast \
  --gpu-integrator=flame --build-amrex --offline
~~~

Verify the actual toolkit and generated binary names rather than assuming the
commands will work in the current environment.

## Closeout discipline

Append all new experiments to `../NOTES.md`; do not rewrite the historical
blocked records under the 20260709 task.  Update `RESULT.md` and `REVIEW.md`
only when evidence changes.  Run `git diff --check` and `benchmark/status.sh`
at each handoff.  Add `results/DONE` only after every checkbox in `../PLAN.md`
is genuinely complete.
