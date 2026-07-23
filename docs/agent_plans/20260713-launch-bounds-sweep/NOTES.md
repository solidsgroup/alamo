# NOTES — launch-bounds-sweep A100 leg (submit phase, 2026-07-13)

Covers PLAN.md Step 3: submit 4-arm A100 sweep on NOVA. **Nothing waited
on** — this is the queued state as of submission time. Followed the
`fapply-322b-a100` NOVA job pattern (docs/agent_plans/20260713-fapply-322b-a100/NOTES.md)
almost exactly; deviations noted below.

## Arms

| Arm | ALAMO_ELASTIC_MIN_BLOCKS | Role |
|-----|-------------------------:|------|
| mb1 | 1 | sanity anchor — resource-identical to knob-off; Fapply wall should land near this morning's MODIFIED arm's 255.7s |
| mb2 | 2 | 128-reg cap |
| mb3 | 3 | 80-reg cap |
| mb4 | 4 | 64-reg cap |

All 4 arms built from the same commit: worktree `/home/jackplum/Projects/alamo-fapply-322b`,
branch `launch-bounds-sweep`, tip `23fc0f3b9` (rebased onto chamber-gpu,
includes 3.2b Fapply/Fsmooth surgery + the launch-bounds knob + the local
ptxas sweep). Only the build-time define differs between arms — no
per-arm source checkout differences (unlike the baseline/modified split
this morning, which was a real source diff).

## Local verification of the flag-routing mechanism (before touching NOVA)

Confirmed the mechanism this task's setup asked me to verify, before
submitting anything:

1. `src/Operator/ElasticLaunch.H` (worktree) gates on
   `#ifndef ALAMO_ELASTIC_MIN_BLOCKS #define ALAMO_ELASTIC_MIN_BLOCKS 0`,
   i.e. any `-DALAMO_ELASTIC_MIN_BLOCKS=N` reaching the Elastic.cpp /
   Operator.cpp translation units controls the knob.
2. The repo's generated `Makefile` compiles every `.cpp.o` via
   `$(COMP_CMD) $< -o $@ ${ALAMO_INCLUDE} ${CXX_COMPILE_FLAGS}`, and
   `configure`'s CUDA branch resets `CXX_COMPILE_FLAGS = ` (plain `=`) once,
   then appends `+=` several more times. GNU Make semantics: a
   command-line-supplied value for `CXX_COMPILE_FLAGS` overrides the `=`
   reset but is *kept* and extended by every subsequent `+=` in the
   makefile — so `make ... CXX_COMPILE_FLAGS="-DALAMO_ELASTIC_MIN_BLOCKS=N"`
   is a safe injection point without touching configure or the Makefile.
3. Verified empirically, not just by reading: ran
   `CUDA_HOME=/home/jackplum/Projects/alamo/.local/cuda-12.6.3-redist
   ./configure --comp=g++ --dim 3 --cuda 80 --profile --offline` in the
   worktree (AMReX sm_80 tree already built there from the Step 2 ptxas
   sweep, so no rebuild needed), then:
   - `make -n bin/alamo_gpu CXX_COMPILE_FLAGS="-DALAMO_ELASTIC_MIN_BLOCKS=2"`
     — confirmed `-DALAMO_ELASTIC_MIN_BLOCKS=2` appears at the tail of the
     real `nvcc -x cu -dc ... src/Operator/Elastic.cpp` command line.
   - Actually compiled (not dry-run) `obj/obj-3d-profile-cuda80-g++/Operator/Elastic.cpp.o`
     for min_blocks=1 and min_blocks=2 with `-Xptxas -v` appended: min_blocks=1
     shows the Fapply body at 254 registers (matches Step 2's ptxas table
     exactly); min_blocks=2 shows it clamped to 128 registers (also matches
     the table's mb2 row), while unrelated reduce-op kernels in the same TU
     stay at 254 — confirming the knob only affects the 3 wired call sites,
     not a blanket recompile-everything effect.
   - `git status --porcelain` in the worktree after this test: clean (obj/
     is gitignored, nothing else touched).

This is the verification the task asked for before trusting the build
mechanism on NOVA. **Flag routing confirmed working — did not need to
stop.**

## NOVA paths

- Bundle: `/work/brunnels/jackplum/alamo-lbsweep.bundle` (single ref
  `_bundle_lbsweep` -> `23fc0f3b9`; temp ref created locally via
  `git update-ref`, bundled, then deleted locally post-bundle — same
  pattern as `alamo322b.bundle` this morning, ~275MB, full history).
- Arm checkouts: `/work/brunnels/jackplum/alamo-lbsweep/{mb1,mb2,mb3,mb4}/`,
  each `git clone`d from the bundle and checked out to `_bundle_lbsweep`.
  Verified on NOVA: all 4 `git rev-parse HEAD` = `23fc0f3b95f13929f122eff4bdf22a149beb8977` ✓
- sbatch scripts (written directly on NOVA via heredoc, templated from
  this morning's `build_322b.sbatch` / `wall_322b.sbatch` / `ncu_322b.sbatch`
  — not added under the local repo's `benchmark/`, live only on NOVA per
  the same convention as the 322b-a100 precedent):
  - `<arm>/build_lb.sbatch` — `ARCHES=80` only (`--cuda 80`), same as
    precedent; `make ... bin/alamo_gpu CXX_COMPILE_FLAGS="-DALAMO_ELASTIC_MIN_BLOCKS=<n>"`
    for the arm's `n`; appends a `cuobjdump --dump-resource-usage` register
    check for both Elastic.cpp.o (Fapply) and Operator.cpp.o (Fsmooth)
    straight into the build log, so the register-verification evidence
    lands automatically without a separate step.
  - `<arm>/wall_lb.sbatch` — identical to `wall_322b.sbatch`: TinyProfiler
    wall run, deck `input_3d_centre_bore_256_a2`,
    `stop_time=0.06 elastic.solver.verbose=4 amr.plot_int=30`,
    `--time=01:00:00` (task said this is enough; this morning's ran 16 min).
  - `<arm>/ncu_lb.sbatch` — same NVTX-gated pattern
    (`--nvtx-include "Operator::Elastic::Fapply()/"`, `--set basic`,
    `--launch-count 40`) PLUS
    `--metrics smsp__inst_executed.sum,sm__cycles_active.sum` added per
    this task's instruction (this morning's basic set lacked raw
    instruction/cycle counts); the closing per-kernel summary table also
    picks up the two new metrics.

## Job IDs + dependency graph

| Step  | Arm | min_blocks | Job ID   | Depends on |
|-------|-----|-----------:|---------:|------------|
| build | mb1 | 1 | 11640750 | — |
| build | mb2 | 2 | 11640751 | — |
| build | mb3 | 3 | 11640752 | — |
| build | mb4 | 4 | 11640753 | — |
| wall  | mb1 | 1 | 11640754 | afterok:11640750 |
| ncu   | mb1 | 1 | 11640755 | afterok:11640750 |
| wall  | mb2 | 2 | 11640756 | afterok:11640751 |
| ncu   | mb2 | 2 | 11640757 | afterok:11640751 |
| wall  | mb3 | 3 | 11640758 | afterok:11640752 |
| ncu   | mb3 | 3 | 11640759 | afterok:11640752 |
| wall  | mb4 | 4 | 11640760 | afterok:11640753 |
| ncu   | mb4 | 4 | 11640761 | afterok:11640753 |

Confirmed via `squeue -u jackplum` immediately after submission: 4 build
jobs `PD` (queued, Priority reason — cluster busy at submit time, unlike
this morning where they started running immediately), all 8 dependent
jobs `PD` with `Dependency` reason and the correct unfulfilled
`afterok:<jobid>` per arm.

## Register-verification plan per arm (Step 3 CHECK)

Each `build_lb.sbatch` prints `cuobjdump --dump-resource-usage` output for
the Fapply and Fsmooth symbols directly into `build_lb.<jobid>.out`, right
after the build. Expected register counts, from the Step 2 local ptxas
table (`results/ptxas_table.md`, worktree):

| Arm | Fapply<Sym::Major> expected regs | Diagonal expected regs | Fsmooth expected regs |
|-----|----------------------------------:|------------------------:|------------------------:|
| mb1 | 254 | 184 | 96 |
| mb2 | 128 | 128 | 96 |
| mb3 | 80  | 80  | 79 |
| mb4 | 64  | 64  | 64 |

Collection-phase check: grep each `build_lb.<jobid>.out` for the Fapply/
Fsmooth register lines and diff against this table. A mismatch = misrouted
flag for that arm — kill that arm's wall + ncu jobs
(`scancel <wall_jobid> <ncu_jobid>`) and report before trusting any A100
numbers from it. (`cuobjdump --dump-resource-usage` labels differently
than `ptxas -v` — collection phase should also cross check against the raw
`nvcc -Xptxas -v` output if `cuobjdump`'s symbol demangling makes the
Fapply/Diagonal/Fsmooth entries ambiguous; the local verification above
used raw `-Xptxas -v` and it was unambiguous there.)

## Output file locations (once jobs complete)

- Build logs: `<arm>/build_lb.<jobid>.out` / `.err` (includes the register
  check); binary `<arm>/bin/alamo_gpu-3d-profile-cuda80-g++`.
- Wall logs: `<arm>/wall_lb.<jobid>.out` / `.err` and
  `<arm>/wall_lb_<arm>.log` (TinyProfiler table); plotfiles under
  `<arm>/out_lb_<arm>/plot*` at steps 0/30/... for fcompare parity vs the
  mb1 anchor.
- ncu: `<arm>/profiles_lb/<arm>_nvtx.ncu-rep`,
  `<arm>/profiles_lb/<arm>_nvtx.log`, plus the per-kernel CSV summary
  (now including `smsp__inst_executed.sum,sm__cycles_active.sum`) printed
  at the end of `<arm>/ncu_lb.<jobid>.out`.

## Exact commands used

```bash
# Local: verify flag routing (worktree /home/jackplum/Projects/alamo-fapply-322b)
CUDA_HOME=/home/jackplum/Projects/alamo/.local/cuda-12.6.3-redist \
  ./configure --comp=g++ --dim 3 --cuda 80 --profile --offline
make -n bin/alamo_gpu CXX_COMPILE_FLAGS="-DALAMO_ELASTIC_MIN_BLOCKS=2"   # dry-run check
rm -f obj/obj-3d-profile-cuda80-g++/Operator/Elastic.cpp.o
make CXX_COMPILE_FLAGS="-DALAMO_ELASTIC_MIN_BLOCKS=1 -Xptxas -v" \
  obj/obj-3d-profile-cuda80-g++/Operator/Elastic.cpp.o   # -> 254 regs, matches table
rm -f obj/obj-3d-profile-cuda80-g++/Operator/Elastic.cpp.o
make CXX_COMPILE_FLAGS="-DALAMO_ELASTIC_MIN_BLOCKS=2 -Xptxas -v" \
  obj/obj-3d-profile-cuda80-g++/Operator/Elastic.cpp.o   # -> 128 regs, matches table

# Local: bundle the launch-bounds-sweep tip
git update-ref refs/heads/_bundle_lbsweep 23fc0f3b9
git bundle create /tmp/alamo-lbsweep.bundle refs/heads/_bundle_lbsweep
git update-ref -d refs/heads/_bundle_lbsweep
scp /tmp/alamo-lbsweep.bundle nova:/work/brunnels/jackplum/alamo-lbsweep.bundle

# NOVA: clone 4 arms
mkdir -p /work/brunnels/jackplum/alamo-lbsweep
cd /work/brunnels/jackplum/alamo-lbsweep
for arm in mb1 mb2 mb3 mb4; do
  git clone /work/brunnels/jackplum/alamo-lbsweep.bundle "$arm"
  git -C "$arm" checkout _bundle_lbsweep
done

# NOVA: submit build jobs, then dependent wall+ncu jobs
cd /work/brunnels/jackplum/alamo-lbsweep
sbatch --parsable mb1/build_lb.sbatch   # -> 11640750
sbatch --parsable mb2/build_lb.sbatch   # -> 11640751
sbatch --parsable mb3/build_lb.sbatch   # -> 11640752
sbatch --parsable mb4/build_lb.sbatch   # -> 11640753
sbatch --parsable --dependency=afterok:11640750 mb1/wall_lb.sbatch   # -> 11640754
sbatch --parsable --dependency=afterok:11640750 mb1/ncu_lb.sbatch    # -> 11640755
sbatch --parsable --dependency=afterok:11640751 mb2/wall_lb.sbatch   # -> 11640756
sbatch --parsable --dependency=afterok:11640751 mb2/ncu_lb.sbatch    # -> 11640757
sbatch --parsable --dependency=afterok:11640752 mb3/wall_lb.sbatch   # -> 11640758
sbatch --parsable --dependency=afterok:11640752 mb3/ncu_lb.sbatch    # -> 11640759
sbatch --parsable --dependency=afterok:11640753 mb4/wall_lb.sbatch   # -> 11640760
sbatch --parsable --dependency=afterok:11640753 mb4/ncu_lb.sbatch    # -> 11640761
```

## Deviations from plan / precedent

- Used 4 separate arm dirs (`mb1..mb4`) each cloned from the same commit
  (`23fc0f3b9`) rather than one checkout with a swapped file, since the
  difference between arms is a build-time `CXX_COMPILE_FLAGS` define, not
  a source diff — this is actually *simpler* than the 322b-a100 baseline/
  modified split (no per-arm source difference to keep straight), but
  matches its directory-per-arm structure for consistency with the
  existing wall/ncu sbatch templates.
- Register verification lives inside each `build_lb.sbatch` itself (a
  `cuobjdump --dump-resource-usage` check appended right after the `make`
  call) rather than as a separate collection-phase step, so the evidence
  is captured automatically in the build log even if a later session
  forgets to re-derive it. Collection phase should still cross-check
  against `results/ptxas_table.md` (worktree) as the source of truth,
  and fall back to raw `nvcc -Xptxas -v` diffing if `cuobjdump`'s output
  is ambiguous for any arm (untested here — the local verification used
  raw ptxas -v, not cuobjdump, so cuobjdump's exact output format for this
  binary is unconfirmed).
- ncu `--metrics` list extended with `smsp__inst_executed.sum,sm__cycles_active.sum`
  per explicit task instruction; everything else in the ncu job (NVTX
  range, `--set basic`, `--launch-count 40`, resource footprint) is
  unchanged from the 322b-a100 precedent.
- Did not create a wrapper script under `docs/agent_plans/.../` in the
  main repo (sbatch content lives only on NOVA, matching precedent); this
  NOTES.md's "Exact commands used" section is the reproduction record.
- NOVA cluster was busier at submission time than this morning (build
  jobs queued `PD`/Priority instead of starting immediately as `R`) — not
  a routing problem, just queue state; flagging so the collection-phase
  agent doesn't read early `PD` status as a submission failure.

## CORRECTION 1 (2026-07-13) — flag routing was BROKEN, first round scancelled

The first submission (builds 11640750-53, dependents 11640754-61) **failed**.
The coordinator caught it: passing `-DALAMO_ELASTIC_MIN_BLOCKS=n` as a
**command-line** `make CXX_COMPILE_FLAGS="..."` argument OVERRODE every
`CXX_COMPILE_FLAGS +=` in the makefile (GNU make: command-line vars beat all
makefile assignments including `+=`). That stripped `-std=c++20`, `-Winline`,
and the `METADATA_*` defines from every TU; `WriteMetaData.cpp` died on
undefined `METADATA_GITHASH` (build_lb.11640750.err:9252). My original local
verification missed this because I only compiled `Elastic.cpp` in isolation,
which doesn't reference the METADATA defines — the single-file test was not a
faithful proxy for the full `bin/alamo_gpu` link.

**The coordinator's proposed fix (env var instead of command-line var) does
NOT work either — I verified this before resubmitting.** `configure`
regenerates `.make/Makefile.pre.conf`, whose **line 10 is a plain
`CXX_COMPILE_FLAGS = ` (empty) reset**, and that fragment is `-include`d at the
very top of the Makefile. GNU make: a makefile `=` assignment overrides an
environment variable. So an exported `CXX_COMPILE_FLAGS` gets wiped by line 10
just like it does the makefile's own `+=` chain. Confirmed empirically with
`make -n` (worktree): under the env-var approach the `Elastic.cpp` compile line
contains `-std=c++20` and `-Winline` but is **missing**
`-DALAMO_ELASTIC_MIN_BLOCKS` entirely — i.e. it would have silently produced
knob-OFF (min_blocks=0) binaries for all 4 arms. The register check would have
caught it (all arms 254 regs), but only after another wasted build+queue cycle.

**Working fix (verified end-to-end locally):** after `./configure`, append our
own line to the TAIL of the regenerated fragment:

```bash
echo 'CXX_COMPILE_FLAGS += -DALAMO_ELASTIC_MIN_BLOCKS=<n>' >> .make/Makefile.pre.conf
```

This lands AFTER the `=` reset (pre.conf:10) and accumulates with every later
`+=` (GIT_DIFF on pre.conf:16, and the main Makefile:41 `-std=c++20`/`-Winline`/
`METADATA_*` line). Verified in-worktree via `make -n` that the `Elastic.cpp`
compile line then contains ALL of `-DALAMO_ELASTIC_MIN_BLOCKS`, `-std=c++20`,
`-Winline`, AND `-DMETADATA_GITHASH`; and by actually compiling `Elastic.cpp.o`
with `-Xptxas -v` that the Fapply body clamps to 128 regs for n=2 (matches the
ptxas table) while unrelated reduce-op siblings stay at 254.

Each `build_lb.sbatch` now:
1. appends the `+=` line to `.make/Makefile.pre.conf` after configure (and
   `grep`s it back to prove it landed);
2. runs `make ... bin/alamo_gpu QUIET=` (empty `QUIET` disables the Makefile's
   `QUIET ?= @` recipe-echo suppression — `VERBOSE=` is not honored by this
   Makefile) so the real compile command lines are captured in
   `build_lb_make.log`;
3. asserts the `Elastic.cpp` compile line contains BOTH
   `-DALAMO_ELASTIC_MIN_BLOCKS=<n>` AND `-std=c++20` (prints `FLAG_APPEND_OK`
   or `FLAG_APPEND_FAIL`) — a direct guard against the override class of bug;
4. keeps the `cuobjdump --dump-resource-usage` register check.

**Cleanup done:** `scancel`led orphaned dependents 11640754-11640761, the
doomed mb4 build 11640753, and the piggyback devarena job 11640783. Queue
confirmed empty before resubmit.

### Corrected job IDs + dependency graph (round 2)

| Step  | Arm | min_blocks | Job ID   | Depends on |
|-------|-----|-----------:|---------:|------------|
| build | mb1 | 1 | 11641025 | — |
| build | mb2 | 2 | 11641026 | — |
| build | mb3 | 3 | 11641027 | — |
| build | mb4 | 4 | 11641028 | — |
| wall  | mb1 | 1 | 11641029 | afterok:11641025 |
| ncu   | mb1 | 1 | 11641030 | afterok:11641025 |
| wall  | mb2 | 2 | 11641031 | afterok:11641026 |
| ncu   | mb2 | 2 | 11641032 | afterok:11641026 |
| wall  | mb3 | 3 | 11641033 | afterok:11641027 |
| ncu   | mb3 | 3 | 11641034 | afterok:11641027 |
| wall  | mb4 | 4 | 11641035 | afterok:11641028 |
| ncu   | mb4 | 4 | 11641036 | afterok:11641028 |
| wall (devarena, mb1) | mb1 | 1 | 11641037 | afterok:11641025 |

Confirmed via `squeue`: 4 builds `PD`/Priority, all 9 dependents `PD`/Dependency
with correct unfulfilled `afterok`. The round-1 IDs (11640750-61, 11640783) are
dead — do not reference them for results.

The register-verification table (expected Fapply/Diagonal/Fsmooth regs per arm)
in the earlier section is unchanged and still the source of truth for the
collection phase; each build log now also carries the `FLAG_APPEND_OK`/`FAIL`
line as a faster first-line guard.

## Not done (next phase)

- Waiting for all 12 jobs to complete, grepping build logs for the
  register-verification table above, extracting TinyProfiler Fapply
  exclusive-wall numbers per arm, running fcompare mb1 vs mb2/mb3/mb4
  plotfiles for stress parity, ncu CSV analysis (including the new
  instruction/cycle metrics), and writing `results/RESULT.md` with the
  verdict — that is PLAN.md Step 4, out of scope for this submit-only
  session.

## Orchestrator addition (2026-07-13): arena-policy A/B piggyback

PLAN.md backlog item "Arena policy A/B (Phase 5)" — zero-src-change runtime
flag. Job 11640783 = mb1 wall job cloned with `amrex.the_arena_is_managed=0`
(device arena), afterok:11640750 (mb1 build). Compare its TinyProfiler tables
against mb1's managed run (job 11640754). Outputs:
`mb1/wall_lb_devarena.<jobid>.out`, `mb1/wall_lb_mb1_devarena.log`. If the
device-arena run OOMs on the 40GB A100, that is itself the answer (256^3 deck
does not fit unmanaged) — record and move on.
