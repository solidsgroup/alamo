#!/usr/bin/env python3
"""F6: scripted MLMG convergence benchmark over contrast x resolution x
coarsening depth, sweeping tests/LowMachElasticPressure.

Part of the diagnostics-first pass on MLMG convergence at high modulus
contrast (~/.claude/plans/improve-MLMG-solver-ideas.md). The only prior
evidence was a hand-run table in
tests/LowMachElasticPressure/Readme.rst:145-252, produced outside
scripts/runtests.py, and it never varied resolution -- so it cannot show
whether iteration count is h-independent, which is the actual acceptance
test for any coarse-grid fix (D1a/D1b) and the one axis
elastic.max_coarsening_level capping can never recover on its own.

Design notes (see the plan for the full reasoning):

- Invokes ./bin/lowmach-<dim>d-<comp> directly, NOT scripts/runtests.py.
  runtests.py is built for pass/fail against a reference and cleans up
  its output directories; this wants a parameter grid and a table, and
  the full cross-product here is too large to express as `#@` sections
  in the input file without slowing down the regular test suite.
- Reuses the input file's own `#@` header (scripts/solverlib.read_test_config)
  for the base args of a chosen section (default: 2d-pressure-1Pa) rather
  than duplicating pressure.ic.constant.value / elastic.apply_fluid_pressure
  here -- if that input file changes, this stays in sync.
- Reuses scripts/solverlib.lowmach_pressure_errors for the disp_y/sigma_yy
  error columns, which is the same ray sample and closed-form comparison
  tests/LowMachElasticPressure/test itself uses.
- The h-refine axis's error columns are NOT an h-independence readout:
  the interface (LowMach's diffuse eta band) has a fixed physical width,
  so refining amr.n_cell changes the interface width in cells and moves
  the discretization error for entirely legitimate reasons. The
  h-independence readout is n_iter (and convfactor), full stop.
- Every run gets its own directory and is invoked with an absolute
  executable path, absolute input path, and cwd=<rundir>, so parallel
  runs do not collide on the Backtrace.N file alamo writes into cwd on
  abort.

Usage:
    ./scripts/solver_benchmark.py --preset contrast --csv out.csv
    ./scripts/solver_benchmark.py --preset cap-bisect --csv out.csv
    ./scripts/solver_benchmark.py --preset h-refine --csv out.csv
    ./scripts/solver_benchmark.py --contrast 1e3 1e4 --cap 2 -1 --ncell 32x64 --csv out.csv
"""
import argparse
import csv
import math
import os
import subprocess
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import solverlib

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TESTDIR = os.path.join(REPO_ROOT, "tests", "LowMachElasticPressure")
INPUT = os.path.join(TESTDIR, "input")

# The Readme's 11-point contrast sweep (tests/LowMachElasticPressure/Readme.rst:157-219).
README_CONTRASTS = [1, 3, 10, 33, 100, 333, 1000, 3333, 10000, 33333, 100000]

PRESETS = {
    # Reproduces the existing Readme table. This is the harness's own
    # acceptance test -- if this slice does not come back matching
    # Readme.rst (1x -> ~17 it.; 1000x uncapped -> diverges; 1e4x capped
    # -> ~12-18 it.), the harness is wrong, not the solver.
    #
    # One documented exception: at cap=2, 33,333x is a razor's-edge case
    # (the Readme itself calls out 52 iterations there vs. a flat 12-18
    # everywhere else, "near the edge of that range"). A ~1e-5 relative
    # change in elastic.void.model.mu/kappa (e.g. from rounding 140/33333
    # to 4 vs. 17 significant figures) has been observed to flip that one
    # cell between converging in ~51 iterations and diverging in 2 -- this
    # is genuine sensitivity of the operator at that specific contrast, not
    # a harness bug. Every other cell in the table (1x-10,000x, and both
    # branches at 100,000x) reproduces the Readme numbers exactly.
    "contrast": dict(contrasts=README_CONTRASTS, caps=[-1, 2], ncells=["32x64"]),
    # A2: bisect the coarsening depth at a few representative contrasts.
    "cap-bisect": dict(contrasts=[1000, 10000, 33333],
                        caps=[0, 1, 2, 3, 4, 5, 6, -1], ncells=["32x64"]),
    # The missing axis: does iteration count stay flat under refinement?
    "h-refine": dict(contrasts=[1000, 10000, 33333],
                      caps=[-1, 2], ncells=["32x64", "64x128", "128x256"]),
}


def parse_ncell(s):
    nx, ny = s.lower().split("x")
    return int(nx), int(ny)


def exe_path(exe, dim, comp):
    return os.path.join(REPO_ROOT, "bin", f"{exe}-{dim}d-{comp}")


def build_args(base_args, contrast, cap, ncell, verbose):
    nx, ny = parse_ncell(ncell)
    mu_solid, kappa_solid = 140.0, 150.0
    args = list(base_args.split())
    args += [
        f"elastic.void.model.mu={mu_solid / contrast}",
        f"elastic.void.model.kappa={kappa_solid / contrast}",
        f"amr.n_cell={nx} {ny}",
        f"elastic.solver.verbose={verbose}",
        "elastic.solver.abort_on_fail=0",
        "amrex.signal_handling=0",
        "amrex.throw_exception=1",
        # One solve: plot_int (not plot_dt) guarantees a plotfile after the
        # single step regardless of stop_time/plot_dt (Integrator.cpp's
        # "always plot the final step" fallback only fires for plot_int).
        "max_step=1",
        "amr.plot_int=1",
    ]
    if cap is not None and cap >= 0:
        args.append(f"elastic.max_coarsening_level={cap}")
    return args


def run_one(exestr, input_path, args, rundir, timeout):
    os.makedirs(rundir, exist_ok=True)
    plot_file = os.path.join(rundir, "output")
    cmd = [exestr, input_path] + args + [f"plot_file={plot_file}"]
    timed_out = False
    t0 = time.time()
    try:
        proc = subprocess.run(cmd, cwd=rundir, capture_output=True, text=True,
                               timeout=timeout)
        returncode, stdout = proc.returncode, proc.stdout + proc.stderr
    except subprocess.TimeoutExpired as e:
        timed_out = True
        returncode = -1
        def _s(x):
            return x.decode(errors="replace") if isinstance(x, bytes) else (x or "")
        stdout = _s(e.stdout) + _s(e.stderr)
    wall_s = time.time() - t0

    with open(os.path.join(rundir, "run.log"), "w") as f:
        f.write(stdout)

    return returncode, stdout, timed_out, wall_s, plot_file


def final_plotfile(plot_file_base):
    # plot_file names a directory; numbered plotfiles (e.g. 00000cell,
    # 00001cell) live inside it (see tests/LowMachElasticPressure/test's
    # equivalent glob "{outdir}/*cell/").
    import glob
    candidates = sorted(glob.glob(os.path.join(plot_file_base, "*cell")))
    return candidates[-1] if candidates else None


def sweep_row(exestr, input_path, base_args, contrast, cap, ncell, outdir, timeout,
              keep_plotfiles, verbose):
    tag = f"c{contrast:g}_cap{cap}_n{ncell}".replace(".", "p")
    rundir = os.path.join(outdir, tag)
    args = build_args(base_args, contrast, cap, ncell, verbose)
    returncode, stdout, timed_out, wall_s, plot_file_base = run_one(
        exestr, input_path, args, rundir, timeout)

    status = solverlib.classify_run(returncode, stdout, timed_out)
    records = solverlib.parse_solver_stats(stdout) or solverlib.parse_solver_stats_fallback(stdout)
    last = records[-1] if records else None

    disp_y_err = sigma_yy_err = math.nan
    plotfile = final_plotfile(plot_file_base)
    if status == "OK" and plotfile:
        try:
            pressure = 1.0  # matches base section's pressure.ic.constant.value=1.0
            disp_y_err, sigma_yy_err, _ = solverlib.lowmach_pressure_errors(plotfile, pressure)
        except Exception as e:  # noqa: broad -- error columns are best-effort
            print(f"  (warning: error-column computation failed for {tag}: {e})",
                  file=sys.stderr)

    if not keep_plotfiles and os.path.isdir(rundir):
        import shutil
        for entry in os.listdir(rundir):
            if entry.endswith("cell"):
                shutil.rmtree(os.path.join(rundir, entry), ignore_errors=True)

    return dict(
        contrast=contrast, ncell=ncell, cap_requested=cap,
        status=status,
        n_iter=last.iters if last else -1,
        resid0=last.resid0 if last else math.nan,
        residN=last.resid if last else math.nan,
        convfactor=last.convfactor if last else math.nan,
        disp_y_err=disp_y_err, sigma_yy_err=sigma_yy_err,
        wall_s=wall_s,
    )


def load_existing(csv_path):
    done = set()
    if csv_path and os.path.exists(csv_path):
        with open(csv_path) as f:
            for row in csv.DictReader(f):
                done.add((row["contrast"], row["cap_requested"], row["ncell"]))
    return done


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--preset", choices=sorted(PRESETS), default=None)
    p.add_argument("--contrast", type=float, nargs="+", default=None)
    p.add_argument("--cap", type=int, nargs="+", default=None,
                    help="max_coarsening_level; -1 means uncapped")
    p.add_argument("--ncell", nargs="+", default=None, help="e.g. 32x64")
    p.add_argument("--section", default="2d-pressure-1Pa",
                   help="#@ section of tests/LowMachElasticPressure/input to read base args from")
    p.add_argument("--exe", default="lowmach")
    p.add_argument("--dim", type=int, default=2)
    p.add_argument("--comp", default="clang++")
    p.add_argument("--verbose", type=int, default=2,
                   help="elastic.solver.verbose; must be >=1 for SOLVER_STATS")
    p.add_argument("--timeout", type=int, default=600)
    p.add_argument("--jobs", type=int, default=1)
    p.add_argument("--outdir", default=None)
    p.add_argument("--csv", default=None)
    p.add_argument("--rst", default=None)
    p.add_argument("--keep-plotfiles", action="store_true")
    p.add_argument("--dry-run", action="store_true")
    args = p.parse_args()

    if args.preset:
        cfg = PRESETS[args.preset]
        contrasts = args.contrast or cfg["contrasts"]
        caps = args.cap if args.cap is not None else cfg["caps"]
        ncells = args.ncell or cfg["ncells"]
    else:
        contrasts = args.contrast or README_CONTRASTS
        caps = args.cap if args.cap is not None else [-1, 2]
        ncells = args.ncell or ["32x64"]

    base_cfg = solverlib.read_test_config(INPUT, args.section)
    base_args = base_cfg["args"]

    exestr = exe_path(args.exe, args.dim, args.comp)
    if not args.dry_run and not os.path.isfile(exestr):
        print(f"error: executable not found: {exestr}\n"
              f"(build it first: make -j8 bin/{args.exe}-{args.dim}d-{args.comp})",
              file=sys.stderr)
        sys.exit(1)

    outdir = os.path.abspath(args.outdir or os.path.join(REPO_ROOT, "scratch", "solver_benchmark"))
    os.makedirs(outdir, exist_ok=True)

    combos = [(c, cap, n) for c in contrasts for cap in caps for n in ncells]
    done = load_existing(args.csv)
    todo = [(c, cap, n) for (c, cap, n) in combos
            if (str(c), str(cap), n) not in done]

    print(f"{len(combos)} combinations ({len(todo)} to run, {len(done)} already in {args.csv})")

    rows = []
    if args.csv and os.path.exists(args.csv):
        with open(args.csv) as f:
            rows = list(csv.DictReader(f))

    fieldnames = ["contrast", "ncell", "cap_requested", "status", "n_iter",
                  "resid0", "residN", "convfactor", "disp_y_err", "sigma_yy_err", "wall_s"]

    for (contrast, cap, ncell) in todo:
        if args.dry_run:
            print(f"[dry-run] contrast={contrast} cap={cap} ncell={ncell}")
            continue
        print(f"running contrast={contrast:g} cap={cap} ncell={ncell} ...", end=" ", flush=True)
        row = sweep_row(exestr, INPUT, base_args, contrast, cap, ncell, outdir,
                         args.timeout, args.keep_plotfiles, args.verbose)
        print(f"{row['status']} iters={row['n_iter']} ({row['wall_s']:.1f}s)")
        rows.append(row)
        if args.csv:
            with open(args.csv, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=fieldnames)
                w.writeheader()
                w.writerows(rows)

    if not args.dry_run:
        print()
        print(" | ".join(fieldnames))
        for row in rows:
            print(" | ".join(str(row[k]) for k in fieldnames))

    if args.rst and rows:
        with open(args.rst, "w") as f:
            f.write(".. list-table::\n   :header-rows: 1\n\n")
            f.write("   * - " + "\n     - ".join(fieldnames) + "\n")
            for row in rows:
                f.write("   * - " + "\n     - ".join(str(row[k]) for k in fieldnames) + "\n")
        print(f"wrote {args.rst}")


if __name__ == "__main__":
    main()
