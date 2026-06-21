#!/usr/bin/env python3
"""Phase-3 memory-budget calculator (roadmap 3.2).

Given a per-node memory cost (bytes/node) and a ghost/AMR overhead model, print
the largest wide-shallow 3D grid that fits on each target GPU. This is a *planning*
tool: all numbers are first-order estimates and are labelled approximate. No GPU,
no CUDA, no SLURM -- pure stdlib Python.

bytes/node model (from Phase-3 recon, PLAN.md "Bytes/node figure"):
  - Phase-3 GPU runs DISABLE elastic (decision D1), so the dominant resident cost is
    the phase-field/thermal fields plus AMReX bookkeeping, NOT the elastic model fab.
  - elastic-DISABLED phase-field default ~= 256 B/node, broken down roughly as:
        eta + phi + temp                 ~ 3 doubles  = 24 B
        old/new state copies + scratch   ~ several doubles
        AMReX MultiFab / metadata overhead, alignment, BC scratch
    rounded conservatively to 256 B/node.
  - elastic-ENABLED (for comparison) ~= 512 B/node:
        Matrix4<3> elastic model fab     ~ 360 B/node
        displacement vector (3 doubles)  = 24 B
        RHS / residual / stress scratch  + the phase-field fields above
    rounded conservatively to 512 B/node.
  Pass --bytes-per-node to override; pass 512 to model the elastic-enabled config.

Ghost / AMR overhead:
  - --ghost-factor multiplies the base-level node count to account for ghost cells
    (a halo of ghost cells around every box; default 1.3, i.e. ~30% inflation).
  - shallow AMR (--max-level, default 1, per the wide-shallow strategy) adds finer
    levels that cover only a --refine-fill fraction of the domain. Each finer level
    has 8x the node density (2x per dimension in 3D) over the fraction it covers.

Grid aspect: wide-shallow (roadmap 3.3) -> Nx = Ny = 2 * Nz. We solve for the
largest such grid whose estimated total device memory stays under a safety headroom
(default 85%) of the device RAM.
"""

import argparse
import sys

# Built-in device table: name -> memory in GB (GiB-style binary GB).
# These are nominal device capacities; usable memory is lower, hence --headroom.
DEVICES = [
    ("A1000", 8),
    ("A100-40", 40),
    ("A100-80", 80),
    ("H100-80", 80),
    ("H200", 141),
]

GB = 1024 ** 3  # bytes per GiB


def amr_node_multiplier(max_level, refine_fill):
    """Approx total-node multiplier over the base level for shallow AMR.

    Base level contributes factor 1.0. Each finer level l (1..max_level) covers
    `refine_fill` of the domain at 8x the node density of the level above it,
    so it adds refine_fill * 8**l nodes relative to the base level.
    This is a deliberately simple, conservative first-order model.
    """
    mult = 1.0
    for level in range(1, max_level + 1):
        mult += refine_fill * (8.0 ** level)
    return mult


def total_bytes_for_base(nz, bytes_per_node, ghost_factor, amr_mult):
    """Estimate total device bytes for a wide-shallow base grid with Nz = nz.

    Wide-shallow aspect: Nx = Ny = 2*Nz. base nodes = (2*nz)*(2*nz)*nz = 4*nz**3.
    """
    base_nodes = 4 * (nz ** 3)
    effective_nodes = base_nodes * ghost_factor * amr_mult
    return effective_nodes * bytes_per_node, base_nodes, effective_nodes


def largest_grid(mem_bytes, bytes_per_node, ghost_factor, amr_mult, headroom):
    """Largest Nz (wide-shallow) whose estimated memory fits under headroom*mem."""
    budget = mem_bytes * headroom
    nz = 0
    last = None
    # Grow Nz until we exceed the budget; step in multiples of 32 to respect a
    # blocking_factor=32 wide-shallow strategy (round-ish, planning granularity).
    step = 32
    candidate = step
    while True:
        est_bytes, base_nodes, eff_nodes = total_bytes_for_base(
            candidate, bytes_per_node, ghost_factor, amr_mult
        )
        if est_bytes > budget:
            break
        nz = candidate
        last = (nz, base_nodes, eff_nodes, est_bytes)
        candidate += step
    if last is None:
        # Even the smallest step does not fit; report the smallest step anyway.
        est_bytes, base_nodes, eff_nodes = total_bytes_for_base(
            step, bytes_per_node, ghost_factor, amr_mult
        )
        return step, base_nodes, eff_nodes, est_bytes, False
    return last + (True,)


def human_bytes(n):
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.2f} {unit}"
        n /= 1024.0


def build_parser():
    p = argparse.ArgumentParser(
        prog="phase3_memory_budget.py",
        description="Phase-3 memory-budget calculator: largest 3D wide-shallow grid "
        "that fits per GPU (planning estimate, roadmap 3.2).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="All figures are FIRST-ORDER ESTIMATES for planning only. "
        "Default 256 B/node models the elastic-DISABLED phase-field config (D1); "
        "pass --bytes-per-node 512 for the elastic-enabled comparison.",
    )
    p.add_argument(
        "--bytes-per-node",
        type=float,
        default=256.0,
        help="Resident device bytes per grid node (default: 256, elastic-disabled "
        "phase-field config). Use 512 for the elastic-enabled figure.",
    )
    p.add_argument(
        "--ghost-factor",
        type=float,
        default=1.3,
        help="Multiplier for ghost-cell halo overhead (default: 1.3).",
    )
    p.add_argument(
        "--max-level",
        type=int,
        default=1,
        help="AMR max_level for the wide-shallow strategy (default: 1).",
    )
    p.add_argument(
        "--refine-fill",
        type=float,
        default=0.1,
        help="Fraction of the domain covered by each finer AMR level (default: 0.1).",
    )
    p.add_argument(
        "--headroom",
        type=float,
        default=0.85,
        help="Usable fraction of device RAM to target (default: 0.85).",
    )
    p.add_argument(
        "--device",
        type=str,
        default=None,
        help="Restrict output to one device by name (e.g. A100-80). "
        "Default: print all devices.",
    )
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)

    if args.bytes_per_node <= 0:
        print("error: --bytes-per-node must be positive", file=sys.stderr)
        return 2
    if not (0 < args.headroom <= 1.0):
        print("error: --headroom must be in (0, 1]", file=sys.stderr)
        return 2

    devices = DEVICES
    if args.device is not None:
        wanted = args.device.strip().lower()
        devices = [d for d in DEVICES if d[0].lower() == wanted]
        if not devices:
            names = ", ".join(d[0] for d in DEVICES)
            print(
                f"error: unknown device '{args.device}'. Known: {names}",
                file=sys.stderr,
            )
            return 2

    amr_mult = amr_node_multiplier(args.max_level, args.refine_fill)

    print("Phase-3 memory budget -- largest wide-shallow 3D grid per GPU "
          "(APPROXIMATE / planning only)")
    print()
    print(f"  bytes/node     : {args.bytes_per_node:g}  "
          f"({'elastic-disabled (D1)' if args.bytes_per_node <= 300 else 'elastic-enabled regime'})")
    print(f"  ghost factor   : {args.ghost_factor:g}")
    print(f"  AMR max_level  : {args.max_level}   refine-fill: {args.refine_fill:g}  "
          f"-> AMR node multiplier ~ {amr_mult:.3f}x")
    print(f"  headroom       : {args.headroom:g} of device RAM")
    print(f"  grid aspect    : wide-shallow, Nx = Ny = 2*Nz (blocking_factor=32 steps)")
    print()

    header = (
        f"{'device':<10} {'mem':>8} {'max Nx x Ny x Nz':>22} "
        f"{'#nodes(base)':>14} {'est. bytes':>12} {'fit':>4}"
    )
    print(header)
    print("-" * len(header))

    for name, gb in devices:
        mem_bytes = gb * GB
        nz, base_nodes, eff_nodes, est_bytes, fit = largest_grid(
            mem_bytes,
            args.bytes_per_node,
            args.ghost_factor,
            amr_mult,
            args.headroom,
        )
        nx = ny = 2 * nz
        grid_str = f"{nx} x {ny} x {nz}"
        fit_str = "ok" if fit else "TIGHT"
        print(
            f"{name:<10} {str(gb) + ' GB':>8} {grid_str:>22} "
            f"{base_nodes:>14,d} {human_bytes(est_bytes):>12} {fit_str:>4}"
        )

    print()
    print("Notes (APPROXIMATE):")
    print("  * #nodes(base) is the base-level node count; est. bytes includes the")
    print("    ghost factor and the AMR multiplier across finer levels.")
    print("  * Nz is rounded to a multiple of 32 (wide-shallow blocking_factor).")
    print("  * 'TIGHT' means even the smallest 32-node grid exceeds the headroom budget.")
    print("  * Override the model with --bytes-per-node / --ghost-factor / "
          "--max-level / --refine-fill.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
