# Port validation contract

This is a port-supplied contract, not a fixed Flame test deck. Every port copies
`templates/VALIDATION.csv`, fills one row for every category below, and records
the exact executable, input, hardware, result, evidence path, and tolerance
rationale. A category may be `not-applicable` only with a physics or scope
reason and an owner. `not-run` and `blocked` never count as pass.

## Required categories

| Category | Contract |
|----------|----------|
| strict-build | Compile the declared closure with strict device diagnostics in every supported dimension and preserve the CPU build. |
| analytic-exact | Compare with an analytic or exact solution when one exists; otherwise record why none applies and name the stronger substitute. |
| conservation | Check every quantity the declared physics conserves; mark only non-conserved quantities not-applicable. |
| golden-regression | Compare the same physics, initial data, timestep policy, and output quantities against an approved CPU or pre-port reference. |
| restart-parity | Compare uninterrupted and restart/resume trajectories at the same final time. |
| multi-box | Compare a representative one-box and multi-box/decomposition run, including asynchronous and reduction paths. |
| sanitizer | Run the representative multi-box path under the selected backend's memory and synchronization checker. |

A correctness pattern's `Verify` field selects categories from this table and
adds its pattern-specific observation. It never substitutes a branch-specific
command for the contract.

## Tolerance rationale

Every numerical comparison records the quantity and units, norm, absolute and
relative tolerances, reference scale, precision/build mode, timestep and mesh,
reason for the bound, and evidence that the bound detects a known bad result.
An unexplained coarse field tolerance is not a correctness oracle. Exact
integer, enum, topology, and restart-metadata comparisons use exact equality.

## Gate

The validation gate passes only when every applicable row is `pass`, every
`not-applicable` row has an owner-approved rationale, and the scope and closure
identifiers match the tested binary. Analytic/exact, conservation, and golden
checks complement one another; a coarse golden comparison cannot silently
waive a failed invariant. Compiler diagnostics and recorded inspection outrank
scanner output.

## Pilot status

The historical PFC and HeatConduction closed-book runs in
`build/phase6/GATE_RUNS.csv` demonstrate pattern application outside Flame, but
they are not completed port-validation instantiations. Because an
integrator-specific oracle deck is outside this task, the required first
non-Flame filled contract remains `pending-pilot` and must be supplied by the
next authorized port.

