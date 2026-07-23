# Current-source advisory scan

Port ID: `hostile-current`
Source revision: `13342c7cdf3141d86011fa672ee26c658fb23a7f`
Output: `current-source-coverage.csv`

The deterministic repository scan contains 1,593 rows across 154 application
source files: 1,559 candidates and 34 lexically converted sites. This is task
evidence, not a completed-port ledger; its candidates have deliberately not
been dispositioned and it proves neither device reachability nor correctness.

Hydro.cpp has 50 advisory rows. They include GPU-002 declaration shapes,
GPU-004 member/indirect-call shapes, and GPU-013 rows at the `dt_max_handle`
writes on lines 466–469. Fracture.H has 55 rows, including GPU-002 container
declarations, GPU-004 member/indexed calls, and GPU-013 at the
`driving_force_norm +=` write on line 696. The scanner cannot establish that
these sites are defects; compiler-first closure and ledger disposition remain
authoritative.

GPU-031 is intentionally manual because a regex cannot prove numerical call
closure, termination, receiver ownership, or device suitability. Its current
Hydro/Riemann and Fracture evidence is recorded in
`docs/gpu_manual/evidence/host-only-numerical-kernels.md`.
