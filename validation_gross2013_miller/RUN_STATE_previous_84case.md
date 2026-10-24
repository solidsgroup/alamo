# Continuation record

The user explicitly selected local execution on this machine. Do not ask about
a cluster again. Cheap launch agents manage simulations; root owns scientific
postprocessing and interpretation. Frozen solid properties remain identical to
the selected Chen calibration. The reactive fixed-coefficient gas-model
approximation and missing Gross composition-dependent inputs are documented in
README.md.

Current executable SHA256:
`843796b9b7aac079874a3b2e5f49fd5ff4ccb1cc908eec94d5410a6f674243d4`.
An immutable copy is in `binaries/`. Use `launch_case.py` for every new launch.

## Active jobs at 2026-09-11 18:44 UTC

| Case suffix | Agent | Agent-owned exec session | Physical stop |
|---|---|---:|---:|
| M03 dx1_w4_v5_restart | launch_q500 | 70917 | 3 ms |
| M03 dx0.5_w8_v5_restart | launch_q200 | 45232 | 3 ms |
| M03 dx0.5_w4_v5_restart | launch_q1000 | 35724 | 3 ms |
| M24_p67.3863_seed101_pilot_ignition | launch_q500 | 91573 | 3 ms |

The first three use the original prelaunch executable hash
`4fff085fa108ca4c66f577e18dff59880dd9f20425742db9b08f1dbe2556cabb`.
Their old manual launches need actual exit receipts from the owning agents.
Never substitute the newly built executable hash. The M24 launcher writes its
receipt automatically and uses the current executable.

Root-owned session 9654 runs `watch_pilots.py`. It refreshes `STATUS.md`,
per-case analysis histories, and `analysis/ignition_resolution_pilots.png` every
minute when new snapshots appear. It does not launch production. Avoid running
another analyzer concurrently against the same cases/cache. If the watcher
has stopped, inspect its log before restarting; its file lock prevents a
duplicate watcher.

## Completed in this continuation

- Fixed generic restart copies to communicate across differing MPI rank maps.
- Initialized restart cell storage consistently with fresh/nodal allocation.
- Verified exact initial saved fields on 1/4/8 ranks; the final allocation fix
  passed 100 steps on four ranks plus two restart initialization repeats.
- Compared evolved serial/parallel AMR fields. Maximum temperature difference
  is 1.6e-6 K; maximum solid-fraction difference is about 1.4e-9.
- Started the full M24 disk-pack ignition pilot on eight ranks.
- Preserved the failed fine-v3 receipt with executed binary hash explicitly
  unverified: the original digest was measured only after launch/rebuild.
- Added restart ancestry to recession histories and checked a smaller surface
  interpolation strip against the full covering grid (exact agreement).
- Regenerated only the 84 unlaunched production decks: effectively unlimited
  step ceiling and final plot output enabled. Input hashes all match metadata.
- Measured preliminary packed cost: 100 steps/49.51 microseconds cost 113 s on
  one rank, 56 s on four ranks, and 64 s on eight ranks. These separate runs
  include startup/output and shared-machine load; developed-flame cost is
  unknown. The user was told the full sweep may take months at this resolution.

## Outstanding scientific work

1. Finish M03 post-ignition and resolution checks. Coarse flame is still active
   after heating ends at 1 ms; fine pilots have not yet reached that cutoff.
   No validated rate or resolution result exists yet.
2. Finish packed ignition/AMR robustness pilot. Its 3 ms duration is not enough
   particle recession for a statistically adequate Figure 10 comparison.
3. Choose production grid/interface width from actual convergence evidence;
   verify developed-flame MPI behavior and use local batches with measured
   throughput/memory/disk needs. Approximately 105 GiB was free at restart.
4. Launch the production sweep through cheap agents only after the pilot
   findings support it. `production_cases.txt` lists 84 prepared cases; zero
   are launched or accepted as of this record.
5. Root must check sustained burning, time-window drift, sufficient recession
   through packed layers, and three-seed variability before final aggregation.
   `summarize_sweep.py` never treats a pilot rate as a comparison point. The
   final Figure 10 errors and plot are still outstanding.
