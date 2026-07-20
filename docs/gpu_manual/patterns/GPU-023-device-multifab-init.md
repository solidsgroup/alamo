# GPU-023: Initialize device-backed MultiFabs on device
Status: draft
Class: correctness
Recognizer: regex: `\b[A-Za-z_]\w*(?:\[[^]]+\])?\.(?:setVal|mult|plus)\s*\(`
Applies: MultiFab operations whose data may reside in the device arena during GPU kernels.
Transform:
  Before:
    `fab.setVal(0.0);` or algebra using implicit execution mode.
  After:
    `fab.setVal<amrex::RunOn::Device>(0.0);` and explicit device execution for dependent `mult`/`plus` operations.
Constraints: Keep scope narrow. Mandatory when a recognizer hit is confirmed in a device-backed GPU path; receiver type and reachability require inspection. Keep operation ordering, component ranges, and synchronization unchanged.
Verify: `make -j4`; run elastic/operator GPU golden and sanitizer checks; expect no invalid access and unchanged numerical output. Repeat in the configured 2-D and 3-D modes where applicable, and retain the command/output in the build log.
Failure modes: Host execution on device memory can fail or silently leave stale values; wrong execution mode creates race-dependent residuals. Review the failing signature before changing the pattern; do not broaden its scope to silence an unrelated failure.
Evidence: chamber-gpu commit d522e1ac08729b306a215ad896d8a304983d55de; `src/Operator/Elastic.cpp:417-425`, `src/Operator/Operator.cpp:50-85`.
