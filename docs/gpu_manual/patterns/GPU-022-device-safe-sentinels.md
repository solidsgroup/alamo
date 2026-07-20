# GPU-022: Replace host globals with device-safe sentinel access
Status: draft
Class: correctness
Recognizer: regex: `return\s+Set::Garbage\b|static\s+.*Garbage\s*=|extern.*Garbage`
Applies: Device-callable math paths returning a host-only global sentinel.
Transform:
  Before:
    `return Set::Garbage;` from a function reachable by a GPU kernel.
  After:
    Use a host/device `Garbage()` accessor backed by device and host storage selected at compile time.
Constraints: Keep scope narrow. Mandatory when recognizer matches. Sentinel remains diagnostic only; preserve call-site behavior and avoid dereferencing host storage in device code.
Verify: `./configure --dim=2 --debug && make -j4 && scripts/runtests.py --dim=2 --serial tests/Unit`; repeat with `--dim=3`; expect no device-link/illegal-address error and unchanged matrix checks.
Failure modes: Host global capture fails nvcc/device linking or faults at runtime; changing sentinel initialization masks invalid-state diagnostics. Review the failing signature before changing the pattern; do not broaden its scope to silence an unrelated failure.
Evidence: chamber-gpu commit f296eae5d31d39d2bf0bbfefdfddc5ea8142b0da; `src/Set/Base.H:114-120`, `src/Set/Matrix4_Major.H`, `src/Set/Matrix4_MajorMinor.H`.
