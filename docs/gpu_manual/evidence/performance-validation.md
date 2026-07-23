# Performance validation evidence

Classification: chamber-gpu corpus evidence. The records illustrate acceptable
evidence shape; they do not establish another port's baseline or optimization.
The normative, physics-agnostic contract is `../PERFORMANCE.md`.

- `docs/agent_plans/20260709-fapply-kernel-surgery/results/RESULT.md`: performance edits passed `GOLDEN_MODE=cpu benchmark/ci_golden_compare.sh`; sanitizer coverage was device-safety evidence only where the deck later aborted for an unrelated solver condition.
- `docs/agent_plans/20260709-fsmooth-launch-fusion/results/RESULT.md`: component-launch fusion retained the strict CPU golden result, supporting GPU-015's correctness gate.
- `docs/agent_plans/20260713-fapply-322b-a100/results/RESULT.md`: Fapply changes require profiler evidence plus regression/golden checks, supporting GPU-018 through GPU-020 constraints.
- `benchmark/GPU_TEST_PERF_TRACKING.md`: compare the same binary mode and case; strict/no-fast-math cases are correctness signals, while fast builds are throughput signals.

No timing claim is transferable without the exact command, binary, GPU, dimensionality, box layout, and golden result. A sanitizer pass does not establish convergence or performance.
