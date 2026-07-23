# GPU-025: Declare and enforce the supported GPU closure
Transform status: draft
Class: scaffolding
Detection: Advisory manual/build inspection in `recognizers/table.csv`; record selected integrator sources and object closure.
Invariant: A GPU target must contain every source and device-safe call edge required by the selected execution path; exclusion cannot prove correctness.
Port contract: Declare integrator scope, source/object closure, unsupported paths, closure record, and retirement trigger. Record assumptions and dispositions in the port ledger.
Transform:
  Before: A launcher compiles every integrator and relies on accidental compiler compatibility.
  After: Select a narrow supported closure and fail unsupported paths explicitly.
Corpus example: chamber-gpu Flame policy files evidence one closure instance, not the universal supported integrator.
Constraints: Follow the scaffolding-lifecycle policy in `ARCHITECTURE_POLICIES.md`: shrink quarantine/closure, never grow it silently; add sources only after correctness evidence.
Verify: Pass the `ONBOARDING.md` closure gate, then require `VALIDATION.md` `strict-build`, `multi-box`, and `sanitizer` rows with a written tolerance rationale.
Failure modes: Over-broad closure reintroduces errors; over-pruning causes links or missing physics. Review diagnostics before changing scope.
Evidence: Primary: `evidence/primary-sources.md` (CUDA compilation/link rules). Corpus: chamber-gpu commit `54a941433b7582578cb5d56db794b1d648fb03cc`; policy/source and guide paths.
