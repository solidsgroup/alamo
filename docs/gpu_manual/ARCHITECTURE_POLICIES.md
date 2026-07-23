# Cross-cutting architecture policies

These policies are the single normative homes for recurring decisions. Tier 1
patterns point here and do not restate competing policies. The port owner records
each decision in the scope or inspection artifact.

## Device-side error propagation

Owner decision: the port owner defines a compact device-safe error code and the
host observation boundary. A leaf or deep templated callee records the first
error in device storage and returns device-safe status; intermediate callees
propagate status without logging, throwing, allocating, or dereferencing host
state. The launch owner retains the flag through all contributing streams,
synchronizes at the host observation/lifetime boundary, translates the code and
captured indices to one host diagnostic, then resets or destroys the flag.

Worked example: `model.eval()` detects non-finite state and atomically records
`{nonfinite, level, i, j, k}`; `flux()` and a templated stencil return early;
after the MFIter launch group, the host checks once and emits the full message.
GPU-012 is the direct-kernel instance of this policy, not its limit.

## Scaffolding lifecycle

Owner decision: every quarantine or narrowed closure records owner, reason,
introduced date, affected scope, observable fallback behavior, removal trigger,
and retirement gate. Scaffolding may bound incremental work but must not hide a
required production path or silently turn hot work into a CPU fallback. New
work cannot enlarge an existing quarantine without an explicit scope decision.

Worked example: an unsupported diagnostic-only plugin is excluded from the
first closure with a test that fails if selected; its removal trigger is a
device-callable replacement plus the full validation contract. At port close it
is removed or carried as a named limitation with a new owner and review date.
GPU-010 and GPU-025 implement this policy at loop and build-closure boundaries.

## FEATURE and [NUM] surfacing

Owner decision: when conversion exposes a capability or numerical choice, stop
the affected ledger row, preserve existing behavior, and add a decision record
using `templates/FEATURE_DECISIONS.csv`. A FEATURE requires task-level scope
opt-in. `[NUM]` additionally requires the validation quantities and tolerance
rationale to be approved before implementation. Rejected and deferred choices
stay visible; neither may be smuggled in as device compatibility work.

Worked example: an FFT wrapper exposes a tempting change to which AMR levels
advance. The device-safe wrapper replacement proceeds with existing scheduling;
the AMR scheduling change is recorded `[NUM]`, owner-deferred, and excluded from
the conversion oracle. `FEATURES.md` and `ONE_OFFS.md` remain the chamber-gpu
corpus ledgers, not universal backlogs.

