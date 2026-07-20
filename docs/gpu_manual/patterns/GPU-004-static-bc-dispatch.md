# GPU-004: Remove virtual boundary-condition dispatch from kernels
Status: draft
Class: correctness
Recognizer: regex: `(?:\(\s*\*\s*m_bc\s*\)|m_bc\s*->\s*GetBC\s*\(\s*\))\s*\(`
Applies: A device lambda invokes a virtual or base-pointer boundary-condition accessor.
Transform:
  Before:
    `(*m_bc)(u, gradu, sigma, i, j, k, domain)` executes virtual BC dispatch in a device body.
  After:
    On host, copy `m_bc->GetBcTypeArray()` into POD state; in the GPU branch call static device-safe `Elastic::eval(bc_type, ...)`, while the CPU branch retains `(*m_bc)(...)`.
Constraints: Mandatory when recognizer matches. Keep host-only BC construction/parsing and virtual calls outside kernels; preserve CPU dispatch plus component, face, and ghost-cell semantics. This recognizer is exclusive to BC invocation, not Gas-family pointers or parser forwarding.
Verify: `make -j4`; compare CPU/GPU outputs at identical timestep and dimensions; expected nvcc emits no virtual-function/device-call error and boundary values remain unchanged in a Neumann/Dirichlet smoke case.
Failure modes: Any mismatch or compile diagnostic is a failed conversion. device compilation rejects vtable use, or a copied BC contains host-owned state and faults at runtime. A broad replacement can alter boundary ordering or dispatch the wrong concrete condition. Validation must include the smallest representative input and a CPU reference; do not waive a failure as numerical noise.
Evidence: commits 1da57132e730ebaf97014ff63bc74f0b72cbf870, 54a941433b7582578cb5d56db794b1d648fb03cc
