# Host-only numerical-kernel candidate evidence

Classification: current-target inspection evidence, not a successful port.
Source revision: `13342c7cdf3141d86011fa672ee26c658fb23a7f`.

## Hydro/Riemann candidate

`src/Integrator/Hydro.cpp` has a device launch that calls
`riemannsolver->Solve(...)`, uses Gas methods, and contains a host-style
try/catch diagnostic path. `src/Integrator/Hydro.H` owns the solver through a
base pointer. `src/Solver/Local/Riemann/Riemann.H` declares virtual `Solve`, and
Roe/HLLC/HLLE implementations are unannotated, branch-heavy numerical
algorithms that call Gas thermodynamic methods. Gas and its nested models own
variable-sized host containers.

This evidence supports a draft work class: the solver is a numerical sub-port,
not an annotation hunk. It does not select the eventual Riemann algorithm,
dispatch representation, kernel boundary, or tolerances.

## Fracture candidate

`src/Integrator/Fracture.H` device launches index vector-backed crack/material
collections and call crack-model methods; crack interfaces include virtual
methods. Another launch compound-writes `crack.driving_force_norm`, which needs
reduction/ownership adjudication.

This evidence demonstrates the same generic risks—container ownership, runtime
dispatch, deep callability, branching, and aggregate writes—without claiming
that the Hydro solution transfers unchanged.

## Required proof

A consuming port records the full compiler call chain, mathematical invariants,
iteration/termination behavior, branch map, device-safe state representation,
error policy, reduction identities, and `GPU_NATIVE_SHAPE.md` resource evidence.
Until strict builds and port oracles pass, both examples remain blocked source
findings.

