# GPU-002: Replace runtime model pointers with static value dispatch
Status: draft
Class: correctness
Recognizer: regex: `(?:Thermo::Thermo|Transport::Transport|EOS::EOS)\s*\*\s*[A-Za-z_]\w*`
Applies: Gas/model state stores polymorphic base pointers used from a kernel.
Transform:
  Before:
    `Thermo::Thermo* thermo; thermo->cp_mol(...)`
  After:
    `Thermo::Thermo<CpConstant> thermo; thermo.cp_mol(...)` with concrete tuple-held selection where alternatives are required.
Constraints: Mandatory when recognizer matches. Preserve model selection, species data, and ownership; convert only runtime Gas-family base-pointer members. Do not apply to genuinely host-only plugin lifetimes or unrelated raw pointers.
Verify: `make -j4`; compare CPU/GPU outputs at identical timestep and dimensions; expected GPU target compiles without virtual-call or pointer-capture errors and model values match CPU smoke output.
Failure modes: Any mismatch or compile diagnostic is a failed conversion. nvcc rejects vtable use, a device copy retains a host address, or selected model state is lost. A missing concrete alternative can silently select the wrong tuple member. Validation must include the smallest representative input and a CPU reference; do not waive a failure as numerical noise.
Evidence: commits 3d382a5655c715c7483a15d53c3f8e5c908330af
