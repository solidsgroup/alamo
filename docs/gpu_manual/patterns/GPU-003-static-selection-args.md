# GPU-003: Forward context through static selection
Status: draft
Class: correctness
Recognizer: regex: `static_polymorphism_parser<[^>\n]*>\s*\(\s*value\s*\)`
Applies: A tuple-based parser calls `Parse` without forwarding constructor/context arguments.
Transform:
  Before:
    `select(name, value)` -> `static_polymorphism_parser(value)` -> `Parse(obj, pp)`
  After:
    `select(name, value, args...)` forwards `std::forward<Args>(args)...` through recursion to `Parse`.
Constraints: Mandatory when recognizer matches. This recognizes only the old exact no-args call. Keep argument order and value lifetime when adding forwarding; do not introduce dynamic dispatch or host pointers into parsed values.
Verify: `make -j4`; compare CPU/GPU outputs at identical timestep and dimensions; expected parser compiles and selected models receive required species/geometry context.
Failure modes: Any mismatch or compile diagnostic is a failed conversion. template recursion drops arguments, selected model has default/uninitialized context, or nvcc reports an unresolved overload. The recognizer must stop matching once `std::forward<Args>` is present. Validation must include the smallest representative input and a CPU reference; do not waive a failure as numerical noise. Include launch-region behavior, ownership, and dimensional assumptions in review; the transform is complete only when those remain explicit.
Evidence: commits 3d382a5655c715c7483a15d53c3f8e5c908330af
