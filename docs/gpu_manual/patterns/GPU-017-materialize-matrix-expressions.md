# GPU-017: Materialize chained matrix expressions
Status: draft
Class: correctness
Recognizer: regex: `\.inverse\(\)\.transpose\(\)|\.cross\(\)\.normalized\(\)`
Applies: Eigen expression-template chains evaluated inside GPU device code.
Transform:
  Before:
    `Eigen::Matrix3d FinvT = F.inverse().transpose();`
  After:
    `Eigen::Matrix3d Finv = F.inverse();` then `Eigen::Matrix3d FinvT = Finv.transpose();`
Constraints: Keep scope narrow. Mandatory when the recognizer matches. Preserve operation order and scalar type; do not replace a chain with a mathematically different expression.
Verify: `make -j4`; run the 3-D elastic solve/regression; expect nvcc compilation and no CUDA launch failure or changed golden values. Repeat in the configured 2-D and 3-D modes where applicable, and retain the command/output in the build log.
Failure modes: Leaving nested expressions can compile on CPU but fault on device (often CUDA error 719); changing evaluation order can alter constitutive results. Review the failing signature before changing the pattern; do not broaden its scope to silence an unrelated failure.
Evidence: chamber-gpu commit a5c1b2ddfd63e09848d69980d83496654d490ceb; `src/Model/Solid/Finite/NeoHookean.H:36-90#2`.
