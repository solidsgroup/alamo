# GPU-015: Fuse component launches
Status: draft
Class: performance
Recognizer: regex: `for\s*\([^)]*\bn\s*<\s*ncomp[^)]*\)[\s\S]{0,160}ParallelFor\s*\([^,\n]+,\s*\[`
Applies: Component-wise AMReX work where each component launch repeats identical index work.
Transform:
  Before:
    `for (int n=0; n<ncomp; ++n) ParallelFor(box, [=](i,j,k) { f(i,j,k,n); });`
  After:
    `ParallelFor(box, ncomp, [=](i,j,k,n) { f(i,j,k,n); });`
Constraints: Keep scope narrow. Apply only with profiling justification; preserve component ordering, synchronization, and CPU/golden results. Do not fuse kernels with different bounds, dependencies, or reductions.
Verify: `make -j4`; run the affected MLMG/operator regression and compare golden output; expect build success, unchanged numerical checks, and measured launch/time improvement. Repeat in the configured 2-D and 3-D modes where applicable, and retain the command/output in the build log.
Failure modes: Wrong component indexing can corrupt neighboring components; altered launch geometry can change races or golden values. Review the failing signature before changing the pattern; do not broaden its scope to silence an unrelated failure.
Evidence: chamber-gpu commit dc02baf1779b6668a6ea5725ade2e79b22118dfb; `src/Operator/Operator.cpp:116-150`.
