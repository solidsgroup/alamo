#!/bin/sh

if command -v llvm-cov >/dev/null 2>&1; then
	exec llvm-cov gcov "$@"
fi

clang_major=$(clang++ -dumpversion | cut -d. -f1)
if command -v "llvm-cov-$clang_major" >/dev/null 2>&1; then
	exec "llvm-cov-$clang_major" gcov "$@"
fi

echo "llvm-cov is required for Clang coverage" >&2
exit 127
