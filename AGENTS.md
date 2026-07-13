# Repository Guidelines

## Project Structure & Module Organization

Alamo is a C++20 AMReX-based scientific code. Core implementation lives in `src/`, organized by modules such as `Integrator/`, `Model/`, `Operator/`, `Solver/`, `IC/`, `IO/`, and `GPU/`. Main executables are built from top-level `src/*.cc` files into `bin/`. Unit-test code is under `src/Test/`, `src/Unit/`, and `src/test.cc`. Regression tests live in `tests/<TestName>/` with an `input`, executable `test` script, and optional `reference/` data. Documentation sources are in `docs/source/`; process and agent planning files are in `docs/llm/` and `docs/agent_plans/`. Utility scripts are in `scripts/`, benchmark harnesses in `benchmark/`, and profiling/report tooling in `analysis/`.

## Build, Test, and Development Commands

- `./configure`: configure a default 3D production build and fetch/configure AMReX as needed.
- `./configure --dim=2 --debug`: configure a 2D debug build.
- `make -j4`: compile the configured target; executables appear in `bin/`.
- `make test`: run tab checks, build docs, then run the regression suite.
- `scripts/runtests.py --dim=3 --serial tests/Unit`: run a focused serial unit test.
- `make check`: run documentation checks, tab checks, and `eclint check src`.
- `make docs`: build Doxygen/Sphinx documentation.
- `make clean` or `make realclean`: remove build outputs; use `realclean` after changing MPI.

## Coding Style & Naming Conventions

Follow `.editorconfig` and `.clang-format`: C++ uses 4-space indentation, no tabs, C++20, GNU-derived formatting, non-indented namespaces, and wrapped braces. Python uses 4 spaces; Makefiles use tabs. C++ headers generally use `.H`; implementations use `.cpp` or `.cc`. Match existing module naming: directories, classes, and test folders use descriptive PascalCase or domain names, while scripts use lowercase snake_case where already established.

## Testing Guidelines

Add regression tests as `tests/<Name>/input` plus a `test` checker and `reference/` files when output comparisons are needed. Add low-level C++ tests near related code in `src/Test/` or `src/Unit/` and wire them through `src/test.cc`. Before merging to `development`, run both 2D and 3D builds followed by `make test`.

## Commit & Pull Request Guidelines

Recent history uses scoped subjects such as `benchmark: ...` and `process: ...`, often with campaign tags in parentheses. Keep commits imperative, scoped, and specific. Pull requests should describe the change, list build/test commands run, link issues or task plans, and include plots, screenshots, or report paths for benchmark-visible changes.

## Agent-Specific Instructions

For fresh clones, enable hooks with `git config core.hooksPath .githooks`. Agent work should respect `CLAUDE.md`: start assigned tasks from `docs/llm/PLAN.md`, create templated plans under `docs/agent_plans/`, and do not edit archived docs or changelog entries in place.
