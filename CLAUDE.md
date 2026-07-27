# alamo — working guidelines

## Pull requests
- **Never open pull requests.** Commit and push branches only; leave PR creation to me.

## Building
- This is a C++ codebase. After **any** code change, recompile with `make -j8`.
- A change is not done until it compiles.

## Branching
- When planning a project, **ask which branch** the changes should go on, or whether
  to create a new branch — do not assume the current branch.

## Commits
- When planning a project, place **commits at logical checkpoints** in the plan so
  work lands in reviewable, self-contained increments.
