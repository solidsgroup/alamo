#!/usr/bin/env bash
set -eu

postfix=ci
patch=".diff/${postfix}.diff.patch"
development_patch=".diff/${postfix}.development-diff.patch"

if [[ -n ${TMPDIR:-} && ${TMPDIR} != /* ]]; then
    export TMPDIR="${PWD}/${TMPDIR}"
fi
scratch=$(mktemp -d)
trap 'rm -rf -- "${scratch}"' EXIT
mkdir -p "${scratch}/source/"{src,docs/source}
cp Makefile "${scratch}/source/"
cd "${scratch}/source"

# Keep Git operations inside the test repository, even when called from a hook.
git_local_vars=$(git rev-parse --local-env-vars)
while IFS= read -r name; do
    unset "${name}"
done <<< "${git_local_vars}"
unset GIT_EXTERNAL_DIFF GIT_DIFF_OPTS
export GIT_CONFIG_GLOBAL=/dev/null GIT_CONFIG_NOSYSTEM=1 GIT_ATTR_NOSYSTEM=1

git init -q --template=
git config user.name 'CI patch check'
git config user.email 'ci@example.invalid'
git config core.attributesFile /dev/null
printf 'base\n' > staged.txt
printf 'base\n' > working.txt
printf '\000base' > data.bin
git add staged.txt working.txt data.bin
git commit -qm 'Development base'
git update-ref refs/remotes/origin/development HEAD

printf 'feature\n' > feature.txt
git add feature.txt
git commit -qm 'Feature commit'
printf 'staged edit\n' > staged.txt
git add staged.txt
staged_tree=$(git write-tree)
printf 'unstaged edit\n' > working.txt
printf '\000edited' > data.bin

"${MAKE:-make}" POSTFIX="${postfix}" "${patch}" "${development_patch}"
git add --update
expected_tree=$(git write-tree)

# Apply each patch to its starting tree and compare the tracked files.
replay="${scratch}/replay"
git worktree add -q --detach "${replay}" HEAD
git -C "${replay}" read-tree --reset -u "${staged_tree}"
git -C "${replay}" apply --index "${PWD}/${patch}"
replayed_tree=$(git -C "${replay}" write-tree)
test "${replayed_tree}" = "${expected_tree}"
git -C "${replay}" reset -q --hard origin/development
git -C "${replay}" apply --index "${PWD}/${development_patch}"
replayed_tree=$(git -C "${replay}" write-tree)
test "${replayed_tree}" = "${expected_tree}"

# A missing development ref must also clear any stale patch.
git read-tree "${staged_tree}"
git update-ref -d refs/remotes/origin/development
"${MAKE:-make}" POSTFIX="${postfix}" "${patch}" "${development_patch}"
test -s "${patch}"
test ! -e "${development_patch}"
echo 'Git patch checks passed'
