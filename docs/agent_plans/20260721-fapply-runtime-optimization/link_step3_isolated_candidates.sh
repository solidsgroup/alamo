#!/usr/bin/env bash
set -euo pipefail
set -x

REPO=/home/jackplum/Projects/alamo
TASK="$REPO/docs/agent_plans/20260721-fapply-runtime-optimization"
ART="$TASK/artifacts/a1000-sm86-2e6a8f8f-20260721"
BASE_OBJ="$ART/baseline/obj"
CAND="$ART/candidate-step3"
ISO="$CAND/isolated"

cd "$REPO"
test ! -e "$ISO"
mkdir -p "$ISO/bin" "$ISO/obj"

link_one () {
  local dim="$1" fp="$2" profile="$3" postfix="$4" output="$5"
  local frozen="$BASE_OBJ/obj-$postfix"
  local candidate="obj/obj-$postfix"
  local isolated="$ISO/obj/obj-$postfix"
  test -d "$frozen"
  test -f "$candidate/Operator/Elastic.cpp.o"
  test ! -e "$isolated"
  test ! -e "$output"

  cp -a "$frozen" "$isolated"
  chmod -R u+w "$isolated"
  cp "$candidate/Operator/Elastic.cpp.o" "$isolated/Operator/Elastic.cpp.o"

  local config=(--comp=g++ --dim "$dim" --cuda 86 --cuda-fp "$fp")
  if [[ "$profile" == 1 ]]; then config+=(--profile); fi
  ./configure "${config[@]}"

  local objects
  objects="$(find "$isolated" -type f -name '*.o' -print | sort | tr '\n' ' ')"
  make -s \
    --eval 'isolated-link:;$(LINK_CMD) -o $(ISO_OUT) $(ISO_OBJ) $(LIB) $(MPI_LIB) $(LINKER_FLAGS)' \
    ISO_OUT="$output" ISO_OBJ="$objects" isolated-link
  chmod +x "$output"
  sha256sum "$output" "$candidate/Operator/Elastic.cpp.o" \
    "$frozen/Operator/Elastic.cpp.o"
}

link_one 2 strict 0 2d-nofast-cuda86-g++ \
  "$ISO/bin/alamo_gpu-2d-strict-sm86-step3-isolated"
link_one 2 fast 0 2d-cuda86-g++ \
  "$ISO/bin/alamo_gpu-2d-fast-sm86-step3-isolated"
link_one 2 fast 1 2d-profile-cuda86-g++ \
  "$ISO/bin/alamo_gpu-2d-profile-fast-sm86-step3-isolated"
link_one 3 strict 0 3d-nofast-cuda86-g++ \
  "$ISO/bin/alamo_gpu-3d-strict-sm86-step3-isolated"
link_one 3 fast 0 3d-cuda86-g++ \
  "$ISO/bin/alamo_gpu-3d-fast-sm86-step3-isolated"
link_one 3 fast 1 3d-profile-cuda86-g++ \
  "$ISO/bin/alamo_gpu-3d-profile-fast-sm86-step3-isolated"

find "$ISO/bin" -maxdepth 1 -type f -print0 | sort -z | xargs -0 sha256sum \
  >"$ISO/BINARY_SHA256SUMS"
set +x
echo STEP3_ISOLATED_LINK_PASS
