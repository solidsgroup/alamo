#!/usr/bin/env bash
#
# package_phase3_nova_bundle.sh
#
# Bundle the NOVA Phase 3 GPU/CPU comparison artifacts into a tarball and
# transfer it to a local machine via scp.
#
# Usage:
#   SCP_DEST=user@host:/path/to/dir bash analysis/package_phase3_nova_bundle.sh
#
# Optional environment variables:
#   BUNDLE_DIR   bundle directory to package. If unset, the script picks the
#                newest analysis/phase3_nova_bundle_* directory under the
#                repo root; if none exist, it falls back to the repo root.
#   OUT_DIR      where to write the tarball locally before scp (default: /tmp)
#   TAR_NAME     tarball filename (default: phase3_nova_bundle_<timestamp>.tar.gz)
#   SCP_DEST     required scp destination, e.g. jackplum@kermit:/home/jackplum/Downloads/
#   DRY_RUN      set to 1 to create the tarball but skip scp
#
set -euo pipefail

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"
ts="$(date +%Y%m%d_%H%M%S)"
out_dir="${OUT_DIR:-/tmp}"
tar_name="${TAR_NAME:-phase3_nova_bundle_${ts}.tar.gz}"
tar_path="${out_dir%/}/${tar_name}"
scp_dest="${SCP_DEST:-}"
dry_run="${DRY_RUN:-0}"

if [[ ! -w "$out_dir" ]]; then
  echo "error: OUT_DIR is not writable: $out_dir" >&2
  exit 1
fi

if [[ -n "${BUNDLE_DIR:-}" ]]; then
  bundle_dir="$BUNDLE_DIR"
else
  mapfile -t bundle_candidates < <(
    find "$repo_root/analysis" -maxdepth 1 -mindepth 1 -type d -name 'phase3_nova_bundle_*' -printf '%T@ %p\n' 2>/dev/null | sort -nr | awk '{print $2}'
  )
  if [[ "${#bundle_candidates[@]}" -gt 0 ]]; then
    bundle_dir="${bundle_candidates[0]}"
  else
    bundle_dir="$repo_root"
  fi
fi

if [[ ! -d "$bundle_dir" ]]; then
  echo "error: bundle directory does not exist: $bundle_dir" >&2
  exit 1
fi

stage_dir="$(mktemp -d "${out_dir%/}/phase3_nova_stage.XXXXXX")"
trap 'rm -rf "$stage_dir"' EXIT

declare -a files=()

stage_file() {
  local rel="$1"
  local src="$bundle_dir/$rel"
  local dst="$stage_dir/$rel"
  if [[ -e "$src" ]]; then
    mkdir -p "$(dirname "$dst")"
    cp -a "$src" "$dst"
    files+=("$rel")
  fi
}

stage_manifest() {
  {
    echo "bundle_dir=$bundle_dir"
    echo "timestamp=$ts"
    echo "hostname=$(hostname)"
    echo
    if git -C "$bundle_dir" rev-parse --is-inside-work-tree >/dev/null 2>&1; then
      echo "git rev-parse HEAD:"
      git -C "$bundle_dir" rev-parse HEAD
      echo
      echo "git status --short:"
      git -C "$bundle_dir" status --short || true
    else
      echo "git rev-parse HEAD:"
      echo "(not a git checkout)"
    fi
    echo
    echo "captured files:"
    for f in "${files[@]}"; do
      echo "$f"
    done
  } > "$stage_dir/manifest.txt"
  files+=("manifest.txt")
}

for rel in \
  commands.txt \
  quick_summary.txt \
  RESULTS_ANALYSIS.md \
  build_alamo_gpu_3d.sbatch \
  rq.sh \
  benchmark/build_alamo_nova_3d.sh \
  benchmark/phase3_nova_oneshot.sh \
  benchmark/phase3_scaling_sweep.sh \
  benchmark/nova_flame_gpu_3d.slurm \
  benchmark/nova_flame_gpu_3d_multi.slurm \
  benchmark/nova_flame_cpu_3d.slurm \
  benchmark/PHASE3_3D_READINESS.md \
  benchmark/PHASE3_R3_crossover.md
do
  stage_file "$rel"
done

while IFS= read -r -d '' path; do
  path="${path#./}"
  stage_file "$path"
done < <(find "$bundle_dir" -maxdepth 1 -type f \( -name '*.err' -o -name '*.out' \) -print0 | sort -z)

if [[ "${#files[@]}" -eq 0 ]]; then
  echo "error: no files found to package in $bundle_dir" >&2
  exit 1
fi

stage_manifest
tar -czf "$tar_path" -C "$stage_dir" .

echo "bundle_dir: $bundle_dir"
echo "created: $tar_path"
echo "contents:"
printf '  %s\n' "${files[@]}"

if [[ "$dry_run" == "1" ]]; then
  echo "DRY_RUN=1, skipping scp"
  exit 0
fi

if [[ -z "$scp_dest" ]]; then
  echo "error: SCP_DEST is required, for example:" >&2
  echo "  SCP_DEST=jackplum@10.24.220.162:/home/jackplum/Downloads/" >&2
  exit 1
fi

scp "$tar_path" "$scp_dest"
echo "copied to: $scp_dest"
