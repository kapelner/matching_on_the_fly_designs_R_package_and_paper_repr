#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
pkg_dir="$repo_root/EDI"
r_bin="${R_BIN:-R}"

cleanup_dirs=()
cleanup() {
  local dir
  for dir in "${cleanup_dirs[@]}"; do
    rm -rf "$dir"
  done
}
trap cleanup EXIT

run_build() {
  local label="$1"
  shift

  local elapsed start end build_dir
  build_dir="$(mktemp -d "${TMPDIR:-/tmp}/edi-build.${label//[^A-Za-z0-9_-]/_}.XXXXXX")"
  cleanup_dirs+=("$build_dir")

  printf '\n== %s ==\n' "$label"
  start="$(date +%s)"
  env "$@" "$r_bin" CMD INSTALL --preclean --no-multiarch -l "$build_dir" "$pkg_dir"
  end="$(date +%s)"
  elapsed=$((end - start))
  printf '%s build time: %ss\n' "$label" "$elapsed"
}

run_build "portable"
run_build "native" EDI_NATIVE_SPEED=1
run_build "native+lto" EDI_NATIVE_SPEED=1 EDI_NATIVE_LTO=1
