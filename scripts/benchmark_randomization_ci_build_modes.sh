#!/usr/bin/env bash
set -euo pipefail

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
pkg_dir="$repo_root/EDI"
bench_r="$repo_root/scripts/benchmark_randomization_ci_ordinal_ppo.R"
r_bin="${R_BIN:-R}"

cleanup_dirs=()
cleanup() {
  local dir
  for dir in "${cleanup_dirs[@]}"; do
    rm -rf "$dir"
  done
}
trap cleanup EXIT

run_mode() {
  local label="$1"
  shift

  local build_dir build_start build_end bench_start bench_end
  build_dir="$(mktemp -d "${TMPDIR:-/tmp}/edi-ci.${label//[^A-Za-z0-9_-]/_}.XXXXXX")"
  cleanup_dirs+=("$build_dir")

  printf '\n== %s ==\n' "$label"
  build_start="$(date +%s)"
  env "$@" "$r_bin" CMD INSTALL --preclean --no-multiarch -l "$build_dir" "$pkg_dir"
  build_end="$(date +%s)"
  printf '%s install time: %ss\n' "$label" "$((build_end - build_start))"

  bench_start="$(date +%s)"
  EDI_LIB="$build_dir" \
  EDI_LABEL="$label" \
  EDI_NUM_CORES=3 \
  EDI_R=201 \
  EDI_REPS=3 \
  "$r_bin" --vanilla "$bench_r"
  bench_end="$(date +%s)"
  printf '%s benchmark wall time: %ss\n' "$label" "$((bench_end - bench_start))"
}

run_mode "portable"
run_mode "native" EDI_NATIVE_SPEED=1
run_mode "native+lto" EDI_NATIVE_SPEED=1 EDI_NATIVE_LTO=1
