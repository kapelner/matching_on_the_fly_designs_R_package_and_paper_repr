#!/usr/bin/env bash
set -euo pipefail

REPS="${REPS:-5}"
OUT_DIR="${OUT_DIR:-benchmark/simd_matrix}"

Rscript scripts/benchmark_simd_matrix.R --reps="${REPS}" --out_dir="${OUT_DIR}"
