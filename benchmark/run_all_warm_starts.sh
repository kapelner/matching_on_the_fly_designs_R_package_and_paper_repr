#!/bin/bash
# Run all warm-start benchmarks for n=100,200,500,1000 with 100 timing reps,
# including all paths that are disabled by default.

set -euo pipefail
cd "$(dirname "$0")/.."

export WARM_START_BENCH_REPS=100
export WARM_START_BENCH_P=5
export WARM_START_BENCH_FIXED_N=1

for N in 100 200 500 1000; do
    echo "=== n=$N: run_benchmark_final (Rand, Boot, JK_unused, PB) ==="
    WARM_START_BENCH_N=$N \
    WARM_START_BENCH_B=10 \
    WARM_START_BENCH_R=10 \
    WARM_START_BENCH_J=5 \
    WARM_START_BENCH_PB=2 \
    WARM_START_BENCH_RESULTS="warm_starts_final_results_n${N}.csv" \
    Rscript benchmark/run_benchmark_final.R

    echo "=== n=$N: benchmark_bayes_jackknife (Bayesian Bootstrap, JK) ==="
    WARM_START_BENCH_N=$N \
    WARM_START_BENCH_BB=10 \
    WARM_START_BENCH_J=5 \
    WARM_START_BENCH_RESULTS="warm_starts_bayes_jackknife_results_n${N}.csv" \
    Rscript benchmark/benchmark_bayes_jackknife_warm_starts.R
done

echo "=== Assembling tables into warm_starts.md ==="
Rscript benchmark/assemble_warm_starts.R

echo "=== All done ==="
