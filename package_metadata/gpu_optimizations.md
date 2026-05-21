# Potential GPU Optimization Targets

Date: 2026-05-21

## Scope

This report identifies places in the EDI package where GPU acceleration could plausibly speed up computations. It is based on the current native C++ and R code paths, the existing OpenMP/SIMD build setup, and the kernel timing notes in `package_metadata/perf_experiments_final.md`.

The package already has meaningful CPU optimization infrastructure:

- `EDI/src/Makevars` enables OpenMP and native SIMD by default, with optional `-O3` and LTO controls.
- Many resampling kernels are already parallelized across permutations or bootstrap replicates with OpenMP.
- The main fitted-model kernels use Eigen, BLAS/LAPACK, hand-written likelihood loops, and several specialized C++ implementations.

GPU work should therefore be selective. The best candidates are workloads with many independent replicates, large dense matrix operations, or embarrassingly parallel pairwise computations. Small `n`, small `p`, branch-heavy exact tests, R callback loops, and short optimizer paths are poor GPU targets because transfer overhead and kernel-launch overhead will dominate.

## Highest-Priority Candidates

### 1. Resampling and randomization distributions

Relevant files:

- `EDI/src/ols_distr_parallel.cpp`
- `EDI/src/fast_wilcox_parallel.cpp`
- `EDI/src/fast_kk_wilcox_parallel.cpp`
- `EDI/src/kk_compound_distr_parallel.cpp`
- `EDI/src/ridit_distr_parallel.cpp`
- R callers through `compute_rand_two_sided_pval()`, `approximate_randomization_distribution_beta_hat_T()`, and `approximate_bootstrap_distribution_beta_hat_T()` across the inference classes.

Why this is attractive:

- These kernels already expose the natural GPU dimension: one permutation, bootstrap replicate, or simulated assignment per output column.
- The CPU code parallelizes over replicates, which is also the GPU mapping: one block or warp group per replicate, with reductions over subjects.
- Large `B` / `r` values are common in bootstrap and randomization inference, so fixed data transfer costs can be amortized.

Concrete opportunities:

- `compute_ols_distr_parallel_cpp()` computes a treatment coefficient for every assignment column in `w_mat`. A GPU version could batch the per-replicate sufficient statistics (`sum_w`, `X' w`, `X' y_sim`, `w' y_sim`) and solve many small linear systems. For small `p`, the best GPU approach is likely a custom batched small-matrix solver or cuSOLVER batched QR/Cholesky.
- `compute_ols_bootstrap_parallel_cpp()` currently materializes `X_b`, `y_b`, and `w_b` per bootstrap replicate. GPU code could avoid materializing full bootstrap matrices by gathering rows from device-resident `X`, `y`, and `w`, then accumulating `X'X` and `X'y` per replicate.
- `compute_wilcox_distr_parallel_cpp()` has a very GPU-friendly fast path when `delta == 0`: ranks are precomputed once and each replicate only needs reductions of `ranks[i] * w[i]`. This is a strong first prototype because the computation is simple and numerically low-risk.
- `compute_matching_compound_distr_parallel_cpp()` is less regular but still replicate-parallel. It should only be considered after the simpler OLS and Wilcoxon GPU prototypes.

Expected payoff:

- High when `n` is at least several thousand or when `B` / `r` is in the thousands to tens of thousands.
- Low when users run the current default-scale `B = 501` or smaller with modest `n`, because OpenMP may already be faster once GPU transfer overhead is included.

Implementation notes:

- Add GPU as an optional backend, not a replacement for the current C++ path.
- Keep inputs column-major to match R/Eigen memory and minimize reshaping.
- Return exactly the same vector shape and `NA_REAL` handling as the CPU kernels.
- Expose a dispatch policy similar to the existing parallel/optimization dispatch controls rather than silently using a GPU.

### 2. Pair-distance matrix construction

Relevant file:

- `EDI/src/pair_dist_helpers.cpp`

Why this is attractive:

- `compute_pair_distance_matrix_cpp()` computes all pairwise weighted squared distances among `m` pair averages.
- The work is `O(m^2 p)` and each `(i, j)` output is independent.
- This maps cleanly to CUDA/OpenCL: one thread block per tile of the distance matrix, with the covariate dimension reduced inside the block.

Current performance context:

- The performance report identified `pair_distance_matrix` as a nontrivial native kernel in broader sweeps.
- A previous CPU rewrite already removed Rcpp accessor overhead and gave a large win, so the remaining opportunity is mostly from scaling to large `m`.

Expected payoff:

- High for large matching problems where `m` is large enough that the full distance matrix is expensive and reused.
- Low for small `m` or when memory bandwidth and transfer dominate.

Implementation notes:

- The full `m x m` matrix can be large. GPU acceleration should be paired with memory-aware tiling.
- A more ambitious variant would avoid materializing the full matrix and compute only nearest-neighbor or candidate-pair distances if the downstream matching logic permits it. That algorithmic change could beat both CPU and GPU full-matrix builds.

### 3. GLMM and frailty likelihood evaluation

Relevant files:

- `EDI/src/_glmm_engine.h`
- `EDI/src/fast_logistic_glmm.cpp`
- `EDI/src/fast_poisson_glmm.cpp`
- `EDI/src/fast_ordinal_clmm.cpp`
- `EDI/src/fast_ordinal_glmm.cpp`
- `EDI/src/fast_weibull_frailty.cpp`
- `EDI/src/fast_clogit_plus_glmm.cpp`

Why this is plausible:

- The hot work is repeated likelihood, gradient, and Hessian evaluation over groups and Gauss-Hermite quadrature nodes.
- The performance report shows these kernels are dominated by objective loops, `exp`/`log`/`log1p`, Eigen matrix-vector products, and grouped reductions.
- The structure is parallel over groups and quadrature nodes for a fixed parameter vector.

Why this is risky:

- Optimizers call the objective repeatedly with small-to-moderate matrices, so GPU launch and transfer costs can dominate.
- The current implementations have many numerically sensitive branches, soft barriers, log-sum-exp reductions, and model-specific derivative paths.
- If the parameter vector, design matrix, responses, group layout, and quadrature rule are not kept resident on the GPU across optimizer iterations, the GPU path will probably lose.

Best use case:

- Large KK/GLMM fits with many groups, many observations per fit, and repeated resampling where the same design matrix is used many times with different weights or simulated responses.

Implementation notes:

- Do not start with the generic `_glmm_engine.h` path. Prototype on one specialized class first, probably `LogisticGLMMObjective`, because the math is simpler than ordinal CLMM or Weibull frailty.
- Keep the optimizer on the CPU initially, but keep data on the GPU and evaluate objective/gradient on-device.
- Return only objective value and gradient per iteration, not large intermediate matrices.
- Consider GPU acceleration only for estimate-only paths first. Hessian paths are more complex and often not the dominant cost after recent early-return work.

## Medium-Priority Candidates

### 4. D-optimal and A-optimal search

Relevant file:

- `EDI/src/optimal_design_search.cpp`

Why this is tempting:

- The performance report shows `d_optimal_search` as the largest kernel in a broader non-GLMM sweep.
- The expensive inner step scans all treatment/control swaps, which is a large independent candidate-evaluation problem.

Why this is not the first GPU target:

- The outer loop is an iterative greedy search with data-dependent stopping and state updates after each accepted swap.
- Each accepted swap changes `Pw` and optionally `Hw`, so the GPU would repeatedly evaluate a candidate grid, reduce to the best candidate, transfer or synchronize, and update state.
- Existing notes already identify algorithmic pruning as the real next opportunity.

Potential GPU strategy:

- For each greedy iteration, compute all swap deltas on the GPU and reduce to the best swap.
- Keep `P`, `H`, `Pw`, `Hw`, treatment indices, and control indices resident on the device.
- Use the GPU only when `n_T * n_C` is large enough to amortize repeated reductions.

Recommendation:

- Try algorithmic pruning or candidate-set restriction first.
- Revisit GPU if real workloads have large `n` and many `nsim`, and if the full swap scan remains necessary.

### 5. Zero-inflated / zero-augmented likelihood Hessians

Relevant files:

- `EDI/src/fast_zinb.cpp`
- `EDI/src/fast_zero_augmented_poisson.cpp`
- `EDI/src/fast_hurdle_negbin.cpp`
- `EDI/src/fast_hurdle_poisson_glmm.cpp`

Why this is plausible:

- ZINB and ZAP were large in the broader performance sweep.
- Hessian construction is a per-observation accumulation of small outer products, which can be parallelized on GPU.

Why this is not an easy win:

- ZINB spends substantial time in special functions such as digamma, trigamma, `lgamma`, powers, logs, and exponentials.
- GPU special-function implementations may differ in accuracy, speed, and availability.
- Reducing many per-observation small outer products into one Hessian requires careful reduction strategy to avoid atomics or nondeterministic floating-point noise.

Recommendation:

- Keep CPU as default.
- Consider GPU only for large `n`, fixed `p`, full-path variance computations, and benchmarked workloads where Hessian time dominates.
- A CPU-side improvement such as preallocated special-function tables may be lower-risk and should be exhausted first.

## Lower-Priority or Poor GPU Fits

### R callback bootstrap loops

Relevant files:

- `EDI/src/base_bootstrap_loop.cpp`
- `EDI/src/kk_bootstrap_loop.cpp`

These loops call back into R functions and manipulate R objects. They are not GPU candidates as written. GPU acceleration would require lifting the full estimator into native C++/GPU code so each replicate can run without R callbacks.

### Small exact-test and scalar CI kernels

Relevant files include:

- `EDI/src/miettinen_nurminen_speedups.cpp`
- `EDI/src/newcombe_speedups.cpp`
- `EDI/src/zhang_exact_speedups.cpp`
- `EDI/src/bisection_ci.cpp`
- `EDI/src/bisection_ci_loop.cpp`
- `EDI/src/lrt_ci_newton.cpp`

These are often scalar, branch-heavy, or involve short root-finding loops. They may benefit from CPU vectorization or batching, but they are unlikely to justify GPU complexity unless a future API batches thousands of independent intervals at once.

### Rank/sort-heavy Wilcoxon paths with nonzero shift

Relevant files:

- `EDI/src/fast_wilcox_parallel.cpp`
- `EDI/src/fast_kk_wilcox_parallel.cpp`

When `delta != 0`, each replicate may require sorting shifted values. GPU sorting is possible, but it adds complexity and is not the easiest first target. The `delta == 0` pre-ranked path should be the first Wilcoxon GPU experiment.

## Suggested GPU Backend Design

### Recommended backend

I would use `ggml` as the first GPU backend layer for this package, with CPU/OpenMP kept as the default and fallback.

The reason is portability. `ggml` exposes a common C/C++ tensor backend interface and already has backends for the hardware families we would want to support: CPU, CUDA for NVIDIA GPUs, Metal for Apple Silicon, Vulkan for cross-vendor GPUs, SYCL for Intel-oriented deployments, and HIP/other backends through the broader `ggml` / `llama.cpp` ecosystem. The official ggml Vulkan docs describe Vulkan as cross-vendor acceleration for NVIDIA, AMD, Intel, Qualcomm, and ARM Mali GPUs on Linux, Windows, and Android. The llama.cpp backend docs also list Metal, CUDA, HIP, Vulkan, SYCL, BLAS, and other targets. This makes `ggml` a better first choice than a CUDA-only path if the goal is reusable package infrastructure rather than a one-off NVIDIA benchmark.

The second reason is deployment. This is an R package, so adding a hard dependency on CUDA, ROCm, or oneAPI would make installation much harder for many users. A `ggml`-based optional backend can be built only when requested, while the existing CPU/OpenMP implementation remains the portable baseline.

The third reason is fallback behavior. `ggml` supports backend scheduling patterns where unsupported operations can fall back to CPU. That is useful here because EDI has mixed workloads: some are clean dense tensor reductions, some are rank/sort-heavy, and some call into R or use branch-heavy exact-test logic.

Important caveat: `ggml` should not be treated as a magic replacement for all kernels. It is strongest when the computation can be represented as tensor operations, matrix products, elementwise transforms, reductions, and repeated graph execution. For custom rank sorting, greedy swap search, batched tiny QR solves, or model-specific special functions, a hand-written CUDA/HIP/SYCL kernel or the current CPU/OpenMP code may still be faster and simpler. The backend plan should therefore be:

1. Use `ggml` for portable tensor/reduction prototypes.
2. Use the `ggml` Vulkan backend as the broadest cross-vendor first GPU target.
3. Use `ggml` Metal on Apple Silicon, CUDA on NVIDIA, and SYCL/HIP where those backends are available and benchmark better.
4. Add native vendor-library kernels only for proven bottlenecks that `ggml` cannot express efficiently, such as batched OLS solves through cuSOLVER/cuBLAS on NVIDIA.

References:

- `ggml` Vulkan backend: <https://ggml-org-ggml.mintlify.app/backends/vulkan>
- `llama.cpp` backend overview: <https://www.mintlify.com/ggml-org/llama.cpp/concepts/backends>
- `llama.cpp` supported backend list: <https://github.com/ggml-org/llama.cpp>

### Optional backend

Keep GPU support optional. R packages with CUDA dependencies are harder to build, test, distribute, and install. A practical architecture is:

- CPU/OpenMP remains the default.
- GPU support is enabled by a compile-time flag or optional companion package.
- Runtime dispatch checks data size, requested backend, and GPU availability.
- CPU and GPU outputs are compared in tests with explicit tolerances.

### Candidate libraries

Possible implementation routes:

- `ggml` as the preferred first backend layer for portable tensor graphs and reductions.
- CUDA C++ only for NVIDIA-specific kernels where `ggml` is not expressive enough or where cuBLAS/cuSOLVER batched linear algebra is clearly faster.
- OpenCL, SYCL, HIP, or raw Vulkan only if `ggml` blocks a needed operation and the workload justifies maintaining backend-specific code.
- R-level packages such as `torch` or GPU matrix packages are probably less suitable for the current native C++ package internals, but could be useful for prototypes.

### Dispatch threshold

GPU dispatch should require enough work to amortize overhead. Useful heuristics:

- For replicate-parallel kernels: `n * B` or `n * r` must exceed a benchmarked threshold.
- For pair distances: `m * m * p` must exceed a benchmarked threshold, and the output matrix must fit comfortably in GPU memory.
- For GLMM: number of groups times quadrature nodes times observations per optimizer evaluation must be large, and data must stay resident across iterations.

## Recommended Implementation Order

1. Prototype GPU `compute_wilcox_distr_parallel_cpp()` for the `delta == 0` fast path.
   - Minimal math, simple reductions, clear parity checks.
   - Good way to establish build, dispatch, device-memory, and test infrastructure.

2. Prototype GPU `compute_pair_distance_matrix_cpp()`.
   - Simple independent output cells and a clear large-problem scaling story.
   - Benchmark against the current raw-pointer CPU version.

3. Prototype GPU OLS randomization/bootstrap distributions.
   - Start with `compute_ols_distr_parallel_cpp()`, then consider `compute_ols_bootstrap_parallel_cpp()`.
   - Use batched sufficient-statistic accumulation and small-system solves.

4. Reassess GLMM GPU work only after the first three prototypes.
   - Start with logistic GLMM estimate-only objective/gradient.
   - Only proceed if large real workloads show optimizer objective evaluations still dominate after CPU improvements.

5. Treat D/A-optimal search and ZINB/ZAP Hessians as second-wave experiments.
   - They may have high upside but require more custom algorithm work and more careful numerical validation.

## Benchmarking Plan

For every GPU candidate, compare:

- CPU single-thread, CPU OpenMP, and GPU elapsed time.
- Transfer-included and transfer-excluded timing.
- Small, medium, and large workload sizes.
- Numerical parity against current CPU results.
- Determinism or acceptable floating-point tolerance.

Suggested initial benchmark matrix:

| Kernel | Small | Medium | Large |
|---|---:|---:|---:|
| Wilcoxon `delta == 0` | `n=500, r=501` | `n=5,000, r=5,001` | `n=50,000, r=10,001` |
| Pair distances | `m=500, p=20` | `m=5,000, p=30` | memory-limited |
| OLS randomization | `n=800, p=30, r=501` | `n=5,000, p=50, r=5,001` | `n=20,000, p=100, r=10,001` |
| Logistic GLMM | current benchmark scale | 10x groups | 100x groups |

## Bottom Line

The most practical GPU wins are not the individual scalar model fits. They are the batched workloads already present in the package: randomization distributions, bootstrap distributions, Wilcoxon rank-sum reductions, OLS replicate fits, and pair-distance matrix construction.

The GLMM and frailty likelihoods are possible GPU targets, but only after data-resident objective/gradient evaluation is designed carefully. D/A-optimal search and ZINB/ZAP Hessians have visible CPU cost, but their first next step is probably algorithmic or CPU-side cleanup rather than immediate GPU porting.
