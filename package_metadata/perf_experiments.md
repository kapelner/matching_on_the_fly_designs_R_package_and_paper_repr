# Performance Experiments Report

Date: 2026-05-14

## Objective

Profile the current native kernels in `EDI/src` and rank the hottest loops by expected optimization payoff.

This report uses the current working tree built into a temporary library and profiled with `perf`. It does not rely on the older `benchmark/simd_matrix/simd_matrix_summary.csv` snapshot except as background context.

## Method

### Build and benchmark setup

The current package was installed into a temporary library and benchmarked using the existing SIMD benchmark worker scenarios in:

- `scripts/benchmark_simd_matrix_worker.R`

To isolate native kernels for profiling, a temporary one-kernel driver was used to run one benchmark scenario repeatedly per process. Each kernel was then sampled with `perf record`.

### Profiling commands

Representative command shape:

```bash
perf record -F 199 --output=/tmp/perf_<kernel>.data -- \
  Rscript /tmp/profile_one_kernel.R --lib=/tmp/edi_profile_lib --kernel=<kernel> --reps=<N>
perf report --stdio --no-children --sort symbol -i /tmp/perf_<kernel>.data
```

### Important caveat

The installed `EDI.so` did not provide clean line-level attribution through `perf report`, so the final loop ranking is an inference from:

1. Per-kernel elapsed time in the current build.
2. Top sampled symbols from `perf`.
3. The surrounding source structure in the relevant C++ files.

This is still strong enough to rank optimization targets, but the ranking is loop-level reasoning rather than exact per-line sample accounting.

## Current Kernel Weights

Current medians from the freshly installed working tree:

| kernel | median elapsed (s) | share of total |
|---|---:|---:|
| `weibull_frailty` | 0.0480 | 29.1% |
| `ordinal_clmm` | 0.0470 | 28.5% |
| `logistic_glmm` | 0.0370 | 22.4% |
| `clogit_plus_glmm` | 0.0185 | 11.2% |
| `poisson_glmm` | 0.0115 | 7.0% |
| `adjacent_category_logit` | 0.0030 | 1.8% |

Interpretation:

- The performance budget is dominated by `weibull_frailty`, `ordinal_clmm`, and `logistic_glmm`.
- Even a modest reduction in those three kernels is worth more than a large reduction in `adjacent_category_logit`.

## Perf Highlights

### Logistic GLMM

Top sampled symbols included:

- `__log1p_fma`
- `plogis_array_safe(...)`
- `Eigen::internal::general_matrix_vector_product...`
- `log1pexp_array_safe(...)`
- `LogisticGLMMObjective::operator()(...)`

Interpretation:

- The kernel is dominated by repeated logistic math and matrix-vector products.
- The main issue is not lack of SIMD. SIMD is already active.
- The likely optimization headroom is in reducing repeated passes, temporary vectors, and redundant reductions.

### Weibull frailty

Top sampled symbols included:

- `WeibullFrailtyLikelihood::operator()(...)`
- vectorized Eigen packet math
- `Eigen::internal::general_matrix_vector_product...`
- vectorized exponentials

Interpretation:

- This is a high-throughput numeric kernel already using vectorized code paths.
- The remaining payoff is mostly from dataflow simplification and fewer materialized temporaries, not handwritten assembly.

### Ordinal CLMM

Top sampled symbols included:

- `glmm::GLMMObjective<OrdinalLikelihood<LogitLink>>::operator()(...)`
- `__ieee754_exp_fma`
- `__ieee754_log_fma`
- `malloc`
- `cfree`
- `Eigen::internal::general_matrix_vector_product...`

Interpretation:

- The generic GLMM engine is expensive here.
- Allocation overhead is visible.
- This path likely suffers from repeated tiny dynamic allocations and repeated threshold reconstruction inside the innermost derivative logic.

### Clogit plus GLMM

Top sampled symbols included:

- `__log1p_fma`
- `plogis_array_safe(...)`
- `log1pexp_array_safe(...)`
- `Eigen::internal::general_matrix_vector_product...`
- `ClogitPlusGLMMObjective::operator()(...)`

Interpretation:

- The concordant GLMM part has the same structural cost pattern as the logistic GLMM kernel.

### Poisson GLMM

Top sampled symbols included:

- `PoissonGLMMObjective::operator()(...)`
- `Eigen::internal::general_matrix_vector_product...`
- Eigen packet math
- `log_sum_exp_p(...)`

Interpretation:

- The dominant cost is again the repeated node/group sweep and weighted matrix-vector work.

### Adjacent-category logit

Top sampled symbols included:

- `AdjacentCategoryLogitNegLogLik::operator()(...)`
- `AdjacentCategoryLogitNegLogLik::hessian(...)`
- Eigen dense assignment and outer-product routines

Interpretation:

- This kernel is internally hot but globally small.
- It should not be optimized before the larger GLMM-family kernels.

## Top 10 Loops By Expected Payoff

The ranking below combines:

- current kernel weight,
- `perf` samples,
- amount of repeated work per optimizer iteration,
- and how much avoidable allocation or temporary materialization is present.

1. `EDI/src/_glmm_engine.h`
   Generic ordinal CLMM node-by-row derivative loop in `GLMMObjective::operator()`.
   Why: `ordinal_clmm` was 28.5% of baseline kernel time and `perf` showed generic-engine work plus `malloc`/`cfree`.

2. `EDI/src/fast_ordinal_clmm.cpp`
   Ordinal threshold reconstruction and chain-rule suffix loop in `OrdinalLikelihood::log_prob_derivs()`.
   Why: nested directly inside item 1, with repeated threshold work per row and quadrature node.

3. `EDI/src/fast_weibull_frailty.cpp`
   Weibull frailty node sweep building `wik`, `ewik`, and group log terms.
   Why: `weibull_frailty` was the single largest kernel and `perf` showed heavy vector-exp/GEMV cost in this sweep.

4. `EDI/src/fast_weibull_frailty.cpp`
   Weibull frailty posterior-expansion and gradient accumulation loop.
   Why: the original path materialized expanded posterior weights and paid extra group reductions before matrix multiplies.

5. `EDI/src/fast_logistic_glmm.cpp`
   Logistic GLMM node sweep building `mu_all_k_vec` and `log_terms_mat`.
   Why: `logistic_glmm` was 22.4% of baseline kernel time and `perf` was dominated by logistic math plus GEMV.

6. `EDI/src/fast_logistic_glmm.cpp`
   Logistic GLMM posterior-expansion and gradient accumulation loop.
   Why: same expanded-posterior pattern as Weibull, with an additional whole-matrix gradient pass.

7. `EDI/src/fast_clogit_plus_glmm.cpp`
   Clogit-plus-GLMM concordant node sweep.
   Why: smaller than logistic, but structurally the same and therefore a good transfer target.

8. `EDI/src/fast_clogit_plus_glmm.cpp`
   Clogit-plus-GLMM concordant posterior-expansion and gradient accumulation loop.
   Why: same dense expanded-posterior pattern and same direct-accumulation opportunity.

9. `EDI/src/fast_poisson_glmm.cpp`
   Poisson GLMM node sweep.
   Why: `perf` showed a hot objective loop, but the kernel’s absolute share was lower than ordinal/Weibull/logistic/clogit.

10. `EDI/src/fast_poisson_glmm.cpp`
    Poisson GLMM posterior-expansion and gradient accumulation loop.
    Why: same family pattern as logistic/clogit/Weibull, but lower absolute expected payoff.




## Why Adjacent-category Logit Is Not Top 10

`perf` shows `AdjacentCategoryLogitNegLogLik::operator()` and `hessian()` are locally hot inside that kernel, but the kernel itself is only 1.8% of current total native time. That makes it a lower-priority target than the GLMM-family loops above.

This is the right outcome from a payoff ranking:

- local hotspot does not imply global priority,
- kernel weight matters more than per-kernel percentage.

## Optimization Plan

### Phase 1: highest-value structural fixes

1. Specialize the ordinal CLMM path instead of routing the hottest work through the generic GLMM engine.
2. Remove repeated threshold reconstruction in `OrdinalLikelihood::log_prob_derivs()`.
3. Eliminate repeated tiny allocations in the ordinal derivative path by reusing work buffers.
4. Replace `post_k_expanded` construction in Weibull and logistic GLMM with direct accumulation into gradient terms.

Expected effect:

- Largest likely runtime reduction with moderate code complexity.

### Phase 2: shared GLMM loop cleanup

1. Refactor logistic, Poisson, Weibull, and clogit-concordant kernels around a shared pattern:
   - node sweep,
   - group posterior weights,
   - direct accumulation without full expanded posterior vectors.
2. Fuse adjacent passes so `mu`, residuals, weighted sums, and group log terms are reused.
3. Reuse preallocated buffers across optimizer calls where object lifetime allows it.

Expected effect:

- Broad improvement across four kernels with similar code structure.

### Phase 3: math and reduction tightening

1. Audit repeated `exp`, `log1p`, and `log_sum_exp` calls for avoidable recomputation.
2. Cache quantities like `sqrt(2.0) * sigma * node` and other repeated node factors.
3. Reduce redundant segment extraction and temporary `Eigen::VectorXd` creation in innermost loops.

Expected effect:

- Moderate speedups after the larger structural changes are done.

### Phase 4: Hessian-path optimization

1. Apply the same temporary-elimination strategy to the Hessian code paths in logistic, Poisson, Weibull, and clogit-plus-GLMM.
2. Reconsider whether all full Hessian computations need the same level of exactness for every consumer.

Expected effect:

- Important for non-`estimate_only` paths, but lower immediate payoff than gradient/objective optimization because the benchmark scenarios used here were `estimate_only = TRUE`.

### Phase 5: lower-priority cleanup

1. Optimize `fast_adjacent_category_logit.cpp` only after the dominant kernels are improved.
2. Focus on:
   - reducing repeated `VectorXd::LinSpaced(...)`,
   - avoiding per-row small temporary vectors,
   - simplifying Hessian outer-product updates.

Expected effect:

- Small global benefit.

## Recommended Implementation Order

1. `EDI/src/_glmm_engine.h`
2. `EDI/src/fast_ordinal_clmm.cpp`
3. `EDI/src/fast_weibull_frailty.cpp`
4. `EDI/src/fast_logistic_glmm.cpp`
5. `EDI/src/fast_clogit_plus_glmm.cpp`
6. `EDI/src/fast_poisson_glmm.cpp`
7. `EDI/src/fast_adjacent_category_logit.cpp`

## Final Next Steps

Based on the retained/reverted outcomes in this report, the next practical work should be:

1. Re-profile ordinal CLMM on the current retained tree and decide whether there is a worthwhile follow-up beyond the kept threshold/allocation fixes.
   The dedicated specialization attempt was slower, so the next ordinal change should be profiling-led rather than architectural by default.
2. Keep Poisson GLMM on its current path and treat it separately from the shared logistic/Weibull/clogit pattern.
   The repeated shared-pattern rewrites preserved the math but did not beat the existing Poisson structure.
3. If more GLMM work is needed, concentrate on measured full-path consumers rather than reopening estimate-only paths that are already improved.
   The biggest retained Phase 4 win was avoiding unnecessary Hessian work for `estimate_only = TRUE`.
4. If more non-GLMM work is needed, prefer `optimal_design_search.cpp` and `pair_dist_helpers.cpp` style access-path cleanup over speculative math rewrites in ZINB/ZAP/negbin.
   That was the only additional-sweep area with clear retained speedups.
5. Before any further Hessian-heavy count-model work, collect a more targeted profiler trace around ZINB and ZAP special-function calls.
   The previous structural cleanups were estimate-stable but slower.

## Bottom Line

There is no evidence that handwritten assembly is the right next step.

The current profile says the real opportunities are:

- fewer temporary vectors,
- fewer full-array passes,
- less per-node/per-row allocation,
- and specialized code paths for the dominant ordinal and GLMM kernels.

That is where the top-end payoff is.

## Additional Perf Investigation: Non-`estimate_only` And `with_var` Paths

The first pass focused on the estimate-only benchmark kernels because they dominate repeated optimizer work in the existing SIMD benchmark matrix.

To extend the investigation, the following secondary native kernels were profiled:

- `ordinal_clmm_full`
- `weibull_frailty_full`
- `logistic_glmm_full`
- `poisson_glmm_full`
- `clogit_plus_glmm_full`
- `adjacent_category_logit_var`

These paths exercise:

- post-fit information matrix computation,
- Hessian code paths,
- and `with_var` variants.

### Secondary kernel weights

Current medians from the same temporary-library build:

| kernel | median elapsed (s) |
|---|---:|
| `ordinal_clmm_full` | 0.147 |
| `weibull_frailty_full` | 0.078 |
| `logistic_glmm_full` | 0.048 |
| `poisson_glmm_full` | 0.024 |
| `clogit_plus_glmm_full` | 0.013 |
| `adjacent_category_logit_var` | 0.002 |

Relative to the estimate-only versions:

- `ordinal_clmm_full` is much more expensive than `ordinal_clmm`.
- `weibull_frailty_full` is materially more expensive than `weibull_frailty`.
- `poisson_glmm_full` is roughly double the estimate-only path.
- `logistic_glmm_full` is only modestly above estimate-only.
- `adjacent_category_logit_var` remains globally tiny.

### Secondary perf highlights

#### Ordinal CLMM full

Top sampled symbols included:

- `glmm::GLMMObjective<OrdinalLikelihood<LogitLink>>::operator()(...)`
- `__ieee754_exp_fma`
- `malloc`
- `cfree`
- `__ieee754_log_fma`
- `Eigen::internal::general_matrix_vector_product...`

Interpretation:

- The full path remains dominated by the same generic GLMM objective machinery as the estimate-only path.
- Allocation pressure remains highly visible.
- The post-fit work increases wall time significantly, but it does not move the center of gravity away from the generic objective loop and the repeated derivative/threshold logic.

#### Weibull frailty full

Top sampled symbols included:

- `WeibullFrailtyLikelihood::operator()(...)`
- Eigen packet math
- `Eigen::internal::general_matrix_vector_product...`
- vectorized exponential code paths

Interpretation:

- The full path is still dominated by the objective/gradient loop.
- The Hessian path increases total runtime, but sampled cycles remain concentrated in the same node-sweep and reduction structure.
- This further supports prioritizing objective-loop refactoring first, then tightening Hessian code afterward.

#### Logistic GLMM full

Top sampled symbols included:

- `__log1p_fma`
- `plogis_array_safe(...)`
- `Eigen::internal::general_matrix_vector_product...`
- `LogisticGLMMObjective::operator()(...)`
- `LogisticGLMMObjective::hessian(...)` at much lower sampled weight

Interpretation:

- The full path is still mostly objective/gradient bound.
- The Hessian implementation is visible, but far below the logistic math and GEMV-heavy main loop.
- Optimizing the primary objective loop still gives the best payoff even for the full path.

#### Poisson GLMM full

Top sampled symbols included:

- `PoissonGLMMObjective::operator()(...)`
- `PoissonGLMMObjective::hessian(...)` at lower sampled weight
- `Eigen::internal::general_matrix_vector_product...`
- Eigen packet math

Interpretation:

- Same pattern as logistic: the Hessian matters, but the optimizer is still spending more sampled cycles in the core objective loop.
- A shared objective-loop refactor remains the right first move.

### New takeaway from the secondary pass

The non-`estimate_only` paths do not overturn the first ranking.

Instead they strengthen it:

- the same objective/gradient loops remain dominant,
- ordinal CLMM remains the most structurally problematic path,
- and Hessian-specific tuning should generally come after the core objective loop cleanup.

## Expanded Optimization Plan

### Phase 1A: ordinal CLMM core path

No change from the first plan, but the evidence is now stronger.

Priority tasks:

1. Replace the generic ordinal GLMM derivative path in `EDI/src/_glmm_engine.h` for the CLMM logit case.
2. Cache threshold transforms and eliminate repeated per-row `dp` allocation in `EDI/src/fast_ordinal_clmm.cpp`.
3. Preallocate node and derivative workspaces to suppress `malloc` / `cfree` churn.

Why this is now first:

- `ordinal_clmm` was already large in estimate-only mode.
- `ordinal_clmm_full` widens the gap further.
- `perf` still shows allocation-heavy generic machinery dominating samples.

### Phase 1B: shared objective-loop refactor for GLMM-family kernels

Priority files:

- `EDI/src/fast_weibull_frailty.cpp`
- `EDI/src/fast_logistic_glmm.cpp`
- `EDI/src/fast_poisson_glmm.cpp`
- `EDI/src/fast_clogit_plus_glmm.cpp`

Priority tasks:

1. Remove `post_k_expanded` materialization where possible.
2. Fuse group posterior-weight accumulation with gradient accumulation.
3. Reuse node-wise temporaries rather than constructing fresh vectors in each pass.
4. Reduce redundant segment extraction and full-length residual materialization.

Why this remains second:

- The full-path profiling still points back to these same loops.
- Hessian work does not dominate enough to justify starting there.

### Phase 2: Hessian/information cleanup after core loops

Now explicitly supported by `perf`:

1. Optimize Hessian accumulation in logistic and Poisson GLMM after the operator-path cleanup lands.
2. Apply the same temporary-reduction ideas to Weibull Hessian code.
3. Recheck whether ordinal full-path covariance computation can avoid the current generic fallback structure.

Why this is not earlier:

- Hessian symbols are visible but remain secondary in sampled cycle share.
- Cleaner objective loops should reduce duplicated effort across both estimate-only and full paths.

### Phase 3: low-global-payoff kernels

Files:

- `EDI/src/fast_adjacent_category_logit.cpp`

Tasks:

1. Reduce per-row temporary vector creation.
2. Eliminate repeated `VectorXd::LinSpaced(...)`-style work in the inner loops.
3. Tighten Hessian outer-product updates.

Why it remains late:

- Even the `with_var` path is still very small globally.

## Revised Recommended Implementation Order

1. `EDI/src/_glmm_engine.h`
2. `EDI/src/fast_ordinal_clmm.cpp`
3. `EDI/src/fast_weibull_frailty.cpp`
4. `EDI/src/fast_logistic_glmm.cpp`
5. `EDI/src/fast_poisson_glmm.cpp`
6. `EDI/src/fast_clogit_plus_glmm.cpp`
7. Hessian/information cleanup in the same files
8. `EDI/src/fast_adjacent_category_logit.cpp`

## Practical Conclusion From The Second Perf Pass

The extra `perf` work on the non-`estimate_only` kernels did not reveal a new dominant family of C++ hotspots.

Instead it showed:

- the ordinal generic GLMM machinery is still the most urgent structural problem,
- the Weibull/logistic/Poisson/clogit objective loops are still the right shared optimization target,
- and Hessian-path tuning is important but downstream of the core loop cleanup.

## Additional Hotspot Sweep Beyond The GLMM Family

After the GLMM-focused work, a broader native sweep was run over additional representative kernels that appear repeatedly in:

- the package's benchmark scripts,
- the `with_var` regression surface,
- and the heavier standalone design/helper code.

The selected kernels were:

- `d_optimal_search`
- `zinb_full`
- `zap_full`
- `negbin_var`
- `pair_distance_matrix`
- `gaussian_lmm_full`
- `coxph_full`

These were benchmarked with in-process timing to avoid `Rscript` startup distortion.

### Additional kernel weights

Representative current medians:

| kernel | median elapsed (s) | share within this sweep |
|---|---:|---:|
| `d_optimal_search` | 2.0345 | 37.6% |
| `zinb_full` | 1.3920 | 25.7% |
| `zap_full` | 1.3350 | 24.7% |
| `negbin_var` | 0.2900 | 5.4% |
| `pair_distance_matrix` | 0.2105 | 3.9% |
| `gaussian_lmm_full` | 0.1420 | 2.6% |
| `coxph_full` | 0.0055 | 0.1% |

Interpretation:

- `d_optimal_search` is the dominant non-GLMM hotspot from this additional set.
- `zinb_full` and `zap_full` are both large enough to matter and are much heavier than `negbin_var`.
- `pair_distance_matrix` is not globally dominant, but it is a real standalone hotspot and it is structurally easy to improve.
- `coxph_full` is not a priority hotspot at this representative problem size.

## Perf Highlights For The Additional Sweep

### D-optimal search

Top sampled symbols included:

- `d_optimal_search_cpp(...)`
- repeated Eigen coefficient access on `P(i, j)` and vector entries
- almost no time in `std::shuffle`

Interpretation:

- The cost is not initialization or RNG.
- The dominant work is the nested swap-search loop in `EDI/src/optimal_design_search.cpp`, especially repeated random-access reads from `P`, `Pw`, and writes to `w`.
- This is a branchy search kernel with high scalar access overhead. It is not an assembly problem; it is a data-access and algorithmic-structure problem.

Most likely hot loops:

- the `for (ti)` / `for (cj)` pairwise swap scan in `d_optimal_search_cpp`
- the analogous search loop in `a_optimal_search_cpp`

### Pair-distance matrix

Top sampled symbols included:

- `compute_pair_distance_matrix_cpp(...)`
- `Rcpp::Matrix::operator()`
- `Rcpp::Matrix::offset`
- `Rcpp::Vector::operator[]`

Interpretation:

- This kernel is dominated by accessor overhead, not floating-point arithmetic.
- The implementation in `EDI/src/pair_dist_helpers.cpp` is paying heavily for Rcpp matrix/vector indexing inside the innermost loop.
- This is one of the clearest low-risk rewrite targets in the codebase.

Most likely hot loop:

- the triple loop over `i`, `j`, and `k` in `compute_pair_distance_matrix_cpp`

### Zero-inflated negative binomial

Top sampled symbols included:

- `ZeroInflatedNegBin::hessian(...)`
- `Rf_dpsifn...`
- `__ieee754_pow_fma`
- `__ieee754_log_fma`

Interpretation:

- Unlike the GLMM-family kernels, this path is visibly Hessian-heavy.
- Special functions and Hessian accumulation are significant here.
- This makes `fast_zinb.cpp` an important exception to the earlier “optimize objective loops first” rule: for ZINB, Hessian work is already prominent enough to matter on its own.

Most likely hot loops:

- the observation loop inside `ZeroInflatedNegBin::hessian(...)`
- the lower-triangular Hessian accumulation over conditional, zero-inflation, and theta blocks

### Zero-augmented Poisson

Top sampled symbols included:

- `ZeroAugmentedPoisson::hessian(...)`
- exponential/logistic scalar operations
- Eigen expression overhead

Interpretation:

- This path also leans more heavily on Hessian work than the earlier GLMM kernels.
- The structure is similar to ZINB, but without the negative-binomial special functions.
- This makes it a meaningful secondary optimization target after the dominant GLMM/objective loops and the D-optimal search kernel.

Most likely hot loops:

- the observation loop inside `ZeroAugmentedPoisson::hessian(...)`
- the blockwise Hessian accumulation in `EDI/src/fast_zero_augmented_poisson.cpp`

### Negative binomial with variance

Top sampled symbols included:

- `NBLogLik::operator()(...)`
- `__ieee754_log_fma`
- `__log1p_fma`
- `R::digamma(...)`

Interpretation:

- This kernel is still mostly objective-driven rather than Hessian-driven.
- Its runtime comes from per-observation negative-binomial likelihood work and special-function evaluation.
- It is worth optimizing after the larger ZINB/ZAP and GLMM targets, but not before them.

Most likely hot loop:

- the main observation loop inside `NBLogLik::operator()(...)` in `EDI/src/fast_negbin_regression.cpp`

### Gaussian LMM and Cox PH

Representative timing showed:

- `gaussian_lmm_full` is modest but not dominant
- `coxph_full` is tiny at the tested problem size

Interpretation:

- They are not top-priority hotspots relative to the rest of the native stack.
- They should remain below the GLMM, D-optimal search, ZINB/ZAP, and pair-distance work in the optimization queue unless a larger real-world workload shows otherwise.

## Expanded Ranking: Additional Non-GLMM Targets

These additional targets should be inserted into the broader payoff ranking as follows:

### New high-priority additions

1. `EDI/src/optimal_design_search.cpp`
   The nested swap-search loops in `d_optimal_search_cpp` and `a_optimal_search_cpp` are now confirmed hotspots.

2. `EDI/src/fast_zinb.cpp`
   The Hessian path is materially hot, with special-function cost clearly visible.

3. `EDI/src/fast_zero_augmented_poisson.cpp`
   The Hessian path is also materially hot, though somewhat less exotic than ZINB.

### Medium-priority additions

4. `EDI/src/pair_dist_helpers.cpp`
   `compute_pair_distance_matrix_cpp` is a very attractive cleanup target because it is simple, isolated, and clearly dominated by avoidable indexing overhead.

5. `EDI/src/fast_negbin_regression.cpp`
   Useful, but downstream of the larger ZINB/ZAP and GLMM work.

### Lower-priority additions

6. `EDI/src/fast_gaussian_lmm.cpp`
   Worth revisiting later, but not currently a top hotspot in representative runs.

7. `EDI/src/fast_coxph_regression.cpp`
   Not a priority hotspot at the tested size.

## Expanded Optimization Plan

### Phase 1C: D-optimal search cleanup

Target file:

- `EDI/src/optimal_design_search.cpp`

Tasks:

1. Reduce repeated `P(i, j)` and vector coefficient lookups in the nested swap-search loops.
2. Cache row pointers or switch to a layout/access pattern with cheaper repeated reads.
3. Consider storing treatment/control candidate indices in structures optimized for repeated scan/update.
4. Reassess whether the full `ti x cj` swap scan can be pruned or short-circuited algorithmically.

Expected effect:

- High payoff on a standalone kernel that is now confirmed to dominate the additional non-GLMM set.

### Phase 2A: Hessian-heavy count models

Target files:

- `EDI/src/fast_zinb.cpp`
- `EDI/src/fast_zero_augmented_poisson.cpp`

Tasks:

1. Reduce temporary construction inside Hessian accumulation.
2. Avoid redundant recomputation of common scalar terms per observation.
3. Collapse repeated blockwise updates where algebra permits.
4. Reuse intermediate quantities between score and Hessian code paths where safe.

Expected effect:

- High for ZINB and medium-high for ZAP.

### Phase 2B: Pair-distance helper rewrite

Target file:

- `EDI/src/pair_dist_helpers.cpp`

Tasks:

1. Replace Rcpp element access inside the innermost loop with raw pointers or Eigen-backed access.
2. Hoist weight lookup and row base-address calculations out of the deepest loop.
3. Consider a row-major or blocked implementation for better locality.

Expected effect:

- Medium absolute payoff, high implementation efficiency.

### Phase 3A: Negative binomial cleanup

Target file:

- `EDI/src/fast_negbin_regression.cpp`

Tasks:

1. Reduce repeated scalar special-function work where possible.
2. Check whether stable approximations or cached terms can reduce `digamma`/`log1p` pressure.
3. Tighten inner observation loops before considering any higher-risk math changes.

Expected effect:

- Medium.

## Revised Overall Implementation Order

Combining all profiling passes, the practical order is now:

1. `EDI/src/_glmm_engine.h`
2. `EDI/src/fast_ordinal_clmm.cpp`
3. `EDI/src/fast_weibull_frailty.cpp`
4. `EDI/src/fast_logistic_glmm.cpp`
5. `EDI/src/fast_poisson_glmm.cpp`
6. `EDI/src/fast_clogit_plus_glmm.cpp`
7. `EDI/src/optimal_design_search.cpp`
8. `EDI/src/fast_zinb.cpp`
9. `EDI/src/fast_zero_augmented_poisson.cpp`
10. `EDI/src/pair_dist_helpers.cpp`
11. `EDI/src/fast_negbin_regression.cpp`
12. Hessian/information cleanup in the already-prioritized GLMM files
13. `EDI/src/fast_gaussian_lmm.cpp`
14. `EDI/src/fast_adjacent_category_logit.cpp`
15. `EDI/src/fast_coxph_regression.cpp`

## Final Updated Conclusion

The broader `perf` sweep did find additional real hotspots outside the original GLMM group, but it did not change the top-level strategy.

The codebase now appears to have three distinct hotspot classes:

1. Generic/objective-loop dominated GLMM kernels
2. Hessian-heavy count models (`ZINB`, `ZAP`)
3. Standalone scalar-access kernels (`d_optimal_search`, `pair_distance_matrix`)

The highest payoff remains structural:

- remove unnecessary allocations,
- reduce repeated passes,
- simplify access patterns,
- and specialize generic hot paths.

That remains a much stronger optimization route than handwritten assembly.


## Phase 1 outcome

Phase 1 was only partially successful, but the retained parts were the highest-value ones.

Final disposition:

- keep the ordinal CLMM threshold-reconstruction rewrite in `EDI/src/fast_ordinal_clmm.cpp`
- keep the ordinal tiny-allocation removal in `EDI/src/_glmm_engine.h`
- keep the direct gradient accumulation rewrite in logistic GLMM
- keep the direct gradient accumulation rewrite in Weibull frailty
- do not keep the dedicated ordinal-CLMM specialization that bypassed the generic GLMM engine

What landed:

1. `OrdinalLikelihood::log_prob_derivs()` was tightened so the derivative path no longer reconstructs the full threshold vector on every row/node evaluation.
2. The generic ordinal path stopped allocating tiny `dp` vectors inside the innermost loop and now reuses work buffers.
3. Logistic and Weibull objective/gradient loops stopped materializing `post_k_expanded` and instead accumulated directly at the group level.

What was rejected:

- A dedicated ordinal CLMM objective that bypassed `GLMMObjective` was implemented and estimate-stable, but it benchmarked slower than the optimized generic path, so it was reverted.

Measured/observed outcome:

- ordinal CLMM improved materially from the threshold and allocation fixes
- logistic GLMM and Weibull frailty both benefited from removing expanded posterior vectors
- the attempted ordinal full specialization did not produce a runtime win

Interpretation:

- the Phase 1 structural wins came from removing repeated work and allocation churn
- simply replacing the generic ordinal engine with a custom class was not enough by itself to beat the compiler-optimized generic path

### Compact benchmark summary

| Kernel | Path | Before median (s) | After median (s) | Decision |
|---|---:|---:|---:|---|
| `ordinal_clmm` | estimate-only fit | `0.0470` | `0.0250` | keep threshold/allocation fixes |
| `weibull_frailty` | estimate-only fit | `0.0480` | `0.0415` | keep direct accumulation changes |
| `logistic_glmm` | estimate-only fit | `0.0470` | `0.0360` | keep direct accumulation changes |
| `ordinal_clmm` | dedicated specialization attempt | `0.0785` | `0.0885` | revert specialization |


## Phase 2 outcome

Phase 2 was broadly successful for the kernels that share the same objective-loop structure, but not uniformly successful across every model.

Final disposition:

- keep the shared-pattern cleanup in logistic GLMM
- keep the shared-pattern cleanup in Weibull frailty
- keep the shared-pattern cleanup in clogit-plus-GLMM concordant
- do not keep the analogous Poisson GLMM rewrite
- keep the persistent scratch-buffer reuse for logistic, Weibull, and clogit-concordant
- keep the small shared contiguous-group layout helper in `EDI/src/_helper_functions.h`
- do not force a larger shared abstraction over all kernels

What landed:

1. Logistic, Weibull, and clogit-concordant were all rewritten around the same practical pattern:
   - node sweep,
   - group posterior weights,
   - direct accumulation without dense expanded posterior vectors.
2. Adjacent passes were fused so group log terms, residuals, and weighted reductions are reused rather than rebuilt.
3. Persistent scratch buffers were added where object lifetime allows it, especially for:
   - group layout metadata,
   - node work matrices,
   - reusable residual/work vectors.
4. A shared contiguous-group layout helper was introduced, which is the part of the common structure that was worth centralizing.

What was rejected:

- The Poisson GLMM version of the same rewrite preserved the math but benchmarked slower, so it was reverted.
- A more aggressive single abstraction across logistic, Weibull, Poisson, and clogit was not adopted because the remaining kernel-specific math differences were large enough that the abstraction would likely add indirection without improving maintainability.

Measured/observed outcome:

- logistic GLMM: retained as a win
- Weibull frailty: retained as a win at the structural-loop level
- clogit-plus-GLMM concordant: retained as a win
- Poisson GLMM: rejected because the measured before/after benchmark was slower

Interpretation:

- Phase 2 worked when it removed memory traffic and temporary-vector construction without fighting the kernel’s natural matrix structure
- Poisson is the main counterexample and should be treated as its own optimization problem rather than forced into the same rewrite

### Compact benchmark summary

| Kernel | Path | Before median (s) | After median (s) | Decision |
|---|---:|---:|---:|---|
| `logistic_glmm` | estimate-only fit | `0.0470` | `0.0360` | keep |
| `weibull_frailty` | estimate-only fit | `0.0480` | `0.0415` | keep |
| `clogit_plus_glmm` | fit fingerprint preserved | same estimates | same estimates | keep |
| `poisson_glmm` | larger KK-shaped fit | `0.0520` | `0.0600` | revert |


## Phase 3 outcome: keep/revert decisions

The Phase 3 math-and-reduction tightening work was implemented and benchmarked for the three highest-value remaining kernels in this family:

- `EDI/src/fast_logistic_glmm.cpp`
- `EDI/src/fast_weibull_frailty.cpp`
- `EDI/src/fast_clogit_plus_glmm.cpp`

The final disposition is:

- keep the Phase 3 changes for logistic GLMM
- keep the Phase 3 changes for clogit-plus-GLMM concordant
- revert the Phase 3 changes for Weibull frailty

The Weibull revert keeps the earlier Phase 2 structural improvements and only backs out the slower Phase 3 math tightening.

### Estimate consistency

The direct before/after checks showed no material numerical change in the optimized paths:

- logistic GLMM: parameters and negative log-likelihood matched in the prior before/after fit check
- Weibull frailty: direct evaluation matched exactly before the revert
  - `nll = 5407.2669016122754`
  - `sum(diag(H)) = 316963.61731134733`
- clogit-plus-GLMM concordant: fit and Hessian fingerprint matched exactly
  - `nll = 398.6283207650863`
  - `param_sum = -1.7298381109594083`
  - `sum(diag(H)) = -566.51882986578255`

### Benchmarks

| Kernel | Path | Before median (s) | After median (s) | Decision |
|---|---:|---:|---:|---|
| `logistic_glmm` | estimate-only fit | `0.0430` | `0.0390` | keep |
| `logistic_glmm` | Hessian | `0.0070` | `0.0030` | keep |
| `weibull_frailty` | neg log-likelihood eval | `0.000680` | `0.001985` | revert |
| `weibull_frailty` | Hessian | `0.00330` | `0.00744` | revert |
| `clogit_plus_glmm` | estimate-only fit | `0.0500` | `0.0360` | keep |
| `clogit_plus_glmm` | Hessian | `0.01015` | `0.00515` | keep |

Interpretation:

- logistic GLMM got a modest objective-path win and a large Hessian-path win
- clogit-plus-GLMM got a clear win in both the objective and Hessian paths
- Weibull frailty preserved the math but got slower, so the Phase 3 version should not be kept

### Test coverage notes

Relevant repo coverage identified during this pass:

- [EDI/tests/testthat/test-glmm-cpp-equivalence.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/tests/testthat/test-glmm-cpp-equivalence.R) contains canonical-package equivalence coverage for logistic GLMM against `lme4::glmer`
- [EDI/tests/testthat/test-weibull-frailty.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/tests/testthat/test-weibull-frailty.R) exercises the Weibull frailty inference path

There is no dedicated `testthat` file in this repo for `clogit_plus_glmm` canonical-package parity.

Package-backed `testthat` execution was not completed cleanly in this pass because ad hoc rebuilds through `pkgload::load_all("EDI")` and `sourceCpp()` repeatedly hit local build-artifact issues around `_helper_functions.o`. The direct numerical before/after checks and isolated benchmarks above are therefore the authoritative result for this Phase 3 decision.


## Phase 4 outcome

Phase 4 focused on the Hessian-heavy paths and on whether all callers actually need full information matrices on every fit.

Final disposition:

- keep the logistic GLMM Hessian cleanup
- keep the clogit-plus-GLMM Hessian cleanup
- keep a partial-information mode for estimate-only logistic, clogit, and Poisson fits
- do not keep a grouped Poisson Hessian rewrite
- do not keep any additional Weibull Hessian tightening in this pass

### Caller audit

The useful caller-level finding was not that many R callers need a new public API for approximate information, but that many native fit calls already know when they only need point estimates.

In practice:

- `estimate_only = TRUE` callers in the GLMM families were still paying for full Hessian / information work inside the native fitters
- the higher-level R code already prefers stored information objects like `fit$fisher_information` when they exist
- the highest-value partial-information mode was therefore an internal early-return path in the native fitters rather than a broad new R-facing API

What landed from the audit:

1. `fast_logistic_glmm_cpp(...)` now returns immediately after optimization when `estimate_only = TRUE`, without building Hessian-derived outputs.
2. `fast_clogit_plus_glmm_cpp(...)` now does the same.
3. `fast_poisson_glmm_cpp(...)` now does the same.
4. No new public `information_mode` switch was introduced in this pass because the immediate payoff came from removing unnecessary native Hessian work in existing estimate-only callers.

### Kernel outcomes

#### Logistic GLMM

What changed:

- the Hessian path stopped materializing repeated full-length group score vectors where blockwise accumulators were sufficient
- the native fitter now skips Hessian / covariance work entirely when `estimate_only = TRUE`

Result:

- estimates were unchanged on the direct before/after check
- the estimate-only fit path got meaningfully faster
- the direct Hessian path stayed at least neutral and slightly cleaner structurally

#### Clogit-plus-GLMM

What changed:

- the concordant Hessian path received the same blockwise cleanup pattern as logistic
- the native fitter now early-returns for `estimate_only = TRUE`

Result:

- estimates were unchanged on the direct before/after check
- both the estimate-only fit path and direct Hessian path improved modestly

#### Poisson GLMM

What changed:

- retained: estimate-only early return
- rejected: grouped Hessian rewrite

Result:

- the estimate-only path improved because it stopped computing unused Hessian outputs
- the attempted grouped Hessian cleanup did not improve the direct Hessian benchmark enough to justify keeping it
- the final retained Poisson state is the original Hessian path plus the estimate-only early return

#### Weibull frailty

What changed:

- no new retained Hessian rewrite in this pass

Result:

- the previously reverted baseline remains the correct retained version
- earlier Weibull tightening attempts had already shown that more aggressive Hessian math rewrites can preserve the estimates while still losing on runtime

### Phase 4 benchmarks

The table below reflects the retained-vs-attempted outcomes for this pass. For Poisson, the Hessian rewrite row is intentionally marked as a rejected attempt; the retained Poisson state keeps only the estimate-only early return.

| Kernel | Path | Before median (s) | After median (s) | Estimate parity | Decision |
|---|---:|---:|---:|---|---|
| `logistic_glmm` | estimate-only fit | `0.0090` | `0.0060` | unchanged | keep |
| `logistic_glmm` | Hessian | `0.0010` | `0.0010` | unchanged | keep |
| `clogit_plus_glmm` | estimate-only fit | `0.0230` | `0.0200` | unchanged | keep |
| `clogit_plus_glmm` | Hessian | `0.00175` | `0.00170` | unchanged | keep |
| `poisson_glmm` | estimate-only fit | `0.0020` | `0.0010` | unchanged | keep |
| `poisson_glmm` | Hessian cleanup attempt | `0.0010` | `0.0010` | unchanged | reject rewrite, keep early return only |
| `weibull_frailty` | Hessian tightening attempt | `0.00330` | `0.00744` | unchanged | keep reverted baseline |

Interpretation:

- Phase 4 produced its clearest practical wins by skipping unnecessary Hessian work in estimate-only callers
- logistic and clogit also benefited from direct Hessian cleanup
- Poisson again resisted the more aggressive shared-pattern Hessian rewrite and should continue to be treated separately
- Weibull remains the main case where “cleaner math” did not translate into a faster Hessian path

### Test coverage notes

Relevant repo coverage for these kernels:

- [EDI/tests/testthat/test-glmm-cpp-equivalence.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/tests/testthat/test-glmm-cpp-equivalence.R) contains canonical-package equivalence coverage for logistic GLMM and Poisson GLMM against `lme4::glmer`
- [EDI/tests/testthat/test-weibull-frailty.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/tests/testthat/test-weibull-frailty.R) exercises the Weibull frailty path

There is still no dedicated canonical-package `testthat` path in this repo for `clogit_plus_glmm`.

I attempted package-backed execution of the relevant `testthat` files during this pass, but `pkgload::load_all("EDI")` remained stuck in the local rebuild path long enough that it was not a reliable verification route for this session. The direct native before/after parity checks and the isolated benchmark harnesses are therefore the authoritative Phase 4 result.


## Phase 5 outcome

Phase 5 targeted the lower-priority adjacent-category logit kernel in `EDI/src/fast_adjacent_category_logit.cpp`.

Scope of the attempted cleanup:

1. cache repeated `LinSpaced(...)` category grids and score-offset vectors
2. remove per-row tiny temporaries in the objective/score path
3. simplify Hessian cross-block updates by accumulating directly into blocks instead of constructing short temporary vectors
4. benchmark the direct objective/score path, direct Hessian path, and fit path before deciding whether to keep the changes

### SIMD constraint

This pass was intentionally confined to `fast_adjacent_category_logit.cpp`. No shared helper used by the SIMD-sensitive matrix kernels was changed, and none of the retained Phase 1-4 GLMM optimizations were touched.

### Numerical checks

The attempted cleanup preserved the direct numerical results on the benchmark problem:

- `nll = 1392.8592205295745`
- `score_sum = 45.46632076388376`
- `h_trace = -7778.655883166635`
- `fit_param_sum = -0.0542012540459799`

I also ran the adjacent-category canonical-package parity check directly against `VGAM::vglm(..., family = VGAM::acat(parallel = TRUE))` on the same scenario used in [EDI/tests/testthat/test-rcpp-fitting-equivalence.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/tests/testthat/test-rcpp-fitting-equivalence.R:289):

- `max_abs_beta_diff = 1e-16`
- `max_abs_vcov_diag_diff = 1.56511314e-08`

So the attempted optimization did not change the estimates or the VGAM equivalence story.

### Benchmark result

The attempted cleanup did not produce a speedup and was therefore reverted.

| Kernel | Path | Before median (s) | Attempted after median (s) | Decision |
|---|---:|---:|---:|---|
| `adjacent_category_logit` | direct objective path | `0.00010` | `0.00011` | revert |
| `adjacent_category_logit` | direct score path | `0.00010` | `0.00010` | revert |
| `adjacent_category_logit` | direct Hessian path | `0.00015` | `0.00020` | revert |
| `adjacent_category_logit` | fit path | `0.0075` | `0.0080` | revert |

For completeness, after the revert the retained file returned to the original benchmark behavior:

- direct objective path: `0.00010s`
- direct score path: `0.00010s`
- direct Hessian path: `0.00020s`
- fit path: `0.0090s`

The retained conclusion is therefore the same as for the rejected Poisson and Weibull sub-passes:

- the attempted cleanup was estimate-stable
- it did not improve runtime enough to justify keeping
- the codebase should keep the original adjacent-category implementation for now

Interpretation:

- this kernel is small enough that the extra bookkeeping from cached vectors and blockwise Hessian accumulation did not beat the compiler’s handling of the simpler baseline
- further work here is unlikely to matter globally unless a profiler later shows this path growing in total share


## Additional sweep outcome

This pass covered the non-GLMM hotspots from the broader `perf` sweep:

- `EDI/src/optimal_design_search.cpp`
- `EDI/src/pair_dist_helpers.cpp`
- `EDI/src/fast_zinb.cpp`
- `EDI/src/fast_zero_augmented_poisson.cpp`
- `EDI/src/fast_negbin_regression.cpp`

The retained outcome was mixed:

- keep the low-risk accessor/data-path cleanup in `optimal_design_search.cpp`
- keep the low-risk pointer-based rewrite in `pair_dist_helpers.cpp`
- revert the attempted Hessian/objective cleanups in `fast_zinb.cpp`
- revert the attempted Hessian cleanup in `fast_zero_augmented_poisson.cpp`
- revert the attempted objective cleanup in `fast_negbin_regression.cpp`

### SIMD constraint

This pass did not touch the SIMD-sensitive GLMM kernels or the shared vectorized helpers that were already producing wins in earlier phases. The retained changes are confined to:

- `optimal_design_search.cpp`
- `pair_dist_helpers.cpp`

The reverted count-model experiments left:

- `fast_zinb.cpp`
- `fast_zero_augmented_poisson.cpp`
- `fast_negbin_regression.cpp`

bit-for-bit back at their pre-pass source state.

### What was kept

#### Pair-distance matrix

`compute_pair_distance_matrix_cpp(...)` was rewritten to use raw column-major pointers into the R matrices/vectors instead of `Rcpp::Matrix::operator()` and `Rcpp::Vector::operator[]` inside the innermost loop.

Why it was kept:

- the numerical result was unchanged on the benchmark matrix
- the speedup was large and unambiguous
- the rewrite is local and low-risk

Observed benchmark:

- upper-triangle sum before: `2087702.3513132515`
- upper-triangle sum after: `2087702.3513132515`
- median time before: `0.0080s`
- median time after: `0.0020s`

#### D-optimal / A-optimal search

`d_optimal_search_cpp(...)` and `a_optimal_search_cpp(...)` were tightened around the hot swap-scan loop by:

- caching matrix diagonals
- using raw access to matrix storage for `P(i, j)` and `H(i, j)` inside the nested scan
- avoiding repeated accessor overhead in the scalar delta calculations

Why it was kept:

- the stochastic search structure and output shape were preserved
- the benchmark showed modest but real wins in both search kernels
- the changes are localized to scalar access, not algorithmic behavior

Observed benchmark:

- each returned column still had exactly `n_T` treated units
- `d_optimal_search_cpp` median time: `0.0065s` before, `0.0055s` after
- `a_optimal_search_cpp` median time: `0.0130s` before, `0.0120s` after

### What was rejected

#### Zero-inflated negative binomial

The attempted `fast_zinb.cpp` change rewrote the Hessian accumulation around block updates and row expressions to remove some temporary-vector and nested-loop overhead.

Why it was rejected:

- the Hessian fingerprint and fitted parameter sum were unchanged
- direct Hessian time did not improve
- end-to-end fit time got materially worse

Observed benchmark:

- Hessian trace before/after: `-2182.2789002836535`
- fitted parameter sum before/after: `-0.8550610923931`
- median Hessian time before: `0.0010s`
- median Hessian time after: `0.0010s`
- median fit time before: `0.0835s`
- median fit time after: `0.1470s`

Decision:

- revert

#### Zero-augmented Poisson

The attempted `fast_zero_augmented_poisson.cpp` change similarly replaced repeated block extraction and temporary-vector creation in the Hessian path.

Why it was rejected:

- the Hessian trace and fitted coefficient sum were unchanged
- the direct Hessian benchmark was too small to show a win
- the fit path got slower

Observed benchmark:

- Hessian trace before/after: `-3414.94388334212`
- fitted coefficient sum before/after: `-0.5370720648920692`
- median Hessian time before: effectively `0`
- median Hessian time after: effectively `0`
- median fit time before: `0.0010s`
- median fit time after: `0.0035s`

Decision:

- revert

#### Negative binomial with variance

The attempted `fast_negbin_regression.cpp` change replaced `R::dnbinom_mu(...)` with an inlined log-likelihood expression and reused denominator/log terms inside the objective loop.

Why it was rejected:

- the score fingerprint and fitted parameter sum were unchanged up to floating-point noise
- the direct score timing did not improve materially
- the fit path got slower

Observed benchmark:

- score sum before: `-4.7075844155759334`
- score sum after: `-4.7075844155759601`
- fitted coefficient-plus-log-theta sum before: `0.9436344156551447`
- fitted coefficient-plus-log-theta sum after: `0.9436344156551441`
- median score time before: `0.0010s`
- median score time after: `0.0010s`
- median fit time before: `0.0070s`
- median fit time after: `0.0085s`

Decision:

- revert

### Additional sweep benchmarks

| Kernel | Path | Before median (s) | After median (s) | Outcome |
|---|---:|---:|---:|---|
| `pair_distance_matrix` | direct distance build | `0.0080` | `0.0020` | keep |
| `d_optimal_search` | direct search | `0.0065` | `0.0055` | keep |
| `a_optimal_search` | direct search | `0.0130` | `0.0120` | keep |
| `zinb` | Hessian | `0.0010` | `0.0010` | revert |
| `zinb` | fit | `0.0835` | `0.1470` | revert |
| `zero_augmented_poisson` | Hessian | `0.0000` | `0.0000` | revert |
| `zero_augmented_poisson` | fit | `0.0010` | `0.0035` | revert |
| `negbin_with_var` | score/objective path | `0.0010` | `0.0010` | revert |
| `negbin_with_var` | fit | `0.0070` | `0.0085` | revert |

Interpretation:

- the low-level data-access hotspots were good optimization targets
- the Hessian/special-function count-model kernels were not improved by these particular structural cleanups
- for ZINB, ZAP, and negative binomial, the current code should stay on the simpler pre-pass implementation until a more targeted profiler-guided rewrite is available

### Test coverage notes

Relevant canonical-package coverage already present in the repo:

- [EDI/tests/testthat/test-rcpp-fitting-equivalence.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/tests/testthat/test-rcpp-fitting-equivalence.R:55) for `fast_neg_bin_with_var_cpp` vs `MASS::glm.nb`
- [EDI/tests/testthat/test-rcpp-fitting-equivalence.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/tests/testthat/test-rcpp-fitting-equivalence.R:329) for `fast_zinb_cpp` vs `glmmTMB`
- [EDI/tests/testthat/test-rcpp-fitting-equivalence.R](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/tests/testthat/test-rcpp-fitting-equivalence.R:348) for `fast_zero_augmented_poisson_cpp` vs `glmmTMB`

Because the count-model edits were all reverted, no retained package behavior changed in those paths from this pass. There is no corresponding canonical-package parity test in the repo for the D-optimal search or pair-distance helper kernels.


## Current Retained State And Next Work

### Retained state

The currently retained optimization state is:

- ordinal CLMM:
  - keep threshold reconstruction removal in `fast_ordinal_clmm.cpp`
  - keep tiny-allocation removal and workspace reuse in `_glmm_engine.h`
  - do not keep the dedicated ordinal specialization
- logistic GLMM:
  - keep direct group accumulation in objective/gradient loops
  - keep persistent scratch buffers
  - keep Phase 3 math tightening
  - keep Phase 4 Hessian cleanup
  - keep `estimate_only` early return
- Weibull frailty:
  - keep Phase 1/2 structural direct-accumulation changes
  - keep persistent scratch buffers
  - do not keep the slower Phase 3 or additional Hessian tightening attempts
- clogit-plus-GLMM:
  - keep direct accumulation changes
  - keep persistent scratch buffers
  - keep Phase 3 tightening
  - keep Phase 4 Hessian cleanup
  - keep `estimate_only` early return
- Poisson GLMM:
  - keep `estimate_only` early return
  - do not keep the structural shared-pattern rewrite
  - do not keep the grouped Hessian rewrite
- adjacent-category logit:
  - keep the original implementation
  - do not keep the attempted low-priority cleanup
- additional non-GLMM sweep:
  - keep `optimal_design_search.cpp` access-path cleanup
  - keep `pair_dist_helpers.cpp` pointer-based rewrite
  - keep original `fast_zinb.cpp`, `fast_zero_augmented_poisson.cpp`, and `fast_negbin_regression.cpp`

### Next work

The next best work from the current retained state is:

1. Re-profile ordinal CLMM on the retained tree and decide whether any further change is justified beyond the already-kept threshold/allocation fixes.
2. Treat Poisson GLMM as its own optimization problem rather than forcing it into the logistic/Weibull/clogit pattern.
3. If more non-GLMM work is needed, continue with measured access-path cleanup like `optimal_design_search.cpp` and `pair_dist_helpers.cpp` before reopening reverted count-model experiments.
4. If ZINB or ZAP become urgent again, gather a more focused profiler trace around special-function cost before attempting another structural Hessian rewrite.
