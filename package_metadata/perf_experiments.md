# Performance Experiments Report

Date: 2026-05-13

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

## Concrete Next Steps

1. Implement a dedicated optimized ordinal CLMM objective path to bypass the generic per-row dynamic derivative machinery.
2. Refactor Weibull and logistic GLMM objective/gradient loops to remove `post_k_expanded`.
3. Re-run the same `perf` procedure after each change set and compare:
   - kernel median elapsed time,
   - top `perf` symbols,
   - and allocation-heavy symbols like `malloc`, `cfree`, and `R_gc_internal`.

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


DONE: 


### 4. Weibull frailty posterior-expansion and gradient accumulation loop

Location:

- `EDI/src/fast_weibull_frailty.cpp` lines around the second `for (int k = 0; k < n_nodes; ++k)` loop in `operator()`

Why it ranks fourth:

- It materializes `post_k_expanded` for the whole dataset per node.
- It performs additional group reductions and then a weighted matrix-vector multiplication.
- This is a classic case where fusing reductions can remove memory traffic and temporaries.

Expected payoff:

- Very high.

### 5. Logistic GLMM node sweep building `mu_all_k_vec` and `log_terms_mat`

Location:

- `EDI/src/fast_logistic_glmm.cpp` around the first `for (int k = 0; k < n_nodes; ++k)` loop in `operator()`

Why it ranks fifth:

- `logistic_glmm` is 22.4% of total current kernel time.
- `perf` shows dominant time in `__log1p_fma`, `plogis_array_safe`, and GEMV.
- This loop drives those calls.

Expected payoff:

- High.

### 6. Logistic GLMM posterior-expansion and gradient accumulation loop

Location:

- `EDI/src/fast_logistic_glmm.cpp` around the second `for (int k = 0; k < n_nodes; ++k)` loop in `operator()`

Why it ranks sixth:

- Same structural issue as Weibull: `post_k_expanded` is built across all observations per node.
- The gradient then does another whole-matrix pass.

Expected payoff:

- High.


### 7. Clogit-plus-GLMM concordant node sweep

Location:

- `EDI/src/fast_clogit_plus_glmm.cpp` around the first concordant `for (int k = 0; k < n_nodes; ++k)` loop in `operator()`

Why it ranks seventh:

- The kernel is smaller than logistic GLMM but structurally very similar.
- Any pattern learned while optimizing logistic GLMM should transfer here.

Expected payoff:

- Medium-high.

### 8. Clogit-plus-GLMM concordant posterior-expansion and gradient accumulation loop

Location:

- `EDI/src/fast_clogit_plus_glmm.cpp` around the second concordant `for (int k = 0; k < n_nodes; ++k)` loop in `operator()`

Why it ranks eighth:

- Same full-length `post_k_expanded` materialization pattern.
- Same opportunity to accumulate group contributions without constructing a dense expanded vector.

Expected payoff:

- Medium-high.



### 9. Poisson GLMM node sweep

Location:

- `EDI/src/fast_poisson_glmm.cpp` around the first `for (int k = 0; k < n_nodes; ++k)` loop in `operator()`

Why it ranks ninth:

- `perf` shows this kernel is dominated by the objective loop and GEMV.
- The kernel is smaller than the three leaders, so its absolute payoff is lower even if the local loop is hot.

Expected payoff:

- Medium.

### 10. Poisson GLMM posterior-expansion and gradient accumulation loop

Location:

- `EDI/src/fast_poisson_glmm.cpp` around the second `for (int k = 0; k < n_nodes; ++k)` loop in `operator()`

Why it ranks tenth:

- Same overall pattern as logistic/clogit/Weibull.
- Worth fixing after the larger kernels unless the code changes can be shared.

Expected payoff:

- Medium.