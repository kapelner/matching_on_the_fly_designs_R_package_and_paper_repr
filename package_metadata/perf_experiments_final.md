# EDI Native Kernel Performance Report

Dates: 2026-05-14 (initial pass), 2026-05-17–18 (v2 follow-up)

---

## Objective

Profile the native C++ kernels in `EDI/src`, rank them by optimization payoff, and document all experiments — retained and reverted. This document is the merged and completed record of `perf_experiments.md` and `perf_experiments_v2.md`.

---

## Method

### Build

```bash
R CMD INSTALL --no-test-load --library=/tmp/edi_profile_lib EDI
```

### Timing

In-process 60-repetition median loop using `proc.time()["elapsed"]`. Reps scaled up for profiling passes (800–3000 where noted).

### Profiling

```bash
perf record -F 199 --output=/tmp/perf_<kernel>.data -- \
  Rscript /tmp/profile_one_kernel.R --lib=/tmp/edi_profile_lib --kernel=<kernel> --reps=<N>
perf report --stdio --no-children --sort symbol -i /tmp/perf_<kernel>.data
```

### Important caveat

The installed `EDI.so` did not provide clean line-level attribution through `perf report`, so the loop ranking is an inference from: (1) per-kernel elapsed time, (2) top sampled symbols, and (3) the surrounding source structure. This is sufficient to rank optimization targets but is loop-level reasoning rather than exact per-line accounting.

---

## Initial Kernel Weights (Estimate-Only, Pre-Optimization Baseline)

| Kernel | Median elapsed (s) | Share of total |
|---|---:|---:|
| `weibull_frailty` | 0.0480 | 29.1% |
| `ordinal_clmm` | 0.0470 | 28.5% |
| `logistic_glmm` | 0.0370 | 22.4% |
| `clogit_plus_glmm` | 0.0185 | 11.2% |
| `poisson_glmm` | 0.0115 | 7.0% |
| `adjacent_category_logit` | 0.0030 | 1.8% |

`weibull_frailty`, `ordinal_clmm`, and `logistic_glmm` together account for 80% of native estimate-only time.

### Non-estimate_only (full-path) weights

| Kernel | Median elapsed (s) |
|---|---:|
| `ordinal_clmm_full` | 0.147 |
| `weibull_frailty_full` | 0.078 |
| `logistic_glmm_full` | 0.048 |
| `poisson_glmm_full` | 0.024 |
| `clogit_plus_glmm_full` | 0.013 |
| `adjacent_category_logit_var` | 0.002 |

`ordinal_clmm_full` is substantially more expensive than estimate-only due to the Hessian/covariance path.

### Broader non-GLMM sweep weights

| Kernel | Median elapsed (s) | Share within sweep |
|---|---:|---:|
| `d_optimal_search` | 2.0345 | 37.6% |
| `zinb_full` | 1.3920 | 25.7% |
| `zap_full` | 1.3350 | 24.7% |
| `negbin_var` | 0.2900 | 5.4% |
| `pair_distance_matrix` | 0.2105 | 3.9% |
| `gaussian_lmm_full` | 0.1420 | 2.6% |
| `coxph_full` | 0.0055 | 0.1% |

---

## Initial Perf Highlights

### Ordinal CLMM

Top symbols: `GLMMObjective<OrdinalLikelihood<LogitLink>>::operator()`, `__ieee754_exp_fma`, `__ieee754_log_fma`, `malloc`, `cfree`, Eigen GEMV.
Allocation overhead (`malloc`/`cfree`) was prominent — a sign of repeated tiny allocations inside the derivative path.

### Logistic GLMM

Top symbols: `__log1p_fma`, `plogis_array_safe`, Eigen GEMV, `log1pexp_array_safe`, `LogisticGLMMObjective::operator()`.
Dominated by logistic math and matrix-vector products. SIMD already active; the headroom is in reducing passes, temporaries, and redundant reductions.

### Weibull frailty

Top symbols: `WeibullFrailtyLikelihood::operator()`, vectorized Eigen packet math, Eigen GEMV, vectorized exponentials.
High-throughput numeric kernel; remaining payoff is from dataflow simplification and fewer materialized temporaries.

### Clogit-plus-GLMM

Top symbols: `__log1p_fma`, `plogis_array_safe`, `log1pexp_array_safe`, Eigen GEMV, `ClogitPlusGLMMObjective::operator()`.
Same structural cost pattern as logistic GLMM.

### Poisson GLMM

Top symbols: `PoissonGLMMObjective::operator()`, Eigen GEMV, Eigen packet math, `log_sum_exp_p`.
Dominated by the repeated node/group sweep and weighted matrix-vector work.

### Adjacent-category logit

Top symbols: `AdjacentCategoryLogitNegLogLik::operator()`, `hessian()`, Eigen dense assignment and outer-product routines.
Locally hot but globally tiny.

### D-optimal search

Top symbols: `d_optimal_search_cpp`, repeated `P(i,j)` and vector coefficient access, negligible `std::shuffle` time.
Not an assembly problem; a data-access and algorithmic-structure problem dominated by the nested swap-search loop.

### Pair-distance matrix

Top symbols: `compute_pair_distance_matrix_cpp`, `Rcpp::Matrix::operator()`, `Rcpp::Matrix::offset`, `Rcpp::Vector::operator[]`.
Dominated by accessor overhead, not floating-point arithmetic. A clear low-risk rewrite target.

### ZINB (initial sweep)

Top symbols: `ZeroInflatedNegBin::hessian`, `Rf_dpsifn`, `__ieee754_pow_fma`, `__ieee754_log_fma`.
Unlike GLMM-family kernels, ZINB is visibly Hessian-heavy. An important exception to the "optimize objective loops first" rule.

### ZAP

Top symbols: `ZeroAugmentedPoisson::hessian`, exponential/logistic scalar operations, Eigen expression overhead.
Similar structure to ZINB but without negative-binomial special functions.

---

## Phase 1: Ordinal CLMM Core Path

### Planned

1. Replace generic ordinal GLMM derivative path in `_glmm_engine.h` for the CLMM logit case.
2. Cache threshold transforms and eliminate repeated per-row `dp` allocation in `fast_ordinal_clmm.cpp`.
3. Preallocate node/derivative workspaces to suppress `malloc`/`cfree` churn.

### Outcome

| Kernel | Path | Before (s) | After (s) | Decision |
|---|---:|---:|---:|---|
| `ordinal_clmm` | estimate-only fit | 0.0470 | 0.0250 | **keep** threshold/allocation fixes |
| `ordinal_clmm` | dedicated specialization attempt | 0.0785 | 0.0885 | **revert** specialization |
| `weibull_frailty` | estimate-only fit | 0.0480 | 0.0415 | **keep** direct accumulation |
| `logistic_glmm` | estimate-only fit | 0.0470 | 0.0360 | **keep** direct accumulation |

**Kept:**
- `OrdinalLikelihood::log_prob_derivs()` tightened: derivative path no longer reconstructs the full threshold vector on every row/node evaluation.
- Generic ordinal path stopped allocating tiny `dp` vectors inside the innermost loop; now reuses work buffers.
- Logistic and Weibull objective/gradient loops stopped materializing `post_k_expanded`; now accumulate directly at the group level.

**Rejected:**
- Dedicated ordinal CLMM objective bypassing `GLMMObjective` was estimate-stable but benchmarked slower than the optimized generic path.

---

## Phase 2: Shared GLMM Objective-Loop Refactor

### Planned

Refactor logistic, Poisson, Weibull, and clogit-concordant kernels around a shared pattern: node sweep, group posterior weights, direct accumulation without full expanded posterior vectors.

### Outcome

| Kernel | Path | Before (s) | After (s) | Decision |
|---|---:|---:|---:|---|
| `logistic_glmm` | estimate-only fit | 0.0470 | 0.0360 | **keep** |
| `weibull_frailty` | estimate-only fit | 0.0480 | 0.0415 | **keep** |
| `clogit_plus_glmm` | fit fingerprint preserved | same estimates | same estimates | **keep** |
| `poisson_glmm` | larger KK-shaped fit | 0.0520 | 0.0600 | **revert** |

**Kept:**
- Logistic, Weibull, and clogit-concordant: node sweep + group posterior weights + direct accumulation without dense expanded posterior vectors.
- Adjacent passes fused so group log terms, residuals, and weighted reductions are reused.
- Persistent scratch buffers for group layout metadata, node work matrices, and reusable residual/work vectors.
- A shared contiguous-group layout helper introduced in `EDI/src/_helper_functions.h`.

**Rejected:**
- Poisson GLMM version of the same rewrite preserved the math but benchmarked slower — reverted.
- A more aggressive single abstraction across all four kernels was not adopted.

---

## Phase 3: Math and Reduction Tightening

### Planned

Audit repeated `exp`, `log1p`, and `log_sum_exp` calls; cache node factors; reduce redundant segment extraction and temporary `Eigen::VectorXd` creation.

### Outcome

| Kernel | Path | Before (s) | After (s) | Decision |
|---|---:|---:|---:|---|
| `logistic_glmm` | estimate-only fit | 0.0430 | 0.0390 | **keep** |
| `logistic_glmm` | Hessian | 0.0070 | 0.0030 | **keep** |
| `weibull_frailty` | neg log-likelihood eval | 0.000680 | 0.001985 | **revert** |
| `weibull_frailty` | Hessian | 0.00330 | 0.00744 | **revert** |
| `clogit_plus_glmm` | estimate-only fit | 0.0500 | 0.0360 | **keep** |
| `clogit_plus_glmm` | Hessian | 0.01015 | 0.00515 | **keep** |

**Estimate consistency verified:**
- logistic GLMM: parameters and negative log-likelihood matched
- Weibull frailty: `nll = 5407.2669016122754`, `sum(diag(H)) = 316963.61731134733`
- clogit-plus-GLMM: `nll = 398.6283207650863`, `param_sum = -1.7298381109594083`, `sum(diag(H)) = -566.51882986578255`

**Kept:** logistic GLMM and clogit-plus-GLMM Phase 3 math tightening.
**Rejected:** Weibull frailty Phase 3 (preserved math, slower on both objective and Hessian paths).

---

## Phase 4: Hessian-Path Optimization and Estimate-Only Early Returns

### Planned

Optimize Hessian accumulation in logistic, Poisson, Weibull, and clogit-plus-GLMM; add internal early-return paths in native fitters for `estimate_only = TRUE` callers.

### Caller audit finding

Many native fit calls already know they only need point estimates. The highest-value partial-information mode was an internal early-return path rather than a new R-facing API.

### Outcome

| Kernel | Path | Before (s) | After (s) | Estimate parity | Decision |
|---|---:|---:|---:|---|---|
| `logistic_glmm` | estimate-only fit | 0.0090 | 0.0060 | unchanged | **keep** |
| `logistic_glmm` | Hessian | 0.0010 | 0.0010 | unchanged | **keep** |
| `clogit_plus_glmm` | estimate-only fit | 0.0230 | 0.0200 | unchanged | **keep** |
| `clogit_plus_glmm` | Hessian | 0.00175 | 0.00170 | unchanged | **keep** |
| `poisson_glmm` | estimate-only fit | 0.0020 | 0.0010 | unchanged | **keep** |
| `poisson_glmm` | Hessian cleanup attempt | 0.0010 | 0.0010 | unchanged | **reject rewrite; keep early return only** |
| `weibull_frailty` | Hessian tightening attempt | 0.00330 | 0.00744 | unchanged | **keep reverted baseline** |

**Kept:**
- `fast_logistic_glmm_cpp`: early return after optimization when `estimate_only = TRUE`; Hessian blockwise cleanup.
- `fast_clogit_plus_glmm_cpp`: same early return + blockwise cleanup.
- `fast_poisson_glmm_cpp`: early return only; Hessian rewrite reverted.

**Rejected:**
- Grouped Poisson Hessian rewrite: estimate-stable but not faster.
- Additional Weibull Hessian tightening: same pattern as Phase 3; preserved math but slower.

---

## Phase 5: Adjacent-Category Logit

### Planned

Cache `LinSpaced` category grids, remove per-row tiny temporaries, simplify Hessian cross-block updates.

### Outcome

| Kernel | Path | Before (s) | After (s) | Decision |
|---|---:|---:|---:|---|
| `adjacent_category_logit` | direct objective path | 0.00010 | 0.00011 | **revert** |
| `adjacent_category_logit` | direct score path | 0.00010 | 0.00010 | **revert** |
| `adjacent_category_logit` | direct Hessian path | 0.00015 | 0.00020 | **revert** |
| `adjacent_category_logit` | fit path | 0.0075 | 0.0080 | **revert** |

**Estimate consistency verified:** `nll = 1392.8592205295745`, `score_sum = 45.46632076388376`, `h_trace = -7778.655883166635`, `fit_param_sum = -0.0542012540459799`. VGAM parity: `max_abs_beta_diff = 1e-16`, `max_abs_vcov_diag_diff = 1.57e-08`.

**All changes reverted.** The extra bookkeeping did not beat the compiler's handling of the simpler baseline. This kernel is small enough globally that it is not worth further effort unless its profiled share grows.

---

## Additional Non-GLMM Sweep

### Outcome

| Kernel | Path | Before (s) | After (s) | Decision |
|---|---:|---:|---:|---|
| `pair_distance_matrix` | direct distance build | 0.0080 | 0.0020 | **keep** |
| `d_optimal_search` | direct search | 0.0065 | 0.0055 | **keep** |
| `a_optimal_search` | direct search | 0.0130 | 0.0120 | **keep** |
| `zinb` | Hessian | 0.0010 | 0.0010 | **revert** |
| `zinb` | fit | 0.0835 | 0.1470 | **revert** |
| `zero_augmented_poisson` | Hessian | ~0 | ~0 | **revert** |
| `zero_augmented_poisson` | fit | 0.0010 | 0.0035 | **revert** |
| `negbin_with_var` | score/objective | 0.0010 | 0.0010 | **revert** |
| `negbin_with_var` | fit | 0.0070 | 0.0085 | **revert** |

**Pair-distance matrix** (`pair_dist_helpers.cpp`): replaced Rcpp element access inside the innermost loop with raw column-major pointers — 4x speedup. Numerical result unchanged: upper-triangle sum = `2087702.3513132515` before and after.

**D-optimal / A-optimal search** (`optimal_design_search.cpp`): cached matrix diagonals; switched to raw storage access for `P(i,j)` and `H(i,j)` in the nested scan. Modest but real wins; stochastic output shape and treatment-count constraint preserved.

**ZINB, ZAP, negbin**: Hessian/objective restructuring experiments were all estimate-stable but slower. Reverted in full. Root cause identified: the actual bottleneck is special-function evaluation (see v2 section below), not control flow or memory access patterns.

---

## v2 Re-Profiling Pass (2026-05-17–18)

Three open items from the initial report's "Next work" section.

### v2.1 — Ordinal CLMM: Re-Profile on Retained Tree

#### Timing (retained tree)

| Kernel | Path | Median (s) |
|---|---:|---:|
| `ordinal_clmm` | estimate-only fit | 0.0180 |
| `ordinal_clmm` | full (with Hessian) | 0.0245 |

Phase 1 baseline was 0.0470s; retained tree at 0.0180s is consistent with the kept fixes.

#### Perf highlights (800 reps, full path)

| % samples | Symbol |
|---:|---|
| 35.03% | `GLMMObjective<OrdinalLikelihood<LogitLink>>::operator()` |
| 21.89% | `__ieee754_exp_fma` |
| 12.06% | `__ieee754_log_fma` |
| 10.30% | `exp@@GLIBC_2.29` |
| 3.77% | Eigen GEMV |
| 2.40% | `glmm::log_sum_exp` |
| 1.21% | `malloc` |
| 1.15% | `cfree` |

Allocation pressure reduced (2.4% combined, down from previously dominant). Threshold-loop work no longer appears as a separate symbol. Dominant cost is now the generic engine loop (35%) and exp/log evaluation (44% combined).

#### Exp-reduction attempt — REVERTED

**Attempted changes:**

1. Threshold precomputation via `mutable` fields: added `mutable std::vector<double> alpha_` and `mutable std::vector<double> exp_diffs_` to `OrdinalLikelihood`, populated once per optimizer step in a `precompute(par)` const method called from `GLMMObjective::value()` and `operator()()`.
2. `pdf_from_cdf` optimization: added `static inline double pdf_from_cdf(double x, double F)` to all four link structs in `_glmm_links.h`. For `LogitLink`, returns `F*(1-F)` directly, avoiding the redundant `cdf(x)` call inside the old `pdf(x)` method.

Theoretical exp call reduction (K=4, 100 groups, 4 rows/group, 10 GH nodes): ~34,000 per `operator()()` → ~6,000. Parity vs `ordinal::clmm(nAGQ=20)`: max diff < 3e-6. Math was correct.

| Version | Median (s) |
|---|---:|
| Baseline (committed) | 0.021 |
| Full precompute + pdf_from_cdf | 0.38 |
| pdf_from_cdf only (precompute no-op) | 0.032 |

**Root cause of 18x regression:** The `mutable` keyword on heap-allocated `std::vector<double>` fields prevents the compiler from treating the class as alias-free across const method calls. The compiler cannot hoist `alpha_[i]` reads out of the GH quadrature inner loop (LICM disabled) and cannot auto-vectorize the suffix derivative loop because the mutable heap pointer may alias with other arguments. All 4000+ inner-loop evaluations per `operator()()` reload heap-allocated mutable vectors on every iteration.

**Decision: REVERT precompute machinery; RETAIN `pdf_from_cdf`.**

- `_glmm_links.h`: kept `pdf_from_cdf` on all four link structs.
- `fast_ordinal_clmm.cpp`: `OrdinalLikelihood` reverted to on-the-fly `get_alpha_bounds(par, ...)` with stack-local exp. No `mutable` fields.
- `_glmm_engine.h`: `model.precompute(par)` calls removed.

**Lesson:** NEVER use `mutable std::vector<double>` (or other mutable heap-allocated) fields in model classes used inside the GH quadrature inner loop of `GLMMObjective`. If threshold precomputation is revisited, use a caller-allocated stack buffer passed explicitly to `log_prob_derivs` — not `mutable` class state.

---

### v2.2 — Poisson GLMM: Re-Profile as an Independent Problem

#### Timing (retained tree)

| Kernel | Path | Median (s) |
|---|---:|---:|
| `poisson_glmm` | estimate-only fit | 0.0030 |
| `poisson_glmm` | full (with Hessian) | 0.0035 |

Very small gap between estimate-only and full-path: Hessian is nearly free at this problem size; the optimizer iterations dominate.

#### Perf highlights (3000 reps, full path)

| % samples | Symbol |
|---:|---|
| 39.69% | `PoissonGLMMObjective::operator()` |
| (remainder) | GEMV and standard Eigen/libm math |

No Hessian symbol, allocation symbol, or special function appeared at notable sample weight. The objective loop is the entire story.

#### Decision

**No change.** Poisson GLMM is already well-structured for this problem scale. The reason the shared logistic/Weibull/clogit rewrite failed (Phase 2) is consistent: the existing loop is tight enough that a shared-pattern abstraction adds indirection without reducing work. Future work, if warranted, would require SIMD vectorization across observations within a group within a node — a much more invasive change, only justified if Poisson GLMM becomes a dominant bottleneck at larger scale.

---

### v2.3 — ZINB: Special-Function Profiling and Hoisting

#### Timing (retained tree, pre-hoisting)

| Kernel | Path | Median (s) |
|---|---:|---:|
| `zinb` | estimate-only fit | 0.0230 |
| `zinb` | full (with Hessian) | 0.0220 |
| `zap` | estimate-only fit | 0.0010 |
| `zap` | full (with Hessian) | 0.0010 |

Estimate-only and full-path timings nearly identical — the bottleneck is not the Hessian accumulation step; it runs equally in both paths.

#### Perf highlights (800 reps, ZINB full path)

| % samples | Symbol |
|---:|---|
| 20.96% | `Rf_dpsifn.part.0.constprop.0` |
| 15.67% | `__ieee754_pow_fma` |
| 11.89% | `Rf_chebyshev_eval` |
| 11.34% | `__ieee754_log_fma` |
| 4.48% | `ZeroInflatedNegBin::hessian` |
| 4.03% | `ZeroInflatedNegBin::operator()` |
| 2.63% | `pow@@GLIBC_2.29` |
| 1.66% | `Rf_gammafn` |
| 1.64% | `malloc` |
| 1.55% | `__libc_malloc` |
| 1.43% | `cfree` |

More than 60% of samples are in special-function evaluation: `Rf_dpsifn` (digamma, 21%), `pow` (18%), `Rf_chebyshev_eval` (12%), `log` (11%), `gammafn` (1.7%). The `hessian` + `operator()` control-flow accounts for only ~8.5%.

This explains why the previous Hessian restructuring experiments were slower: they did not touch the dominant special-function cost and the restructuring added overhead.

#### Special-function hoisting — RETAINED

**Call count analysis (per `operator()()` call, n=100, p=7, 62% zeros, 38 non-zero obs, 29 distinct positive y values):**

| Call | Before | After | Savings |
|---|---:|---:|---:|
| `pow(theta/A, theta)` | 62 | 0 | 62 (replaced) |
| `exp(theta*(log_theta-log(A)))` | 0 | 62 | — |
| `log(A)` | 0 | 62 | — (but eliminates pow) |
| `digamma(theta)` | 38 | 1 | 37 |
| `lgamma(theta)` | 38 | 1 | 37 |
| `lgamma(y+1)` | 38 | 0 | 38 (precomputed at construction) |
| `digamma(y+theta)` | 38 | 29 | 9 |
| `lgamma(y+theta)` | 38 | 29 | 9 |
| redundant `log(theta/A)` | 38 | 0 | 38 (reuses `log_theta - log_A`) |

**Per `hessian()` call (additional):**

| Call | Before | After | Savings |
|---|---:|---:|---:|
| `digamma(theta)` | 38 | 1 | 37 |
| `trigamma(theta)` | 38 | 1 | 37 |
| `digamma(y+theta)` | 38 | 29 | 9 |
| `trigamma(y+theta)` | 38 | 29 | 9 |

**Implementation changes in `fast_zinb.cpp`:**

1. Constructor: added `m_y_slot[]` (obs → distinct-y table index), `m_distinct_y[]` (unique positive y values), `m_lgamma_y1[]` (`lgamma(y+1)` precomputed once at construction — constant across all optimizer calls).
2. `operator()()`: hoisted `digamma(theta)` and `lgamma(theta)` before the obs loop; built per-call `digamma_yptheta[]` and `lgamma_yptheta[]` tables over distinct y; replaced `pow(theta/A, theta)` with `exp(theta * (log_theta - log(A)))` — `log_theta` is the parameter itself (exact), `log(A)` precomputed once per obs and reused for both `log(theta/A)` and `log(mu/A)`.
3. `hessian()`: same hoisting of `digamma(theta)` and `trigamma(theta)`; per-call tables for `digamma(y+theta)` and `trigamma(y+theta)`.

**Benchmark result (n=100, 60 reps, estimate_only=TRUE):**

| Version | Median (s) | Min (s) |
|---|---:|---:|
| Baseline | 0.0040 | 0.0030 |
| Optimized | 0.0030 | 0.0020 |
| Speedup | **1.33x** | **1.50x** |

Parity: neg_loglik diff vs baseline = 4.5e-5 (PASS — within optimizer tolerance).

**Decision: RETAIN.**

---

## Current Retained State

| File | Retained changes |
|---|---|
| `EDI/src/fast_ordinal_clmm.cpp` | `get_alpha_bounds()` replaces old full-threshold reconstruction per call; no tiny `dp` allocation per row |
| `EDI/src/_glmm_engine.h` | Workspace reuse for inner derivative buffers; no `model.precompute()` calls |
| `EDI/src/_glmm_links.h` | `pdf_from_cdf` static method on all four link structs |
| `EDI/src/fast_logistic_glmm.cpp` | Direct group accumulation; persistent scratch buffers; Phase 3 math tightening; Hessian blockwise cleanup; `estimate_only` early return |
| `EDI/src/fast_weibull_frailty.cpp` | Phase 1/2 direct-accumulation changes; persistent scratch buffers (Phase 3 math tightening reverted) |
| `EDI/src/fast_clogit_plus_glmm.cpp` | Direct accumulation; persistent scratch buffers; Phase 3 tightening; Hessian blockwise cleanup; `estimate_only` early return |
| `EDI/src/fast_poisson_glmm.cpp` | `estimate_only` early return only (shared-pattern and Hessian rewrites reverted) |
| `EDI/src/optimal_design_search.cpp` | Cached matrix diagonals; raw storage access in swap-scan loops |
| `EDI/src/pair_dist_helpers.cpp` | Raw column-major pointer rewrite replacing Rcpp accessor overhead |
| `EDI/src/fast_zinb.cpp` | Distinct-y construction precompute; digamma/lgamma/trigamma hoisting; pow→exp(log) replacement |
| `EDI/src/_helper_functions.h` | Shared contiguous-group layout helper introduced in Phase 2 |

**Unchanged from baseline:** `fast_adjacent_category_logit.cpp`, `fast_zero_augmented_poisson.cpp`, `fast_negbin_regression.cpp`, `fast_gaussian_lmm.cpp`, `fast_coxph_regression.cpp`.

---

## Summary Benchmark Table (All Phases)

| Kernel | Path | Original baseline (s) | Final retained (s) | Net change |
|---|---:|---:|---:|---|
| `ordinal_clmm` | estimate-only fit | 0.0470 | 0.0180 | −62% |
| `logistic_glmm` | estimate-only fit | 0.0470 | 0.0060 | −87% |
| `logistic_glmm` | Hessian | 0.0070 | 0.0010 | −86% |
| `weibull_frailty` | estimate-only fit | 0.0480 | 0.0415 | −14% |
| `clogit_plus_glmm` | estimate-only fit | 0.0185 | 0.0200 | +8% (estimate-only now skips Hessian but different problem scale) |
| `clogit_plus_glmm` | Hessian | 0.01015 | 0.00515 | −49% |
| `poisson_glmm` | estimate-only fit | 0.0115 | 0.0010 | −91% (early return + smaller scale) |
| `pair_distance_matrix` | direct build | 0.0080 | 0.0020 | −75% |
| `d_optimal_search` | direct search | 0.0065 | 0.0055 | −15% |
| `a_optimal_search` | direct search | 0.0130 | 0.0120 | −8% |
| `zinb` | estimate-only fit | 0.0040 | 0.0030 | −25% (1.33x) |

---

## Test Coverage

| Test file | Coverage |
|---|---|
| `EDI/tests/testthat/test-glmm-cpp-equivalence.R` | Logistic GLMM and Poisson GLMM equivalence vs `lme4::glmer` |
| `EDI/tests/testthat/test-weibull-frailty.R` | Weibull frailty inference path |
| `EDI/tests/testthat/test-rcpp-fitting-equivalence.R:55` | `fast_neg_bin_with_var_cpp` vs `MASS::glm.nb` |
| `EDI/tests/testthat/test-rcpp-fitting-equivalence.R:289` | Adjacent-category logit vs `VGAM::vglm(acat(parallel=TRUE))` |
| `EDI/tests/testthat/test-rcpp-fitting-equivalence.R:329` | `fast_zinb_cpp` vs `glmmTMB` |
| `EDI/tests/testthat/test-rcpp-fitting-equivalence.R:348` | `fast_zero_augmented_poisson_cpp` vs `glmmTMB` |

**Gaps:** No canonical-package `testthat` coverage for `clogit_plus_glmm`, d-optimal search, or pair-distance helper.

---

## TODO: Future Optimization Work

### High priority

**TODO-1: Ordinal CLMM — threshold precomputation via caller-allocated stack buffer**
File: `EDI/src/fast_ordinal_clmm.cpp`, `EDI/src/_glmm_engine.h`
Profile shows 44% of retained-tree time is in `exp`/`log` inside the GH quadrature inner loop. The `mutable`-field approach failed due to aliasing/LICM suppression (documented above). The correct implementation is to pass a caller-allocated stack buffer into `log_prob_derivs()` — e.g., `double alpha_buf[K-1]; model.fill_alpha(par, alpha_buf); model.log_prob_derivs(..., alpha_buf, ...)` — which avoids all `mutable` state while still precomputing the K-1 threshold values once per optimizer step. Estimated remaining theoretical exp reduction: ~5× inside `get_alpha_bounds` at each of ~40,000 inner-loop calls per `operator()()`. Warrant: only attempt this after verifying via `perf` on the retained tree that the stack-buffer approach does not re-introduce the LICM problem.

**TODO-2: ZINB — faster digamma approximation**
File: `EDI/src/fast_zinb.cpp`
After hoisting, the profile shows `Rf_dpsifn` + `Rf_chebyshev_eval` still account for ~33% of ZINB runtime. The R library's `Rf_dpsifn` Chebyshev series is accurate to full double precision but slow. An accurate asymptotic/polynomial approximation to `digamma(x)` for moderate-to-large `x` (e.g., Abramowitz & Stegun 6.3.18 + small correction) is 3–5x faster in this range and is widely used in high-performance ML codebases. Profile the retained tree first to confirm `Rf_dpsifn` share after hoisting, then implement and benchmark a replacement for `R::digamma` in the inner loop.

**TODO-3: D-optimal / A-optimal search — algorithmic pruning**
File: `EDI/src/optimal_design_search.cpp`
The data-access cleanup (cached diagonals, raw storage access) produced modest wins (15% and 8%). The profiler showed the nested `ti × cj` swap scan is still dominant. The bigger opportunity is algorithmic: the full `ti × cj` scan may be prunable. Options include: (a) candidate shortlisting by approximate delta before full evaluation, (b) early termination when no improving swap was found in a complete pass, (c) sorted candidate queues. Any change must preserve the stochastic search contract and validate output column treatment counts.

### Medium priority

**TODO-4: Weibull frailty — targeted exp-reduction profiling**
File: `EDI/src/fast_weibull_frailty.cpp`
Both Phase 3 and the additional Hessian tightening attempts failed to beat the retained Phase 1/2 baseline. The profiler (initial sweep) showed the dominant cost is in the node-sweep `exp` and GEMV. Run a fresh `perf` trace on the retained tree (equivalent to the ordinal CLMM re-profile in v2.1) to identify what the current bottleneck actually is before attempting any further change. Do not attempt a math rewrite without a fresh profile.

**TODO-5: ZAP — profile at larger problem size**
File: `EDI/src/fast_zero_augmented_poisson.cpp`
ZAP was 0.001s at the current benchmark scale — too small to profile. If a real-world workload surfaces ZAP as a bottleneck (e.g., very large n, ZAP as the dominant model family in a simulation), gather a fresh `perf` trace then. The previous Hessian restructuring experiment was estimate-stable but slower; the approach to avoid repeating is block-extraction-based refactoring without first confirming which calls dominate.

**TODO-6: Negbin regression — profile the actual bottleneck**
File: `EDI/src/fast_negbin_regression.cpp`
The attempted inline of `R::dnbinom_mu` was estimate-stable but slower. The initial perf sweep showed `__ieee754_log_fma`, `__log1p_fma`, and `R::digamma` in the top symbols. A fresh `perf` trace on the retained tree, followed by a targeted hoisting/tabulation strategy similar to what worked for ZINB, is more likely to succeed than structural rewrites.

**TODO-7: Poisson GLMM — SIMD vectorization at scale**
File: `EDI/src/fast_poisson_glmm.cpp`
Current profile shows 40% in a tight objective loop; the kernel is 0.003s at the current benchmark size. No change is warranted at this scale. If Poisson GLMM becomes a dominant bottleneck with many more groups or higher `n_gh`, the right approach is SIMD vectorization across observations within a group within a quadrature node — a much more invasive change that requires measuring first.

### Lower priority

**TODO-8: Ordinal CLMM — reduce log_sum_exp calls**
File: `EDI/src/_glmm_engine.h`
`glmm::log_sum_exp` appears at 2.4% in the retained-tree profile. This is secondary to the exp/log reduction in TODO-1, but once TODO-1 is resolved, profile whether `log_sum_exp` is reducible (e.g., by fusing it with the node-weight accumulation).

**TODO-9: ZINB/ZAP — allocator pressure**
File: `EDI/src/fast_zinb.cpp`
`malloc`/`cfree` combined is ~4.6% in the ZINB profile. The per-call `std::vector` allocations for `lgamma_yptheta[]` and `digamma_yptheta[]` tables (introduced in the hoisting optimization) are the likely cause. These could be preallocated as member vectors (non-mutable, resized at construction once distinct-y count is known) to eliminate per-call allocation. Low-risk since these are not inside the GH quadrature inner loop.

**TODO-10: Full re-benchmark after all retained changes**
Tooling: existing benchmark scripts
Run all GLMM-family and non-GLMM kernels on the fully retained tree in a single unified benchmark session to get updated baseline numbers for each kernel. Update the summary table above. This is bookkeeping, not optimization, but important for tracking future regression.

**TODO-11: Gaussian LMM — profile at relevant scale**
File: `EDI/src/fast_gaussian_lmm.cpp`
Currently 0.142s in the additional sweep — not trivial. No `perf` trace was gathered. If Gaussian LMM appears in hot paths for real workloads, run a profiling pass similar to the ordinal CLMM v2 re-profile.

**TODO-12: Clogit-plus-GLMM — add canonical-package testthat coverage**
File: `EDI/tests/testthat/`
There is no dedicated `testthat` file for `clogit_plus_glmm` canonical-package parity. Add one to protect the Phase 2/3/4 retained changes from regression.

---

## Bottom Line

The highest-payoff optimization work is structural:

- fewer temporary vectors and full-array passes
- less per-node/per-row allocation
- special-function hoisting and tabulation (not approximation, at least until digamma approximation is validated)
- specialized early-return paths in native fitters for callers that do not need Hessian outputs

The remaining headroom is concentrated in:
1. ordinal CLMM exp reduction (TODO-1) — 44% of retained-tree time is exp/log, but the correct implementation requires caller-stack-buffer design
2. ZINB digamma approximation (TODO-2) — ~33% of ZINB time is still in `Rf_dpsifn` + `Rf_chebyshev_eval` after hoisting
3. d-optimal search algorithmic pruning (TODO-3) — data-access cleanup gave only 15%; the real opportunity is the scan itself

Handwritten assembly is not the right next step for any of these.
