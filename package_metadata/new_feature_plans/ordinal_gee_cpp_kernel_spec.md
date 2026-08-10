# Ordinal GEE C++ Kernel Spec

Generated: 2026-07-28

## Scope

This spec defines an implementation plan for a C++ kernel that fits the
KK-design ordinal GEE model currently served only by `multgee::ordLORgee`,
so that `InferenceOrdinalKKGEE` gets a `use_rcpp` fast path like every other
KK-GEE inference class in this codebase.

Related documents:

- [../EDI/R/inference_mixin_kk_gee_shared.R](../EDI/R/inference_mixin_kk_gee_shared.R) —
  shared GEE dispatch/fallback machinery every KK-GEE daughter class uses.
- [../EDI/R/inference_ordinal_KK_combined.R](../EDI/R/inference_ordinal_KK_combined.R) —
  `InferenceOrdinalKKGEE`, the class this spec modifies.
- [../EDI/src/fast_gee.cpp](../EDI/src/fast_gee.cpp) — existing
  pairs/singletons GEE kernel for gaussian/binomial/poisson families; the
  scaffolding this spec's kernel mirrors.
- [../EDI/src/fast_ordinal_regression.cpp](../EDI/src/fast_ordinal_regression.cpp) and
  [../EDI/src/ordinal_fixed_link_helpers.h](../EDI/src/ordinal_fixed_link_helpers.h) —
  existing cumulative-link mean-model code this spec reuses.
- [path_audits_source.R](../package_tests/path_audits_source.R) line 102 —
  `InferenceOrdinalKKGEE` audit row, to be updated once `use_rcpp` exists.

## Motivation

`InferenceOrdinalKKGEE` hardcodes `use_rcpp = FALSE` when calling
`init_kk_gee_shared` (`inference_ordinal_KK_combined.R:37`) because no C++
kernel exists for ordinal GEE with a local-odds-ratio (LOR) association
structure. Every model fit — the point estimate, every asymptotic CI/p-value,
every randomization-inference replicate, and every bootstrap-weighted refit —
goes through `multgee::ordLORgee`, which is markedly slower than this
codebase's C++ kernels for the other four GEE families (gaussian, binomial,
poisson via `fast_gee.cpp`). This is the last KK-GEE response family without
a fast path.

## Non-Goals

- Do not implement multgee's general `LORstr` options (`"category.exch"`,
  `"RC"`, etc.). Only `LORstr = "uniform"` — a single scalar log-odds-ratio
  parameter shared across all cutpoint pairs — is ever passed by this
  codebase (`inference_ordinal_KK_combined.R:147`), and it is the only
  structure this spec implements.
- Do not support general cluster sizes. KK designs only ever produce
  singleton or paired clusters, the same constraint
  `gee_pairs_singletons_cpp_impl` already hard-asserts
  (`fast_gee.cpp:146`: `"only singletons and pairs are supported"`).
- Do not remove `multgee` as a dependency or drop the `use_rcpp = FALSE`
  fallback path. `multgee::ordLORgee` remains the parity-test oracle and the
  fallback, exactly like `geepack` for the other GEE families.
- Do not change `compute_estimate_with_bootstrap_weights`'s current
  `fast_ordinal_regression_weighted_cpp`-based approximation for anything
  other than `InferenceOrdinalKKGEE` — other ordinal classes that call that
  function for unrelated (non-GEE) reasons are out of scope.

## Statistical Background

The mean model is a proportional-odds cumulative-logit model shared with
every other ordinal class in this codebase: K categories, K−1 ordered
cutpoints `alpha_1 < ... < alpha_{K-1}`, linear predictor `eta = X * beta`,
and cumulative probability `gamma_k(x) = P(Y <= k | x) = expit(alpha_k - eta)`.
This is exactly what `edi_ordinal::FixedOrdinalRegression`
(`ordinal_fixed_link_helpers.h`) and `OrdinalRegression`
(`fast_ordinal_regression.cpp`) already compute for the independence case —
reused as-is here, with `Link::Logit`.

For a matched pair, define K−1 binary indicators per subject,
`Z^{(k)} = 1(Y <= k)`. The **local odds ratio** at cutpoint pair `(k, l)`
between the two cluster members is

```text
theta_{k,l} = [P(Z1<=k, Z2<=l) * P(Z1>k, Z2>l)] / [P(Z1<=k, Z2>l) * P(Z1>k, Z2<=l)]
```

Under `LORstr = "uniform"`, `log(theta_{k,l}) = alpha` (one scalar) for every
`(k, l)` in every pair. This is the association parameter multgee estimates
alongside `beta`; it plays the role `fast_gee.cpp`'s scalar exchangeable
correlation (`gee_estimate_exchangeable_alpha`, `fast_gee.cpp:92`) plays for
the other three families, but the mapping from `alpha` to a working
covariance is different because the response is multi-category.

## Algorithm

This implements Touloumis, Agresti & Kateri (2013, JASA), "GEE for
Multinomial Responses Using a Local Odds Ratios Parameterization" — the
algorithm `multgee::ordLORgee` implements — restricted to the case every
cluster has size 1 or 2. This restriction is the same simplification that
makes `fast_gee.cpp` compact relative to a general-cluster-size GEE library:
a pair's local-odds-ratio table only ever has one 2x2 cell to reconcile per
cutpoint pair, computed via a single closed-form Plackett-quadratic
evaluation, rather than the larger correlation-structure bookkeeping
multgee's general code needs for bigger clusters.

Per outer iteration:

1. **Beta update.** Given the current `alpha` (and hence the LOR-implied
   joint cell probabilities for every pair, via the Plackett quadratic
   applied to each of the `(K-1)^2` cutpoint cells), take one GEE
   Fisher-scoring step: `beta_new = beta_old + (D' V^-1 D)^-1 D' V^-1 (Z - gamma)`,
   where `D` is the Jacobian of the cumulative probabilities with respect to
   `(alpha_1, ..., alpha_{K-1}, beta)` (available from
   `edi_ordinal::FixedOrdinalRegression`'s existing score/Hessian machinery)
   and `V` is the pair's LOR-implied working covariance of `Z`. Singleton
   (reservoir) subjects contribute independently with no LOR term, mirroring
   how `gee_pairs_singletons_cpp_impl` handles `grp_size == 1`.
2. **Alpha update.** Given the current `beta` (hence marginal `gamma`'s),
   re-estimate the scalar `alpha` via a 1-D Newton/moment solve that pools
   the `(K-1)^2` cutpoint cells across all pairs — a single scalar root-find,
   not a matrix estimating equation, precisely because `LORstr = "uniform"`
   collapses the association to one number.
3. Iterate 1–2 until both `beta` and `alpha` change by less than `tol`
   (default `1e-8`, matching `fast_gee.cpp`'s convention) or `maxit` is
   reached.
4. **Variance.** Standard GEE sandwich estimator: bread
   `(D' V^-1 D)^-1` aggregated over clusters, meat from the empirical
   outer product of each cluster's residual score contribution
   `D_i' V_i^-1 (Z_i - gamma_i)`. This is what `stats::vcov(mod)` returns
   from `ordLORgee` today (`inference_ordinal_KK_combined.R:164`,
   `:325` post-fix) and what the kernel's `vcov` output must match.

## C++ Interface

New file `EDI/src/fast_ordinal_gee.cpp`, mirroring `fast_gee.cpp`'s exported
signatures:

```cpp
// [[Rcpp::export]]
List ordinal_gee_pairs_singletons_cpp(
    SEXP X_r, SEXP y_r, SEXP group_id_r,
    Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
    int maxit = 100, double tol = 1e-8
);

// [[Rcpp::export]]
List ordinal_gee_pairs_singletons_weighted_cpp(
    SEXP X_r, SEXP y_r, SEXP group_id_r, SEXP weights_r,
    Rcpp::Nullable<Rcpp::NumericVector> warm_start_beta = R_NilValue,
    Rcpp::Nullable<Rcpp::NumericMatrix> warm_start_fisher_info = R_NilValue,
    int maxit = 100, double tol = 1e-8
);
```

Both return `List(beta, alpha, vcov, fisher_information, converged, niter, K)`
— the same field names as `gee_pairs_singletons_cpp`'s result plus `K`
(category count), so `InferenceOrdinalKKGEE`'s extraction code
(`gee_treatment_index`, `extract_gee_treatment_estimate`, etc.) needs no
family-specific branching beyond what already exists. There is no
`family_str` argument (unlike `fast_gee.cpp`) since this kernel only ever
fits the ordinal cumulative-logit family.

Internals reuse `edi_ordinal::FixedOrdinalRegression` and the `Link::Logit`
primitives from `ordinal_fixed_link_helpers.h` for the mean-model
probabilities and Jacobian, rather than duplicating cumulative-logit math.
The pairwise LOR/Plackett-quadratic logic and the GEE1 outer loop are new
code, factored into a small internal header (`_ordinal_gee_helpers.h`,
following the existing `_helper_functions.h` naming convention) shared
between the two exported entry points.

## R Integration

1. Add `use_rcpp = TRUE` to `InferenceOrdinalKKGEE$initialize` and thread it
   into `init_kk_gee_shared(des_obj, use_rcpp = use_rcpp)`
   (`inference_ordinal_KK_combined.R:37`), matching the pattern every other
   KK-GEE daughter already uses.
2. `fit_ordinal_gee_mod()` (added in the 2026-07-28 randomization-inference
   fix) gains a `private$use_rcpp` branch: `TRUE` calls
   `ordinal_gee_pairs_singletons_cpp`, `FALSE` keeps calling
   `multgee::ordLORgee` exactly as today. `shared_gee_dispatch` and
   `compute_treatment_estimate_during_randomization_inference` both already
   route through `fit_ordinal_gee_mod()`, so no further changes are needed
   there.
3. `compute_estimate_with_bootstrap_weights` switches from
   `fast_ordinal_regression_weighted_cpp` (independence approximation) to
   `ordinal_gee_pairs_singletons_weighted_cpp` when `use_rcpp = TRUE`, giving
   Bayesian-bootstrap/BRT paths the real GEE-LOR estimator instead of an
   independence-only stand-in. When `use_rcpp = FALSE`, weighted refits fall
   back to a weighted `multgee::ordLORgee` call (multgee supports case
   weights natively).
4. `init_kk_gee_shared`'s existing `!use_rcpp` branch already requires
   `geepack` for other daughters; `InferenceOrdinalKKGEE` should instead
   assert `multgee` is installed when `use_rcpp = FALSE` (it already is,
   unconditionally, via the class's own `initialize`), and should NOT
   trigger the generic `geepack` check, since this class never calls
   `geepack`.

## Testing

Extend `EDI/tests/testthat/test-kk-gee-parity.R` with a
`compare_kk_gee_wrapper_paths`-style test for `InferenceOrdinalKKGEE`,
parameterized over:

- `K in {3, 4, 5}` categories
- pair-heavy, singleton-heavy, and mixed designs
- a sparse/near-boundary case (one category with very low prevalence) to
  stress the Plackett-quadratic solve

Tolerances should follow the existing convention (tight for well-conditioned
cases, loosened for the iterative solves the way `InferenceCountPoissonKKGEE`
already uses `2e-2`/`2e-3` rather than the `1e-4`/`1e-5` used for binomial).

Keep the `"ordinal KK GEE randomization p-value is estimable for >2-level
responses"` regression test (added 2026-07-28) passing under both
`use_rcpp = TRUE` and `use_rcpp = FALSE`.

## Benchmarking

Add `benchmark/fast_ordinal_gee_correctness.{R,cpp}` and
`benchmark/fast_ordinal_gee_profile.{R,cpp}`, following the existing
correctness/profile pairing convention in this repo (e.g.
`benchmark/fast_trigamma_correctness.{R,cpp}` and
`benchmark/fast_trigamma_profile.{R,cpp}`). The correctness script
cross-checks the kernel against `multgee::ordLORgee` over many random
datasets varying `K`, `n`, and pair/singleton mix; the profile script times
the C++ kernel against `multgee` across an `n`/`K` grid to quantify the
actual speedup this spec is chasing.

## Risks / Open Questions

- **Plackett-quadratic degeneracy.** Near-zero cutpoint cell counts (sparse
  categories, small pairs-only designs) can make the quadratic solve
  ill-conditioned. Needs the same kind of guarded fallback
  (`tryCatch`-equivalent in C++, returning `converged = false`) that
  `gee_pairs_singletons_cpp_impl` already uses for its own non-convergence
  cases.
- **Convergence-tolerance matching.** multgee's default `tol`/`maxiter` need
  to be matched closely enough that `use_rcpp = TRUE` and `FALSE` agree
  within the parity tolerances above; exact algorithmic equivalence is not
  guaranteed by matching the published estimating equations alone.
- **Complete/quasi-complete separation.** Small-n cumulative-logit fits can
  separate, same as the existing independence-model ordinal kernels; no new
  handling beyond what `edi_ordinal::FixedOrdinalRegression` already does is
  assumed, but this should be explicitly verified for the paired-LOR case
  during implementation.

## Implementation Phases

### Phase 1: Core kernel
- [ ] TODO-1: Implement `ordinal_gee_pairs_singletons_cpp` (unweighted) in
`fast_ordinal_gee.cpp`, reusing `edi_ordinal::FixedOrdinalRegression` for the
mean model. Validate against `multgee::ordLORgee` directly (no R wrapper
class yet), matching the style of
`EDI/tests/testthat/test-kk-gee-parity.R`'s existing direct-solver test
("fast KK GEE direct solver matches geepack for binomial and Poisson").

### Phase 2: Weighted variant
- [ ] TODO-2: Implement `ordinal_gee_pairs_singletons_weighted_cpp`, validated against a
weighted `multgee::ordLORgee` call.

### Phase 3: R integration
- [ ] TODO-3: Wire `use_rcpp` through `InferenceOrdinalKKGEE` as described above. Update
`compute_estimate`, CI/p-value paths, `compute_treatment_estimate_during_randomization_inference`,
and `compute_estimate_with_bootstrap_weights`.

### Phase 4: Test & benchmark suite
- [ ] TODO-4: Extend `test-kk-gee-parity.R` per the Testing section; add the
correctness/profile benchmark pair.

### Phase 5: Documentation & audit update
- [ ] TODO-5: Update `path_audits_source.R:102`'s `InferenceOrdinalKKGEE` notes to record
`use_rcpp` support and drop the now-obsolete "InferenceAsymp" framing if it
no longer reflects the class's fast/slow duality. Regenerate
`path_audits.html`.

## Acceptance Criteria

- `ordinal_gee_pairs_singletons_cpp` and `_weighted_cpp` match
  `multgee::ordLORgee` within the tolerances specified in Testing, across
  the full parameter grid described there.
- `InferenceOrdinalKKGEE$new(des, use_rcpp = TRUE)` produces
  `compute_estimate()`, `compute_asymp_confidence_interval()`,
  `compute_asymp_two_sided_pval()`, `compute_rand_two_sided_pval()`, and
  `compute_rand_confidence_interval()` results matching the
  `use_rcpp = FALSE` path within test tolerance.
- The 2026-07-28 randomization-inference regression test passes under both
  `use_rcpp` values.
- Profile benchmark demonstrates a measurable speedup over `multgee` at
  representative `n`/`K` combinations (no fixed target multiplier — report
  actual numbers).
- `multgee` remains a `Suggests` dependency; no existing `use_rcpp = FALSE`
  behavior changes.
