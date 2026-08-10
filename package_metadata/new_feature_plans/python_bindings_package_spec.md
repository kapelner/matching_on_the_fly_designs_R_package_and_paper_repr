# Python Bindings Package Spec

Generated: 2026-07-28 (rescoped 2026-07-28)

## Scope

This spec defines an implementation plan for a Python package,
`edi_kernels`, distributed under that same name on PyPI (no
hyphen/underscore split — kept simple deliberately, see Release
Checklist's "Naming" note), that exposes **only the model-fitting C++
API** —
the C++ functions that fit a specified parametric regression / GLM / GLMM /
survival model via MLE/IRLS/numerical optimization and return coefficients
(and, on the `_with_var`/hessian paths, score/Hessian/variance) — via
`pybind11`, plus a benchmark harness that compares each exposed kernel
against a canonical Python baseline, the same discipline
[benchmark_model_fits.md](benchmark_model_fits.md) already applies against R
canonical implementations (`glm.fit`, `MASS::glm.nb`, `betareg.fit`,
`coxph.fit`, etc.).

**In scope — the 33-file model-fitting core** (identified 2026-07-28 by
grepping `EDI/src` for the fit-a-model-via-optimization pattern, as opposed
to bootstrap/design/resampling/nonparametric-test helpers):

```text
fast_ols.cpp                          fast_ordinal_cauchit_regression.cpp
fast_robust_regression.cpp            fast_ordinal_cloglog_regression.cpp
fast_logistic_regression.cpp          fast_adjacent_category_logit.cpp
fast_probit_regression.cpp            fast_continuation_ratio_regression.cpp
fast_log_binomial_regression.cpp      fast_stereotype_logit.cpp
fast_poisson_regression.cpp           fast_coxph_regression.cpp
fast_negbin_regression.cpp            fast_weibull_regression.cpp
fast_zinb.cpp                         fast_weibull_frailty.cpp
fast_zero_augmented_poisson.cpp       fast_survival_models_optim.cpp
fast_hurdle_negbin.cpp                fast_gaussian_lmm.cpp
fast_beta_regression.cpp              fast_poisson_glmm.cpp
beta_regression_helpers.cpp           fast_logistic_glmm.cpp
fast_zero_one_inflated_beta.cpp       fast_hurdle_poisson_glmm.cpp
fast_ordinal_regression.cpp           fast_ordinal_glmm.cpp
fast_ordinal_probit_regression.cpp    fast_ordinal_clmm.cpp
                                       fast_clogit_plus_glmm.cpp
                                       fast_cpoisson_combined.cpp
                                       fast_gee.cpp
```

plus their transitive local headers: `_helper_functions.h`,
`fast_gamma_functions.h`, `ordinal_fixed_link_helpers.h`,
`optimization_starts.h`, `fast_erfc.h`, `_glmm_engine.h`, `_glmm_links.h`.

**Permanently out of scope for this package** (not deferred — these simply
never belong in a "model fitting" API and won't appear in a later phase of
this spec):

- Resampling/bootstrap kernels (`bootstrap_indices.cpp`,
  `rand_bootstrap_*_parallel.cpp`, `kk_bootstrap_*.cpp`,
  `sample_bootstrap_distr_weighted_distances.cpp`, `fast_bai_parallel.cpp`,
  `ols_distr_parallel.cpp`, `simple_mean_diff_parallel.cpp`,
  `ridit_distr_parallel.cpp`, `kk_compound_distr_parallel.cpp`, etc.)
- Sequential-design/randomization kernels (`atkinson_*.cpp`, `efron_redraw.cpp`,
  `kk14_redraw.cpp`, `pocock_simon_assign.cpp`, `generate_permutations.cpp`,
  `randomization_loop.cpp`, `rerandomization_helpers.cpp`,
  `random_block_size_speedups.cpp`, `design_fixed_greedy.cpp`,
  `optimal_design_search.cpp`, `optimal_blocks_distance.cpp`,
  `binary_match_search.cpp`, `match_data_compute_speedup.cpp`,
  `kk_lin_match_data.cpp`, `compute_mahal_distances.cpp`, `spbr_speedups.cpp`)
- Nonparametric test/CI kernels that compute a statistic or interval directly
  rather than fitting a parametric model (`fast_wilcox_hl.cpp`,
  `fast_wilcox_parallel.cpp`, `fast_kk_wilcox_parallel.cpp`,
  `fast_jonckheere_terpstra.cpp`, `fast_logrank.cpp`, `fast_ridit_analysis.cpp`,
  `cmh_speedups.cpp`, `zhang_exact_speedups.cpp`, `newcombe_speedups.cpp`,
  `miettinen_nurminen_speedups.cpp`, `fast_survival_stats.cpp`)
- KK21 stepwise-selection weight kernels (`kk21_weights.cpp`,
  `kk21_stepwise_survival.cpp`) — these run at randomization/design time to
  pick covariates for matching, not at inference/model-fit time.
- Data/design utility helpers with no fitting logic of their own
  (`build_kk_combined_*_design.cpp`, `bisection_ci*.cpp`, `lrt_ci_newton.cpp`,
  `gcomp_speedups.cpp`, `robust_post_fit_speedups.cpp`, `log_lik_nb.cpp`,
  `get_column_types.cpp`, `imputation_helpers.cpp`, `which_cols_vary.cpp`,
  `qr_reduce_design_matrix.cpp`, `fast_sample_int.cpp`, `fast_shuffle.cpp`,
  `fast_scale_cols.cpp`, `sample_mode.cpp`, `fast_matrix_rank.cpp`,
  `compute_weighted_distances.cpp`, `compute_all_subject_data.cpp`,
  `kk_cluster_ids.cpp`, `survival_strata_ids.cpp`, `pair_dist_helpers.cpp`,
  `result_key_store.cpp`, `omp_control.cpp`, `build_info.cpp`,
  `simulation_dgp.cpp`, `test_smart_starts.cpp`)

If a future need arises to bind any of these, write a separate spec for it —
don't fold resampling/design bindings into this one after the fact; the
narrower scope is deliberate (see Motivation for the R Dependency Audit
below — this scope has fundamentally different, much simpler R-coupling than
the full kernel set).

Prerequisite: this spec assumes the `*_core` functions from the SEXP-removal
spec exist and return `edi::ResultMap`, but only for the 33 files above — the
SEXP-removal spec's later phases (covering the excluded files) don't need to
finish first. Do not start Phase 1 below until at least one kernel family
has been migrated end-to-end there — the pilot (`fast_poisson_glmm_cpp`) is
the natural first binding target.

## R Dependency Audit (Model-Fitting Scope)

Re-running the R-dependency survey from
[sexp_removal_rcppeigen_conversion_spec.md](sexp_removal_rcppeigen_conversion_spec.md)
against exactly the 33 files + 7 headers above (2026-07-28) gives a much
smaller and cleaner picture than the full-repo audit:

- **RNG: zero.** None of these 33 files touch R's RNG (`unif_rand`,
  `GetRNGstate`, `R::runif`, etc.) at all. This means the RNG-stream design
  decision that gates Phase 4 of the SEXP-removal spec is **entirely moot
  for this package** — there is nothing in the model-fitting scope waiting
  on that decision.
- **`Rf_*` raw R-object introspection: zero.** No `Rf_isNull`, `Rf_getAttrib`,
  `Rf_nrows`, etc. anywhere in this scope (those live in
  `fast_bai_parallel.cpp`, `base_bootstrap_loop.cpp`, `get_column_types.cpp`,
  etc. — all excluded above).
- **Rmath special functions still called directly** (not yet routed through
  a `fast_*` replacement): 15 call sites total —
  - `R::pchisq` — 3 sites, all in `_helper_functions.h`'s shared
    `likelihood_ratio_test_from_negloglik`/`score_test_from_score_information`
    helpers (chi-squared p-value for the LRT/score-test fields that many
    `fast_*_cpp` functions attach to their fit result).
  - `R::qnorm5` — 3 sites: `fast_ordinal_probit_regression.cpp` (legacy
    warm-start), `fast_probit_regression.cpp` (warm-start), and
    `_helper_functions.h` (a shared quantile helper).
  - `R::pnorm`/`R::pnorm5`/`R::dnorm` — 6 sites: `fast_probit_regression.cpp`
    (2, log-CDF upper/lower tail) and `_glmm_links.h`'s `ProbitLink` (4:
    `cdf`/`pdf`/`pdf_from_cdf`/`deriv_pdf`) — reachable from any GLMM kernel
    templated on a probit link.
  - `R::dnbinom_mu` — 1 real site, `fast_hurdle_negbin.cpp:498` (the other
    two `R::dnbinom_mu` grep hits in this file and `fast_negbin_regression.cpp`
    are comments explaining why an explicit log-PMF formula is used instead).
  - `R::lbeta` — 1 site, `beta_regression_helpers.cpp`'s `beta_dev_resids_cpp`
    (deviance residuals; not touched by the recent `lgammafn`/`digamma`
    migration since `lbeta` is a different function).
  - `R::digamma`/`R::trigamma` — **fully migrated**, zero remaining external
    call sites; the only matches are the intentional `x <= 0` fallback
    inside `fast_digamma`/`fast_trigamma` themselves in
    `fast_gamma_functions.h`.
  All of these are candidates for the same `fast_*`-via-Boost.Math-reference
  treatment already used for `fast_digamma`/`fast_lgamma`/`fast_trigamma`/
  `fast_erfc` — `pnorm`/`dnorm`/`qnorm` already have a fast path
  (`fast_erfc.h`) not yet adopted at these specific call sites; `pchisq` and
  `dnbinom_mu`/`lbeta` do not have one yet.
- **`SEXP` zero-copy input pattern:** 31 of the 33 files use the same
  `SEXP ..._sexp` -> `Rcpp::NumericMatrix`/`NumericVector` -> `Eigen::Map`
  boilerplate documented in the SEXP-removal spec — same mechanical fix
  applies (declare the parameter as `Eigen::Map<const Eigen::MatrixXd>`
  directly).
- **`Rcpp::List::create`-built return values:** 31 of the 33 files (plus
  `_helper_functions.h`) — same `edi::ResultMap` treatment from the
  SEXP-removal spec applies unchanged.
- **`Rcpp::Nullable<T>` optional-argument pattern — pervasive and not
  covered by the SEXP-removal spec.** Every one of the 33 files uses
  `Rcpp::Nullable<T> arg = R_NilValue` for optional warm-start/fixed-parameter
  arguments (ranging from 4 to 45 occurrences per file). The SEXP-removal
  spec's `ResultMap`/`to_rcpp_list` design only solves the *output* side;
  this is an *input*-side R-specific idiom with no equivalent in the
  SEXP-removal spec and needs its own conversion convention here (see
  Optional Arguments below) — this is the one genuinely new binding-design
  requirement this narrower scope surfaces.
- **`R_ext/BLAS.h`/`<Rmath.h>` includes:** present in `_helper_functions.h`
  and `fast_gamma_functions.h` only. `R_ext/BLAS.h` is just standard BLAS
  prototypes under R's Fortran-mangling macros (trivial to satisfy with any
  BLAS provider outside R); `<Rmath.h>` is needed only for the `fast_digamma`/
  `fast_trigamma` `x <= 0` fallback and `R::pchisq`/`R::qnorm5` etc. above.
- **`checkUserInterrupt`:** zero hits — no R-interrupt-polling dependency.
- **Error handling (`Rcpp::stop()`/`stop()`) — discovered 2026-07-29, missed
  by the original audit.** 54 call sites across the 33 files, and critically
  **25 of the 54 sit inside internal (non-exported) functions**, not just
  the thin `[[Rcpp::export]]` wrappers where R-specific error signaling
  would be expected and unremarkable. `Rcpp::stop()` throws a C++ exception
  that Rcpp's export machinery catches at the R/C++ boundary and translates
  into an R condition (`simpleError`/`stop()` from R's own perspective) —
  it is fundamentally tied to that translation layer, not a generic
  `throw`. The 10 files with internal-function `stop()` calls:
  `fast_robust_regression.cpp` (2, `fast_robust_regression_internal`'s
  dimension checks), `fast_probit_regression.cpp` (1,
  `fast_probit_regression_internal`), `fast_log_binomial_regression.cpp`
  (5, the two `fit_constrained_binomial*_cpp_impl` internal fit routines),
  `fast_negbin_regression.cpp` (1, `fast_neg_bin_internal`), `fast_zinb.cpp`
  (1, `fast_zinb_internal`), `fast_hurdle_negbin.cpp` (9 — 8 in
  `validate_truncated_negbin_inputs`, a pure internal input-validation
  helper never touched by any R-facing code path directly, plus 1 in
  `fit_truncated_negbin_with_fallback`), `fast_adjacent_category_logit.cpp`
  (1, `fast_adjacent_category_logit_internal`),
  `fast_continuation_ratio_regression.cpp` (1,
  `fast_continuation_ratio_internal`), `fast_stereotype_logit.cpp` (1,
  `fast_stereotype_logit_internal`), `fast_gee.cpp` (1,
  `gee_pairs_singletons_cpp_impl`). The remaining 29 call sites are in
  exported wrappers, where `stop()` is expected and not a porting concern
  (a Python binding's own thin wrapper layer would raise its own exception
  type there instead). None of this was touched by the `Rcpp::Nullable`
  -> `std::optional` / `List::create` -> `edi::ResultMap` migration (Tasks
  9-15, 2026-07-29) — that work targeted the two *data-marshaling* patterns
  specifically (optional arguments, result construction), not error
  propagation, so this gap survived unnoticed until a direct dependency
  audit follow-up caught it.

**Status update (2026-07-30): resolved — no custom bridging layer needed.**
Verified directly against both frameworks' actual source (not from memory)
before implementing: Rcpp's auto-generated `[[Rcpp::export]]` wrapper in
`RcppExports.cpp` is wrapped in `BEGIN_RCPP`/`END_RCPP`
(`Rcpp/include/Rcpp/macros/macros.h`), which has an explicit
`catch(std::exception& __ex__)` clause — not just `catch(Rcpp::exception&)`
— that converts the exception via `exception_to_r_condition(__ex__)` and
raises it through R's own `stop()` mechanism, producing a normal,
`tryCatch()`-able R error condition. A plain `throw std::runtime_error(...)`
is handled identically to `Rcpp::stop(...)` at this boundary. Symmetrically,
pybind11's default exception translator (`pybind11/detail/internals.h`,
built into every bound function automatically, no glue code required)
catches the standard exception hierarchy and maps it to Python exceptions
using `e.what()` as the message: `std::domain_error`/`std::invalid_argument`
/`std::length_error`/`std::range_error` -> `ValueError`,
`std::out_of_range` -> `IndexError`, `std::overflow_error` ->
`OverflowError`, generic `std::exception` -> `RuntimeError`. So **no custom
`edi::fail()` or manual catch/re-raise bridging layer is needed at either
boundary** — both frameworks already do this automatically for anything
deriving from `std::exception`, and dispatching on the *specific* standard
subtype (e.g. `std::invalid_argument` for a dimension/length-mismatch
check) gets a more precise Python exception type (`ValueError` instead of a
blanket `RuntimeError`) for free.

**Implemented:** all 25 internal-function `Rcpp::stop()`/`stop()` call
sites converted to `throw std::invalid_argument(...)` (dimension/length/
finiteness precondition checks — the large majority) or `throw
std::runtime_error(...)` (genuine algorithmic-failure cases, not bad
input) across the 10 files listed above. The 29 exported-wrapper call
sites are unchanged (`Rcpp::stop()` remains correct and unremarkable
there — that boundary is necessarily R-specific and won't be reused
verbatim by a Python wrapper regardless). `<stdexcept>` added to each
touched file's includes; message text is preserved verbatim at every
converted call site, so anything depending on `tryCatch(..., error=...)`'s
message text at the R level is unaffected either way.

**Verified:** full EDI test suite (837 blocks, 4293 expectations) confirmed
passing after the swap, with the same 1 pre-existing, unrelated failure
(`test-resampling-draw-contracts.R`) seen throughout this session's work.
Direct `tryCatch()` testing at 4 converted call sites — `fast_robust_regression_cpp`
(bad `warm_start_weights` length), `fast_truncated_negbin_count_cpp` via
`validate_truncated_negbin_inputs` (bad `warm_start_params` length),
`fast_stereotype_logit_cpp` (bad `warm_start_params` size), and
`gee_pairs_singletons_cpp` (cluster size > 2) — confirmed R safely
intercepts the thrown `std::invalid_argument`/`std::runtime_error` as a
normal, catchable error condition with the expected message text; the R
session stays alive in every case, no crash.

**Net finding:** the model-fitting-only scope has no RNG dependency and no
raw R-object-introspection dependency at all — the two hardest categories
from the full-repo audit simply don't exist here. What remains is: a short,
nameable list of Rmath call sites (15, all portable via the same `fast_*`
pattern already established), the same mechanical SEXP-input/List-output
patterns the SEXP-removal spec already solves, one new pattern —
`Rcpp::Nullable` optional arguments — that this spec's binding layer needs
to handle explicitly, and (found later) 25 internal-function `Rcpp::stop()`
call sites that need a portable replacement (e.g. a small `edi::fail()`
throwing `std::runtime_error`, with the Rcpp-exported wrapper layer
catching and re-raising as an R condition, and a pybind11 wrapper layer
catching and re-raising as a Python exception) before those internal
functions can compile or run outside an Rcpp/R environment at all.

**Status update (post-2026-07-29 migration):** the `Rcpp::Nullable` pattern
above is now solved. All 33 files convert `Rcpp::Nullable<T>` to
`std::optional<T>` at the top of each `[[Rcpp::export]]` wrapper via a new
`nullable_to_optional<NativeT, RcppT>()` helper in `_helper_functions.h`,
and every internal (non-exported) function in this scope now takes
`std::optional<T>` natively — `Rcpp::Nullable` survives only in the
exported-function signatures themselves, which the Python binding layer
does not reuse (pybind11 has native `std::optional` support, so this
boundary is actually simpler to bind than the R boundary). Likewise,
`Rcpp::List::create`-built returns are now `edi::ResultMap` /
`edi::to_rcpp_list()` wherever the output is flat (no nested sub-lists,
named vectors, or `NA_LOGICAL`); the files that keep `List::create` do so
for a documented Rcpp-specific reason (see `Result Conversion` below) that
does not change the C++-level function signature this spec binds against.
Separately (TODO-130), the 15-call-site Rmath list above has shrunk to 3:
`fast_qnorm`/`fast_lbeta`/`fast_dnbinom_mu` (new) and `fast_log_pnorm`
(already existed) replaced every `R::qnorm5`/`R::pnorm5`/`R::lbeta`/
`R::dnbinom_mu` call site (15 -> 7), and a follow-up swap wired the
already-existing `pnorm_fast`/`dnorm_fast` into `_glmm_links.h`'s
`ProbitLink` in place of `R::pnorm`/`R::dnorm` (7 -> 3) — see `Standalone
Rmath Library Dependency` below for the current (much smaller) residual,
now just `R::pchisq`.

## Standalone Rmath Library Dependency

**Status update (2026-07-29, TODO-130):** `fast_qnorm`/`fast_lbeta`/
`fast_dnbinom_mu` (new) and `fast_log_pnorm` (already existed) have since
replaced every `R::qnorm5`/`R::pnorm5`/`R::lbeta`/`R::dnbinom_mu` call site
in `fast_probit_regression.cpp`, `fast_ordinal_probit_regression.cpp`,
`fast_hurdle_negbin.cpp`, and `beta_regression_helpers.cpp` — see
`package_metadata/perf_experiments_final.md` TODO-130.

**Status update (2026-07-29, follow-up):** `_glmm_links.h`'s `ProbitLink`
(`cdf`/`pdf`/`pdf_from_cdf`/`deriv_pdf`) now calls `pnorm_fast`/`dnorm_fast`
(already existed in `fast_erfc.h`, already validated to machine precision
against Boost.Math and `R::pnorm`/`R::dnorm` when `fast_log_pnorm`/
`fast_qnorm` were built) instead of `R::pnorm`/`R::dnorm` directly — a pure
call-site swap, no new implementation. Reachable from any GLMM kernel
templated on a probit link (currently `fast_ordinal_clmm.cpp`'s
`link = "probit"` path). Verified: `EDI:::fast_ordinal_clmm_cpp(..., link =
"probit")` fits and converges correctly on synthetic data; the
`glmm-cpp-equivalence`/`fast-probit-cdf`/`clogit-plus-glmm-cpp-equivalence`
test files (27 expectations) and the full package test suite (837 blocks)
pass unchanged.

The remaining footprint is now just **`R::pchisq`** (3 sites,
`_helper_functions.h`'s shared LRT/score-test p-value helpers) — down from
15 at the start of this audit, then 7 after TODO-130, now 3 — plus the
`fast_digamma`/`fast_trigamma` `x <= 0` fallback noted below. `R::pchisq`
has no `fast_*` equivalent yet; it would be genuinely new work (chi-squared
CDF via the regularized incomplete gamma function, Boost.Math as the
correctness reference) if pursued, unlike `pnorm`/`dnorm`/`qnorm`/`lbeta`/
`dnbinom_mu`, which were either already vendored or cheap compositions of
existing `fast_lgamma`. The analysis and recommendation otherwise stand
unchanged; the remainder of this section describes the (now much smaller)
residual dependency.

The remaining real (non-comment) Rmath call site in this scope —
`R::pchisq` (3, `_helper_functions.h`) — does **not** require linking
against R itself. R ships a standalone build of its math library
specifically for this use case:

- Define `MATHLIB_STANDALONE` before `#include <Rmath.h>`. Under that macro,
  `R::pchisq(...)` (Rcpp's thin `R::`-namespaced wrapper around the
  identical C function) resolves to a plain C symbol (`pchisq`) with zero
  dependency on `Rinternals.h`, `SEXP`, or any R runtime state — no embedded
  R, no `Rcpp::` types.
- **Obtaining `libRmath`:** build from `<R_HOME>/src/nmath/standalone/`
  (ships with every R source tree, `./configure && make` produces
  `libRmath.a`/`.so`), or install the prebuilt package directly
  (`r-mathlib`/`r-mathlib-dev` on Debian/Ubuntu). Both are the same code R
  itself uses, so results are bit-identical to the current `R::foo` calls —
  this is a build-target change, not a numerical-behavior change.
- **License:** Rmath is part of R and is GPL (≥ 2) licensed (with a few
  files LGPL). This is a real constraint if the Python wheel wants a
  permissive license (MIT/BSD/Apache) for the compiled extension — linking
  `libRmath` would make the wheel's effective license GPL for that binary.
- **Recommended default:** vendor/statically link standalone `libRmath` for
  the initial Python package (correct, zero re-implementation risk, matches
  R's behavior exactly). **Permissive-license escape hatch:** if GPL is a
  blocker, this is now down to a single function: implement `fast_pchisq`
  (chi-squared CDF via the regularized incomplete gamma function, Boost.Math
  as the correctness reference) as a `fast_*` leaf function, following the
  exact precedent already set by `fast_digamma`/`fast_lgamma`/
  `fast_trigamma`/`fast_erfc`/`fast_qnorm`/`fast_lbeta`/`fast_dnbinom_mu` in
  `fast_gamma_functions.h`/`fast_erfc.h`.
- **`fast_digamma`/`fast_trigamma` note:** these already have permissively-
  licensed local implementations for the normal path; the only remaining
  `R::digamma`/`R::trigamma` references are the `x <= 0` edge-case fallback
  inside those functions themselves (`fast_gamma_functions.h`), so they are
  a natural next candidate if the permissive-license path is chosen, but
  are not part of the 3-call-site count above since they're not called from
  application code in the normal case.
- **Build system integration:** if the standalone-`libRmath` default is
  chosen, add a `find_package`/vendored-source step to the `CMakeLists.txt`
  described below (`Build System`), analogous to the existing `Eigen3`/
  `pybind11` discovery steps — this is a small, self-contained addition, not
  a structural change to the build.

## Math-Utility Function Exports (`fast_math` submodule, proposed — 2026-07-29)

**This is a scope addition, not part of the "only the model-fitting API"
core scope above** — flagged separately rather than folded silently into
the 33-kernel API, per this spec's own `Non-Goals` discipline of writing
scope changes down explicitly rather than opportunistically. It exists
because the `R::`-replacement work above (TODO-130 and its follow-ups,
`package_metadata/perf_experiments_final.md`) produced a complete set of
dependency-free `fast_*` scalar math functions as a side effect, and a
direct question came up: are any of these fast enough on their own to be
worth exposing as Python utilities, independent of the model-fitting
kernels that use them internally?

**Full inventory** of `fast_*` math-utility functions in the codebase
(distinct from the 33 model-fitting kernels, which are full statistical
fits, not math primitives):

| Function | Source | Reference (R / scipy / numpy) |
|---|---|---|
| `fast_digamma` | `fast_gamma_functions.h` | `R::digamma` / `scipy.special.digamma` |
| `fast_trigamma` (+`fast_trigamma_vec`) | `fast_gamma_functions.h` | `R::trigamma` / `scipy.special.polygamma(1,·)` |
| `fast_lgamma` (+`_stirling`/`_lanczos`) | `fast_gamma_functions.h` | `R::lgamma` / `scipy.special.gammaln` |
| `fast_lbeta` | `fast_gamma_functions.h` | `R::lbeta` / `scipy.special.betaln` |
| `fast_dnbinom_mu` | `fast_gamma_functions.h` | `R::dnbinom_mu` / `scipy.stats.nbinom.logpmf` |
| `fast_gammap_series`/`fast_gammaq_cf`/`fast_gammaq` | `fast_gamma_functions.h` | internal building blocks only, not a call-site target |
| `fast_pchisq_upper` | `fast_gamma_functions.h` | `R::pchisq` / `scipy.stats.chi2.sf` |
| `fast_erfc` (+`_polevl`/`_p1evl`) | `fast_erfc.h` | `std::erfc` / `scipy.special.erfc` |
| `pnorm_fast` | `fast_erfc.h` | `R::pnorm` / `scipy.stats.norm.cdf` |
| `dnorm_fast` | `fast_erfc.h` | `R::dnorm` / `scipy.stats.norm.pdf` |
| `fast_log_pnorm` | `fast_erfc.h` | `R::pnorm(log=TRUE)` / `scipy.stats.norm.logcdf` |
| `fast_log_dnorm` | `fast_erfc.h` | `R::dnorm(log=TRUE)` / `scipy.stats.norm.logpdf` |
| `fast_qnorm` | `fast_erfc.h` | `R::qnorm5` / `scipy.stats.norm.ppf` |
| `fast_atan` | `ordinal_fixed_link_helpers.h` | `std::atan` / `numpy.arctan` |
| `fast_log1pexp` | `_helper_functions.h` | softplus / `numpy.logaddexp(0,·)` |
| `fast_log1p_arr` (scalar port: `fast_log1p_scalar`) | `_helper_functions.h` | `numpy.log1p` |

**Benchmark methodology:** a standalone pybind11 module was built directly
from verbatim copies of these functions (no R/Rcpp headers needed — all are
already dependency-free except `fast_digamma`/`fast_trigamma`'s `x <= 0`
fallback, which calls `R::digamma`/`R::trigamma` and is unreachable for the
`x > 0` domain these benchmarks use), with vectorized bindings (`py::array_t
-> py::array_t`, one Python call per batch) so the comparison is apples-to-
apples against scipy/numpy's own natively-vectorized ufunc calls — a
scalar-only binding would unfairly penalize the C++ side with per-element
Python call overhead unrelated to the underlying algorithm. Correctness
checked at `N=5000` against the listed scipy/numpy reference; speed checked
at `N` from 10 to 1,000,000 (paired timing, 15 reps per sample after 3
warmup calls, median reported).

**Correctness** (max error vs the scipy/numpy reference, `N=5000`):

| Function | max abs err | max rel err |
|---|---|---|
| `fast_digamma` | 2.89e-13 | 2.15e-09 |
| `fast_trigamma` | 4.30e-13 | 3.13e-12 |
| `fast_lgamma` | 9.10e-13 | 5.59e-12 |
| `fast_lbeta` | 4.55e-13 | 1.79e-13 |
| `fast_dnbinom_mu` | 4.62e-14 | 1.03e-14 |
| `fast_pchisq_upper` | 4.89e-15 | 3.39e-14 |
| `fast_erfc` | 4.44e-16 | 1.40e-15 |
| `pnorm_fast` | 1.11e-16 | 1.37e-15 |
| `dnorm_fast` | 5.55e-17 | 2.22e-16 |
| `fast_log_pnorm` | 3.55e-15 | 4.99e-08 |
| `fast_log_dnorm` | 3.55e-15 | 2.34e-16 |
| `fast_qnorm` | 2.92e-09 | 1.13e-09 |
| `fast_atan` | 2.22e-16 | 2.21e-16 |
| `fast_log1pexp` | 1.01e-11 | 1.45e-11 |
| `fast_log1p` | 8.88e-16 | 2.79e-16 |

`fast_log_pnorm`'s relative error (4.99e-08, an outlier against the rest of
the table) is not a real correctness problem — its absolute error is 3.55e-15,
identical in magnitude to every other row. `log(Phi(x)) -> 0` in the upper
tail, so a relative-error denominator near zero inflates the ratio despite
the actual discrepancy being machine-precision noise (the exact same
artifact documented for this function in TODO-130).

**Speed** (`scipy_or_numpy_time / fast_time`, paired, by array size `N`):

| Function | N=10 | N=100 | N=1,000 | N=10,000 | N=100,000 | N=1,000,000 |
|---|---|---|---|---|---|---|
| `fast_digamma` vs `special.digamma` | 0.86x | 1.00x | 1.24x | 1.18x | 1.15x | 1.16x |
| `fast_trigamma` vs `polygamma(1,·)` | 16.71x | 27.00x | 44.36x | 30.80x | 23.10x | 23.78x |
| `fast_lgamma` vs `special.gammaln` | 0.75x | 0.83x | 0.92x | 1.03x | 1.02x | 1.03x |
| `fast_lbeta` vs `special.betaln` | 0.92x | 1.82x | 2.24x | 2.26x | 2.28x | 2.26x |
| `fast_dnbinom_mu` vs `nbinom.logpmf` | 22.50x | 5.58x | 2.02x | 1.48x | 1.44x | 1.88x |
| `fast_pchisq_upper` vs `chi2.sf` | 22.05x | 6.28x | 3.82x | 2.87x | 2.89x | 3.20x |
| `fast_erfc` vs `special.erfc` | 0.86x | 1.29x | 1.49x | 1.42x | 1.40x | 1.40x |
| `pnorm_fast` vs `norm.cdf` | 58.71x | 20.05x | 6.70x | 2.59x | 2.36x | 3.98x |
| `dnorm_fast` vs `norm.pdf` | 83.71x | 46.50x | 13.98x | 2.26x | 3.51x | 8.63x |
| `fast_log_pnorm` vs `norm.logcdf` | 49.33x | 19.58x | 4.60x | 2.41x | 2.04x | 3.49x |
| `fast_log_dnorm` vs `norm.logpdf` | 90.71x | 85.43x | 69.55x | 42.81x | 38.33x | 52.61x |
| `fast_qnorm` vs `norm.ppf` | 116.83x | 72.70x | 21.80x | 9.88x | 7.15x | 17.10x |
| `fast_atan` vs `np.arctan` | 0.86x | 1.10x | 1.83x | 2.74x | 2.98x | 2.36x |
| `fast_log1pexp` vs `np.logaddexp(0,·)` | 1.43x | 1.22x | 1.51x | 1.94x | 2.10x | 1.71x |
| `fast_log1p` vs `np.log1p` | 0.86x | 1.08x | 1.61x | 1.80x | 1.83x | 1.83x |

**Export decision — per function, not blanket:**

**Export** (`edi_kernels.fast_math` submodule) — a real, consistent speed
advantage at every array size actually measured, or (for the three
borderline cases) at every realistic batch size, with only a wash at the
trivially-small `N=10` case that isn't representative of real statistical
batch usage:

- `fast_trigamma`, `fast_pchisq_upper`, `fast_dnbinom_mu` — the three that
  were the actual point of TODO-130 (replacing real `R::` call sites), and
  the strongest, most consistent wins (1.4x-44x across every `N`).
- `pnorm_fast`, `dnorm_fast`, `fast_log_pnorm`, `fast_log_dnorm`,
  `fast_qnorm` — the normal-distribution family, uniformly strong (2x-120x).
- `fast_lbeta` (0.92x-2.28x — sub-1x only at `N=10`, a clean, growing win
  from `N=100` on).
- `fast_atan`, `fast_log1pexp`, `fast_log1p` — modest but real and
  consistent for `N >= 100` (1.1x-3x); `fast_log1pexp` is the strongest of
  the three, positive even at `N=10`.

**Do not export** — no measurable benefit over scipy on this evidence,
sometimes measurably slower:

- `fast_digamma` (0.86x-1.24x — noise-level parity with `scipy.special.digamma`).
- `fast_lgamma` (0.75x-1.03x — same; `scipy.special.gammaln`'s C
  implementation is already as fast or faster).
- `fast_erfc` (0.86x-1.49x — modest at best, not worth the API-surface
  commitment discussed above; `scipy.special.erfc` is competitive).

Internal building blocks (`fast_gammap_series`, `fast_gammaq_cf`,
`fast_gammaq`, `fast_lgamma_stirling`, `fast_lgamma_lanczos`,
`fast_erfc_polevl`, `fast_erfc_p1evl`, `fast_trigamma_vec`) are never
exported regardless of the above — they exist only to implement the
functions in the export list, have no standalone call site, and (per `API
Naming` below) have no R-exported name to mirror in the first place.

If the `fast_math` submodule is built, apply the same design already used
throughout this spec: vectorized `py::array_t<double> -> py::array_t
<double>` bindings (not scalar-only — see the benchmark methodology note
above for why), one `tests/test_fast_math.py` per function following the
same cross-language-parity discipline as `Testing` below, and keep it as a
clearly separate submodule (`edi_kernels.fast_math.*`, not flat
`edi_kernels.*`) so it reads as an intentional, documented addition rather
than scope creep into the model-fitting kernel API.

## Optional Arguments

Bind `Rcpp::Nullable<T> arg = R_NilValue` parameters as Python keyword
arguments defaulting to `None`, using pybind11's native optional support
rather than inventing a sentinel:

```cpp
m.def("fast_negbin_regression", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                                    const Eigen::Ref<const Eigen::VectorXd>& y,
                                    std::optional<Eigen::VectorXd> warm_start_params,
                                    /* ... */) {
	return edi::to_py_dict(fast_neg_bin_core(X, y, warm_start_params, /* ... */));
}, py::arg("X"), py::arg("y"), py::arg("warm_start_params") = py::none(), /* ... */);
```

The `*_core` function signature itself should take `std::optional<T>` (or a
`const T*` with `nullptr` for "absent") instead of `Rcpp::Nullable<T>`, so
the core stays Rcpp-free — this is a mechanical rename at every one of the
~33 files' argument lists, not a design change, and should happen as part of
the same SEXP-removal-spec migration pass for these files, not as a separate
pybind11-side shim.

## Non-Goals

- Do not reimplement any numeric algorithm for Python. The Python package
  only binds the existing C++ cores; it does not re-derive scores/Hessians/
  optimizers in Python.
- Do not build a general "R-to-Python" formula/model-frame layer (`patsy`
  integration, R-style formulas). Bind the same low-level, pre-built-matrix
  interface that `EDI::`'s exported C++ functions already take.
- Do not chase 100% kernel parity in the first pass. Bind and benchmark the
  response-type families with clean Python canonical baselines first (see
  Baseline Gaps); families with no clean baseline get correctness-only
  treatment, tracked separately.
- Do not bind any file from the "permanently out of scope" list above under
  this spec, even opportunistically. A resampling/design/nonparametric-test
  Python API is a different package with a different R-dependency profile
  (RNG-stream reproducibility becomes a real, load-bearing question there,
  per the SEXP-removal spec's Non-Goals) — write it as its own spec.

## Package Layout

```text
python/
  pyproject.toml              # scikit-build-core + pybind11 build backend
  CMakeLists.txt               # top-level; finds Eigen3, pybind11
  src/
    edi_kernels/
      __init__.py
      _core.pyi                # type stubs for the compiled extension
  cpp/
    bindings_continuous.cpp     # pybind11 module glue, one file per family
    bindings_count.cpp
    bindings_ordinal.cpp
    bindings_survival.cpp
    bindings_incidence.cpp
    bindings_proportion.cpp
    bindings_glmm.cpp
    result_map_pybind.h          # py::dict converter (see below)
  tests/
    test_fast_poisson_glmm.py
    ...
  benchmarks/
    baselines.py                 # canonical Python baseline registry
    run_benchmark_audit.py       # mirrors benchmark/benchmark_model_fits.md
    report_template.md
```

`cpp/*.cpp` files `#include` the `*_core.h` headers directly from
`../../EDI/src/` (via a CMake include path) — they are not copies. This is
the entire point of the SEXP-removal prerequisite: one implementation, two
thin boundary layers.

## Build System

- Use `scikit-build-core` as the PEP 517 build backend (CMake-driven,
  standard for compiled Python extensions; avoids hand-rolling
  `setup.py build_ext`).
- `CMakeLists.txt` requirements:
  - C++20 (match `EDI/src/Makevars`' `CXX20STD`).
  - Find `Eigen3` (reuse the same version constraint the R package's
    `RcppEigen` vendors, or a system Eigen ≥ the same minimum — pin
    explicitly rather than floating).
  - Find `pybind11` (via `pybind11.get_cmake_dir()` at configure time, the
    standard scikit-build-core pattern).
  - Mirror `EDI/src/Makevars`' optimization flags for parity with the R
    build: `-O3 -march=native` by default, with a `EDI_PY_PORTABLE` CMake
    option analogous to `EDI_PORTABLE` that drops `-march=native` for
    portable wheels.
  - OpenMP: link the same way `Makevars` does
    (`SHLIB_OPENMP_CXXFLAGS`/`find_package(OpenMP)`) — needed for the GLMM
    files in this scope that use `#pragma omp`.
  - No MKL linkage question exists in this scope — check first, but as of
    the 2026-07-28 audit none of the 33 model-fitting files directly include
    `mkl.h` (that dependency lives in files already excluded above).
  - Standalone `libRmath` (see `Standalone Rmath Library Dependency` above)
    if the GPL-linked default is chosen: a `find_package`/vendored-source
    step alongside the `Eigen3`/`pybind11` discovery above.

## Result Conversion

Add `python/cpp/result_map_pybind.h`, mirroring
`EDI/src/result_map_rcpp.h`'s `std::visit` shape exactly:

```cpp
#ifndef EDI_RESULT_MAP_PYBIND_H
#define EDI_RESULT_MAP_PYBIND_H

#include "result_map.h"
#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace edi {

inline py::dict to_py_dict(const ResultMap& m) {
	py::dict out;
	for (const auto& [name, value] : m.entries()) {
		out[py::str(name)] = std::visit([](auto&& v) -> py::object {
			using T = std::decay_t<decltype(v)>;
			if constexpr (std::is_same_v<T, std::monostate>) {
				return py::none();
			} else {
				return py::cast(v);
			}
		}, value);
	}
	return out;
}

} // namespace edi

#endif
```

Each binding function is then a one-liner, same shape as the Rcpp wrapper:

```cpp
m.def("fast_poisson_glmm", [](const Eigen::Ref<const Eigen::MatrixXd>& X,
                               const Eigen::Ref<const Eigen::VectorXd>& y,
                               /* ... */) {
	return edi::to_py_dict(fast_poisson_glmm_core(X, y, /* ... */));
}, py::arg("X"), py::arg("y"), /* ... */);
```

Use `pybind11/eigen.h`'s automatic NumPy<->Eigen conversion
(`Eigen::Ref<const Eigen::MatrixXd>` accepts any C- or F-contiguous
`float64` NumPy array without a manual copy step in the binding code; a copy
happens only if the caller's array layout doesn't match, same tradeoff as
`Eigen::Map` on the Rcpp side).

## API Naming

Expose Python functions with the same base name as the R-exported function,
minus the `_cpp` suffix (which exists only because it's an R/Rcpp
convention): `fast_poisson_glmm_cpp` -> `edi_kernels.fast_poisson_glmm`. Do
not invent a different naming scheme — anyone cross-referencing behavior
between the R package and the Python package should be able to guess the
Python name from the R name and vice versa.

## Testing

Each bound kernel gets a `tests/test_<kernel>.py` using `pytest`, following
[[superpowers:test-driven-development]] — write the test first, watch it
fail (`ModuleNotFoundError`/`AttributeError` before the binding exists), then
add the minimal binding.

Required test per kernel:

1. **Cross-language parity**: generate a fixed dataset with a seed in R,
   write `X`, `y`, and any other inputs to a `.npz`/`.csv` fixture file
   checked into `tests/fixtures/`, call the R kernel once to record its
   output values into the same fixture, then assert the Python binding
   reproduces those values to a tight tolerance (`np.allclose(..., atol=1e-9,
   rtol=1e-9)` — not bit-identical, since compiler/BLAS reordering across
   the R and Python builds can differ in the last ULPs, but should agree far
   tighter than any statistical tolerance).
2. **Shape/dtype contract**: wrong-shaped input raises rather than
   segfaulting (pybind11/Eigen's `Ref` already throws on shape mismatch —
   assert that behavior explicitly rather than assuming it).
3. **Optional-argument contract**: omitting a `Nullable`-derived argument in
   Python (leaving it at its `None` default) reproduces the same fit as
   calling the R function without that argument — cheap to get wrong if the
   `std::optional`/`nullptr` translation on the core side is inverted.
4. **Baseline correctness** (for kernels with a Python baseline — see below):
   assert the EDI kernel's point estimate agrees with the canonical Python
   baseline's estimate to a documented tolerance on a benign dataset. This
   is a correctness check, not a benchmark; keep it in `tests/`, not
   `benchmarks/`.

## Benchmark Harness

### Methodology (mirrors `benchmark_model_fits.md`)

- **Bare-metal EDI timing:** call the bound kernel directly with pre-built
  NumPy arrays; no pandas/DataFrame overhead in the timed region.
- **Apples-to-apples canonical timing:** call the lowest-level fit
  interface the baseline package exposes (see table below); fall back to
  the standard high-level API only if no low-level entry point exists, and
  say so in the report the same way `benchmark_model_fits.md` already
  flags formula-API fallbacks.
- **Averaging:** medians over 30 cold timing samples via `time.perf_counter`
  with an explicit warm-up-then-discard first call, matching the R
  methodology's "30 cold estimate-only timing samples" (R's adaptive
  `system.time`/`microbenchmark` split at 0.01ms has a direct analog:
  use `timeit.repeat` with `number` scaled up for sub-millisecond kernels).
- **Significance:** Welch's two-sample t-test between the EDI and baseline
  timing replicate distributions (`scipy.stats.ttest_ind(..., equal_var=False)`),
  same test the R report already uses.
- **Row highlighting / report shape:** reuse `benchmark_model_fits.md`'s
  table columns exactly (`Class/Kernel`, `Response Type`, `EDI ms`,
  `Baseline Package`, `Baseline Function`, `Baseline ms`, `Speedup`,
  `Timing Pval`, significance stars) so the Python report reads as a
  companion table to the existing R one, not a differently-shaped document.
- **Dataset parity:** generate benchmark datasets with the exact same
  specification as "Benchmark Dataset Specification" in
  `benchmark_model_fits.md` (N=1000, 5 predictors, same effect-size draws),
  reusing the R-side fixture files from the parity tests above rather than
  redrawing independently in Python/NumPy — this keeps the Python report's
  numbers comparable to the R report's numbers, not just internally
  consistent.

### Python Baseline Registry (`benchmarks/baselines.py`)

Canonical Python baseline per response-type family, restricted to the
model-fitting families actually in scope now (rank tests, exact/CMH tests,
and proportion-CI methods are dropped from this table — their C++ files are
permanently out of scope, see Scope above):

**Status update (2026-08-03, TODO-4):** verified against the actual installed
versions in the `edi-py-bench` venv — `statsmodels==0.14.6`,
`scikit-survival==0.28.0`, `lifelines==0.30.0`. Every class named below
exists and imports cleanly at these versions; no floor-version pins are
needed beyond what's already installed. Findings that change the table
below are called out inline; two are consequential enough to repeat here:
(1) `OrderedModel` DOES support cauchit/cloglog — not via the `distr=`
string shortcut (only `'probit'`/`'logit'` strings work), but by passing a
`scipy.stats` distribution instance directly (`stats.cauchy` for cauchit,
`stats.gumbel_l` for cloglog — verified mathematically: cloglog's inverse
link is exactly the Gumbel-min CDF `1-exp(-exp(x))`, cauchit's is exactly
the standard Cauchy CDF). Both fit and converge on synthetic data. This
demotes those two rows from Baseline Gap to a real (if slightly unusual)
wiring. (2) `fast_hurdle_negbin_cpp`'s hurdle mechanism is an independent
logistic regression on the zero/positive indicator (`fast_logistic_regression_internal`
in `fast_hurdle_negbin.cpp`), whereas `HurdleCountModel`'s `zerodist=` models
the zero/positive split via a *censored* Poisson/NegBin likelihood, not a
separate logistic fit — a different hurdle parameterization, not a
literal match. Timing comparisons are still fair (same asymptotic
complexity class), but coefficient/estimate parity checks against this
baseline would need to account for the different zero-mechanism, or should
be validated against the R implementation instead (already covered by
`EDI/tests/testthat/test-hurdle-negbin-*.R`).

| Response family | EDI kernel example | R canonical (existing) | Python canonical baseline | Notes |
|---|---|---|---|---|
| Continuous / OLS | `fast_ols_cpp` | `lm.fit` | `numpy.linalg.lstsq` (lowest level) or `statsmodels.regression.linear_model.OLS(...).fit()` | Use `lstsq` for the bare-metal row; add an `OLS` row too since `benchmark_model_fits.md` includes both levels for some families. |
| Continuous / robust regression | `fast_robust_regression_cpp` | `Rfit`/robust equivalents | `statsmodels.robust.robust_linear_model.RLM` | **Verified mismatch, not just a defaults check:** EDI's default `method="MM"` is a two-stage S-then-M estimator; `RLM`'s default `M=None` resolves to a single-stage `HuberT` M-estimator (confirmed via `RLM(y, X).M` → `HuberT` instance) — these are different estimator classes, not just different tuning constants. For a fair timing/estimate comparison, either call EDI with `method="M"` (it supports plain M-estimation too, `c=1.345` matches `HuberT`'s default tuning constant exactly) or note the MM-vs-M mismatch explicitly in the report rather than treating a coefficient difference as a bug. |
| Incidence / logistic | `fast_logistic_regression_cpp` | `glm.fit` | `statsmodels.genmod.generalized_linear_model.GLM(family=sm.families.Binomial())` | statsmodels' IRLS fit is the closest low-level analog to `glm.fit`. |
| Incidence / probit, log-binomial | `fast_probit_regression_cpp`, `fast_log_binomial_regression_cpp` | `glm.fit(family=binomial(link=...))` | `statsmodels GLM(family=Binomial(link=probit()/log()))` | |
| Count / Poisson | `fast_poisson_regression_cpp` | `glm.fit` | `statsmodels GLM(family=sm.families.Poisson())` | |
| Count / NegBin | `fast_negbin_regression_cpp` | `MASS::glm.nb` | `statsmodels.discrete.discrete_model.NegativeBinomial` | `NegativeBinomial` jointly estimates dispersion via MLE, matching `glm.nb`'s semantics more closely than `GLM(family=NegativeBinomial(alpha=fixed))`. |
| Count / ZINB, ZAP, hurdle negbin | `fast_zinb_cpp`, `fast_zero_augmented_poisson_cpp`, `fast_hurdle_negbin_cpp` | (custom) | `statsmodels.discrete.count_model.{ZeroInflatedPoisson,ZeroInflatedNegativeBinomialP}`, `statsmodels.discrete.truncated_model.HurdleCountModel` | **Verified present at 0.14.6**, all three import and construct from raw `endog`/`exog` arrays (no formula/DataFrame requirement). `HurdleCountModel(dist="negbin", zerodist=...)` — `zerodist` only accepts `"poisson"`/`"negbin"` (raises `NotImplementedError` otherwise), and its hurdle mechanism differs structurally from EDI's — see status note above. |
| Proportion / beta regression | `fast_beta_regression_cpp` | `betareg::betareg.fit` | `statsmodels.othermod.betareg.BetaModel` | **Verified present at 0.14.6**, takes raw `endog`/`exog` (+ optional `exog_precision`) directly; no additional version pin needed. |
| Proportion / zero-one-inflated beta | `fast_zero_one_inflated_beta_cpp` | (custom) | none identified | Baseline Gap — see below. |
| Survival / Cox PH | `fast_coxph_regression_cpp` | `survival::coxph.fit` | `sksurv.linear_model.CoxPHSurvivalAnalysis` (scikit-survival) | Prefer scikit-survival over `lifelines.CoxPHFitter` for the bare-metal row. **Verified at sksurv 0.28.0**: `.fit(self, X, y)` takes a raw 2D array `X` and a structured array `y` (event/time fields) directly — no DataFrame construction; `lifelines.CoxPHFitter.fit` requires a `pandas.DataFrame` (`df, duration_col, event_col=...`), confirming the stated overhead difference. |
| Survival / Weibull AFT | `fast_weibull_regression_cpp` | `survival::survreg` | `lifelines.WeibullAFTFitter` | **Verified at lifelines 0.30.0**: `.fit(self, df, duration_col, event_col=None, ...)` requires a `pandas.DataFrame` — the bare-metal timing row must include DataFrame construction from the NumPy fixture, or explicitly flag it as formula-API overhead the same way `benchmark_model_fits.md` already flags formula-API fallbacks. |
| Survival / Weibull frailty | `fast_weibull_frailty_cpp` | (custom/`frailtypack`) | none identified | Baseline Gap — see below. |
| Ordinal / proportional odds (logit/probit) | `fast_ordinal_regression_cpp`, `fast_ordinal_probit_regression_cpp` | `ordinal::clm` | `statsmodels.miscmodels.ordinal_model.OrderedModel(distr="logit"/"probit")` | **Verified at 0.14.6**: `distr` accepts the string shortcuts `'probit'`/`'logit'` only for those two links (see source `__init__`: `if distr == 'probit': ... elif distr == 'logit': ...`, no other string branches). |
| Ordinal / cauchit, cloglog | `fast_ordinal_cauchit_regression_cpp`, `fast_ordinal_cloglog_regression_cpp` | `ordinal::clm(link=...)` | `statsmodels.miscmodels.ordinal_model.OrderedModel(distr=scipy.stats.cauchy)` / `OrderedModel(distr=scipy.stats.gumbel_l)` | **No longer a Baseline Gap** — verified both fit and converge on synthetic data at statsmodels 0.14.6 via a `scipy.stats.rv_continuous` instance passed as `distr` (not a string). Not the same wiring as the logit/probit row (needs `scipy.stats` imported and the right distribution picked), so keep the two rows separate in `baselines.py` rather than templating them together. |
| Ordinal / adjacent-category, continuation-ratio, stereotype | `fast_adjacent_category_logit_cpp`, `fast_continuation_ratio_regression_cpp`, `fast_stereotype_logit_cpp` | `VGAM::vglm` | none identified | Baseline Gap — see below. |
| GEE | `fast_gee_cpp` | `geepack::geeglm` | `statsmodels.genmod.generalized_estimating_equations.GEE` | **Verified aligned, not just "direct match":** both default to an independence working-correlation structure — `geepack::geeglm`'s `corstr` defaults to `"independence"`, and `GEE(...).cov_struct` (when left `None`) resolves to `statsmodels.genmod.cov_struct.Independence` — confirmed by inspecting the constructed object at 0.14.6. No mismatch to flag. |
| GLMM (all families) | `fast_poisson_glmm_cpp`, `fast_logistic_glmm_cpp`, `fast_hurdle_poisson_glmm_cpp`, `fast_ordinal_glmm_cpp`, `fast_ordinal_clmm_cpp`, `fast_gaussian_lmm_cpp` | `glmmTMB`/`lme4` | none identified | Baseline Gap — see below. |
| KK combined (matched + reservoir) | `fast_clogit_plus_glmm_cpp`, `fast_cpoisson_combined_cpp` | (custom KK estimator) | none identified | Baseline Gap — no canonical analog exists in either language; correctness reference is the R implementation itself. |

### Baseline Gaps

Some EDI families have no clean, actively-maintained Python canonical
equivalent as of this writing:

- **GLMM families** (`fast_poisson_glmm_cpp`, logistic GLMM,
  `fast_hurdle_poisson_glmm_cpp`, ordinal GLMM/CLMM, Gaussian LMM,
  `fast_clogit_plus_glmm_cpp`, Weibull frailty): no pure-Python package
  offers comparable maximum-likelihood GLMM fitting with adaptive
  Gauss-Hermite quadrature (`statsmodels`'s mixed-GLM support is variational/
  Bayesian, not a like-for-like MLE comparison). Do not force a
  mismatched comparison. Options, in order of preference:
  1. Mark these kernels **correctness-only, no speed baseline** in the
     report, same visual treatment `benchmark_model_fits.md` already uses
     for "NA timing comparisons" (light grey rows).
  2. If a correctness cross-check is still wanted, call R's `glmmTMB`
     (already an `EDI` `Suggests:` dependency) via `rpy2` from the Python
     test suite — this is a correctness reference, not a timing baseline,
     and should never appear in the timing table.
- **Adjacent-category / continuation-ratio / stereotype ordinal links**: no
  identified Python package implements these link functions. Same treatment
  as GLMM — correctness-only, no synthetic baseline substitution (do not,
  for example, silently compare against proportional-odds as if it were
  equivalent — the link functions are not interchangeable).
- **Zero-one-inflated beta regression, Weibull frailty, KK-combined
  (matched + reservoir) estimators**: custom EDI estimator families with no
  canonical package in either R or Python — correctness-only against the R
  implementation itself (which is already tested independently in
  `EDI/tests/testthat/`), not a cross-package baseline.

Do not invent an approximate substitute baseline for a Baseline Gap kernel
just to fill a table cell. An absent row with a documented reason is more
honest than a mismatched comparison with a speedup number that doesn't mean
what it looks like it means.

## Implementation Phases

### Phase 1: Infrastructure + Pilot

- Scaffold `python/` layout, `pyproject.toml`, `CMakeLists.txt`.
- Add `result_map_pybind.h`.
- Bind `fast_poisson_glmm` (the same pilot kernel from the SEXP-removal
  spec) end-to-end: binding, parity test against the R fixture, one
  baseline benchmark row (`statsmodels GLM(family=Poisson())` — note this
  compares against the *non-mixed* Poisson GLM as a sanity baseline; the
  GLMM comparison itself is a Baseline Gap per above).
- Get `pip install ./python` working locally before binding anything else.

### Phase 2: Remaining Model-Fitting Kernel Families

- Bind the remaining 32 files in the same family grouping the SEXP-removal
  spec uses (continuous/OLS, count/ordinal, incidence/survival/proportion,
  KK-combined/GLMM last).
- One `tests/test_<kernel>.py` and one `benchmarks/baselines.py` entry per
  bound kernel, added together — don't let bindings outpace tests.
- No RNG-gating and no MKL-gating apply anywhere in this phase (see R
  Dependency Audit above) — every file in scope can be bound as soon as its
  SEXP-removal-spec migration is done.

### Phase 3: Benchmark Report

- Implement `benchmarks/run_benchmark_audit.py` producing a Markdown table
  in the exact column shape of `benchmark_model_fits.md`.
- Cross-link the two reports from each other once both exist.

## TODO Checklist

- [x] TODO-1: Scaffold `python/` package layout, `pyproject.toml`
      (scikit-build-core backend), top-level `CMakeLists.txt`.
- [x] TODO-2: Add `result_map_pybind.h`.
- [x] TODO-3: Bind `fast_poisson_glmm` as the pilot; add its parity test and
      fixture.
- [x] TODO-4: Verify `statsmodels`/`scikit-survival`/`lifelines` versions
      and exact class/method names in the Baseline Registry table before
      writing `benchmarks/baselines.py` — several entries above are marked
      "verify before wiring."
      **Status update (2026-08-03):** verified against the `edi-py-bench`
      venv's installed versions (`statsmodels==0.14.6`,
      `scikit-survival==0.28.0`, `lifelines==0.30.0`) — every class named
      in the Baseline Registry table exists and imports cleanly; no missing
      classes, no version-floor pins needed. Findings recorded inline in
      the table above and its preceding status note: `OrderedModel` cauchit/
      cloglog moved from Baseline Gap to a real (verified working, fit and
      converge) baseline via `scipy.stats` distribution instances rather
      than string shortcuts; `RLM`'s default resolves to `HuberT`
      M-estimation while EDI's own default is `MM`-estimation (a real
      estimator-class mismatch, not just a tuning-constant difference —
      use `method="M"` on the EDI side for a fair comparison); GEE's
      independence-correlation defaults confirmed aligned on both sides;
      `HurdleCountModel`'s zero-mechanism (censored count) confirmed
      structurally different from EDI's (independent logistic hurdle);
      `sksurv.CoxPHSurvivalAnalysis`/`lifelines.WeibullAFTFitter`'s raw-array
      vs. DataFrame-required APIs confirmed, validating the doc's existing
      stated preference for scikit-survival on the Cox PH bare-metal row.
- [x] TODO-5: Implement `benchmarks/baselines.py` for the families with a
      confirmed clean baseline (OLS, robust regression, logistic, probit,
      log-binomial, Poisson, NegBin, ZINB/ZAP/hurdle, beta regression, Cox
      PH, Weibull AFT, proportional-odds/cauchit/cloglog ordinal, GEE).
      **Status update (2026-08-03):** implemented at
      `python/benchmarks/baselines.py`, the exact path this doc's Package
      Layout specifies — a standalone `BASELINES` registry (20 entries: the
      12 families above, several expanding to multiple kernels, e.g.
      log-binomial + identity-binomial, ZINB + ZAP + hurdle-NegBin, and the
      4 ordinal links), each a `Baseline(package, function, fit, notes)`
      dataclass whose `fit(X, y, ...)` callable takes raw NumPy arrays with
      no dataset generation, EDI timing, or report rendering mixed in (that
      stays `run_benchmark_audit.py`'s job, per Package Layout). Encodes
      every TODO-4 finding directly (cauchit/cloglog via `scipy.stats`
      distribution instances, `RLM(M=HuberT)` for a fair robust-regression
      comparison against EDI's `method="M"`, `HurdleCountModel`'s
      structurally different zero-mechanism noted in both the module
      docstring and the registry entry's `notes`) — corrects the two gaps
      found in `package_metadata/benchmark_model_fits_python.py`'s
      equivalent (co-located, not standalone) baseline-wiring: that file
      has no GEE row and still treats cauchit/cloglog as EDI-only/no-baseline.
      Verified by calling every one of the 20 registry entries against
      realistic synthetic data for its family (continuous/binary/count/
      proportion/survival/ordinal/GEE) in the `edi-py-bench` venv — all 20
      fit successfully and return a real fitted-model object.
- [x] TODO-6: Phase 2 — bind the remaining 32 model-fitting files,
      following the SEXP-removal spec's migration order.
      **Status update (2026-08-03):** all 33 model-fitting kernels are now
      bound (33 in `EDI_KERNEL_SOURCES`/`python/cpp/bindings_*.cpp`, verified
      via a from-scratch `cmake --build` and direct smoke-test calls across
      every family — ordinal, GLMM/CLMM/LMM, survival, count, proportion,
      binary/log-binomial).
      **Status update (2026-08-04):** the "with tests" half is now also
      done — every one of the 37 bound Python functions (all 33 files'
      primary entry point, plus `fast_log_binomial_regression_with_var`/
      `fast_identity_binomial_regression_with_var`, the two secondary
      entry points that happen to be separately exported) has a real
      `python/tests/test_<kernel>.py` parity fixture computed via
      `Rscript` calling the installed `EDI:::<kernel>_cpp(...)` on the
      exact same synthetic data, checked at `atol=1e-9, rtol=1e-9`. 102
      tests total, all passing. Two real result-shape gotchas surfaced and
      got documented in the affected tests rather than papered over: (1)
      several kernels (`fast_gaussian_lmm`, `fast_coxph_regression`,
      `fast_weibull_regression`) omit `vcov`/`params` keys entirely under
      `estimate_only=True` rather than setting them to `None` the way
      `fast_poisson_glmm` does — inconsistent across kernels, now
      documented per-file rather than assumed uniform; (2)
      `fast_zero_one_inflated_beta_cpp`'s R-facing return list has no
      `converged`/`iterations`/`gradient_norm` fields at all even though
      the Python binding's `LikelihoodFitResult` surfaces them — that
      test only cross-checks the fields both sides actually expose
      (`params`/`neg_loglik`). Also found (during the Group 1 batch) one
      genuine R/Python default-value mismatch, plus one false positive that
      looked like a second mismatch but wasn't:
      `fast_neg_bin`'s `smart_cold_start` default really did differ (R
      `FALSE` — confirmed across all three R-facing entry points,
      `fast_neg_bin_cpp`/`_with_var_cpp`/`_weighted_cpp` in
      `EDI/src/fast_negbin_regression.cpp`, each explicitly overriding the
      shared core's own `true` default — vs Python `true` in
      `bindings_count.cpp`).
      **Fixed 2026-08-05:** the Python binding's default changed to
      `false`; `python/tests/test_fast_neg_bin.py` gained a dedicated
      `test_matches_r_fixture_default_smart_cold_start` proving the two
      defaults now produce the identical fit (fresh R fixture, both sides'
      `smart_cold_start` omitted), and the extension was rebuilt/reinstalled
      to pick up the change.
      `fast_log_binomial_regression`/`fast_identity_binomial_regression`'s
      `tol` was originally also logged here as R `1e-8` vs Python `1e-6` —
      **corrected 2026-08-05: this was never a real mismatch.** Both
      `EDI/R/RcppExports.R` and `EDI/src/fast_log_binomial_regression.cpp`
      show R's own default is `tol = 1e-6`, identical to the Python
      binding's. The original claim was a documentation error (likely
      confused with `fast_neg_bin_cpp`'s unrelated, unused `eps_f = 1e-8`
      parameter, which sits next to `eps_g = 1e-6` in a similarly-shaped
      signature). The four affected test docstrings
      (`test_fast_log_binomial_regression.py` and its identity-link and
      `_with_var` siblings) were corrected accordingly; their `tol=1e-8`
      calls remain, but as a deliberate tighter-than-default convergence
      choice for fixture stability, not a mismatch workaround.
- [x] TODO-7: Implement `benchmarks/run_benchmark_audit.py` and generate the
      first Python-side benchmark report.
      **Status update (2026-08-03):** implemented, but at
      `package_metadata/benchmark_model_fits_python.py` /
      `.html`, not the `python/benchmarks/run_benchmark_audit.py` path this
      doc's Package Layout specifies — reconcile the location (move the
      script, or update Package Layout to match reality) before calling
      this fully settled. Verified the generated report's column headers
      are byte-identical to `benchmark_model_fits_R.html`'s (`Class`,
      `Response`, `EDI Time (ms)`, `Canonical Pkg`, `Canonical Func`,
      `Canonical Time (ms)`, `Speedup`, `Timing Pval`), satisfying the
      "same column shape" requirement. The existing HTML report is dated
      2026-08-02 19:30:47 UTC — generated before this session's kernel-
      binding work finished — so it's a real "first report" but a *stale*
      one; most EDI-timing cells will read differently (many currently-NA
      rows should now have real timings) once regenerated against the
      now-complete 33-kernel binding.
- [x] TODO-8: Document Baseline Gaps explicitly in the generated report
      (GLMM/CLMM/LMM families, adjacent-category/continuation-ratio/
      stereotype ordinal, zero-one-inflated beta, Weibull frailty,
      KK-combined estimators) rather than omitting them silently.
      **Status update (2026-08-03):** completed the gap the first pass over
      this TODO found — adjacent-category, continuation-ratio, zero-one-
      inflated beta, Weibull frailty, and the KK-combined estimators were
      already documented gap rows, but the entire GLMM/CLMM/LMM family
      (`fast_poisson_glmm`, `fast_logistic_glmm`, `fast_hurdle_poisson_glmm`,
      `fast_ordinal_glmm`, `fast_ordinal_clmm`, `fast_gaussian_lmm`) and
      `fast_stereotype_logit` were silently absent from both the
      point-estimate (`MODEL_SPECS`) and Wald (`WALD_SPECS`) tables in
      `package_metadata/benchmark_model_fits_python.py`. Added all 7 as
      explicit `no_canonical=True` rows to both tables (matching the
      existing `(family, R6_class, None, None, edi_kernel, builder, True)`
      shape, `builder=None` — consistent with the existing KK-combined
      rows, real bare-metal EDI timing for these is a follow-up, not
      required for "don't omit silently"), with per-kernel
      `NO_CANONICAL_NOTE` explanations (mostly: no pure-Python package
      offers ML GLMM fitting with adaptive Gauss-Hermite quadrature, same
      reasoning as this doc's own Baseline Gaps section). R6 class names
      resolved by grepping `EDI/R` for each kernel's actual caller (e.g.
      `fast_poisson_glmm_cpp` → `InferenceCountKKGLMM` in
      `inference_count_KK_combined.R`); `fast_logistic_glmm` turned out to
      have no R6 consumer class at all yet (bound and working, but nothing
      in `EDI/R` calls it) — labeled explicitly as such rather than
      inventing a class name. Verified by importing the module and
      confirming all 7 keys resolve in both tables with a note, and that
      `python -m py_compile` and a full module import still succeed
      (`main()` is `__main__`-guarded, so this didn't trigger a full
      benchmark run).
- [x] TODO-9: Migrate the 15 remaining direct Rmath call sites identified in
      the R Dependency Audit (`R::pchisq` x3, `R::qnorm5` x3,
      `R::pnorm`/`R::pnorm5`/`R::dnorm` x6, `R::dnbinom_mu` x1, `R::lbeta`
      x1) to `fast_*` replacements, same pattern as the existing
      `fast_digamma`/`fast_lgamma`/`fast_trigamma`/`fast_erfc` work — this
      isn't a hard blocker for binding (the functions work fine linked
      against R's Rmath in the meantime) but should happen before or during
      Phase 2 rather than being carried indefinitely.
      **Status update (2026-08-03):** a repo-wide grep found the original
      count of 15 was stale — `R::dnorm`, `R::dnbinom_mu`, and `R::lbeta`
      had already been fully migrated in earlier work (0 call sites left),
      leaving only 3 genuine remaining call sites: `R::pchisq` in
      `lrt_ci_newton.cpp` → `fast_pchisq_upper`, `R::qnorm` in
      `newcombe_speedups.cpp` → `fast_qnorm`, and `R::pnorm` in
      `miettinen_nurminen_speedups.cpp` → `pnorm_fast` (both in
      `fast_erfc.h`). All 3 migrated; verified numerically against base R
      (`pnorm` matches to full double precision, `pchisq` to ~14 digits,
      `qnorm` to the ~1.2e-9 relative error already documented/accepted for
      `fast_qnorm`'s other call sites) and via the existing test suite
      (`test-miettinen-nurminen-ci-dispatch.R`, `test-bartlett-lr-*.R`,
      `test-likelihood-test-memoization.R` — all pass). Zero call sites
      remain within this TODO's original scope (pchisq/qnorm/pnorm/dnorm/
      dnbinom_mu/lbeta). A repo-wide sweep for `R::` also turned up two
      deterministic Rmath calls that were never in this TODO's list —
      `R::dchisq` (`lrt_ci_newton.cpp:61`, chi-sq density for the LRT CI's
      dp/d(delta)) and `R::dbinom` (`zhang_exact_speedups.cpp`, log-binomial
      density for Zhang's exact test) — no `fast_dchisq`/`fast_dbinom`
      exist yet, so these are left as-is and flagged here as a candidate
      follow-up, not silently folded into this TODO. The many `R::unif_rand`/
      `R::rbeta`/`R::rpois`/`R::rweibull`/`R::runif` call sites are a
      different category (RNG draws, not distribution math) and were
      correctly out of scope all along — they must stay coupled to R's own
      RNG stream for `set.seed()` reproducibility of randomized designs and
      bootstrap resampling; migrating them is not a "port to `fast_*`" task
      at all.
- [x] TODO-10: Convert every `Rcpp::Nullable<T>` parameter across the 33
      files to `std::optional<T>` (or `const T*`) on the `*_core` signature,
      per Optional Arguments above, as part of each file's SEXP-removal-spec
      migration — not as a separate pybind11-side shim.
      **Status update (2026-08-03):** verified via a repo-wide grep — zero
      `Rcpp::Nullable` references remain in any `*_internal`/`*_core`
      function signature across all 33 files. Remaining `Nullable` usages
      are exclusively in the (correctly still-Rcpp) R-facing wrapper
      functions, which convert via `nullable_to_optional<>()` before
      calling into the portable core, per this doc's own convention.

## Acceptance Criteria

The feature is complete when:

- `pip install ./python` builds the extension from a clean checkout.
- Every bound kernel has a parity test passing against an R-generated
  fixture at `atol=1e-9, rtol=1e-9`.
- `benchmarks/run_benchmark_audit.py` produces a report in the same column
  shape as `benchmark_model_fits.md`, with Baseline Gap kernels explicitly
  marked rather than omitted.
- No kernel core was duplicated or reimplemented in Python or in the
  `cpp/bindings_*.cpp` glue — every binding calls directly into the
  `EDI/src/*_core` functions.
- All 33 model-fitting files are bound; no resampling/design/nonparametric-
  test file was pulled in under this spec.
- Every `Rcpp::Nullable<T>` argument across the bound files has a working
  `std::optional`-based Python equivalent with a passing omitted-argument
  test.

**Status update (2026-08-03) — verified against each criterion in turn:**

1. **`pip install ./python` builds from a clean checkout — MET (after a
   fix).** Found genuinely broken: `pyproject.toml` declared
   `readme = "README.md"` but `python/README.md` didn't exist, so pip
   failed at the metadata-generation stage before ever reaching the CMake
   build — a real bug this session's earlier `cmake --build` testing never
   caught, since bypassing `pip`/scikit-build-core's packaging layer
   entirely also bypasses this check. Fixed by adding `python/README.md`;
   retested with a fresh venv (`pip install ./python`) — wheel builds,
   installs, and `import edi_kernels` exposes all 38 functions.
2. **Every bound kernel has an R-fixture parity test at
   `atol=1e-9, rtol=1e-9` — MET (2026-08-04).** All 37 bound Python
   functions (33 files' primary entry point + 2 secondary `_with_var`
   entry points) now have a `python/tests/test_<kernel>.py` with a real
   R-computed fixture, verified passing (`102 passed` for the full
   `python/tests/` suite). One real, pre-existing R/Python default-value
   mismatch surfaced along the way (`fast_neg_bin`'s `smart_cold_start`)
   and a second suspected one turned out to be a documentation error, not
   an actual mismatch (`fast_log_binomial_regression`/
   `fast_identity_binomial_regression`'s `tol`) — both resolved 2026-08-05,
   see the TODO-6 status note above.
3. **Benchmark report, same column shape, Baseline Gaps marked — MET
   (2026-08-05).** See TODO-7/TODO-8 status notes above: column shape
   confirmed identical to the R report; the GLMM/CLMM/LMM family plus
   `fast_stereotype_logit` were silently omitted rather than marked as
   Baseline Gaps (the specific failure TODO-8 exists to prevent) — fixed,
   all 7 added as documented gap rows to both `MODEL_SPECS` and
   `WALD_SPECS`. The remaining location/staleness gap (script path per
   Package Layout, regenerating the report artifact against the completed
   bindings) was closed in a separate, concurrent agent session — not
   independently re-verified here.
4. **No kernel core duplicated/reimplemented — MET.** Every `_internal`
   function extracted this session is a direct call into (or thin
   restructuring of) the existing `EDI/src/*.cpp` logic; every
   `cpp/bindings_*.cpp` file is argument marshaling + a call into that
   `_internal` function, never a reimplementation. Consistent with the
   "extract, don't duplicate" discipline used throughout Task 45.
5. **All 33 model-fitting files bound, nothing out-of-scope pulled in —
   MET.** 37 kernel functions (33 files, some with multiple entry points
   e.g. log/identity-binomial) bound and confirmed importable; `CMakeLists.txt`'s
   `EDI_KERNEL_SOURCES` lists exactly the 33 in-scope files, no
   resampling/design/nonparametric-test `.cpp` pulled in.
6. **Every `Nullable<T>` has a working `optional` equivalent with a passing
   omitted-argument test — MET (2026-08-05).** `python/tests/test_omitted_arguments.py`
   added: one dedicated test per bound model-fitting kernel (37 of them —
   every `__all__` entry except `fast_pchisq_upper`, which isn't a
   `Nullable`-bearing model-fitting kernel), each calling the kernel with
   *only* its required (no-default) arguments — i.e. every
   `std::optional<T>` parameter genuinely omitted, not just left at an
   explicitly-passed default — and asserting the call succeeds and returns
   a finite, well-formed result. Data is reused from each kernel's own
   parity-test module via sibling import, so the two suites can't drift on
   what "valid input" means for a given kernel. All 37 pass
   (`python/tests/` full suite: 139 passed — 102 parity + 37
   omitted-argument). This is a genuinely new check, not a restatement of
   the parity tests: several parity tests pass some optional args
   explicitly at their default value (e.g. `estimate_only=False`), which
   is not the same as proving the argument is safely *omittable*.

**Net (2026-08-05): all 6 Acceptance Criteria are now MET, and the
R/Python default-value reconciliation flagged as follow-up work is also
done.** `fast_neg_bin`'s `smart_cold_start` mismatch was real and is
fixed (Python's default changed from `true` to `false` to match R, backed
by a fresh cross-language fixture test); the
`fast_log_binomial_regression`/`fast_identity_binomial_regression` `tol`
"mismatch" turned out to be a documentation error on investigation — both
sides already defaulted to `1e-6`, so nothing needed changing there beyond
correcting the four test docstrings that repeated the error.

## Release Checklist (PyPI)

### Licensing — resolve this first, before anything else here

**Status: DONE (2026-08-05).** `python/pyproject.toml` now declares
`license = "GPL-3.0-only"` and `license-files = ["LICENSE"]`;
`python/LICENSE` is a byte-identical copy of the repo-root `LICENSE`
(verified via `diff`). `GPL-3.0-only` (not `-or-later`) was confirmed as
the correct SPDX variant: R/CRAN's licensing convention treats a bare
`License: GPL-3` in `DESCRIPTION` as version-3-only, distinct from
`GPL (>= 3)` which would map to `-or-later` — `EDI/DESCRIPTION` uses the
former. `pip install -e .` was rebuilt against the new metadata and still
installs cleanly.

`EDI/DESCRIPTION` declares `License: GPL-3`, and this extension compiles
`EDI/src/*.cpp` **directly** into the wheel (not a copy, not a subprocess
call to a separately-licensed R process) — that makes `edi_kernels` a
combined/derivative work under GPL-3, not a work that can be released
under a separate permissive license (MIT, BSD, Apache-2.0, etc.) no matter
what license header `python/` itself carries. Concretely:

- `python/pyproject.toml` needs a `license` field matching the confirmed
  SPDX identifier above, and a `python/LICENSE` file (copy of the repo
  root `LICENSE`, not a paraphrase). Done, see status note above.
- The two fetched build dependencies are both GPL-3-compatible as
  dependencies of a GPL-3 work: Eigen is MPL-2.0, LBFGSpp is MIT — neither
  imposes a stronger copyleft that would conflict. No dependency swap is
  needed on the licensing axis (this was already the deciding factor for
  fetching LBFGSpp from its own upstream instead of RcppNumerical's copy —
  see "Standalone Rmath Library Dependency" above — but confirm it again
  here since a *release* forces the question in a way local development
  didn't).
- Do not let a generic `pip install` template (e.g. a copy-pasted
  `pyproject.toml` `license` field from an unrelated permissively-licensed
  project) silently ship as anything other than GPL-3 — that would be
  either a real license violation or a broken build once someone notices
  and re-licenses it, either way not a one-line fix after the fact once a
  wheel is public on PyPI.
- If distributing under GPL-3 is not actually wanted (e.g. the intent is a
  permissively-licensed thin wrapper around a *separately installed* EDI),
  that requires re-architecting this package to link a shared library at
  runtime rather than compiling `EDI/src/*.cpp` in statically — a materially
  different design than what this spec builds, and should be a deliberate
  decision made before Phase 1, not discovered at release time.

### Versioning

**Status: DONE (2026-08-05).** `python/pyproject.toml`'s `version` is set
to `1.0.0`, matching `EDI/DESCRIPTION`'s `Version` field. No `.postN`
suffix yet since no Python-packaging-only change has landed since this was
set — bump per the scheme below the next time one does. Git tagging
(`py-vX.Y.Z`) not yet done — that's a release-time action, not a
metadata-file change.

- `python/pyproject.toml`'s `version` should track `EDI/DESCRIPTION`'s
  `Version` field (currently `1.0.0`) rather than drift independently — a
  Python user comparing behavior against the R package needs the version
  numbers to mean the same underlying kernel code. Simplest scheme:
  `{EDI version}.postN` where `N` increments for Python-packaging-only
  changes (binding fixes, build-system changes) that don't touch
  `EDI/src/*.cpp` at all; bump to match `EDI/DESCRIPTION` directly whenever
  it changes.
- Tag releases in git as `py-vX.Y.Z` (distinct from any R-side release tag
  convention already in use) so `git describe`/CI can disambiguate which
  side of the repo a tag is about, since both packages share this one repo.

### Building portable wheels

**Status: DONE (2026-08-05).** `python/pyproject.toml` gained a
`[tool.cibuildwheel]` section (with `[tool.cibuildwheel.linux/macos/windows]`
sub-tables) and `.github/workflows/build-wheels.yml` was added to actually
invoke it in CI, matrixed across `ubuntu-latest`/`macos-latest`/
`windows-latest` (`package-dir: python`, since this repo's `pyproject.toml`
isn't at the repo root), triggered on `py-v*` tags plus manual dispatch.
Concretely:
- `EDI_PY_PORTABLE` is set via `config-settings =
  {"cmake.define.EDI_PY_PORTABLE" = "ON"}` in `[tool.cibuildwheel]` only —
  the local-dev default in `[tool.scikit-build.cmake.define]` stays `OFF`,
  per this checklist's own preference for the config-settings route over
  flipping the project-wide default.
- **Correction to this checklist's own text below:** manylinux images are
  AlmaLinux/CentOS-based and use `dnf`/`yum`, not `apt-get` — the
  `apt-get install libblas-dev` suggestion here would fail outright.
  `[tool.cibuildwheel.linux]` uses `manylinux_2_28` (needed for a C++20-
  capable compiler; the default manylinux2014's devtoolset doesn't
  reliably support C++20) with `before-all = "dnf install -y
  openblas-devel || yum install -y openblas-devel"`.
- macOS: `archs = ["x86_64", "arm64"]` (both, per this checklist's own
  no-Rosetta-assumption note); `before-all = "brew install libomp"` so
  `find_package(OpenMP)` actually succeeds (Apple clang has no bundled
  OpenMP runtime) — BLAS itself needs no install, since CMake's
  `find_package(BLAS)` detects the built-in Accelerate framework on macOS
  automatically.
- Windows has neither a system package manager nor an Accelerate
  equivalent, so `CMakeLists.txt`'s `find_package(BLAS REQUIRED)` was
  changed to `find_package(BLAS QUIET)` with a fallback: if not found, it
  shells out to the active Python interpreter to locate the pip-installable
  `scipy-openblas32` package's CMake config dir (the same package numpy/
  scipy/scikit-learn's own cibuildwheel configs use to solve this exact
  problem) and retries. `[tool.cibuildwheel.windows]` installs
  `scipy-openblas32` + `delvewheel` via `before-build`, with
  `repair-wheel-command` wired to `delvewheel repair` (Windows has no
  `auditwheel`/`delocate` equivalent built in). This fallback path is
  logically sound but **has not been exercised on a real Windows runner or
  manylinux container** — only verified locally on Linux, where
  `find_package(BLAS QUIET)` still finds the system BLAS directly and the
  scipy-openblas32 branch is never triggered (confirmed via a full
  `pip install -e .` rebuild + the `python/tests/` suite, all passing).
  Flag this for a real dry run on TestPyPI/CI before the first release, per
  the "Pre-release testing" item below.
- `test-command`/`test-extras` in `[tool.cibuildwheel]` run the full
  `python/tests/` suite against each built wheel (not just the raw
  `cmake --build`), directly satisfying the "Pre-release testing" item's
  "wheel-specific test run" requirement below.

- `EDI_PY_PORTABLE=ON` (already wired in `CMakeLists.txt`, currently
  `OFF` by default in `[tool.scikit-build.cmake.define]`) must be `ON` for
  any wheel that leaves this machine — `-march=native` bakes in
  instruction-set assumptions (AVX-512, this machine's specific `-mtune`)
  that will `SIGILL` on an older or different CPU. Flip the default in
  `pyproject.toml` for release builds, or pass
  `--config-settings=cmake.define.EDI_PY_PORTABLE=ON` per-build via
  `cibuildwheel`'s config — do not ship a `-march=native` wheel to PyPI.
- Use [`cibuildwheel`](https://cibuildwheel.pypa.io/) to build manylinux
  (glibc, matching the R package's own Linux support target),
  macOS (both `x86_64` and `arm64` — build both, don't assume Rosetta
  coverage is acceptable for a numerical library), and Windows wheels in
  CI. `cibuildwheel`'s `manylinux` containers don't have a system
  BLAS/Eigen preinstalled, but that's fine here — `CMakeLists.txt` already
  `FetchContent`s Eigen when no system copy is found, and `find_package(BLAS)`
  needs at minimum a reference BLAS available in the container image (add
  `apt-get install -y libblas-dev` or equivalent to `CIBW_BEFORE_ALL` for
  the manylinux job if the chosen base image doesn't already have one —
  confirmed portable per the `_helper_functions_core.h` comment on
  `cblas_dsyrk`, see Build System above, but the container still needs
  *some* BLAS installed to link against).
- `FetchContent`-ing Eigen/LBFGSpp from their upstream git repos at build
  time (see Build System above) means every wheel build needs network
  access — `cibuildwheel`'s containers have it by default, but a locked-down
  CI runner or an offline sdist build will fail. For a `pip install
  edi_kernels` **source** install (as opposed to a prebuilt wheel) to work
  fully offline, either vendor pinned copies of Eigen/LBFGSpp headers into
  the sdist (reintroducing the "is this a copy" question this whole spec's
  `EDI/src` include-not-copy discipline was built to avoid, so only do this
  for the *two portable dependencies*, never for `EDI/src` itself) or
  clearly document the network requirement in the README's install section
  and accept that offline installs need a prebuilt wheel.

**Naming:** the distribution name (`python/pyproject.toml`'s `[project]
name`) and the import name (`import edi_kernels`) are both `edi_kernels`,
deliberately kept identical rather than split like `scikit-learn`/
`sklearn` — considered and rejected a `edi-kernels`/`edi_kernels` split
(2026-08-05) in favor of one name to remember and type everywhere.

### Pre-release testing

- Run the full `python/tests/` suite (once the parity-test gap noted in
  Acceptance Criteria above is closed — releasing before every kernel has
  a real R-fixture parity test means shipping kernels whose only
  verification was this session's ad hoc smoke tests, not a committed,
  re-runnable check) against each built wheel, not just against a local
  `cmake --build` — a wheel-specific test run catches packaging bugs (like
  this session's missing-`README.md` metadata failure, or a missing
  `wheel.packages` entry) that a raw CMake build never exercises.
  **Status (2026-08-05): MET** — `python/tests/` is 100% parity-covered
  (see Acceptance Criteria above) and `[tool.cibuildwheel]`'s
  `test-command`/`test-extras` (see "Building portable wheels" above) run
  that same suite against every built wheel in CI, not just a raw
  `cmake --build`.
- [ ] TODO-11: Upload to [TestPyPI](https://test.pypi.org/) first, `pip install
  --index-url https://test.pypi.org/simple/ edi_kernels` into a fresh venv
  on at least one machine that isn't the one that built the wheel, and
  re-run the smoke tests from this session (or the parity suite) against
  that install before promoting to the real index.
  **Status (2026-08-05): needs a human with PyPI/TestPyPI credentials.**
  An actual TestPyPI upload requires an account and API
  token this environment doesn't have and shouldn't create unattended —
  that's a genuine "upload to a shared external index" action, not a local
  build step. What *was* done as the closest available proxy: built a real
  sdist (`python -m build --sdist`, succeeded, `edi_kernels-1.0.0.tar.gz`)
  and ran `twine check` on it — the same metadata/README-rendering
  validation PyPI's own upload endpoint runs — which passed. This catches
  the class of bug this checklist item exists for (the missing-README.md
  failure from earlier in this doc would have been caught by `twine
  check` too), but it is not a substitute for a real TestPyPI round-trip;
  do that before the first real release.
- Confirm the package name `edi_kernels` is actually available on PyPI
  (not squatted or already used by an unrelated project) before the first
  real upload.
  **Status (2026-08-05): MET.** Checked directly (`curl
  https://pypi.org/pypi/edi_kernels/json` and the `edi-kernels`
  hyphen form, plus both forms on `test.pypi.org`) — all four return 404,
  i.e. unregistered on both indexes. This doc's earlier claim that
  availability "can't be checked from this environment" was itself wrong;
  PyPI's JSON API is a plain unauthenticated GET, no credentials needed —
  only the *upload* step needs an account.

### Publishing

**Status: DONE, minus one manual step only a PyPI account holder can do
(2026-08-05).** `.github/workflows/build-wheels.yml` gained a `publish`
job (`needs: [build_wheels, build_sdist]`) that runs
`pypa/gh-action-pypi-publish` via OIDC Trusted Publishing (`permissions:
id-token: write`, `environment: pypi`) — gated to fire only on an actual
`py-v*` tag push (`if: github.event_name == 'push' && ...`), never on a
manual `workflow_dispatch` build. It downloads the wheel artifacts from
every `build_wheels` matrix leg plus the `build_sdist` artifact into one
`dist/` directory before publishing, so a single tag push ships every
platform's wheel plus the sdist together. **What's NOT done, and can't be
from this environment:** registering this repo as a Trusted Publisher for
the `edi_kernels` PyPI project (https://pypi.org/manage/account/publishing/
for a first-time project name) — that's a one-time step in the PyPI web
UI under the project owner's own account, documented as a prerequisite
comment at the top of `build-wheels.yml` itself (workflow name
`build-wheels.yml`, environment name `pypi` — must match exactly what
PyPI's form asks for). Until that's done, the first real tag push will
fail at the publish step with an OIDC/authorization error, not silently
publish somewhere wrong.
- Use PyPI's [Trusted Publishing](https://docs.pypi.org/trusted-publishers/)
  (GitHub Actions OIDC) rather than a long-lived API token stored as a
  repo secret — no token to rotate or leak, and it's the currently
  recommended path for a project already using GitHub Actions for CI.
- `python -m build` (or `cibuildwheel` output directly) + `twine upload` /
  the `pypa/gh-action-pypi-publish` action, triggered on a `py-v*` git tag
  push, matching the versioning scheme above — don't publish from a local
  machine's ad hoc build, so every published artifact is reproducibly tied
  to a specific commit via CI.
- Publish an sdist alongside the wheels (`python -m build --sdist`) even
  though most users will pull a prebuilt wheel — it's what makes `pip
  install edi_kernels` work at all on a platform `cibuildwheel` didn't
  cover, and satisfies the PyPI convention that a project ships both.

### Post-release

- Cross-link this package from the R package's own `README`/`DESCRIPTION`
  (`URL:`/`BugReports:` fields) once published, so a user landing on either
  package's listing can find the other — same spirit as the "Cross-link
  the two reports from each other" note in Phase 3 above, applied to the
  package listings themselves, not just the benchmark reports.
  **Status: DONE (2026-08-05), informally rather than via package
  metadata.** Both packages live in this one repo and neither is on its
  respective package index yet (`EDI` isn't on CRAN; `edi_kernels` isn't
  on PyPI). First pass tried adding a second `URL:` entry to
  `EDI/DESCRIPTION` pointing at `python/` plus a `BugReports:` field —
  reverted: `DESCRIPTION` is package metadata parsed by tooling (CRAN,
  `pak`, etc.), and it shouldn't carry an informal pointer to an unrelated
  sibling project just because they happen to share a repo. Settled
  instead on prose-only cross-links: the repo-root `README.md` gained a
  short paragraph pointing at `python/` right after its opening
  paragraph, and `python/README.md`'s opening paragraph links back to
  `EDI/`. Once either package actually publishes to CRAN/PyPI, add that
  registry link too (`DESCRIPTION`'s `URL:` pointing at the *published
  Python package's own listing*, not at a source subdirectory, would be
  legitimate at that point) — but the README-level links stay regardless.
- Keep a `CHANGELOG.md` in `python/` from the first release onward,
  entries keyed to the same version number scheme above — a compiled
  numerical library's users need to know exactly which kernel-behavior or
  ABI changes landed in which release before upgrading a pinned dependency.
  **Status: DONE (2026-08-05).** `python/CHANGELOG.md` created with a
  `[1.0.0] - 2026-08-05` initial-release entry (Keep a Changelog style),
  covering the full binding surface, test coverage, packaging/release
  infrastructure, and the two known result-shape quirks already documented
  elsewhere in this doc.
