# Should `InferenceParamBootstrap` Grow an `estimate_type` Field?

## Context

Today every inference class has exactly one point estimate: `compute_estimate()`
is unambiguous, because there is only one way to estimate the parameter of
interest for a given model. That stopped being strictly true once
`InferenceParamBootstrap` grew a second way to estimate the same
parameter — the parametric-bootstrap bias-corrected estimate implemented in
this session as `compute_param_bootstrap_estimate()` (with its `_pval()` and
`_confidence_interval()` companions), living in
`R/inference_mixin_param_bootstrap_estimate.R` and spliced privately into
`R/inference_all_abstract_param_boot.R`. That method returns a *different*
number from `compute_estimate()` for the same fitted object: `2*theta_hat -
mean(theta_hat_star)` instead of the raw MLE.

This raised the natural design question: if the package is going to grow more
than one kind of estimate for a class (raw MLE, parametric-bootstrap
bias-corrected, and potentially others — Cox-Snell analytic correction,
jackknife bias correction, Firth's penalized-likelihood estimate), should
there be a single `estimate_type` switch analogous to the existing
`testing_type` switch, rather than one differently-named method per
correction?

## The existing precedent: `testing_type`

`InferenceAsympLik` already solves an analogous problem for *testing*: a
fitted model can be probed with a Wald, score, likelihood-ratio, or gradient
test, and the class exposes this as a single stateful field rather than four
separately-named public methods for every downstream consumer.

`R/inference_all_abstract_asymp_lik.R`:

```r
set_testing_type = function(testing_type = c("wald", "score", "gradient", "lik_ratio", "lik_ratio_bartlett_approx", "lik_ratio_bartlett_exact")){
	testing_type = private$normalize_testing_type(testing_type)
	supported = private$get_supported_testing_types_with_bartlett()
	if (!testing_type %in% supported) {
		stop(
			class(self)[1], " does not support testing_type = \"", testing_type,
			"\". Supported values are: ", paste(supported, collapse = ", "),
			call. = FALSE
		)
	}
	private$testing_type = testing_type
	invisible(self)
},
get_testing_type = function(){
	private$testing_type
},
```

Three properties of this design are worth calling out because they're exactly
what a candidate `estimate_type` would need to replicate:

1. **A validated setter.** `set_testing_type()` normalizes the input and
   checks it against `get_supported_testing_types_with_bartlett()` before
   accepting it, erroring immediately (not silently falling back) for
   unsupported values — the same defensive pattern used throughout this
   codebase's `compute_param_bootstrap_*` gating (`supports_param_bootstrap_estimate()`
   stopping with a descriptive message rather than returning `NA` silently).
2. **A per-class supported-values query**, so downstream code (and the
   `path_audits.html` generator) can introspect what's actually available
   for a given concrete class rather than hardcoding assumptions.
3. **A cache key that embeds the mode.** Testing results are memoized per
   `(testing_type, delta)` pair via `likelihood_test_delta_key()`
   (`R/inference_all_abstract.R`), so switching `testing_type` on a live
   object doesn't silently return stale results computed under a different
   mode.

Crucially, `testing_type` **never re-fits the model.** It only selects which
post-fit statistic is computed against the already-cached fit
(`private$cached_mod`). That constraint is what makes the pattern safe to
expose as mutable object state.

## Which candidate estimate corrections fit that same constraint?

| Candidate | Changes the fit? | Fits the `testing_type` pattern? |
|---|---|---|
| Raw MLE (`compute_estimate()`) | — (baseline) | yes |
| Parametric-bootstrap bias correction (`compute_param_bootstrap_estimate()`, implemented) | No — refits *replicates*, but the object's own `cached_mod` is untouched | yes |
| Cox-Snell analytic correction (not yet implemented; see `package_metadata/bias_correction_cox_snell.md`) | No — additive transform of the same theta_hat, using the same Fisher information already computed for `vcov()` | yes |
| Jackknife bias correction (not yet implemented) | No — built on leave-one-out refits of the same model spec, same shape as the parametric-bootstrap case | yes |
| Firth's penalized-likelihood estimate | **Yes** — changes the estimating equation itself (score + penalty = 0), so it changes what `private$fit_combined()`/the optimizer actually solves | **no** |

The dividing line is exactly the one `testing_type` already draws for tests:
does this option merely pick a different *summary of* the existing fit, or
does it require a materially different *fit*? Cox-Snell, jackknife
bias-correction, and the already-implemented parametric-bootstrap correction
are all pure post-hoc transforms of `private$cached_mod` — cheap to switch
between on a live object, safe to cache per-mode. Firth is not: giving it a
value in the same enum as `"cox_snell"` would let a user flip
`estimate_type` between `"cox_snell"` and `"firth"` on an object that was
never actually re-fit under the Firth penalty, silently returning a wrong
number rather than erroring.

## Recommendation

Add `estimate_type` mirroring `testing_type`'s architecture — `get`/`set`,
a `get_supported_estimate_types()` query, and a cache key generalizing
`likelihood_test_delta_key()` (e.g. `estimate_type_key()`) — but scope its
supported values to **pure post-fit transforms of one underlying fit**:

- `"raw"` (delegates to `compute_estimate()`)
- `"param_bootstrap"` (delegates to the already-implemented
  `compute_param_bootstrap_estimate()`)
- `"cox_snell"` (once implemented)
- `"jackknife_bc"` (once implemented)

Keep Firth off this axis entirely. It belongs as a distinct fitting method —
either its own concrete class, or an orthogonal flag that changes what
`fit_combined()` solves — not a value of `estimate_type`. Conflating "which
transform to apply to this fit" with "fit a different model" is exactly the
kind of architectural blur `testing_type`'s design has so far avoided.

## Current implementation status (as of this session)

`estimate_type` does **not** exist yet. What exists today is:

- `compute_param_bootstrap_estimate()` / `compute_param_bootstrap_pval()` /
  `compute_param_bootstrap_confidence_interval()` as **dedicated,
  independently-named public methods** on `InferenceParamBootstrap`
  (`R/inference_all_abstract_param_boot.R`), not dispatched through any
  `estimate_type`-style switch.
- Gating via `private$supports_param_bootstrap_estimate()`, which defaults to
  `isTRUE(private$supports_lik_ratio_param_bootstrap())` — the same boolean
  already tracked per-class as the `pboot` field in
  `package_tests/path_audits_source.R`.
- The private implementation (`run_param_bootstrap_estimate_batch()`,
  `run_param_bootstrap_estimate_replicates()`,
  `compute_param_bootstrap_estimate_impl()`) lives in a dedicated mixin,
  `R/inference_mixin_param_bootstrap_estimate.R`, registered in
  `R/mixin_contracts.R`'s `EDI_MIXIN_CONTRACTS`/`EDI_MIXIN_COMPOSITIONS`.

This means the package currently has **two live, uncoordinated conventions**
for "more than one way to get a number out of a fitted model": `testing_type`
(a stateful switch) and the ad hoc `compute_param_bootstrap_*` naming
(dedicated methods, no switch). If a second post-fit estimate correction
(Cox-Snell or jackknife bias-correction) is implemented before `estimate_type`
lands, it will most likely be added the same ad hoc way — as
`compute_cox_snell_estimate()` or similar — and retrofitting a unifying
`estimate_type` field afterward would mean deciding whether the dedicated
method names stay as thin wrappers around the new switch, or get deprecated
in favor of it. Better to make that call before a third method along the
same lines gets added.

## Open questions if this is pursued

1. **Does `compute_estimate()` itself dispatch on `estimate_type`, or does
   it stay hardcoded to the raw MLE, with a new `compute_estimate_typed()`
   (or similar) doing the dispatch?** Changing `compute_estimate()`'s
   behavior based on mutable object state is a bigger compatibility surface
   than adding a new method, since `compute_estimate()` is called
   unconditionally throughout the existing test harness
   (`package_tests/comprehensive_tests.R`) and downstream consumers with no
   awareness that its return value could now depend on prior state mutation.
2. **Should the dedicated `compute_param_bootstrap_estimate()` /
   `_pval()` / `_confidence_interval()` methods be kept as-is (thin,
   independently documented, independently testable) even after
   `estimate_type` exists**, with the switch calling into them rather than
   replacing them? This session's roxygen2 investigation (documented in the
   commit history of `R/inference_mixin_param_bootstrap_estimate.R`) is a
   concrete argument for keeping named public methods rather than routing
   everything through generic dispatch: named methods get full per-argument
   documentation for free, while a single dispatch method's behavior varies
   by an opaque string field.
3. **Does the CI/pval-inversion machinery need an `estimate_type`-aware
   cache key**, analogous to how `likelihood_test_delta_key()` embeds
   `testing_type`? If `compute_param_bootstrap_confidence_interval()` and a
   future `compute_cox_snell_confidence_interval()` both get folded under
   one dispatch method, their diagnostics (`get_last_param_bootstrap_estimate_diagnostics()`
   today) would need to be keyed by which correction produced them, not
   overwritten by whichever ran last.
