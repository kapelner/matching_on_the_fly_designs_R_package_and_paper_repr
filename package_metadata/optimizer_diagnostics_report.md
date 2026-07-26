# An Advanced Optimizer Diagnostics Layer

## Bottom Line

The package currently has **no centralized way to know why a fit failed, or
whether a fit that looks like it succeeded actually didn't**. Every failure
signal in use today is an ad hoc, independently-invented heuristic bolted on
after the native optimizer returns, and most of the diagnostic information
needed to do this properly is **already computed inside the native fitters
and then thrown away** rather than surfaced.

That second point is the key finding of this report: a real diagnostics layer
here is not a research project. Most of it is a matter of *returning values
that already exist* rather than deriving new ones. Only two of the eight
failure modes below require genuinely new computation.

## Current State

### What the C++ fitters return today is inconsistent

Spot-checking a representative spread of native fitters:

| Fitter | Returns |
|---|---|
| `fast_logistic_regression_cpp` ([fast_logistic_regression.cpp:316-332](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/src/fast_logistic_regression.cpp)) | `b`, `w`/`XtWX`/`fisher_information`, `score`, `neg_ll`, `converged`, `iterations` |
| `fast_poisson_regression_cpp` (`fast_poisson_regression.cpp:531-539`) | `b`, `ssq_b_j`/`ssq_b_2`, `dispersion`, `mu`, `converged`, `iterations` — **no score or information at all** |
| `fast_logistic_glmm_cpp` (`fast_logistic_glmm.cpp:500-565`) | `params`/`b`/`log_sigma`, `ssq_b_T`, `vcov`, `score`, `information`, `converged`, `neg_loglik` — richest of the sample, but `estimate_only` mode drops score/information/vcov entirely |
| `fast_coxph_regression_cpp` (`fast_coxph_regression.cpp:636-711`) | `coefficients`, `converged`, `neg_ll`, `iterations`, `fisher_information` — the only family that threads `iterations` through *every* return branch |
| `fast_zero_augmented_poisson_cpp` (`:370-390`) | on failure: `List::create(Named("converged")=false)` **only** — no `b`, no reason |

`converged` and the coefficient vector are near-universal. `iterations` is
present in maybe half the families. A gradient/score **norm** at termination
is never returned anywhere (some families return the raw score vector, not
its norm). No condition number or Hessian/information eigenvalue diagnostic
is ever surfaced, by any family.

### `converged` cannot be trusted the way its name suggests

`optimize_likelihood_lbfgs` ([_helper_functions.h:988-1014](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/src/_helper_functions.h)) sets

```cpp
fit.converged = (fit.niter < maxit);
```

Convergence is inferred purely from *not hitting the iteration cap* — nothing
about the terminal gradient norm or line-search health. LBFGSpp internally
uses gradient-based stopping criteria (`epsilon`/`epsilon_rel`) to decide when
to stop *before* the cap, but the gradient norm at that stopping point is
computed inside the library and discarded on the way out. This is the exact
mechanism behind the separation problem already documented in
[firth_penalties_report.md](firth_penalties_report.md) and in the earlier
diagnosis of `glm_fit_helpers.R`'s `SEPARATION_THRESHOLD`: a separated
logistic log-likelihood is asymptotically flat, so the optimizer's own
gradient-norm stopping criterion is satisfied while the coefficients have
walked off to a large-but-finite value. `converged == TRUE` currently means
"the optimizer's internal criterion was satisfied or the iteration cap wasn't
hit" — not "this is a trustworthy estimate."

### The eigendecomposition ingredients for a singularity diagnostic already exist

`symmetric_pseudo_inverse` (`_helper_functions.h:377-391`), used by
`covariance_from_information`, already computes a full eigendecomposition and
a relative-eigenvalue-vs-tolerance test to build a pseudo-inverse when the
information matrix is singular or near-singular. Only the resulting
pseudo-inverse matrix is kept — the eigenvalues themselves, which are exactly
a condition-number / near-singularity diagnostic, are computed and discarded
on every single call.

### There is one real precedent for detect-then-structurally-recover

`fit_with_hardened_qr_column_dropping` ([inference_all_abstract.R:889+](/home/kapelner/workspace/matching_on_the_fly_designs_R_package_and_paper_repr/EDI/R/inference_all_abstract.R), used at 45+ call sites via `InferenceMixinKKGEEShared`/`InferenceMixinKKGLMMShared`) runs a rank-revealing pivoted QR on the design matrix, keeps `qr_X$pivot[seq_len(qr_X$rank)]` plus any required columns, and refits on the reduced matrix. This is a genuine example of the pattern this report proposes generalizing: detect a specific failure mode via a cheap decomposition, then recover structurally rather than just erroring. But it is gated behind a `private$harden` opt-in flag, not automatic, and it only covers rank deficiency.

### Two boundary pathologies are silently swallowed, not diagnosed

GLMM variance-component boundary behavior — a random-effect variance
collapsing toward zero, a well-known and statistically meaningful pathology
distinct from separation — is handled with a **silent penalty**, not a
diagnostic. `log_sigma_penalty(log_sigma, center = 5, scale = 10)`
(`fast_logistic_glmm.cpp:56`) is a soft barrier, and `_glmm_engine.h:132`
hard-clamps: `if (std::abs(log_sigma) > max_abs_log_sigma) return 1e100`.
Neither reports whether the fitted `log_sigma` landed near that boundary; the
raw value is returned for the caller to notice or not.

GLMM quadrature adequacy has no check at all. `n_gh` (Gauss-Hermite node
count) is a hardcoded default (e.g. `= 20` in `fast_logistic_glmm.cpp:363`)
with no refit-with-more-nodes-and-compare pattern anywhere in the codebase.

### A ninth failure mode this session measured directly

The Nelder-Mead benchmark run earlier in this session
([firth_penalties_report.md](firth_penalties_report.md), "Why Not Just Use A
Derivative-Free Optimizer?") found R's own `optim(method = "Nelder-Mead")`
hitting its "degenerate simplex" convergence code (`10`) on one run and
succeeding on an *identical* rerun — same data, same seed — specifically on
the log-det-penalized objective, on GLMM/quadrature paths. R's `optim()`
already returns `$convergence` and `$counts` for free; nothing currently
reads or classifies them.

## Failure Taxonomy

| # | Failure mode | Currently detected? | Detection mechanism today |
|---|---|---|---|
| 1 | Separation / boundary divergence | Partially, in 3-4 call sites | Copy-pasted `max(abs(b)) <= 1e6` magnitude threshold |
| 2 | Iteration-cap non-convergence mislabeled as success | No | `converged = (niter < maxit)` conflates "converged well" with "just under the cap" |
| 3 | Rank deficiency / collinearity | Yes, but opt-in only | `fit_with_hardened_qr_column_dropping`, gated by `private$harden` |
| 4 | Near-singular information matrix | No | Eigenvalues computed in `symmetric_pseudo_inverse` and discarded |
| 5 | GLMM variance-component boundary collapse | No | Silent soft penalty / hard clamp, value not compared to boundary afterward |
| 6 | Quadrature (Gauss-Hermite node count) inadequacy | No | Fixed `n_gh`, no adaptive check |
| 7 | Derivative-free optimizer degeneracy (simplex collapse) | No | R's own `$convergence`/`$counts` unused by any caller |
| 8 | Extreme/non-finite bootstrap or LR-statistic replicates | Partially | Separate copy-pasted `1e6` thresholds in `inference_all_abstract_param_boot.R` and `inference_all_abstract_non_param_boot.R` |

Eight distinct failure modes, five different ad hoc detection mechanisms (three of them the same `1e6` constant copy-pasted independently), and no shared vocabulary for reporting any of them.

## Proposed Architecture

Two layers, mirroring the split already established this session between the
native evaluators and the `InferenceMixinOffOptimumLikelihoodEval` /
`InferenceMixinInformationMatrix` R-level consumers of them.

### Layer 1 (native): surface what's already computed

Extend the `List::create(...)` return of each `fast_*_cpp` fitter with a
standardized `diagnostics` sub-list:

```cpp
List::create(
  ...,
  Named("diagnostics") = List::create(
    Named("iterations")            = niter,
    Named("hit_iteration_cap")     = (niter >= maxit),
    Named("gradient_norm")         = grad_norm,       // from the optimizer's own stopping check
    Named("min_information_eigenvalue") = min_eig,    // from symmetric_pseudo_inverse's existing eigendecomposition
    Named("boundary_hit")          = boundary_flag     // e.g. |log_sigma| near max_abs_log_sigma
  )
)
```

Three of these four fields (`gradient_norm`, `min_information_eigenvalue`,
`boundary_hit`) require **no new computation** — they are values the
optimizer or `symmetric_pseudo_inverse` already has in hand at the moment it
currently discards them. This is a mechanical, low-risk, purely-additive
change: touches many files, but each touch is "keep a value that's already
being computed instead of dropping it."

### Layer 2 (R): interpret, classify, and centralize

A new mixin, `InferenceMixinSolverDiagnostics`, following the existing
Pattern-1 shape:

- Centralizes the magnitude-threshold constant (today `1e6`, copy-pasted five
  times with no guarantee the copies stay in sync) into one shared function.
- Classifies the raw native `diagnostics` fields into a standard taxonomy:
  `severity` (`ok` / `warning` / `failure`) × `category` (`separation`,
  `rank_deficiency`, `near_singular_information`, `variance_boundary`,
  `iteration_cap`, `simplex_degenerate`, `quadrature_inadequate`,
  `extreme_replicate`).
- Feeds classified failures into `cache_nonestimable_estimate(reason = ...)`
  with standardized reason strings instead of the current free-form,
  per-family strings (e.g. `"param_bootstrap_estimate_raw_estimate_nonfinite"`
  vs. whatever a different family happens to have written).
- Exposes a `get_last_fit_diagnostics()` accessor for debugging and for
  simulation-study introspection (e.g. "what fraction of bootstrap replicates
  hit `near_singular_information` on this path").
- For any family that calls `stats::optim()` directly (Nelder-Mead or
  otherwise), wraps `$convergence`/`$counts` through the same classifier, so
  R's own base convergence codes (`10` = degenerate simplex, `1` = iteration
  limit, etc.) land in the same taxonomy as the native-fitter diagnostics
  rather than being silently ignored.

This mixin would **not** be added to the standard `InferenceAsympLik`
composition automatically — like `InferenceMixinOffOptimumLikelihoodEval`, it
should be opt-in per family, since wiring the native side (Layer 1) has to
happen family-by-family regardless.

## What's Cheap vs. What's New Work

### Cheap: thread out already-computed values

1. **Gradient norm at termination** — LBFGSpp already computes it internally to decide when to stop; return it.
2. **Information-matrix condition number / minimum eigenvalue** — free *only* on the path that already calls `symmetric_pseudo_inverse`, i.e. when the normal fast path (`eigen_compute_single_entry_on_diagonal_of_inverse_matrix_cpp`, a targeted single-entry computation) or a plain `solve()` has already failed. That path never catches near-singular-but-not-yet-erroring matrices — the more insidious case. Surfacing the diagnostic universally, on every fit, would mean adding a genuine new full eigendecomposition (`O(p^3)`) to a currently cheaper path. See [Performance Implications](#performance-implications).
3. **GLMM variance-boundary flag** — `log_sigma_penalty`/the hard clamp already know the boundary; compare the final `log_sigma` to it and return a flag.
4. **Centralizing the `1e6` threshold** into one function — pure refactor, removes the risk of the copies silently drifting apart (e.g. `inference_mixin_kk_gee_shared.R`'s `max_abs_reasonable_coef` is a separate, independently-configurable value that may not equal `1e6` in practice).
5. **Classifying `optim()`'s own `$convergence`/`$counts`** for any Nelder-Mead-based path — the values already exist in R, this is a lookup table.

None of these require new derivations, new math, or touching the "hard"
likelihood families from the Firth audit. They fix the trust problem in
`converged`, which is the most consequential gap: it's the field every
existing caller already reads and already (wrongly) trusts.

### More involved: genuinely new work

6. **Iteration-cap vs. true convergence, properly separated.** Once (1) is threaded out, redefine `converged` as `gradient_norm < tol` and report `hit_iteration_cap` as an independent field — a real behavior change to what `converged` means, needing a compatibility pass over every caller that currently branches on it.
7. **Quadrature adequacy.** Requires an actual refit at a higher node count and a comparison of the resulting estimate/log-likelihood delta — roughly doubles the cost of the affected fit, so this should be an explicit opt-in diagnostic mode, not a default-on check.
8. **Rank-deficiency visibility without hardening.** Currently, if `private$harden` is off, a near-collinear design is neither corrected nor flagged. Emitting a diagnostic flag whenever the rank-revealing QR *would* have dropped columns — even when hardening isn't applied — closes a real blind spot cheaply, since the QR decomposition itself is not expensive; only the actual re-fit-on-dropped-columns path is opt-in-worthy for cost reasons, not the detection.

## Performance Implications

This is not uniformly free, and this codebase has already been burned by
assuming "small per-call overhead" doesn't matter (see the mutable-`std::vector`
GLMM objective regression on record, an 18x slowdown from exactly this shape
of mistake). Sorting the proposed additions by actual runtime cost:

**Genuinely free** — one extra scalar into a `List::create(...)` that already
happens once per fit, after optimization is already done: gradient norm at
termination, the GLMM variance-boundary flag, centralizing the `1e6`
threshold into one function, and classifying `optim()`'s own
`$convergence`/`$counts`. None of these touch the per-iteration hot loop.

**Conditionally free, not universally free** — the information-matrix
minimum eigenvalue / condition number. It is free exactly when
`symmetric_pseudo_inverse` is already being called, i.e. on the fallback path
after the fast single-entry inverse or a plain `solve()` has already failed.
That path does not cover the more useful case: a matrix that is near-singular
but not singular enough for `solve()` to error. Making the diagnostic
available on *every* fit, not just the already-failing ones, means adding a
new full eigendecomposition (`O(p^3)`) to a path that currently doesn't pay
for one. This should be treated as an explicit cost/benefit decision, not
bundled into Phase 1 as if it were free.

**Real, acknowledged costs** — quadrature adequacy (~2x the affected fit,
opt-in only) and rank-deficiency-visibility when `harden` is off: running
`qr(X)` there adds a decomposition that today simply does not happen, not a
reuse of an existing one.

**The larger risk is where Layer 2's classification logic runs, not any
single field.** Mixin dispatch happens once per fit — negligible against a
GLMM/quadrature fit, but CI inversion, parametric-bootstrap, and permutation
loops call `fit_null`/refit hundreds to thousands of times per confidence
interval or p-value, and on the cheapest paths a single native fit is already
sub-millisecond (the OLS/logistic native-fit benchmarks earlier in this
session ran at roughly 0.0000-0.001s). R6 dispatch plus list construction on
every one of those calls is small in isolation, but "small overhead times a
hot loop" is exactly the failure shape this codebase has hit before.
`InferenceMixinSolverDiagnostics` should be benchmarked specifically against
the fastest repeated-refit paths (plain logistic/OLS under CI inversion, not
just the slow GLMM paths where it clearly doesn't matter) before being
enabled unconditionally, and likely needs a way to be skipped inside
CI-inversion/bootstrap inner loops by default, with diagnostics collected
only on the final reported fit rather than every intermediate one.

## Recommendation

### Phase 1

Thread out the two genuinely free already-computed diagnostic values
(gradient norm, GLMM variance-boundary flag) through the native
`List::create(...)` returns, starting with the families audited above.
Centralize the `1e6` separation threshold into one shared function and point
the existing five call sites at it. Decide the information-matrix
minimum-eigenvalue question separately (see
[Performance Implications](#performance-implications)) rather than bundling
it in as if it were free — either scope it to the already-failing fallback
path only, or explicitly accept the new `O(p^3)` cost for universal coverage.

### Phase 2

Build `InferenceMixinSolverDiagnostics` to classify Layer-1 output into the
standard taxonomy and route it through `cache_nonestimable_estimate()` with
standardized reason strings. Wrap any direct `stats::optim()` caller (present
or future — this is exactly what a Firth-on-the-"No"-tier implementation
would need) through the same classifier. Benchmark dispatch overhead against
the fastest repeated-refit paths before enabling it unconditionally inside
CI-inversion/bootstrap/permutation inner loops; default to collecting
diagnostics only on the final reported fit in those loops, not every
intermediate refit.

### Phase 3

Redefine `converged` to mean gradient-norm-based convergence, with
`hit_iteration_cap` reported separately; audit and update callers that branch
on `converged` today.

### Phase 4

Quadrature-adequacy and always-on rank-deficiency-visibility checks — treat
as separate, lower-priority projects given the added cost (quadrature) or the
narrower blast radius (rank-deficiency visibility only matters where hardening is currently off).

## Why This Matters Beyond Its Own Scope

This isn't an isolated cleanup. It's a direct prerequisite for two things
already discussed this session:

- **Firth vs. Cox-Snell MC dispatch** ([firth_penalties_report.md](firth_penalties_report.md)): cleanly gating which bias-correction path to use requires a *reliable* separation signal. The current heuristic — a magnitude threshold that varies by file and is absent from most families — cannot support that gate today.
- **Any derivative-free optimization of a Firth-penalized objective**: the degenerate-simplex flip-flopping measured empirically this session is exactly failure mode #7 in the taxonomy above. A diagnostics layer that classifies `optim()`'s own convergence code would have caught it automatically instead of requiring a manual rerun-and-compare to notice.

## Final Assessment

- **Centralized diagnostics: worth building, and cheaper than it looks.** The core trust problem — `converged` doesn't mean what callers assume — is fixable by returning values the optimizers already compute, not by deriving anything new.
- **The taxonomy is the actual deliverable.** Eight distinct failure modes currently share, at best, a copy-pasted magnitude constant and no common vocabulary. A shared classification is what turns scattered heuristics into something a caller (or a future Firth implementation) can actually branch on.
- **Only two of eight failure modes need new computation** (quadrature adequacy, and properly separating iteration-cap-hit from true convergence). Everything else is surfacing values the C++ layer already has in hand and currently throws away.
