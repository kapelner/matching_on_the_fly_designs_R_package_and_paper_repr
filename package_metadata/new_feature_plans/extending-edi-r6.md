# Extending EDI with Custom R6 Inference Classes

EDI is implemented with R6 classes. Advanced users can define their own R6
classes outside the package and reuse EDI's design storage, response handling,
randomization, bootstrap, and summary methods.

The custom-extension base classes are intentionally internal while the extension
contract is experimental. Retrieve them with `getFromNamespace()`:

```r
library(EDI)
library(R6)

InferenceCustomAsymp <- getFromNamespace("InferenceCustomAsymp", "EDI")
```

## Inference Contract

A custom asymptotic inference class should inherit from `InferenceCustomAsymp`
and implement a public `fit(estimate_only = FALSE)` method. `fit()` returns a
named list with:

- `estimate`: required numeric scalar treatment-effect estimate.
- `se`: optional numeric scalar standard error.
- `df`: optional degrees of freedom. Use `NA_real_` for z inference.
- `model`: optional fitted model object retained by `get_mod()`.
- `nonestimable_reason`: optional character scalar used when the estimate or
  standard error is unavailable.

Use public accessors instead of private fields:

- `get_response()`
- `get_treatment()`
- `get_covariates()`
- `get_analysis_data()`
- `get_design_object()`
- `get_response_type()`

## Example

```r
InferenceMedianDiff <- R6Class(
  "InferenceMedianDiff",
  inherit = InferenceCustomAsymp,
  public = list(
    fit = function(estimate_only = FALSE) {
      dat <- self$get_analysis_data()
      y_t <- dat$y[dat$w == 1]
      y_c <- dat$y[dat$w == 0]

      est <- stats::median(y_t) - stats::median(y_c)
      if (estimate_only) {
        return(list(estimate = est))
      }

      list(
        estimate = est,
        se = sqrt(stats::var(y_t) / length(y_t) + stats::var(y_c) / length(y_c)),
        df = length(y_t) + length(y_c) - 2,
        model = NULL
      )
    }
  )
)

des <- DesignFixedBernoulli$new(n = 20, response_type = "continuous")
des$add_all_subjects_to_experiment(data.frame(x = seq_len(20)))
des$overwrite_all_subject_assignments(rep(c(0, 1), each = 10))
des$add_all_subject_responses(rnorm(20))

inf <- InferenceMedianDiff$new(des)
inf$compute_estimate()
inf$compute_asymp_two_sided_pval()
inf$compute_asymp_confidence_interval()
inf$compute_bootstrap_two_sided_pval(B = 501)
```

## Custom Designs

EDI also provides internal bases for user-defined designs:

```r
DesignCustomFixed <- getFromNamespace("DesignCustomFixed", "EDI")
DesignCustomSequential <- getFromNamespace("DesignCustomSequential", "EDI")
```

For fixed designs, implement `draw_assignments(r)` and return an `n x r` 0/1
assignment matrix. For sequential designs, implement `assignment_rule()` and
return a scalar 0/1 assignment for the current subject. EDI handles subject
storage, response recording, and validation.

## Implementation TODOs

### Doc Accuracy

- [ ] TODO-1: Fix the wrong class name in the "Custom Designs" section. The doc
  says `DesignCustomFixed <- getFromNamespace("DesignCustomFixed", "EDI")`, but
  the actual internal base class is named `DesignFixedCustom` (defined at
  `EDI/R/design_custom_extensions.R:9`, `R6::R6Class("DesignFixedCustom", ...)`).
  `getFromNamespace("DesignCustomFixed", "EDI")` throws today because no object
  by that name exists. `EDI/tests/testthat/test-custom-extension-contract.R:7,77`
  already uses the correct `DesignFixedCustom` name, confirming the doc — not
  the code — has the name swapped.

### Missing Coverage

- [ ] TODO-2: Document `InferenceCustomRand` and `InferenceCustomBoot`
  (`EDI/R/inference_custom_extensions.R:125` and `:169`), the two other
  user-facing custom-inference base classes. The doc's "Inference Contract"
  section only covers `InferenceCustomAsymp`, but the opening paragraph
  promises "randomization, bootstrap, and summary methods" reuse, and the
  source's own roxygen comments for all three classes explicitly cross-link
  the other two (e.g. `EDI/R/inference_custom_extensions.R:41-44`). Add a
  subsection per class covering: `InferenceCustomRand` (inherits `Inference`,
  `components = "RandomizationTest"`, `likelihood_tier = "none"`, only a
  `fit()`/`compute_estimate()` contract — no `se`/`df` needed) and
  `InferenceCustomBoot` (inherits `InferenceJackknife`, same minimal
  `fit()`/`compute_estimate()` contract, gains jackknife/bootstrap machinery
  from the parent).

### Test Coverage

- [ ] TODO-3: Add a behavioral subclass test for `InferenceCustomRand` and one
  for `InferenceCustomBoot`, mirroring the working end-to-end example already
  present for `InferenceCustomAsymp` in
  `EDI/tests/testthat/test-custom-extension-contract.R:17-75`. Currently those
  two classes are only checked for existence/retrievability
  (`test-custom-extension-contract.R:12-13`) — there is no test that actually
  subclasses either one, implements `fit()`, and calls `compute_estimate()`
  end-to-end, so a future refactor of `InferenceRand`/`InferenceJackknife`
  could silently break either extension contract without any test failing.
