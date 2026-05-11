#' Parametric-Bootstrap-Capable Likelihood Inference
#'
#' Intermediate abstract base for the subset of likelihood-backed inference
#' families that are plausible future targets for parametric null-bootstrap
#' likelihood-ratio calibration.
#'
#' This class intentionally adds only a narrow abstraction layer today:
#'
#' - it keeps `InferenceAsympLik` as the universal likelihood base
#' - it provides a distinct inheritance node for families that may support
#'   likelihood-based parametric bootstrap
#' - it leaves the actual null-simulation machinery abstract / unimplemented
#'
#' Families with highly bespoke partial-likelihood, quadrature, frailty, copula,
#' or custom combined-likelihood geometry can remain direct children of
#' `InferenceAsympLik`.
#'
#' @keywords internal
InferenceParamBootstrap = R6::R6Class("InferenceParamBootstrap",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = list(
		#' @description Placeholder public capability query for future parametric-bootstrap LR support.
		supports_lik_ratio_param_bootstrap = function(){
			FALSE
		}
	)
)
