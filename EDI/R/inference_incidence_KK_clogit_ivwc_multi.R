#' Multivariate Conditional Logistic Regression Inference for KK Designs with Binary Response
#'
#' Fits a compound estimator for KK matching-on-the-fly designs with binary (incidence)
#' responses using the treatment indicator and all recorded covariates. For matched
#' pairs, a conditional logistic regression model is used (via the internal
#' \code{clogit_helper}, which processes discordant pairs and fits logistic regression
#' on within-pair signed differences). For reservoir subjects, a standard logistic
#' regression is used. The two estimates are combined via a variance-weighted linear
#' combination, analogous to \code{InferenceContinMultOLSKK}. If clogit fails
#' (e.g. no discordant pairs), the estimator falls back to logistic regression on the
#' reservoir only.
#'
#' @export
InferenceIncidMultiKKClogitIVWC = R6::R6Class("InferenceIncidMultiKKClogitIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKClogitIVWC,
	public = list(




	),
	private = list(
		include_covariates = function() TRUE
	)
)
