#' Count Composite Likelihood Inference Base
#'
#' Shared branch for count models that use a companion likelihood for testing
#' while the reported estimator may be robust or quasi-likelihood based.
#'
#' @keywords internal
InferenceCountCompositeLikelihood = R6::R6Class("InferenceCountCompositeLikelihood",
	lock_objects = FALSE,
	inherit = InferenceAsympLikStdModCache,
	public = list()
)
