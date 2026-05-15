#' Abstract mixin: Zhang combined randomisation CI for quantile regression
#'
#' Provides \code{compute_rand_confidence_interval()}
#' via Zhang's combined test-inversion method for both Bernoulli (\eqn{m = 0},
#' all subjects in the reservoir) and KK matching-on-the-fly designs
#' (\eqn{m > 0}).
#'
#' @keywords internal
InferenceAbstractQuantileRandCI = R6::R6Class("InferenceAbstractQuantileRandCI",
	lock_objects = FALSE,
	inherit = InferenceKKPassThroughCompound,
	public = c(InferenceMixinQuantileRandCI$public, list()),
	private = c(InferenceMixinQuantileRandCI$private, list())
)
