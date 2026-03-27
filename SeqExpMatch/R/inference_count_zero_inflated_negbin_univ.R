#' Univariate Zero-Inflated Negative Binomial Regression Inference for Count Responses
#'
#' @description
#' Fits a zero-inflated negative binomial regression for count responses using only
#' the treatment indicator in both the count and zero-inflation components. The
#' reported treatment effect is the treatment coefficient from the conditional
#' negative binomial count component, on the log-rate scale.
#'
#' @inherit InferenceRand methods
#' @inherit InferenceBoot methods
#' @inherit InferenceAsymp methods
#' @inherit InferenceRandCI methods
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- DesignSeqOneByOneBernoulli$
#'   new(
#'   n = nrow(x_dat),
#'   response_type = "count",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(c(0, 1, 1, 2, 2, 3, 3, 4))
#' infer <- InferenceCountUnivZeroInflatedNegBinRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceCountUnivZeroInflatedNegBinRegr = R6::R6Class("InferenceCountUnivZeroInflatedNegBinRegr",
	inherit = InferenceCountZeroAugmentedPoissonAbstract,
	private = list(
		za_family = function() glmmTMB::nbinom2(link = "log"),
		za_description = function() "Zero-inflated negative binomial regression"
	)
)
