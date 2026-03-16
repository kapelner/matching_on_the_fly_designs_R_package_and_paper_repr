#' Univariate Zero-Inflated Poisson Regression Inference for Count Responses
#'
#' @description
#' Fits a zero-inflated Poisson regression for count responses using only the
#' treatment indicator in both the count and zero-inflation components. The
#' reported treatment effect is the treatment coefficient from the conditional
#' Poisson count component, on the log-rate scale.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignCRD$
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
#' infer <- SeqDesignInferenceCountUnivZeroInflatedPoissonRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
SeqDesignInferenceCountUnivZeroInflatedPoissonRegr = R6::R6Class("SeqDesignInferenceCountUnivZeroInflatedPoissonRegr",
	inherit = SeqDesignInferenceCountZeroAugmentedPoissonAbstract,
	private = list(
		za_family = function() poisson(link = "log"),
		za_description = function() "Zero-inflated Poisson regression"
	)
)
