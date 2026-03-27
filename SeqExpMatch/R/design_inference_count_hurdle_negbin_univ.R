#' Univariate Hurdle Negative Binomial Regression Inference for Count Responses
#'
#' @description
#' Fits a hurdle negative binomial regression for count responses using only the
#' treatment indicator in both the hurdle and count components. The hurdle
#' indicator model is fit on all subjects, and the zero-truncated negative
#' binomial count component is optimized in C++ on the positive-count subjects.
#' The reported treatment effect is the treatment coefficient from the conditional
#' count component, on the log-rate scale.
#'
#' @inherit DesignInferenceRand methods
#' @inherit DesignInferenceBoot methods
#' @inherit DesignInferenceAsymp methods
#' @inherit DesignInferenceRandCI methods
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignBernoulli$
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
#' infer <- DesignInferenceCountUnivHurdleNegBinRegr$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
DesignInferenceCountUnivHurdleNegBinRegr = R6::R6Class("DesignInferenceCountUnivHurdleNegBinRegr",
	inherit = DesignInferenceCountHurdleNegBinAbstract,
	private = list(
		hurdle_description = function() "Hurdle negative binomial regression"
	)
)
