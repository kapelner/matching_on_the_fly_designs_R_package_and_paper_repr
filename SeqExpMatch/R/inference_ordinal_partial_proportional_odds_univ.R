#' Partial Proportional-Odds Inference for Ordinal Responses
#'
#' @description
#' Fits a treatment-only partial proportional-odds model for ordinal responses.
#' With no baseline covariates, this reduces to the usual proportional-odds
#' cumulative-logit model and uses the package's fast Rcpp ordinal solver.
#'
#' @export
#' @examples
#' set.seed(1)
#' seq_des <- DesignSeqOneByOneBernoulli$new(n = 8, response_type = "ordinal", verbose = FALSE)
#' for (i in seq_len(8)) {
#'   seq_des$add_subject_to_experiment_and_assign(data.frame(x1 = i))
#' }
#' seq_des$add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- InferenceOrdinalUniPartialProportionalOddsRegr$
#'   new(seq_des, verbose = FALSE)
#' infer
InferenceOrdinalUniPartialProportionalOddsRegr = R6::R6Class(
	"InferenceOrdinalUniPartialProportionalOddsRegr",
	inherit = InferenceOrdinalPartialProportionalOddsAbstract,
	public = list(
	),
	private = list(
		ppo_covariate_matrix = function(){
			matrix(0, nrow = private$n, ncol = 0)
		}
	)
)
