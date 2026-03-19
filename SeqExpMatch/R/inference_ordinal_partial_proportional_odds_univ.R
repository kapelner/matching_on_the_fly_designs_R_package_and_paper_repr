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
#' seq_des <- SeqDesignCRD$new(n = 8, response_type = "ordinal", verbose = FALSE)
#' for (i in seq_len(8)) {
#'   seq_des$add_subject_to_experiment_and_assign(data.frame(x1 = i))
#' }
#' seq_des$add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- SeqDesignInferenceOrdinalUniPartialProportionalOddsRegr$
#'   new(seq_des, verbose = FALSE)
#' infer
SeqDesignInferenceOrdinalUniPartialProportionalOddsRegr = R6::R6Class(
	"SeqDesignInferenceOrdinalUniPartialProportionalOddsRegr",
	inherit = SeqDesignInferenceOrdinalPartialProportionalOddsAbstract,
	public = list(
		#' @description
		#' Initialize a treatment-only partial proportional-odds inference object.
		#' @param seq_des_obj A completed \code{SeqDesign} object with an ordinal
		#'   response.
		#' @param num_cores Number of CPU cores for bootstrap/randomization helpers.
		#' @param verbose Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(
				seq_des_obj = seq_des_obj,
				nonparallel = character(0),
				num_cores = num_cores,
				verbose = verbose
			)
		}
	),
	private = list(
		ppo_covariate_matrix = function(){
			matrix(0, nrow = private$n, ncol = 0)
		}
	)
)
