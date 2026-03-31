#' Multivariate G-Computation Risk-Difference Inference for Binary Responses
#'
#' Fits a logistic working model for an incidence outcome using treatment and all
#' recorded covariates, then estimates the marginal risk difference by
#' standardizing predicted risks under all-treated and all-control assignments
#' over the empirical covariate distribution.
#'
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
#'   response_type = "incidence",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(c(0, 1, 0, 1, 0, 1, 1, 0))
#' infer <- InferenceIncidMultiGCompRiskDiff$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceIncidMultiGCompRiskDiff = R6::R6Class("InferenceIncidMultiGCompRiskDiff",
	lock_objects = FALSE,
	inherit = InferenceIncidUnivGCompRiskDiff,
	public = list(

	),

	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)

#' @export
InferenceIncidMultiGCompRiskRatio = R6::R6Class("InferenceIncidMultiGCompRiskRatio",
	lock_objects = FALSE,
	inherit = InferenceIncidUnivGCompRiskRatio,
	public = list(

	),

	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		}
	)
)
