#' G-Computation Risk-Difference Inference for Binary Responses
#'
#' Fits a logistic working model for an incidence outcome using treatment and,
#' optionally, all recorded covariates, then estimates the marginal risk
#' difference by standardizing predicted risks under all-treated and all-control
#' assignments over the empirical covariate distribution.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidGCompRiskDiff$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceIncidGCompRiskDiff = R6::R6Class("InferenceIncidGCompRiskDiff",
	lock_objects = FALSE,
	inherit = InferenceIncidGCompAbstract,
	public = list(
	),

	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		},

		get_estimand_type = function() "RD"
	)
)

#' G-Computation Risk-Ratio Inference for Binary Responses
#'
#' Fits a logistic working model for an incidence outcome using treatment and,
#' optionally, all recorded covariates, then estimates the marginal risk ratio
#' by standardizing predicted risks under all-treated and all-control assignments
#' over the empirical covariate distribution.
#'
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneBernoulli$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidGCompRiskRatio$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceIncidGCompRiskRatio = R6::R6Class("InferenceIncidGCompRiskRatio",
	lock_objects = FALSE,
	inherit = InferenceIncidGCompAbstract,
	public = list(
	),

	private = list(
		build_design_matrix = function(){
			private$create_design_matrix()
		},

		get_estimand_type = function() "RR"
	)
)
