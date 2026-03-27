#' Univariate G-Computation Risk-Difference Inference for Binary Responses
#'
#' @description
#' Fits a univariate logistic working model for an incidence outcome and then
#' estimates the marginal risk difference by standardizing predicted risks under
#' all-treated and all-control assignments over the empirical covariate
#' distribution.
#'
#' @export
InferenceIncidUnivGCompRiskDiff = R6::R6Class("InferenceIncidUnivGCompRiskDiff",
	inherit = InferenceIncidGCompAbstract,
	public = list(





	),

	private = list(
		build_design_matrix = function(){
			cbind(1, private$w)
		},

		get_estimand_type = function() "RD"
	)
)

#' Univariate G-Computation Risk-Ratio Inference for Binary Responses
#'
#' @description
#' Fits a univariate logistic working model for an incidence outcome and then
#' estimates the marginal risk ratio by standardizing predicted risks under
#' all-treated and all-control assignments over the empirical covariate
#' distribution.
#'
#' @details
#' The point estimate is returned on the risk-ratio scale. Confidence intervals
#' and p-values use the delta method on the log-risk-ratio scale and then map
#' back to the risk-ratio scale.
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
#' infer <- InferenceIncidUnivGCompRiskRatio$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceIncidUnivGCompRiskRatio = R6::R6Class("InferenceIncidUnivGCompRiskRatio",
	inherit = InferenceIncidGCompAbstract,
	public = list(





	),

	private = list(
		build_design_matrix = function(){
			cbind(1, private$w)
		},

		get_estimand_type = function() "RR"
	)
)
