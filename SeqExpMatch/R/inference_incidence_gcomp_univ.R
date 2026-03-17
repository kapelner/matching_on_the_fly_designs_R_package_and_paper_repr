#' Univariate G-Computation Risk-Difference Inference for Binary Responses
#'
#' @description
#' Fits a univariate logistic working model for an incidence outcome and then
#' estimates the marginal risk difference by standardizing predicted risks under
#' all-treated and all-control assignments over the empirical covariate
#' distribution.
#'
#' @export
SeqDesignInferenceIncidUnivGCompRiskDiff = R6::R6Class("SeqDesignInferenceIncidUnivGCompRiskDiff",
	inherit = SeqDesignInferenceIncidGCompAbstract,
	public = list(

		#' @description
		#' Initialize the g-computation RD inference object.
		#' @param seq_des_obj A completed \code{SeqDesign} object with an incidence response.
		#' @param num_cores The number of CPU cores to use for bootstrap and randomization inference.
		#' @param verbose Whether to print progress messages.
		#'
		#' @examples
		#' \dontrun{
		#' seq_des = SeqDesignCRD$new(n = 20, response_type = "incidence")
		#' for (i in 1:20) {
		#' 	x_i = data.frame(x1 = rnorm(1), x2 = rnorm(1))
		#' 	w_i = seq_des$add_subject_to_experiment_and_assign(x_i)
		#' 	p_i = plogis(-0.5 + 0.6 * w_i)
		#' 	seq_des$add_subject_response(i, rbinom(1, 1, p_i))
		#' }
		#' seq_des_inf = SeqDesignInferenceIncidUnivGCompRiskDiff$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' }
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		},

		#' @description
		#' Computes the standardized marginal risk-difference estimate.
		compute_treatment_estimate = function(){
			super$compute_treatment_estimate()
		},

		#' @description
		#' Computes a 1 - \code{alpha} confidence interval for the marginal risk difference.
		#' @param alpha The confidence level in the computed confidence interval is 1 - \code{alpha}.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			super$compute_asymp_confidence_interval(alpha = alpha)
		},

		#' @description
		#' Computes a two-sided p-value for the marginal risk difference.
		#' @param delta The null risk difference to test against. The default is 0.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = NULL){
			super$compute_asymp_two_sided_pval_for_treatment_effect(delta = delta)
		},

		#' @description
		#' Computes a bootstrap two-sided p-value for the marginal risk difference.
		#' @param delta The null risk difference to test against. The default is 0.
		#' @param B Number of bootstrap samples.
		#' @param na.rm Whether to remove non-finite bootstrap replicates.
		compute_bootstrap_two_sided_pval = function(delta = NULL, B = 501, na.rm = FALSE){
			super$compute_bootstrap_two_sided_pval(delta = delta, B = B, na.rm = na.rm)
		}
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
#' seq_des <- SeqDesignCRD$
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
#' infer <- SeqDesignInferenceIncidUnivGCompRiskRatio$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
SeqDesignInferenceIncidUnivGCompRiskRatio = R6::R6Class("SeqDesignInferenceIncidUnivGCompRiskRatio",
	inherit = SeqDesignInferenceIncidGCompAbstract,
	public = list(

		#' @description
		#' Initialize the g-computation RR inference object.
		#' @param seq_des_obj A completed \code{SeqDesign} object with an incidence response.
		#' @param num_cores The number of CPU cores to use for bootstrap and randomization inference.
		#' @param verbose Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		},

		#' @description
		#' Computes the standardized marginal risk-ratio estimate.
		compute_treatment_estimate = function(){
			super$compute_treatment_estimate()
		},

		#' @description
		#' Computes a 1 - \code{alpha} confidence interval for the marginal risk ratio.
		#' @param alpha The confidence level in the computed confidence interval is 1 - \code{alpha}.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			super$compute_asymp_confidence_interval(alpha = alpha)
		},

		#' @description
		#' Computes a two-sided p-value for the marginal risk ratio.
		#' @param delta The null risk ratio to test against. The default is 1.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = NULL){
			super$compute_asymp_two_sided_pval_for_treatment_effect(delta = delta)
		},

		#' @description
		#' Computes a bootstrap two-sided p-value for the marginal risk ratio.
		#' @param delta The null risk ratio to test against. The default is 1.
		#' @param B Number of bootstrap samples.
		#' @param na.rm Whether to remove non-finite bootstrap replicates.
		compute_bootstrap_two_sided_pval = function(delta = NULL, B = 501, na.rm = FALSE){
			super$compute_bootstrap_two_sided_pval(delta = delta, B = B, na.rm = na.rm)
		}
	),

	private = list(
		build_design_matrix = function(){
			cbind(1, private$w)
		},

		get_estimand_type = function() "RR"
	)
)
