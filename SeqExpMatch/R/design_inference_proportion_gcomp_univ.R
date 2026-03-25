#' Univariate G-Computation Mean-Difference Inference for Proportion Responses
#'
#' @description
#' Fits a univariate fractional-logit working model for a proportion outcome and
#' then estimates the marginal mean difference by standardizing predicted mean
#' proportions under all-treated and all-control assignments over the empirical
#' covariate distribution.
#'
#' @export
DesignInferencePropUniGCompMeanDiff = R6::R6Class("DesignInferencePropUniGCompMeanDiff",
	inherit = DesignInferencePropGCompAbstract,
	public = list(

		#' @description
		#' Initialize the g-computation mean-difference inference object.
		#' @param seq_des_obj A completed \code{SeqDesign} object with a proportion response.
		#' @param num_cores The number of CPU cores to use for bootstrap and randomization inference.
		#' @param verbose Whether to print progress messages.
		#'
		#' @examples
		#' \dontrun{
		#' seq_des = SeqDesignBernoulli$new(n = 20, response_type = "proportion")
		#' for (i in 1:20) {
		#' 	x_i = data.frame(x1 = rnorm(1), x2 = rnorm(1))
		#' 	w_i = seq_des$add_subject_to_experiment_and_assign(x_i)
		#' 	mu_i = plogis(-0.5 + 0.6 * w_i)
		#' 	seq_des$add_subject_response(i, rbeta(1, shape1 = 8 * mu_i, shape2 = 8 * (1 - mu_i)))
		#' }
		#' seq_des_inf = DesignInferencePropUniGCompMeanDiff$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' }
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
		},

		#' @description
		#' Computes the standardized marginal mean-difference estimate.
		compute_treatment_estimate = function(){
			super$compute_treatment_estimate()
		},

		#' @description
		#' Computes a 1 - \code{alpha} confidence interval for the marginal mean difference.
		#' @param alpha The confidence level in the computed confidence interval is 1 - \code{alpha}.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			super$compute_asymp_confidence_interval(alpha = alpha)
		},

		#' @description
		#' Computes a two-sided p-value for the marginal mean difference.
		#' @param delta The null mean difference to test against. The default is 0.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			super$compute_asymp_two_sided_pval_for_treatment_effect(delta = delta)
		},

		#' @description
		#' Computes a bootstrap two-sided p-value for the marginal mean difference.
		#' @param delta The null mean difference to test against. The default is 0.
		#' @param B Number of bootstrap samples.
		#' @param na.rm Whether to remove non-finite bootstrap replicates.
		compute_bootstrap_two_sided_pval = function(delta = 0, B = 501, na.rm = FALSE){
			super$compute_bootstrap_two_sided_pval(delta = delta, B = B, na.rm = na.rm)
		}
	),

	private = list(
		build_design_matrix = function(){
			cbind(1, private$w)
		}
	)
)
