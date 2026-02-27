#' Simple Mean Difference Inference based on Maximum Likelihood  
#'
#' @description
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring) sequential experimental design estimation and test object after the sequential design is completed.
#' 
#'
#' @examples
#' \dontrun{
#' # Setup design and add responses
#' seq_des = SeqDesignCRD$new(n = 6, response_type = "survival")
#' X = data.table::data.table(MASS::biopsy[1:6, 2:10])
#' for (i in 1:6) seq_des$add_subject_to_experiment_and_assign(X[i, ])
#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
#' 
#' # Initialize inference object
#' seq_des_inf = SeqDesignInferenceSurvivalRestrictedMeanDiff$new(seq_des)
#' 
#' # compute_treatment_estimate
#' seq_des_inf$compute_treatment_estimate()
#' 
#' # compute_mle_confidence_interval
#' seq_des_inf$compute_mle_confidence_interval()
#' 
#' # compute_mle_two_sided_pval_for_treatment_effect
#' seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect()
#' }
#' 
#' @export
SeqDesignInferenceSurvivalRestrictedMeanDiff = R6::R6Class("SeqDesignInferenceSurvivalRestrictedMeanDiff",
	inherit = SeqDesignInference,
	public = list(
		
		#' @description
		#' Initialize a sequential experimental design estimation and test object after the sequential design is completed.
		#' @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference
		#' 							and bootstrap resampling. The default is 1 for serial computation. For simple estimators (e.g. mean difference 
		#' 							and KK compound), parallelization is achieved with zero-overhead C++ OpenMP. For complex models (e.g. GLMs), 
		#' 							parallelization falls back to R's \code{parallel::mclapply} which incurs session-forking overhead.
		#' @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){		
			super$initialize(seq_des_obj, num_cores, verbose)
			assertResponseType(seq_des_obj$get_response_type(), "survival")
		},
		
		#' @description
		#' Computes the appropriate estimate for mean difference
		#' 
		#' @return 	The setting-appropriate (see description) numeric estimate of the treatment effect
		compute_treatment_estimate = function(){
			if (is.null(private$cached_values$beta_hat_T)){
				private$cached_values$beta_hat_T = get_survival_stat_diff(
					private$y,
					private$dead,
					private$w,
					"restricted_mean"
				)
			}
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes a 1-alpha level frequentist confidence interval differently for all response types, estimate types and test types.
		#'
		#' Here we use the theory that MLE's computed for GLM's are asymptotically normal.
		#' Hence these confidence intervals are asymptotically valid and thus approximate for any sample size.
		#'
		#' @param alpha					The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		#'
		#' @return 	A (1 - alpha)-sized frequentist confidence interval for the treatment effect
		compute_mle_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			if (is.null(private$cached_values$beta_hat_T)){
				self$compute_treatment_estimate()
			}
			if (is.null(private$cached_values$s_beta_hat_T)){
				private$compute_s_beta_hat_T()
			}
			if (is.na(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0) {
				return(self$compute_bootstrap_confidence_interval(alpha = alpha, na.rm = TRUE))
			}
			private$cached_values$is_z = TRUE
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes a 2-sided p-value via the log rank test
		#'
		#' @param delta					The null difference to test against. For any treatment effect at all this is set to zero (the default).
		#'
		#' @return 	The approximate frequentist p-value
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)

			if (delta == 0){
				if (is.null(private$cached_values$s_beta_hat_T)){
					private$compute_s_beta_hat_T()
				}
				if (is.na(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0) {
					return(self$compute_bootstrap_two_sided_pval(delta = delta, na.rm = TRUE))
				}
				z_beta_hat_T = private$cached_values$beta_hat_T / private$cached_values$s_beta_hat_T
				2 * min(stats::pnorm(z_beta_hat_T), 1 - stats::pnorm(z_beta_hat_T))
			} else {
				stop("TO-DO")
			}
		},

		#' @description
		#' Computes a 1-alpha level frequentist confidence interval for the randomization test
		#'
		#' @param alpha					The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		#' @param nsim_exact_test		The number of randomization vectors. The default is 501.
		#' @param pval_epsilon			The bisection algorithm tolerance. The default is 0.005.
		#' @param show_progress		Show a text progress indicator.
		#' @return 	A 1 - alpha sized frequentist confidence interval
		compute_confidence_interval_rand = function(alpha = 0.05, nsim_exact_test = 501, pval_epsilon = 0.005, show_progress = TRUE){
			stop("Randomization confidence intervals are not supported for SeqDesignInferenceSurvivalRestrictedMeanDiff due to inconsistent estimator units on the transformed scale (estimates time difference, but randomization test searches for log-time ratio).")
		}
	),

	private = list(
		compute_s_beta_hat_T = function(){
			se_val = get_restricted_mean_se_diff(
				private$y,
				private$dead,
				private$w
			)
			if (is.na(se_val) || se_val <= 0) {
				warning("Restricted mean SE is non-positive or NA; MLE p-value/CI unavailable.")
				private$cached_values$s_beta_hat_T = NA_real_
				return(invisible(NULL))
			}
			private$cached_values$s_beta_hat_T = se_val
		}
	)
)
