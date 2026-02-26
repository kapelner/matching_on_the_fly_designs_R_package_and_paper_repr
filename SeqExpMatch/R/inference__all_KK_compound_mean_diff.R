#' Inference based on Maximum Likelihood for KK designs  
#'
#' @description
#' Inference for mean difference
#' 
#'
#' @export
SeqDesignInferenceAllKKCompoundMeanDiff = R6::R6Class("SeqDesignInferenceAllKKCompoundMeanDiff",
	inherit = SeqDesignInferenceKKPassThroughCompound,
	public = list(
		
		#' @description
		#' Initialize a sequential experimental design estimation and test object after the sequential design is completed.
		#' @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference
		#' 								(which is very slow). The default is 1 for serial computation. This parameter is ignored
		#' 								for \code{test_type = "MLE-or-KM-based"}.
		#' @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			super$initialize(seq_des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)		
		},
		
		#'	
		#' @description
		#' Computes the appropriate estimate for compound mean difference across pairs and reservoir
		#' 
		#' @return 	The setting-appropriate (see description) numeric estimate of the treatment effect
		#' 
		#' @examples
		#' \dontrun{
		#' seq_des = SeqDesignCRD$new(n = 6, response_type = "continuous")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des_inf = SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' }
		#' 	
		compute_treatment_estimate = function(){			
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
				private$compute_reservoir_and_match_statistics()
			}
			private$cached_values$beta_hat_T = 	if (private$cached_values$KKstats$nRT <= 1 || private$cached_values$KKstats$nRC <= 1){
													private$cached_values$KKstats$d_bar	
												} else if (private$cached_values$KKstats$m == 0){ #sometimes there's no matches
													private$cached_values$KKstats$r_bar			
												} else {
													private$cached_values$KKstats$w_star * private$cached_values$KKstats$d_bar + (1 - private$cached_values$KKstats$w_star) * private$cached_values$KKstats$r_bar #proper weighting
												}
			private$cached_values$beta_hat_T
		},
		
		#' @description

		
		#' Computes a 1-alpha level frequentist confidence interval
		#' 
		#' Here we use the theory that MLE's computed for GLM's are asymptotically normal (except in the case 
		#' of estimat_type "median difference" where a nonparametric bootstrap confidence interval (see the \code{controlTest::quantileControlTest} method)
		#' is employed. Hence these confidence intervals are asymptotically valid and thus approximate for any sample size.
		#' 
		#' @param alpha					The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		#' 
		#' @return 	A (1 - alpha)-sized frequentist confidence interval for the treatment effect
		#' 
		#' @examples
		#' \dontrun{
		#' seq_des = SeqDesignCRD$new(n = 6, response_type = "continuous")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des_inf = SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des)
		#' seq_des_inf$compute_mle_confidence_interval()
		#' }
		#'	
		compute_mle_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			
			if (is.null(private$cached_values$s_beta_hat_T)){
				private$shared()
			}		
			private$cached_values$is_z = TRUE
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		
		#' @description

		
		#' Computes a 2-sided p-value
		#'
		#' @param delta	The null difference to test against. For any treatment effect at all this is set to zero (the default).
		#' 
		#' @return 	The approximate frequentist p-value
		#' 
		#' @examples
		#' \dontrun{
		#' seq_des = SeqDesignCRD$new(n = 6, response_type = "continuous")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des_inf = SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des)
		#' seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect()
		#' }
		#' 		
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			if (is.null(private$cached_values$s_beta_hat_T)){
				private$shared()
			}
			private$cached_values$is_z = TRUE
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
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
			if (private$seq_des_obj_priv_int$response_type %in% c("proportion", "count", "survival")) {
				stop("Randomization confidence intervals are not supported for SeqDesignInferenceAllKKCompoundMeanDiff with proportion, count, or survival response types due to inconsistent estimator units on the transformed scale.")
			}
			super$compute_confidence_interval_rand(alpha = alpha, nsim_exact_test = nsim_exact_test, pval_epsilon = pval_epsilon, show_progress = show_progress)
		}
	),

	private = list(
		shared = function(){
			if (is.null(private$cached_values$beta_hat_T)){
				self$compute_treatment_estimate()
			}

			ssqD = private$cached_values$KKstats$ssqD_bar
			ssqR = private$cached_values$KKstats$ssqR

			private$cached_values$s_beta_hat_T =
				if (private$cached_values$KKstats$nRT <= 1 || private$cached_values$KKstats$nRC <= 1){
					# Only matched pairs are usable; fall back to ssqR if ssqD is degenerate
					if (is.finite(ssqD) && ssqD > 0) sqrt(ssqD) else if (is.finite(ssqR) && ssqR > 0) sqrt(ssqR) else NA_real_
				} else if (private$cached_values$KKstats$m == 0){
					# No matched pairs
					if (is.finite(ssqR) && ssqR > 0) sqrt(ssqR) else NA_real_
				} else {
					# Combined: require both components to be positive and finite.
					# If one collapses (e.g. tiny reservoir with near-identical responses),
					# fall back to the other rather than pulling the combined SE to zero.
					if (!is.finite(ssqD) || ssqD <= 0) {
						if (is.finite(ssqR) && ssqR > 0) sqrt(ssqR) else NA_real_
					} else if (!is.finite(ssqR) || ssqR <= 0) {
						sqrt(ssqD)
					} else {
						sqrt(ssqR * ssqD / (ssqR + ssqD))
					}
				}
		}		
	)		
)
