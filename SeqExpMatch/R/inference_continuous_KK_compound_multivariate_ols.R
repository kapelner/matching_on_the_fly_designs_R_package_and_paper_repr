#' Simple Mean Difference Inference based on Maximum Likelihood  
#'
#' @description
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring) sequential experimental design estimation and test object after the sequential design is completed.
#' 
#' @export
SeqDesignInferenceContinMultOLSKK = R6::R6Class("SeqDesignInferenceContinMultOLSKK",
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
			assertResponseType(seq_des_obj$get_response_type(), "continuous")
			super$initialize(seq_des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		},
		
		#' @description
		#' Computes the appropriate estimate
		#' 
		#' @return 	The setting-appropriate (see description) numeric estimate of the treatment effect
		#' 
		#' @examples
		#' seq_des = SeqDesignCRD$new(n = 6, response_type = "continuous")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des_inf = SeqDesignInferenceContinMultOLS$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' 	
			compute_treatment_estimate = function(){			
				if (is.null(private$cached_values$beta_T_reservoir) & is.null(private$cached_values$beta_T_matched)){
					private$shared_for_compute_estimate()
				}		
				if (is.null(private$cached_values[["beta_hat_T"]])){ 
					beta_hat_T = if (private$only_matches()){
															private$cached_values$beta_T_matched
														} else if (private$only_reservoir()){	 		
															private$cached_values$beta_T_reservoir
														} else {
															w_star = private$cached_values$ssq_beta_T_reservoir / 
																		(private$cached_values$ssq_beta_T_reservoir + private$cached_values$ssq_beta_T_matched)	
															w_star * private$cached_values$beta_T_matched + (1 - w_star) * private$cached_values$beta_T_reservoir
														}
					# Avoid floating-point residue when the estimate is mathematically zero.
					if (is.finite(beta_hat_T) && abs(beta_hat_T) < sqrt(.Machine$double.eps)){
						beta_hat_T = 0
					}
					private$cached_values[["beta_hat_T"]] = beta_hat_T
				}
				private$cached_values[["beta_hat_T"]]
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
		#' 
		#' @examples
		#' \dontrun{
		#' seq_des = SeqDesignCRD$new(n = 6)
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des_inf = SeqDesignInferenceContinMultOLS$new(seq_des)
		#' seq_des_inf$compute_mle_confidence_interval()
		#' }
		#'		
		compute_mle_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)	
			private$shared_for_inference()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		
		#' @description

		
		#' Computes a 2-sided p-value
		#'
		#' @param delta					The null difference to test against. For any treatment effect at all this is set to zero (the default).
		#' 
		#' @return 	The approximate frequentist p-value
		#' 
		#' @examples
		#' \dontrun{
		#' seq_des = SeqDesignCRD$new(n = 6)
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des_inf = SeqDesignInferenceContinMultOLS$new(seq_des)
		#' seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect()
		#' }
		#' 				
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)	
			private$shared_for_inference()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),
	
	private = list(		
		shared_for_compute_estimate = function(){
			if (private$only_matches()){
				private$ols_for_matched_pairs() 
			} else if (private$only_reservoir()){			
				private$ols_for_reservoir()		
			} else {
				private$ols_for_matched_pairs() 		
				private$ols_for_reservoir()
			}
		},
		
		shared_for_inference = function(){
			if (is.null(private$cached_values[["beta_hat_T"]])){ 
				self$compute_treatment_estimate()
			}
			ssq_beta_hat_T = 	if (private$only_matches()){
									private$cached_values$ssq_beta_T_matched
								} else if (private$only_reservoir()){			
									private$cached_values$ssq_beta_T_reservoir
								} else {
									private$cached_values$ssq_beta_T_matched * private$cached_values$ssq_beta_T_reservoir / 
										(private$cached_values$ssq_beta_T_matched + private$cached_values$ssq_beta_T_reservoir) #analagous eq's
								}
			private$cached_values$s_beta_hat_T = sqrt(ssq_beta_hat_T)
			private$cached_values$is_z = TRUE #TO-DO: linear combination of degrees of freedom of t's
		},
		
		only_matches = function(){
			private$cached_values$KKstats$nRT <= 2 || private$cached_values$KKstats$nRC <= 2 || (private$cached_values$KKstats$nRT + private$cached_values$KKstats$nRC <= ncol(private$get_X()) + 2)
		},
		
		only_reservoir = function(){
			private$cached_values$KKstats$m == 0
		},
		
		ols_for_matched_pairs = function(){
			# stats::coef(summary(lm(private$cached_values$KKstats$y_matched_diffs ~ private$cached_values$KKstats$X_matched_diffs)))
			mod = fast_ols_with_var_cpp(cbind(1, private$cached_values$KKstats$X_matched_diffs), private$cached_values$KKstats$y_matched_diffs, j = 1) #the only time you need the intercept's ssq
			private$cached_values$beta_T_matched =     mod$b[1]
			private$cached_values$ssq_beta_T_matched = mod$ssq_b_j
		},

		ols_for_reservoir = function(){
			# stats::coef(summary(lm(private$cached_values$KKstats$y_reservoir ~ cbind(private$cached_values$KKstats$w_reservoir, private$cached_values$KKstats$X_reservoir))))
			mod = fast_ols_with_var_cpp(cbind(1, private$cached_values$KKstats$w_reservoir, private$cached_values$KKstats$X_reservoir), private$cached_values$KKstats$y_reservoir)
			private$cached_values$beta_T_reservoir =     mod$b[2]
			private$cached_values$ssq_beta_T_reservoir = mod$ssq_b_j
		}
	)		
)
