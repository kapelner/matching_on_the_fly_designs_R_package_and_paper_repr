#' Simple Mean Difference Inference based on Maximum Likelihood  
#'
#' @description
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring) sequential experimental design estimation and test object after the sequential design is completed.
#' 
#' @export
SeqDesignInferenceContinMultOLSKK = R6::R6Class("SeqDesignInferenceContinMultOLSKK",
	inherit = SeqDesignInferenceMLEorKMKK,
	public = list(
		
		#' @description
		#' Initialize a sequential experimental design estimation and test object after the sequential design is completed.
        #' @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference 
		#' 							or bootstrap inference.
		#' @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		#'
		initialize = function(seq_des_obj, num_cores = 1, verbose = TRUE){			
			assertResponseType(seq_des_obj$get_response_type(), "continuous")
			super$initialize(seq_des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)				
			private$cached_values = super$get_cached_values()
		},
		
		#' @description
		#' Computes the appropriate estimate for mean difference
		#' 
		#' @return 	The setting-appropriate (see description) numeric estimate of the treatment effect
		#' 
		#' @examples
		#' seq_des = SeqDesign$new(n = 6, p = 10, design = "CRD", response_type = "continuous")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des_inf = SeqDesignInferenceContMultOLS$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' 	
		compute_treatment_estimate = function(){	
			if (is.null(private$cached_values$private$cached_values$coefs_matched) & is.null(private$cached_values$private$cached_values$coefs_reservoir)){
				private$shared()
			}		
			if (is.null(private$cached_values$beta_hat_T)){ 
				#no reservoir
				if (private$KKstats$nRT <= 2 || private$KKstats$nRC <= 2 || (private$KKstats$nRT + private$KKstats$nRC <= ncol(private$get_X()) + 2)){
					private$cached_values$beta_hat_T = private$cached_values$coefs_matched[1, 1]
				#and sometimes there's no matches	
				} else if (private$KKstats$m == 0){	 		
					private$cached_values$beta_hat_T = private$cached_values$coefs_reservoir[2, 1]
				#but most of the time... we have matches and a nice-sized reservoir
				} else { 
					beta_match_regression = private$cached_values$coefs_matched[1, 1]			
					beta_reservoir_regression = private$cached_values$coefs_reservoir[2, 1]
					w_star = private$cached_values$ssqd_reservoir_regression / (private$cached_values$ssqd_reservoir_regression + private$cached_values$ssqd_match_regression) #just a convenience for faster runtime	
					private$cached_values$beta_hat_T = w_star * beta_match_regression + (1 - w_star) * beta_reservoir_regression #proper weighting				}
				}
			}
			private$cached_values$beta_hat_T
		},
		
		
		
		#' Compute confidence interval
		#'
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
		#' seq_des = SeqDesign$new(n = 6, p = 10, design = "CRD")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des_inf = SeqDesignInferenceContMultOLS$new(seq_des, test_type = "MLE-or-KM-based")
		#' seq_des_inf$compute_confidence_interval()
		#'		
		compute_mle_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)	
			if (is.null(private$cached_values$beta_hat_T)){ 
				self$compute_treatment_estimate()
			}
			if (private$KKstats$nRT <= 2 || private$KKstats$nRC <= 2 || (private$KKstats$nRT + private$KKstats$nRC <= ncol(private$get_X()) + 2)){
				private$cached_values$s_beta_hat_T = private$cached_values$coefs_matched[1, 2]	
			} else if (private$KKstats$m == 0){			
				private$cached_values$s_beta_hat_T = coefs_reservoir[2, 2]
			} else {
				private$cached_values$s_beta_hat_T = sqrt(private$cached_values$ssqd_match_regression * private$cached_values$ssqd_reservoir_regression / 
					(private$cached_values$ssqd_match_regression + private$cached_values$ssqd_reservoir_regression)) #analagous eq's
			}
			private$cached_values$is_z = TRUE #TO-DO: linear combination of degrees of freedom of t's
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		
		#' Compute p-value
		#'
		#' @description
		#' Computes a 2-sided p-value
		#'
		#' @param delta					The null difference to test against. For any treatment effect at all this is set to zero (the default).
		#' 
		#' @return 	The approximate frequentist p-value
		#' 
		#' @examples
		#' seq_des = SeqDesign$new(n = 6, p = 10, design = "CRD")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des_inf = SeqDesignInferenceContMultOLS$new(seq_des)
		#' seq_des_inf$compute_two_sided_pval_for_treatment_effect()
		#' 				
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			if (is.null(private$cached_values$beta_hat_T)){ 
				self$compute_treatment_estimate()
			}
			if (private$KKstats$nRT <= 2 || private$KKstats$nRC <= 2 || (private$KKstats$nRT + private$KKstats$nRC <= ncol(private$get_X()) + 2)){
				if (delta == 0){
					stop("TO-DO")
				}
				private$cached_values$coefs_matched[1, 4]
			} else if (private$KKstats$m == 0){			
				if (delta == 0){
					stop("TO-DO")
				}
				private$cached_values$coefs_reservoir[2, 4]
			} else {
				s_beta_hat_T = sqrt(private$cached_values$ssqd_match_regression * private$cached_values$ssqd_reservoir_regression / 
					(private$cached_values$ssqd_match_regression + private$cached_values$ssqd_reservoir_regression)) #analagous eq's
				2 * (pnorm(-abs((private$cached_values$beta_hat_T - delta) / s_beta_hat_T))) #approximate by using N(0, 1) distribution			
			}
		}
	),
	
	private = list(		
		cached_values = list(),
		
		shared = function(){
			if (private$KKstats$nRT <= 2 || private$KKstats$nRC <= 2 || (private$KKstats$nRT + private$KKstats$nRC <= ncol(private$get_X()) + 2)){
				private$cached_values$coefs_matched = coef(summary(lm(private$KKstats$y_matched_diffs ~ private$KKstats$X_matched_diffs)))
			} else if (private$KKstats$m == 0){			
				private$cached_values$coefs_reservoir = coef(summary(lm(private$KKstats$y_reservoir ~ cbind(private$KKstats$w_reservoir, private$KKstats$X_reservoir))))		
			} else {
				#compute estimator from matched pairs by regression
				private$cached_values$coefs_matched = coef(summary(lm(private$KKstats$y_matched_diffs ~ private$KKstats$X_matched_diffs)))
				private$cached_values$ssqd_match_regression = private$cached_values$coefs_matched[1, 2]^2 #lin mod returns SE not VAR, so square it				
				
				#compute estimator reservoir sample std error
				private$cached_values$coefs_reservoir = coef(summary(lm(private$KKstats$y_reservoir ~ cbind(private$KKstats$w_reservoir, private$KKstats$X_reservoir))))
				private$cached_values$ssqd_reservoir_regression = private$cached_values$coefs_reservoir[2, 2]^2 #lin mod returns SE not VAR, so square it
			}
		}		
	)		
)