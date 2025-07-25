#' Inference based on Maximum Likelihood for KK designs  
#'
#' @description
#' Inference for mean difference
#' 
#'
#' @export
SeqDesignInferenceAllKKCompoundMeanDiff = R6::R6Class("SeqDesignInferenceAllKKCompoundMeanDiff",
	inherit = SeqDesignInferenceMLEorKMKK,
	public = list(
		
		#' @description
		#' Initialize a sequential experimental design estimation and test object after the sequential design is completed.
        #' @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference 
		#' 							(which is very slow). The default is 1 for serial computation. This parameter is ignored
		#' 							for \code{test_type = "MLE-or-KM-based"}.
		#' @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		#' @param thin		For internal use only. Do not specify. You can thank R6's single constructor-only for this coding noise.
		#' 
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE, thin = FALSE){	
			if (!thin){	
				super$initialize(seq_des_obj, num_cores, verbose)
				assertNoCensoring(private$any_censoring)	
			}		
		},
		
		#' Compute treatment effect
		#'	
		#' @description
		#' Computes the appropriate estimate for compound mean difference across pairs and reservoir
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
		#' seq_des_inf = SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' 	
		compute_treatment_estimate = function(){			
			private$cached_values$beta_hat_T = 	if (private$KKstats$nRT <= 1 || private$KKstats$nRC <= 1){
													private$KKstats$d_bar	
												} else if (private$KKstats$m == 0){ #sometimes there's no matches
													private$KKstats$r_bar			
												} else {
													private$KKstats$w_star * private$KKstats$d_bar + (1 - private$KKstats$w_star) * private$KKstats$r_bar #proper weighting
												}
			private$cached_values$beta_hat_T
		},
		
		#' Compute confidence interval
		#'
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
		#' seq_des = SeqDesign$new(n = 6, p = 10, design = "CRD")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#' 
		#' seq_des_inf = SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des, test_type = "MLE-or-KM-based")
		#' seq_des_inf$compute_confidence_interval()
		#'	
		compute_mle_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			
			if (is.null(private$cached_values$s_beta_hat_T)){
				private$shared()
			}		
			private$cached_values$is_z = TRUE
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		
		#' Compute p-value
		#'
		#' @description
		#' Computes a 2-sided p-value
		#'
		#' @param delta	The null difference to test against. For any treatment effect at all this is set to zero (the default).
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
		#' seq_des_inf = SeqDesignInferenceAllKKCompoundMeanDiff$new(seq_des)
		#' seq_des_inf$compute_two_sided_pval_for_treatment_effect()
		#' 		
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			if (is.null(private$cached_values$s_beta_hat_T)){
				private$shared()
			}			
			2 * pnorm(
					-abs(private$cached_values$beta_hat_T / private$cached_values$s_beta_hat_T)
				) #approximate by using N(0, 1) distribution

		}
	),
	
	private = list(
					
		shared = function(){
			if (is.null(private$cached_values$beta_hat_T)){
				private$compute_KK_compound_mean_difference_estimate()
			}			
			
			private$cached_values$s_beta_hat_T =
				if (private$KKstats$nRT <= 1 || private$KKstats$nRC <= 1){	
					sqrt(private$KKstats$ssqD_bar)
				} else if (private$KKstats$m == 0){ #sometimes there's no matches
					sqrt(private$KKstats$ssqR)		
				} else {
					sqrt(private$KKstats$ssqR * private$KKstats$ssqD_bar / (private$KKstats$ssqR + private$KKstats$ssqD_bar))
				}			
		}		
	)		
)