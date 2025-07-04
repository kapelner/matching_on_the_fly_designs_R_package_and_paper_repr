#' Simple Mean Difference Inference based on Maximum Likelihood  
#'
#' @description
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring) sequential experimental design estimation and test object after the sequential design is completed.
#' 
#'
#' @export
SeqDesignInferenceAllSimpleMeanDiff = R6::R6Class("SeqDesignInferenceAllSimpleMeanDiff",
	inherit = SeqDesignInferenceMLEorKM,
	public = list(
		
		#' @description
		#' Initialize a sequential experimental design estimation and test object after the sequential design is completed.
        #' @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		#' @param num_cores			The number of CPU cores to use to parallelize the sampling during randomization-based inference 
		#' 							(which is very slow). The default is 1 for serial computation. This parameter is ignored
		#' 							for \code{test_type = "MLE-or-KM-based"}.
		#' @param verbose			A flag indicating whether messages should be displayed to the user. Default is \code{TRUE}
		#'
		initialize = function(seq_des_obj, num_cores = 1, verbose = TRUE){						
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
		#' seq_des_inf = SeqDesignInferenceAllSimpleMeanDiffMLE$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' 	
		compute_treatment_estimate = function(){
			if (is.null(private$cached_values$beta_hat_T)){
				private$cached_values$beta_hat_T = mean(super$get_yTs()) - mean(super$get_yCs())
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
		#' seq_des_inf = SeqDesignInferenceAllSimpleMeanDiffMLE$new(seq_des, test_type = "MLE-or-KM-based")
		#' seq_des_inf$compute_confidence_interval()
		#'		
		compute_mle_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)	
			if (is.null(private$cached_values$s_beta_hat_T)){
				private$shared()
			}
			
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},
		
		#' Compute p-value
		#'
		#' @description
		#' Computes a 2-sided p-value for all types of inferential settings written about in the initializer
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
		#' seq_des_inf = SeqDesignInferenceAllSimpleMeanDiffMLE$new(seq_des)
		#' seq_des_inf$compute_two_sided_pval_for_treatment_effect()
		#' 				
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			if (is.null(private$cached_values$df)){
				private$shared()
			}
			
			2 * pt(
					-abs(private$cached_values$beta_hat_T - delta) / private$cached_values$s_beta_hat_T, 
						private$cached_values$df
				)
		}
	),
	
	private = list(		
		cached_values = list(),
		
		shared = function(){
			if (is.null(private$cached_values$beta_hat_T)){
				self$compute_treatment_estimate()
			}
			
			nT = length(super$get_yTs())
			nC = length(super$get_yCs())
			s_1_sq = var(super$get_yTs()) / nT 
			s_2_sq = var(super$get_yCs()) / nC
			private$cached_values$s_beta_hat_T = sqrt(s_1_sq + s_2_sq)
			private$cached_values$df = (s_1_sq + s_2_sq)^2 / (
											s_1_sq^2 / (nT - 1) + s_2_sq^2 / (nC - 1)
										) #Welch-Satterthwaite formula
			private$cached_values$is_z = FALSE
		}
		
	)		
)