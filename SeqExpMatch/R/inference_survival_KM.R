#' Simple Mean Difference Inference based on Maximum Likelihood  
#'
#' @description
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring) sequential experimental design estimation and test object after the sequential design is completed.
#' 
#'
#' @export
SeqDesignInferenceSurvivalKMDiff = R6::R6Class("SeqDesignInferenceSurvivalKMDiff",
	inherit = SeqDesignInference,
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
			assertResponseType(seq_des_obj$get_response_type(), "survival")	
		},
		
		#' @description
		#' Computes the appropriate estimate for mean difference
		#' 
		#' @return 	The setting-appropriate (see description) numeric estimate of the treatment effect
		#' 
		#' @examples
		#' seq_des = SeqDesignCRD$new(n = 6, response_type = "survival")
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43), deads = rep(1, 6))
		#' 
		#' seq_des_inf = SeqDesignInferenceSurvivalKMDiff$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' 	
		compute_treatment_estimate = function(){	
			get_survival_stat_diff(
				private$y, 
				private$dead,
				private$w,
				"median"
			)
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
		#' seq_des_inf = SeqDesignInferenceSurvivalKMDiff$new(seq_des)
		#' seq_des_inf$compute_mle_confidence_interval()
		#' }
		#'		
		compute_mle_confidence_interval = function(alpha = 0.05){
			stop("not implemented --- use the bootstrap instead")
		},
		
		#' @description		
		#' Computes a 2-sided p-value via the log rank test
		#'
		#' @param delta					The null difference to test against. For any treatment effect at all this is set to zero (the default).
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
		#' seq_des_inf = SeqDesignInferenceSurvivalKMDiff$new(seq_des)
		#' seq_des_inf$compute_mle_two_sided_pval_for_treatment_effect()
		#' }
		#' 				
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)

			if (delta == 0){
				survival_obj = survival::Surv(private$y, private$dead)
				surv_diff = survival::survdiff(survival_obj ~ private$w)
				surv_diff$pvalue
			} else {
				stop("TO-DO")
			}
		}
	),
	
	private = list()		
)
