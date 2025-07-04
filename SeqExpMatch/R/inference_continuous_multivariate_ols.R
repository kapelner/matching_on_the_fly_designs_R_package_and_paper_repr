#' Simple Mean Difference Inference based on Maximum Likelihood  
#'
#' @description
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring) sequential experimental design estimation and test object after the sequential design is completed.
#' 
#'
#' @export
SeqDesignInferenceContinMultOLS = R6::R6Class("SeqDesignInferenceContinMultOLS",
	inherit = SeqDesignInferenceMLEorKMSummaryTable,
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
			assertResponseType(seq_des_obj$get_response_type(), "continuous")			
			super$initialize(seq_des_obj, num_cores, verbose)	
			assertNoCensoring(private$any_censoring)
			private$cached_values = super$get_cached_values()
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
			if (is.null(private$cached_values$summary_table)){
				private$shared()
			}			
			if (is.null(private$cached_values$beta_T)){
				private$cached_values$beta_hat_T = self$compute_treatment_estimate()
			}
			private$cached_values$s_beta_hat_T = private$cached_values$summary_table[2, 2]
			private$cached_values$is_z = FALSE 
			private$cached_values$df = private$n - nrow(private$cached_values$summary_table)			
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		}
	),
	
	private = list(		
		cached_values = list(),
		
		shared = function(){
			private$cached_values$summary_table = 
				coef(summary(lm(private$seq_des_obj_priv_int$y ~ ., data = cbind(data.frame(w = private$seq_des_obj_priv_int$w), private$get_X()))))

		}
		
	)		
)