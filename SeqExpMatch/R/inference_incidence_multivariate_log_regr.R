#' Simple Mean Difference Inference based on Maximum Likelihood  
#'
#' @description
#' The methods that support confidence intervals and testing for the mean difference
#' in all response types (except Weibull with censoring) sequential experimental design estimation and test object after the sequential design is completed.
#' 
#'
#' @export
SeqDesignInferenceIncidMultiLogRegr = R6::R6Class("SeqDesignInferenceIncidMultiLogRegr",
	inherit = SeqDesignInferenceIncidUnivLogRegr,
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
			}
		},
		
		#' @description
		#' Computes the appropriate estimate
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
			fast_logistic_regression(cbind(1, private$seq_des_obj_priv_int$w, private$get_X()), private$seq_des_obj_priv_int$y)$b[2]
		}
	),
	
	private = list(		
		generate_mod = function(){
			fast_logistic_regression_with_var(cbind(1, private$seq_des_obj_priv_int$w, private$get_X()), private$seq_des_obj_priv_int$y)
		}
		
#		compute_summary_table = function(){
##			coef(summary_glm_lean(suppressWarnings(glm(private$seq_des_obj_priv_int$y ~ ., 
##						data = cbind(data.frame(w = private$seq_des_obj_priv_int$w), private$get_X()), family = "binomial"))))
#			summary_table = matrix(NA, nrow = 2, ncol = 4)
#			mod = fast_logistic_regression_with_sd(cbind(private$seq_des_obj_priv_int$w, private$get_X()), private$seq_des_obj_priv_int$y)
#			summary_table[, 1] = mod$b
#			summary_table[2, 2] = mod$s_b_2
#			summary_table
#		}	
	)		
)