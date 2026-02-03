#' Inference for A Sequential Design
#'
#' @description
#' An abstract R6 Class that provides MLE-based tests and intervals for a treatment effect in a sequential design
#' where the common denominator is a summary table from a glm.
#' 
SeqDesignInferenceMLEorKMSummaryTable = R6::R6Class("SeqDesignInferenceMLEorKMSummaryTable",
	inherit = SeqDesignInferenceKKPassThrough,
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
		},
		
		#' @description
		#' Computes the appropriate estimate for mean difference
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
		#' seq_des_inf = SeqDesignInferenceContinMultOLS$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' }
		#' 	
		compute_treatment_estimate = function(){
			private$shared()			
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes a 1-alpha level frequentist confidence interval
		#'
		#' @param alpha					The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		#' 
		#' @return 	A (1 - alpha)-sized frequentist confidence interval for the treatment effect
		compute_mle_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes a 2-sided p-value
		#'
		#' @param delta					The null difference to test against. For any treatment effect at all this is set to zero (the default).
		#' 
		#' @return 	The approximate frequentist p-value
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			if (delta == 0){
				private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			} else {
				stop("TO-DO")
			}
		}		
	),
	private = list(
		shared = function(){
			model_output = private$generate_mod() # Implemented by child classes (Weibull, NegBin)
			
			if (is.null(model_output$coefficients) || is.null(model_output$vcov)){
				stop("Model output (coefficients or vcov) is NULL or invalid from generate_mod().")
			}
			
			# Construct summary table (as expected by this class)
			full_coefficients = model_output$coefficients
			full_vcov = model_output$vcov
			full_std_errs = sqrt(diag(full_vcov))
			
			summary_table = matrix(NA, nrow = length(full_coefficients), ncol = 4)
			rownames(summary_table) = names(full_coefficients)
			colnames(summary_table) = c("Value", "Std. Error", "z value", "Pr(>|z|)")
			
			summary_table[, 1] = full_coefficients
			summary_table[, 2] = full_std_errs
			summary_table[, 3] = full_coefficients / full_std_errs # z value
			summary_table[, 4] = 2 * stats::pnorm(-abs(summary_table[, 3])) # p-value
			
			private$cached_values$summary_table = summary_table

			# Populate full_coefficients and full_vcov for consistency across all inference classes
			private$cached_values$full_coefficients = full_coefficients
			private$cached_values$full_vcov = full_vcov

			# Extract beta_hat_T and s_beta_hat_T
			treatment_coef_name = "treatment"
			if (!treatment_coef_name %in% names(full_coefficients)){
				warning("Treatment coefficient 'treatment' not found by name. Attempting to identify it by position (last coefficient). This may lead to incorrect results.")
				treatment_coef_name = names(full_coefficients)[length(full_coefficients)]
			}
			private$cached_values$beta_hat_T = full_coefficients[treatment_coef_name]
			private$cached_values$s_beta_hat_T = sqrt(full_vcov[treatment_coef_name, treatment_coef_name])
			private$cached_values$is_z = TRUE
		}
	)		
)
