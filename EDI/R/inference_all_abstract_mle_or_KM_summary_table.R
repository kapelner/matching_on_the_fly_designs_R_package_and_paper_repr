#' Inference for A Sequential Design
#'
#' An abstract R6 Class that provides asymptotic tests and intervals for a
#' treatment effect in a sequential design
#' where the common denominator is a summary table from a glm.
#'
#' @keywords internal
InferenceMLEorKMSummaryTable = R6::R6Class("InferenceMLEorKMSummaryTable",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Computes the appropriate estimate for mean difference
		#'
		#' @return 	The setting-appropriate (see description) numeric estimate of the treatment effect
		#'
		#' @examples
		#' \dontrun{
		#' seq_des = DesignSeqOneByOneBernoulli$new(n = 6, response_type = "continuous")
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[1, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[2, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[3, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[4, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[5, 2 : 10])
		#' seq_des$add_one_subject_to_experiment_and_assign(MASS::biopsy[6, 2 : 10])
		#' seq_des$add_all_subject_responses(c(4.71, 1.23, 4.78, 6.11, 5.95, 8.43))
		#'
		#' seq_des_inf = InferenceContinMultOLS$new(seq_des)
		#' seq_des_inf$compute_treatment_estimate()
		#' }
		#'
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes a 1-alpha level frequentist confidence interval
		#'
		#' @param alpha                                   The confidence level in the computed
		#'   confidence interval is 1 - \code{alpha}. The default is 0.05.
		#'
		#' @return 	A (1 - alpha)-sized frequentist confidence interval for the treatment effect
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes a 2-sided p-value
		#'
		#' @param delta                                   The null difference to test against. For
		#'   any treatment effect at all this is set to zero (the default).
		#'
		#' @return 	The approximate frequentist p-value
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared()
			if (should_run_asserts()) {
				if (delta == 0){
					private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
				} else {
					stop("TO-DO")
				}
			}
		}
	),
	private = list(
		generate_mod = function() stop(class(self)[1], " must implement generate_mod()"),

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$summary_table)) return(invisible(NULL))

			model_output = private$generate_mod(estimate_only = estimate_only) # Implemented by child classes (Weibull, NegBin)
			private$cached_mod = model_output

			if (should_run_asserts()) {
				if (is.null(model_output$coefficients) || (!estimate_only && is.null(model_output$vcov))){
					stop("Model output (coefficients or vcov) is NULL or invalid from generate_mod().")
				}
			}

			full_coefficients = model_output$coefficients
			treatment_coef_name = "treatment"
			if (!treatment_coef_name %in% names(full_coefficients)){
				treatment_coef_name = names(full_coefficients)[length(full_coefficients)]
			}
			private$cached_values$beta_hat_T = as.numeric(full_coefficients[treatment_coef_name])

			if (estimate_only) return(invisible(NULL))

			# Construct summary table (as expected by this class)
			full_vcov = model_output$vcov
			diag_vcov = diag(full_vcov)
			# Ensure we don't take sqrt of negative numbers
			full_std_errs = ifelse(diag_vcov > 0, sqrt(diag_vcov), NA_real_)

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

			private$cached_values$s_beta_hat_T = as.numeric(full_std_errs[treatment_coef_name])
			private$cached_values$is_z = TRUE
		}
	)
)
