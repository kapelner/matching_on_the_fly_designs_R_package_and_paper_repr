#' Ridit Analysis for Ordinal Responses
#'
#' @description
#' Performs Ridit analysis (Relative to an Identified Distribution unit) for
#' comparing two groups on an ordinal scale. Ridit analysis provides a
#' distribution-free way to estimate the probability that a randomly selected
#' member of the treatment group has a better outcome than a randomly selected
#' member of the control group.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignCRD$new(n = nrow(x_dat), response_type = "ordinal", verbose = FALSE)
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- SeqDesignInferenceOrdinalRidit$new(seq_des, verbose = FALSE)
#' infer
#'
SeqDesignInferenceOrdinalRidit = R6::R6Class("SeqDesignInferenceOrdinalRidit",
	inherit = SeqDesignInference,
	public = list(

		#' @description
		#' Initialize a Ridit analysis inference object.
		#' @param seq_des_obj A SeqDesign object whose entire n subjects are assigned and
		#'   response y is recorded within.
		#' @param reference The group to use as the "Identified Distribution" (reference).
		#'   Must be one of "control", "treatment", or "pooled". Default is "control".
		#' @param num_cores The number of CPU cores to use.
		#' @param verbose A flag indicating whether messages should be displayed.
		initialize = function(seq_des_obj, reference = "control", num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "ordinal")
			assertChoice(reference, c("control", "treatment", "pooled"))
			super$initialize(seq_des_obj, num_cores, verbose)
			private$reference = reference
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Returns the estimated treatment effect (Mean Ridit - 0.5).
		#' @return The numeric estimate.
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Returns the Mean Ridit for the treatment group.
		#' @return The numeric Mean Ridit.
		get_mean_ridit_treatment = function(){
			private$shared()
			private$cached_values$mean_ridit_t
		},

		#' @description
		#' Returns the ridit scores for all subjects.
		#' @return A numeric vector of scores.
		get_ridit_scores = function(){
			private$shared()
			private$cached_values$scores
		},

		#' @description
		#' Computes the asymptotic confidence interval for the treatment effect.
		#' @param alpha Significance level.
		#' @return A numeric vector of length 2.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			private$shared()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes the p-value for the null hypothesis that Mean Ridit = 0.5.
		#' @param delta The null value (centered at 0, so delta=0 means Ridit=0.5).
		#' @return The p-value.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		reference = NULL,

		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			res = fast_ridit_analysis_cpp(
				y = as.integer(private$y),
				w = as.integer(private$w),
				reference = private$reference
			)

			if (is.null(res) || length(res) == 0){
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			private$cached_values$mean_ridit_t = res$mean_ridit_t
			private$cached_values$mean_ridit_c = res$mean_ridit_c
			private$cached_values$beta_hat_T   = res$estimate
			private$cached_values$s_beta_hat_T = res$se
			private$cached_values$scores       = res$scores
			private$cached_values$is_z         = TRUE
		}
	)
)
