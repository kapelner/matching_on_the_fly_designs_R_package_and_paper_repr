#' Jonckheere-Terpstra Test for Ordinal Responses
#'
#' @description
#' Exact Jonckheere-Terpstra rank test for a two-arm ordered alternative with an
#' ordinal response. For treatment versus control, the test statistic is the
#' pairwise superiority score, counting 1 for treatment greater than control and
#' 0.5 for ties. Under a fixed treatment-arm size, the exact p-value is computed
#' by enumerating the tied-score distribution conditional on the observed
#' category counts.
#'
#' This class is test-first: it provides an exact p-value and a simple centered
#' superiority estimate, but it does not implement MLE-based or bootstrap-based
#' inference.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignCRD$
#'   new(
#'   n = nrow(x_dat),
#'   response_type = "ordinal",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(as.integer(c(1, 2, 2, 3, 3, 4, 4, 5)))
#' infer <- SeqDesignInferenceOrdinalJonckheereTerpstraTest$
#'   new(seq_des, verbose = FALSE)
#' infer$
#'   compute_exact_two_sided_pval_for_treatment_effect()
SeqDesignInferenceOrdinalJonckheereTerpstraTest = R6::R6Class("SeqDesignInferenceOrdinalJonckheereTerpstraTest",
	inherit = SeqDesignInference,
	public = list(

		#' @description
		#' Initialize the exact Jonckheere-Terpstra inference object.
		#' @param seq_des_obj A completed \code{SeqDesign} object with an ordinal response.
		#' @param num_cores Number of CPU cores for inherited generic infrastructure.
		#' @param verbose Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "ordinal")
			super$initialize(seq_des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Returns the centered treatment superiority estimate
		#' \eqn{P(Y_T > Y_C) + 0.5 P(Y_T = Y_C) - 0.5}.
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Not implemented for the Jonckheere-Terpstra exact test.
		#' @param alpha Unused confidence level.
		compute_mle_confidence_interval = function(alpha = 0.05){
			stop(class(self)[1], " only implements compute_exact_two_sided_pval_for_treatment_effect. This method is not implemented.")
		},

		#' @description
		#' Not implemented for the Jonckheere-Terpstra exact test.
		#' @param delta Unused null effect.
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			stop(class(self)[1], " only implements compute_exact_two_sided_pval_for_treatment_effect. This method is not implemented.")
		},

		#' @description
		#' Not implemented for the Jonckheere-Terpstra exact test.
		#' @param ... Unused arguments.
		compute_exact_confidence_interval = function(...){
			stop(class(self)[1], " only implements compute_exact_two_sided_pval_for_treatment_effect. This method is not implemented.")
		},

		#' @description
		#' Not implemented for the Jonckheere-Terpstra exact test.
		#' @param ... Unused arguments.
		compute_confidence_interval_rand = function(...){
			stop(class(self)[1], " only implements compute_exact_two_sided_pval_for_treatment_effect. This method is not implemented.")
		},

		#' @description
		#' Not implemented for the Jonckheere-Terpstra exact test.
		#' @param ... Unused arguments.
		compute_two_sided_pval_for_treatment_effect_rand = function(...){
			stop(class(self)[1], " only implements compute_exact_two_sided_pval_for_treatment_effect. This method is not implemented.")
		},

		#' @description
		#' Not implemented for the Jonckheere-Terpstra exact test.
		#' @param ... Unused arguments.
		approximate_bootstrap_distribution_beta_hat_T = function(...){
			stop(class(self)[1], " only implements compute_exact_two_sided_pval_for_treatment_effect. This method is not implemented.")
		},

		#' @description
		#' Not implemented for the Jonckheere-Terpstra exact test.
		#' @param ... Unused arguments.
		compute_bootstrap_confidence_interval = function(...){
			stop(class(self)[1], " only implements compute_exact_two_sided_pval_for_treatment_effect. This method is not implemented.")
		},

		#' @description
		#' Not implemented for the Jonckheere-Terpstra exact test.
		#' @param ... Unused arguments.
		compute_bootstrap_two_sided_pval = function(...){
			stop(class(self)[1], " only implements compute_exact_two_sided_pval_for_treatment_effect. This method is not implemented.")
		},

		#' @description
		#' Computes the exact two-sided Jonckheere-Terpstra p-value for the null of
		#' no ordered treatment effect.
		#' @param delta The null effect. Only \code{delta = 0} is supported.
		compute_exact_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta, len = 1)
			if (!isTRUE(all.equal(delta, 0))) {
				stop(class(self)[1], " only supports delta = 0.")
			}
			private$shared()
			private$cached_values$p_exact
		}
	),

	private = list(
		shared = function(){
			if (!is.null(private$cached_values$p_exact)) return(invisible(NULL))

			res = exact_jonckheere_terpstra_pval_cpp(
				y = as.integer(private$y),
				w = as.integer(private$w)
			)

			private$cached_values$jt_stat2 = res$stat2
			private$cached_values$superiority = res$superiority
			private$cached_values$beta_hat_T = res$superiority - 0.5
			private$cached_values$p_exact = res$p_exact
			private$cached_values$p_lower = res$p_lower
			private$cached_values$p_upper = res$p_upper
		}
	)
)
