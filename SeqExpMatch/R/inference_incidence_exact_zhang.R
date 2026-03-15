#' Exact incidence inference via Zhang's combined reservoir test
#'
#' @description
#' Provides incidence-specific exact inference based on the Zhang (2026) combined
#' method. For CRD designs this uses the reservoir exact test only.
#'
#' @details
#' This class provides the \code{compute_exact_two_sided_pval_for_treatment_effect}
#' and \code{compute_exact_confidence_interval} methods for uncensored incidence
#' outcomes under a completely randomized design. Standard estimator-based methods
#' such as \code{compute_treatment_estimate()} and
#' \code{compute_mle_confidence_interval()} are intentionally unavailable for this
#' exact Zhang interface.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- SeqDesignCRD$new(n = nrow(x_dat), response_type = "incidence", verbose = FALSE)
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$add_all_subject_responses(c(0, 1, 0, 1, 0, 1, 1, 0))
#' infer <- SeqDesignInferenceIncidExactZhang$new(seq_des, verbose = FALSE)
#' infer
#'
SeqDesignInferenceIncidExactZhang = R6::R6Class("SeqDesignInferenceIncidExactZhang",
	inherit = SeqDesignInferenceIncidExactZhangAbstract,
	public = list(

		#' @description
		#' Initialize the exact incidence inference object.
		#' @param	seq_des_obj		A SeqDesign object with an incidence response.
		#' @param	num_cores			Number of CPU cores for parallel processing.
		#' @param	verbose			Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "incidence")
			private$assert_supported_design(seq_des_obj)
			super$initialize(seq_des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Not implemented for Zhang exact incidence inference.
		compute_treatment_estimate = function(){
			super$compute_treatment_estimate()
		},

		#' @description
		#' Not implemented for Zhang exact incidence inference.
		#' @param alpha Significance level.
		compute_mle_confidence_interval = function(alpha = 0.05){
			super$compute_mle_confidence_interval(alpha = alpha)
		},

		#' @description
		#' Not implemented for Zhang exact incidence inference.
		#' @param delta Null treatment effect to test.
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			super$compute_mle_two_sided_pval_for_treatment_effect(delta = delta)
		},

		#' @description
		#' Not implemented for Zhang exact incidence inference.
		#' @param ... Unused.
		compute_confidence_interval_rand = function(...){
			super$compute_confidence_interval_rand(...)
		},

		#' @description
		#' Not implemented for Zhang exact incidence inference.
		#' @param ... Unused.
		compute_two_sided_pval_for_treatment_effect_rand = function(...){
			super$compute_two_sided_pval_for_treatment_effect_rand(...)
		},

		#' @description
		#' Not implemented for Zhang exact incidence inference.
		#' @param ... Unused.
		approximate_bootstrap_distribution_beta_hat_T = function(...){
			super$approximate_bootstrap_distribution_beta_hat_T(...)
		},

		#' @description
		#' Not implemented for Zhang exact incidence inference.
		#' @param ... Unused.
		compute_bootstrap_confidence_interval = function(...){
			super$compute_bootstrap_confidence_interval(...)
		},

		#' @description
		#' Not implemented for Zhang exact incidence inference.
		#' @param ... Unused.
		compute_bootstrap_two_sided_pval = function(...){
			super$compute_bootstrap_two_sided_pval(...)
		},

		#' @description
		#' Computes the exact Zhang confidence interval for the log-odds treatment
		#' effect under CRD.
		#' @param alpha Significance level; the interval covers \code{1 - alpha}.
		#' @param pval_epsilon Bisection convergence tolerance.
		#' @param combination_method Combination rule for component p-values.
		compute_exact_confidence_interval = function(alpha = 0.05, pval_epsilon = 0.005, combination_method = "Fisher"){
			super$compute_exact_confidence_interval(
				alpha = alpha,
				pval_epsilon = pval_epsilon,
				combination_method = combination_method
			)
		},

		#' @description
		#' Computes the exact Zhang two-sided p-value for the log-odds treatment
		#' effect under CRD.
		#' @param delta Null treatment effect (log-odds ratio) to test.
		#' @param combination_method Combination rule for component p-values.
		compute_exact_two_sided_pval_for_treatment_effect = function(delta = 0, combination_method = "Fisher"){
			super$compute_exact_two_sided_pval_for_treatment_effect(
				delta = delta,
				combination_method = combination_method
			)
		}
	),

	private = list(

		assert_supported_design = function(seq_des_obj){
			if (!is(seq_des_obj, "SeqDesignCRD")){
				stop(class(self)[1], " requires a completely randomized design (SeqDesignCRD). Use SeqDesignInferenceIncidKKExactZhang for KK matching-on-the-fly designs.")
			}
		}
	)
)
