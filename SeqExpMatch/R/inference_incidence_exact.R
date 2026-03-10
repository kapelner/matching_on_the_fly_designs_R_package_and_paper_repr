#' Exact incidence inference via Zhang's combined reservoir test
#'
#' @description
#' Provides incidence-specific exact inference based on the Zhang combined method.
#' For CRD designs this uses the reservoir exact test only. Point estimates and
#' Wald-style anchor intervals use an unadjusted log-odds-ratio anchor and are used
#' only to center and bracket the exact CI inversion. KK-specific exact matched-pair
#' logic lives in \code{SeqDesignInferenceIncidKKExact}.
#'
#' @export
SeqDesignInferenceIncidExact = R6::R6Class("SeqDesignInferenceIncidExact",
	inherit = SeqDesignInferenceIncidExactAbstract,
	public = list(

		#' @description
		#' Initialize the exact incidence inference object.
		#' @param seq_des_obj		A SeqDesign object with an incidence response.
		#' @param num_cores			Number of CPU cores for parallel processing.
		#' @param verbose			Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "incidence")
			private$assert_supported_design(seq_des_obj)
			super$initialize(seq_des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Returns the anchor point estimate on the log-odds-ratio scale.
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes a Wald-style confidence interval used as an anchor/bracket for exact inversion.
		#' @param alpha					The confidence level in the computed confidence interval is 1 - \code{alpha}.
		compute_mle_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes a Wald-style two-sided p-value on the log-odds-ratio scale.
		#' @param delta					The null log-odds-ratio to test against.
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},

		#' @description
		#' Computes the Zhang combined exact p-value for an incidence effect.
		#' @param nsim_exact_test		Ignored; kept for interface compatibility.
		#' @param delta					The null log-odds-ratio to test against.
		#' @param transform_responses	Must be \code{"none"}.
		#' @param na.rm					Ignored; kept for interface compatibility.
		#' @param show_progress		Ignored; kept for interface compatibility.
		#' @param permutations		Ignored; kept for interface compatibility.
		#' @param combination_method	How to combine the matched-pair and reservoir p-values.
		compute_two_sided_pval_for_treatment_effect_rand = function(nsim_exact_test = 501, delta = 0, transform_responses = "none", na.rm = TRUE, show_progress = TRUE, permutations = NULL, combination_method = "Fisher"){
			assertNumeric(delta)
			if (!identical(transform_responses, "none")){
				stop("SeqDesignInferenceIncidExact only supports transform_responses = 'none'.")
			}
			combination_method = match.arg(combination_method, c("Fisher", "Stouffer", "min_p"))
			private$compute_combined_exact_pval(delta, combination_method)
		}
	),

	private = list(

		assert_supported_design = function(seq_des_obj){
			if (!is(seq_des_obj, "SeqDesignCRD")){
				stop(class(self)[1], " requires a completely randomized design (SeqDesignCRD). Use SeqDesignInferenceIncidKKExact for KK matching-on-the-fly designs.")
			}
		},

		compute_anchor_log_odds_ratio = function(){
			n11 = sum(private$y[private$w == 1L] == 1L, na.rm = TRUE)
			n10 = sum(private$y[private$w == 1L] == 0L, na.rm = TRUE)
			n01 = sum(private$y[private$w == 0L] == 1L, na.rm = TRUE)
			n00 = sum(private$y[private$w == 0L] == 0L, na.rm = TRUE)

			cells = c(n11, n10, n01, n00)
			if (any(cells == 0L)){
				cells = cells + 0.5
			}

			list(
				beta = log((cells[1] * cells[4]) / (cells[2] * cells[3])),
				se = sqrt(sum(1 / cells))
			)
		},

		shared = function(){
			if (!is.null(private$cached_values$is_z)) return(invisible(NULL))

			anchor = private$compute_anchor_log_odds_ratio()
			private$cached_values$beta_hat_T = anchor$beta
			private$cached_values$s_beta_hat_T = anchor$se
			private$cached_values$is_z = TRUE
		}
	)
)
