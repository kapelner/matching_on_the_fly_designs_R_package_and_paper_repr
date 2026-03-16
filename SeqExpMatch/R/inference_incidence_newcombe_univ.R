#' Univariate Newcombe Risk-Difference Inference for Binary Responses
#'
#' @description
#' Fits the Newcombe hybrid score method (Method 10) for the risk difference in a
#' two-arm binary trial. This method constructs a confidence interval for the
#' difference between two independent proportions by combining Wilson score intervals
#' for each group.
#'
#' @details
#' This class is unadjusted and assumes independent samples (e.g. from a CRD). It
#' ignores any matched-pair structure if present. For matched data, use
#' \code{SeqDesignInferenceIncidUnivKKNewcombeRiskDiff}.
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
#'   response_type = "incidence",
#'   verbose = FALSE
#' )
#' for (i in seq_len(nrow(x_dat))) {
#'   seq_des$
#'   add_subject_to_experiment_and_assign(x_dat[i, , drop = FALSE])
#' }
#' seq_des$
#'   add_all_subject_responses(c(0, 1, 0, 1, 0, 1, 1, 0))
#' infer <- SeqDesignInferenceIncidUnivNewcombeRiskDiff$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
SeqDesignInferenceIncidUnivNewcombeRiskDiff = R6::R6Class("SeqDesignInferenceIncidUnivNewcombeRiskDiff",
	inherit = SeqDesignInference,
	public = list(

		#' @description
		#' Initialize a Newcombe risk-difference inference object.
		#' @param seq_des_obj A completed \code{SeqDesign} object with an incidence response.
		#' @param num_cores The number of CPU cores to use.
		#' @param verbose Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "incidence")
			super$initialize(seq_des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Computes the observed risk-difference estimate.
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes a 1 - \code{alpha} Newcombe confidence interval.
		#' @param alpha The significance level.
		compute_mle_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			counts = private$cached_values$counts
			if (is.null(counts)) return(c(NA_real_, NA_real_))
			
			ci = newcombe_independent_ci_cpp(counts$x_t, counts$n_t, counts$x_c, counts$n_c, alpha)
			names(ci) = paste0(c(alpha / 2, 1 - alpha / 2) * 100, "%")
			ci
		},

		#' @description
		#' Computes a two-sided p-value by inverting the Newcombe interval.
		#' @param delta The null risk difference.
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta, len = 1)
			private$shared()
			counts = private$cached_values$counts
			if (is.null(counts) || counts$n_t == 0 || counts$n_c == 0) return(NA_real_)
			
			# Invert the Newcombe CI to find the largest alpha such that delta is on the boundary
			p_fn = function(a) {
				ci = newcombe_independent_ci_cpp(counts$x_t, counts$n_t, counts$x_c, counts$n_c, a)
				if (delta < private$cached_values$beta_hat_T) ci[1] - delta else ci[2] - delta
			}
			
			res = tryCatch(stats::uniroot(p_fn, interval = c(1e-10, 1 - 1e-10))$root, error = function(e) NA_real_)
			if (!is.finite(res)) return(1.0) # If root not found, delta likely far inside
			res
		}
	),

	private = list(
		get_counts = function(){
			i_t = private$w == 1
			i_c = private$w == 0
			n_t = sum(i_t)
			n_c = sum(i_c)
			x_t = sum(private$y[i_t])
			x_c = sum(private$y[i_c])

			list(
				n_t = n_t, n_c = n_c,
				x_t = x_t, x_c = x_c,
				p_t = if (n_t > 0) x_t / n_t else NA_real_,
				p_c = if (n_c > 0) x_c / n_c else NA_real_
			)
		},

		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			counts = private$get_counts()
			private$cached_values$counts = counts
			private$cached_values$beta_hat_T = counts$p_t - counts$p_c
			private$cached_values$is_z = TRUE
		}
	)
)
