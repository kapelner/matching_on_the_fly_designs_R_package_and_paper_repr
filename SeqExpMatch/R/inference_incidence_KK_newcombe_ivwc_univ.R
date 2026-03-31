#' Univariate IVWC Newcombe Risk-Difference Inference for KK Designs
#'
#' Implements a compound Newcombe risk-difference estimator for KK designs.
#' This class pools information from matched pairs (using the Paired Newcombe
#' method) and the reservoir (using the Independent Newcombe method) via
#' inverse-variance weighted combination (IVWC).
#'
#' @details
#' The matched-pair component uses the discordant pairs to estimate the treatment
#' effect and its variance. The reservoir component treats unmatched subjects as
#' independent samples. The two estimates are combined using the standard IVWC
#' framework of the package.
#'
#' @export
#' @examples
#' set.seed(1)
#' x_dat <- data.frame(
#'   x1 = c(-1.2, -0.7, -0.2, 0.3, 0.8, 1.3, 1.8, 2.3),
#'   x2 = c(0, 1, 0, 1, 0, 1, 0, 1)
#' )
#' seq_des <- DesignSeqOneByOneKK14$
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
#' infer <- InferenceIncidUnivKKNewcombeRiskDiff$
#'   new(
#'   seq_des,
#'   verbose = FALSE
#' )
#' infer
#'
InferenceIncidUnivKKNewcombeRiskDiff = R6::R6Class("InferenceIncidUnivKKNewcombeRiskDiff",
	lock_objects = FALSE,
	inherit = InferenceKKPassThroughCompound,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj A completed KK \code{DesignSeqOneByOne} object.
		#' @param num_cores CPU cores.
		#' @param verbose Flag for progress messages.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			assertResponseType(des_obj$get_response_type(), "incidence")
			if (!is(des_obj, "DesignSeqOneByOneKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design.")
			}
			super$initialize(des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
		},

		#' @description
		#' Returns the IVWC Newcombe estimate.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared_combined()
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Returns the asymptotic confidence interval based on pooled variance.
		#' @param alpha The confidence level in the computed confidence
		#'   interval is 1 - \code{alpha}. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared_combined()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Returns the MLE p-value.
		#' @param delta The null risk difference to test against. Default is zero.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared_combined()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		compute_basic_match_data = function(){
			if (is.null(private$X)){
				private$X = private$get_X()
			}
			# Use the optimized Zhang helper to get counts
			private$cached_values$KKstats = compute_zhang_match_data_cpp(private$w, private$m, private$y, private$X)
		},

		pool_estimates_ivwc = function(est1, var1, est2, var2){
			ok1 = is.finite(est1) && is.finite(var1) && var1 > 0
			ok2 = is.finite(est2) && is.finite(var2) && var2 > 0

			if (ok1 && ok2){
				w1 = var2 / (var1 + var2)
				return(list(
					estimate = w1 * est1 + (1 - w1) * est2,
					variance = var1 * var2 / (var1 + var2)
				))
			} else if (ok1){
				return(list(estimate = est1, variance = var1))
			} else if (ok2){
				return(list(estimate = est2, variance = var2))
			} else {
				return(list(estimate = NA_real_, variance = NA_real_))
			}
		},

		shared_combined = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (is.null(private$cached_values$KKstats)) private$compute_basic_match_data()
			
			KKstats = private$cached_values$KKstats
			m = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			
			# Matched part
			est_m = NA_real_
			var_m = NA_real_
			if (m > 0){
				# Discordant counts from KKstats: d_plus (1,0) and d_minus (0,1)
				# n11 and n00 are also in KKstats
				n = m
				p10 = KKstats$d_plus / n
				p01 = KKstats$d_minus / n
				est_m = p10 - p01
				var_m = (p10 + p01 - (p10 - p01)^2) / n
			}
			
			# Reservoir part
			est_r = NA_real_
			var_r = NA_real_
			if (nRT > 0 && nRC > 0){
				pRT = KKstats$n11 / nRT
				pRC = KKstats$n01 / nRC
				est_r = pRT - pRC
				var_r = pRT * (1 - pRT) / nRT + pRC * (1 - pRC) / nRC
			}
			
			# IVWC Pooling
			res = private$pool_estimates_ivwc(est_m, var_m, est_r, var_r)
			
			private$cached_values$beta_hat_T = res$estimate
			private$cached_values$s_beta_hat_T = sqrt(res$variance)
			private$cached_values$is_z = TRUE
			private$cached_values$df = NA_real_
		}
	)
)
