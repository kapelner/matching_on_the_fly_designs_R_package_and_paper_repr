#' KK Newcombe Risk-Difference IVWC Inference for Binary Responses
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
#' @examples
#' \donttest{
#' seq_des = DesignSeqOneByOneKK14$new(n = 10, response_type = 'incidence')
#' for (i in 1:10) {
#'   seq_des$add_one_subject_to_experiment_and_assign(data.frame(x1 = rnorm(1), x2 = rnorm(1)))
#' }
#' seq_des$add_all_subject_responses(rbinom(10, 1, 0.5))
#' inf = InferenceIncidKKNewcombeRiskDiff$new(seq_des)
#' inf$compute_estimate()
#' }
#' @export
InferenceIncidKKNewcombeRiskDiff = R6::R6Class("InferenceIncidKKNewcombeRiskDiff",
	lock_objects = FALSE,
	inherit = InferenceKKPassThroughCompoundNoParamBootstrap,
	public = utils::modifyList(list(
		#' @description Initialize the inference object.
		#' @param des_obj A completed KK \code{DesignSeqOneByOne} object.
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose Flag for progress messages.
		#' @param smart_cold_start_default Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL,  verbose = FALSE, smart_cold_start_default = TRUE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
			}
			if (should_run_asserts()) {
				if (!inherits(des_obj, "DesignSeqOneByOneKK14") && !inherits(des_obj, "DesignFixedBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design.")
				}
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},
		#' @description Weighted bootstrap estimate using a direct empirical risk-difference surrogate.
		#' @param subject_or_block_weights bootstrap weights at the subject/block level.
		#' @param estimate_only flag.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE) {
			row_weights = private$expand_subject_or_block_weights_to_row_weights(subject_or_block_weights)
			if (weights_are_effectively_constant(row_weights)) {
				beta_hat_T = as.numeric(private$weighted_empirical_risk_difference(row_weights))[1L]
			} else {
				beta_hat_T = as.numeric(private$weighted_empirical_risk_difference(row_weights))[1L]
			}
			private$cached_values$beta_hat_T = beta_hat_T
			private$cached_values$s_beta_hat_T = NA_real_
			private$cached_values$df = NA_real_
			private$cached_values$beta_hat_T
		}
	), list()), # empty list for modifyList
	private = list(
		compute_basic_match_data = function(){
			# Use the optimized Zhang helper to get counts
			private$cached_values$KKstats = compute_zhang_match_data_cpp(private$get_X(), private$y, private$w, private$m)		},
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
		weighted_empirical_risk_difference = function(row_weights){
			ok_t = is.finite(private$w) & private$w == 1 & is.finite(row_weights) & row_weights > 0
			ok_c = is.finite(private$w) & private$w == 0 & is.finite(row_weights) & row_weights > 0
			if (!any(ok_t) || !any(ok_c)) return(NA_real_)
			p_t = stats::weighted.mean(as.numeric(private$y[ok_t]), as.numeric(row_weights[ok_t]))
			p_c = stats::weighted.mean(as.numeric(private$y[ok_c]), as.numeric(row_weights[ok_c]))
			if (!is.finite(p_t) || !is.finite(p_c)) return(NA_real_)
			p_t - p_c
		},
		shared_combined = function(){
			if (!isTRUE(private$has_match_structure)) {
				private$cache_nonestimable_estimate("kk_design_required")
				return(invisible(NULL))
			}
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (is.null(private$cached_values$KKstats)) private$compute_basic_match_data()
			if (is.null(private$cached_values$KKstats)) return(invisible(NULL))
			
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
			private$cached_values$df = NA_real_
		}
	)
)
