#' Non-parametric Wilcoxon-based Compound Inference for KK Designs
#'
#' @description
#' Fits a non-parametric compound estimator for KK matching-on-the-fly designs.
#' For matched pairs, it uses the Wilcoxon Signed-Rank Hodges-Lehmann estimate.
#' For reservoir subjects, it uses the Wilcoxon Rank-Sum (Mann-Whitney U) Hodges-Lehmann
#' estimate. The two estimates are combined via a variance-weighted linear combination.
#' This method is robust to outliers and does not assume a specific parametric
#' distribution for the response.
#'
#' @export
SeqDesignInferenceAllKKWilcox = R6::R6Class("SeqDesignInferenceAllKKWilcox",
	inherit = SeqDesignInferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param seq_des_obj		A SeqDesign object (must be a KK design).
		#' @param num_cores			Number of CPU cores for parallel processing.
		#' @param verbose			Whether to print progress messages.
		initialize = function(seq_des_obj, num_cores = 1, verbose = FALSE){
			res_type = seq_des_obj$get_response_type()
			if (res_type == "incidence"){
				stop("Rank-based compound inference is not recommended for incidence data; clogit or compound mean difference estimators are preferred.")
			}
			assertResponseType(res_type, c("continuous", "count", "proportion", "survival"))
			if (!is(seq_des_obj, "SeqDesignKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (SeqDesignKK14 or subclass).")
			}
			super$initialize(seq_des_obj, num_cores, verbose)
			if (private$any_censoring){
				stop(class(self)[1], " does not currently support censored survival data. Use restricted mean or Cox-based methods instead.")
			}
		},

		#' @description
		#' Returns the estimated treatment effect (Hodges-Lehmann median shift).
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes the non-parametric confidence interval.
		#' @param alpha					The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		compute_mle_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			# Even though estimates are non-parametric, the combined estimator is asymptotically normal
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes the non-parametric p-value.
		#' @param delta					The null difference to test against. For any treatment effect at all this is set to zero (the default).
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$assert_finite_se()
			if (delta == 0){
				private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			} else {
				stop("Testing non-zero delta is not yet implemented for the combined rank estimator.")
			}
		}
	),

	private = list(

		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			# --- Matched pairs: Wilcoxon Signed-Rank HL Estimate ---
			if (m > 0){
				private$rank_for_matched_pairs()
			}
			beta_m   = private$cached_values$beta_T_matched
			ssq_m    = private$cached_values$ssq_beta_T_matched
			m_ok     = !is.null(beta_m) && is.finite(beta_m) &&
			           !is.null(ssq_m)  && is.finite(ssq_m) && ssq_m > 0

			# --- Reservoir: Wilcoxon Rank-Sum HL Estimate ---
			if (nRT > 0 && nRC > 0){
				private$rank_for_reservoir()
			}
			beta_r   = private$cached_values$beta_T_reservoir
			ssq_r    = private$cached_values$ssq_beta_T_reservoir
			r_ok     = !is.null(beta_r) && is.finite(beta_r) &&
			           !is.null(ssq_r)  && is.finite(ssq_r) && ssq_r > 0

			# --- Variance-weighted combination ---
			if (m_ok && r_ok){
				w_star = ssq_r / (ssq_r + ssq_m)
				private$cached_values$beta_hat_T   = w_star * beta_m + (1 - w_star) * beta_r
				private$cached_values$s_beta_hat_T = sqrt(ssq_m * ssq_r / (ssq_m + ssq_r))
			} else if (m_ok){
				private$cached_values$beta_hat_T   = beta_m
				private$cached_values$s_beta_hat_T = sqrt(ssq_m)
			} else if (r_ok){
				private$cached_values$beta_hat_T   = beta_r
				private$cached_values$s_beta_hat_T = sqrt(ssq_r)
			} else {
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
			}
			private$cached_values$is_z = TRUE
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				stop("Rank compound estimator: could not compute a finite standard error.")
			}
		},

		rank_for_matched_pairs = function(){
			diffs = private$cached_values$KKstats$y_matched_diffs
			# signed-rank test requires at least some non-zero differences
			if (all(diffs == 0)) return(invisible(NULL))
			
			mod = tryCatch({
				stats::wilcox.test(diffs, conf.int = TRUE)
			}, error = function(e) NULL)
			
			if (is.null(mod)) return(invisible(NULL))
			
			beta = as.numeric(mod$estimate)
			# Back-calculate SE from 95% CI width
			se = (as.numeric(mod$conf.int[2]) - as.numeric(mod$conf.int[1])) / (2 * 1.96)
			
			private$cached_values$beta_T_matched     = if (is.finite(beta)) beta else NA_real_
			private$cached_values$ssq_beta_T_matched = if (is.finite(se) && se > 0) se^2 else NA_real_
		},

		rank_for_reservoir = function(){
			y_r = private$cached_values$KKstats$y_reservoir
			w_r = private$cached_values$KKstats$w_reservoir
			yT = y_r[w_r == 1]
			yC = y_r[w_r == 0]
			
			mod = tryCatch({
				stats::wilcox.test(yT, yC, conf.int = TRUE)
			}, error = function(e) NULL)
			
			if (is.null(mod)) return(invisible(NULL))
			
			beta = as.numeric(mod$estimate)
			# Back-calculate SE from 95% CI width
			se = (as.numeric(mod$conf.int[2]) - as.numeric(mod$conf.int[1])) / (2 * 1.96)
			
			private$cached_values$beta_T_reservoir     = if (is.finite(beta)) beta else NA_real_
			private$cached_values$ssq_beta_T_reservoir = if (is.finite(se) && se > 0) se^2 else NA_real_
		}
	)
)
