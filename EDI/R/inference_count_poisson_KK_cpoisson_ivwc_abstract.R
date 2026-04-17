#' Abstract class for Conditional Poisson / Negative Binomial Compound Inference
#'
#' This class implements a compound estimator for KK matching-on-the-fly designs with
#' count responses. For matched pairs, it uses conditional Poisson regression (implemented
#' via binomial logistic regression on the differences of covariates). For reservoir
#' subjects, it uses Negative Binomial regression. The two estimates are combined via a
#' variance-weighted linear combination.
#'
#' @keywords internal
InferenceAbstractKKPoissonCPoissonIVWC = R6::R6Class("InferenceAbstractKKPoissonCPoissonIVWC",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj		A DesignSeqOneByOne object (must be a KK design).
		#' @param verbose			Whether to print progress messages.
		initialize = function(des_obj,  verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
			}
			if (should_run_asserts()) {
				if (!is(des_obj, "DesignSeqOneByOneKK14") && !is(des_obj, "FixedDesignBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
				}
			}
			super$initialize(des_obj, verbose)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},

		#' @description
		#' Returns the estimated treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes the asymptotic confidence interval.
		#' @param alpha                                   The confidence level in the computed
		#'   confidence interval is 1 - \code{alpha}. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes the asymptotic p-value.
		#' @param delta                                   The null difference to test against. For
		#'   any treatment effect at all this is set to zero (the default).
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			if (delta == 0){
				private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			} else {
				if (should_run_asserts()) {
					stop("Testing non-zero delta is not yet implemented for this class.")
				}
				NA_real_
			}
		}
	),

	private = list(
		best_Xmm_colnames_matched = NULL,
		best_Xmm_colnames_reservoir = NULL,
		best_Xmm_j_treat_reservoir = NULL,

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Ensure we have the best design from the original data
			if (is.null(private$best_Xmm_colnames_matched) && is.null(private$best_Xmm_colnames_reservoir)){
				private$shared(estimate_only = TRUE)
			}

			if (is.null(private$cached_values$KKstats)) private$compute_basic_match_data()
			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			# --- Matched pairs component ---
			beta_m = NA_real_
			if (m > 0){
				yT = KKstats$yTs_matched
				yC = KKstats$yCs_matched
				y_total = yT + yC
				valid_idx = which(y_total > 0)
				if (length(valid_idx) > 0){
					y_prop = yT[valid_idx] / y_total[valid_idx]
					weights = y_total[valid_idx]
					
					if (private$include_covariates() && !is.null(private$best_Xmm_colnames_matched)){
						X_diff = KKstats$X_matched_diffs[valid_idx, intersect(private$best_Xmm_colnames_matched, colnames(KKstats$X_matched_diffs)), drop = FALSE]
						Xmm = cbind(1, X_diff)
					} else {
						Xmm = matrix(1, nrow = length(valid_idx), ncol = 1)
					}
					
					mod_m = tryCatch(fast_logistic_regression_weighted_cpp(X = Xmm, y = y_prop, weights = weights), error = function(e) NULL)
					if (!is.null(mod_m) && is.finite(mod_m$b[1])) beta_m = as.numeric(mod_m$b[1])
				}
			}

			# --- Reservoir component ---
			beta_r = NA_real_
			if (nRT > 0 && nRC > 0){
				y_r = KKstats$y_reservoir
				w_r = KKstats$w_reservoir
				
				if (private$include_covariates() && !is.null(private$best_Xmm_colnames_reservoir)){
					X_r = as.matrix(KKstats$X_reservoir)
					X_cov = X_r[, intersect(private$best_Xmm_colnames_reservoir, colnames(X_r)), drop = FALSE]
					Xmm = cbind(1, w_r, X_cov)
					j_treat = private$best_Xmm_j_treat_reservoir %||% 2L
				} else {
					Xmm = cbind(1, w_r)
					j_treat = 2L
				}
				
				mod_r = tryCatch(fast_negbin_regression(Xmm, y_r), error = function(e) NULL)
				if (!is.null(mod_r) && length(mod_r$b) >= j_treat && is.finite(mod_r$b[j_treat])) beta_r = as.numeric(mod_r$b[j_treat])
			}

			# Inverse-variance weighted pooling (using original variances for weights)
			m_ok = is.finite(beta_m)
			r_ok = is.finite(beta_r)

			if (m_ok && r_ok){
				ssq_m_orig = private$cached_values$ssq_beta_T_matched
				ssq_r_orig = private$cached_values$ssq_beta_T_reservoir
				if (is.finite(ssq_m_orig) && is.finite(ssq_r_orig)){
					w_star = ssq_r_orig / (ssq_r_orig + ssq_m_orig)
					return(w_star * beta_m + (1 - w_star) * beta_r)
				}
				return((beta_m + beta_r) / 2)
			} else if (m_ok){
				return(beta_m)
			} else if (r_ok){
				return(beta_r)
			}
			NA_real_
		},

		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			private$compute_fast_randomization_distr_via_reused_worker(y, permutations, delta, transform_responses, zero_one_logit_clamp = zero_one_logit_clamp)
		},

		# Abstract: subclasses return TRUE (multivariate) or FALSE (univariate).
		include_covariates = function() stop(class(self)[1], " must implement include_covariates()"),

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			# Recompute KKstats if cache was cleared (e.g., after y transformation for rand CI)
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}

			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			# --- Matched pairs: Conditional Poisson via Binomial Trick ---
			if (m > 0){
				private$cpoisson_for_matched_pairs(estimate_only = estimate_only)
			}
			beta_m   = private$cached_values$beta_T_matched
			ssq_m    = private$cached_values$ssq_beta_T_matched
			m_ok     = !is.null(beta_m) && is.finite(beta_m) &&
			           (estimate_only || (!is.null(ssq_m) && is.finite(ssq_m) && ssq_m > 0))

			# --- Reservoir: Negative Binomial ---
			if (nRT > 0 && nRC > 0){
				private$negbin_for_reservoir(estimate_only = estimate_only)
			}
			beta_r   = private$cached_values$beta_T_reservoir
			ssq_r    = private$cached_values$ssq_beta_T_reservoir
			r_ok     = !is.null(beta_r) && is.finite(beta_r) &&
			           (estimate_only || (!is.null(ssq_r) && is.finite(ssq_r) && ssq_r > 0))

			# --- Variance-weighted combination ---
			if (m_ok && r_ok){
				if (estimate_only) {
					# Simple mean if SEs not available
					private$cached_values$beta_hat_T = (beta_m + beta_r) / 2
					return(invisible(NULL))
				}
				w_star = ssq_r / (ssq_r + ssq_m)
				private$cached_values$beta_hat_T   = w_star * beta_m + (1 - w_star) * beta_r
				private$cached_values$s_beta_hat_T = sqrt(ssq_m * ssq_r / (ssq_m + ssq_r))
			} else if (m_ok){
				private$cached_values$beta_hat_T   = beta_m
				if (!estimate_only) private$cached_values$s_beta_hat_T = sqrt(ssq_m)
			} else if (r_ok){
				private$cached_values$beta_hat_T   = beta_r
				if (!estimate_only) private$cached_values$s_beta_hat_T = sqrt(ssq_r)
			} else {
				private$cache_nonestimable_estimate("kk_cpoisson_ivwc_no_usable_component")
			}
			private$cached_values$is_z = TRUE
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},

		cpoisson_for_matched_pairs = function(estimate_only = FALSE){
			KKstats = private$cached_values$KKstats
			yT = KKstats$yTs_matched
			yC = KKstats$yCs_matched

			# Filter pairs where total count is zero (provide no information for the conditional likelihood)
			y_total = yT + yC
			valid_idx = which(y_total > 0)
			if (length(valid_idx) == 0) return(invisible(NULL))

			y_prop = yT[valid_idx] / y_total[valid_idx]
			weights = y_total[valid_idx]

			if (private$include_covariates()){
				# Use differences of covariates as predictors
				X_diff = KKstats$X_matched_diffs[valid_idx, , drop = FALSE]
				# Ensure treatment effect is the intercept (since treatment diff is always 1)
				Xmm = cbind(1, X_diff)
			} else {
				Xmm = matrix(1, nrow = length(valid_idx), ncol = 1)
			}

			# Fit binomial logistic regression
			mod = tryCatch({
				if (estimate_only) {
					res = fast_logistic_regression_weighted_cpp(X = Xmm, y = y_prop, weights = weights)
					list(b = res$b, ssq_b_1 = NA_real_, X_fit = Xmm)
				} else {
					res = fast_logistic_regression_weighted_cpp(X = Xmm, y = y_prop, weights = weights)
					vcov = solve(res$XtWX)
					list(b = res$b, ssq_b_1 = vcov[1, 1], X_fit = Xmm)
				}
			}, error = function(e) NULL)

			if (is.null(mod)) return(invisible(NULL))

			private$cached_values$beta_T_matched     = as.numeric(mod$b[1])
			if (private$include_covariates()){
				private$best_Xmm_colnames_matched = setdiff(colnames(mod$X_fit), "(Intercept)")
			}
			if (!estimate_only) private$cached_values$ssq_beta_T_matched = as.numeric(mod$ssq_b_1)
		},

		negbin_for_reservoir = function(estimate_only = FALSE){
			y_r    = private$cached_values$KKstats$y_reservoir
			w_r    = private$cached_values$KKstats$w_reservoir
			X_r    = as.matrix(private$cached_values$KKstats$X_reservoir)
			j_treat = 2L

			if (private$include_covariates()){
				X_full = cbind(1, w_r, X_r)
				if (!estimate_only) {
					qr_full = qr(X_full)
					r_full  = qr_full$rank
					if (r_full < ncol(X_full)){
						keep = qr_full$pivot[seq_len(r_full)]
						if (!(2L %in% keep)) keep[r_full] = 2L
						keep    = sort(keep)
						X_full  = X_full[, keep, drop = FALSE]
						j_treat = which(keep == 2L)
					}
				}
			} else {
				X_full = cbind(1, w_r)
			}

			mod = tryCatch({
				if (estimate_only) {
					list(b = fast_negbin_regression(X_full, y_r)$b, ssq_b_j = NA_real_, X_fit = X_full)
				} else {
					res = fast_negbin_regression_with_var(X_full, y_r, j = j_treat)
					res$X_fit = X_full
					res
				}
			}, error = function(e) NULL)
			if (is.null(mod)) return(invisible(NULL))

			private$cached_values$beta_T_reservoir     = as.numeric(mod$b[j_treat])
			if (private$include_covariates()){
				private$best_Xmm_colnames_reservoir = setdiff(colnames(mod$X_fit), c("(Intercept)", "w_r"))
				private$best_Xmm_j_treat_reservoir = j_treat
			}
			if (!estimate_only) private$cached_values$ssq_beta_T_reservoir = as.numeric(mod$ssq_b_j)
		}
	)
)
