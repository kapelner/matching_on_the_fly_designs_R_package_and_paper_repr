#' Abstract class for Conditional Poisson Combined-Likelihood Compound Inference
#'
#' Fits a single joint likelihood over all KK design data for count responses.
#' The matched-pair component uses the conditional Poisson likelihood; conditioning
#' on each pair's total count n_k = yT_k + yC_k reduces it to a weighted binomial
#' logistic form: yT_k | n_k ~ Binomial(n_k, p_k) with logit(p_k) = beta_T + X_diff_k' beta_xs.
#' The reservoir component uses the marginal Poisson likelihood.
#' Both components share beta_T and beta_xs; beta_0 (intercept) is only identified
#' from reservoir rows and is profiled out in the pairs-only case.
#'
#' Because the two components use different link functions (logit vs log), they
#' cannot be stacked into a single GLM.  Instead, the combined log-likelihood
#' L_total = L_cond_Poisson(pairs) + L_Poisson(reservoir)
#' is maximized jointly via BFGS with an analytic gradient.
#' The standard error for beta_T is extracted from the numerical Hessian at the MLE.
#'
#' Degenerate cases:
#' Pairs only     : reduces to weighted logistic (glm.fit), no BFGS needed.
#' Reservoir only : reduces to Poisson regression (glm.fit), no BFGS needed.
#'
#' @keywords internal
InferenceAbstractKKPoissonCPoissonCombinedLikelihood = R6::R6Class("InferenceAbstractKKPoissonCPoissonCombinedLikelihood",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize combined-likelihood conditional-Poisson inference for KK count data.
		#' @param des_obj A completed KK design object.
		#' @param verbose Whether to print progress messages.
		#' @return A new inference object.
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
		#' Compute the treatment estimate.
		#' @param estimate_only Whether to skip standard-error calculations.
		#' @return The treatment estimate.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared_combined_likelihood(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Compute an asymptotic confidence interval.
		#' @param alpha Significance level.
		#' @return A confidence interval.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared_combined_likelihood(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Compute an asymptotic two-sided p-value for the treatment effect.
		#' @param delta Null treatment effect value.
		#' @return A two-sided p-value.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared_combined_likelihood(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		compute_fast_randomization_distr = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = .Machine$double.eps){
			private$compute_fast_randomization_distr_via_reused_worker(y, permutations, delta, transform_responses, zero_one_logit_clamp = zero_one_logit_clamp)
		},

		# Abstract: subclasses return TRUE (multivariate) or FALSE (univariate).
		include_covariates = function() stop(class(self)[1], " must implement include_covariates()"),

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},

		set_failed_combined_cache = function(){
			private$cache_nonestimable_estimate("kk_cpoisson_combined_fit_unavailable")
			invisible(FALSE)
		},

			reduce_combined_covariates = function(X_diff_v, X_r_v, w_r_v){
				cached_keep = private$cached_values$combined_cov_keep
				p = ncol(X_r_v)
				if (is.null(p) || p == 0L) return(integer(0))

				build_combined = function(){
					nd = nrow(X_diff_v)
					nR = nrow(X_r_v)
					X_pairs = cbind(
						matrix(0, nrow = nd, ncol = 2L),
						X_diff_v
					)
					if (nR == 0L){
						X_res = matrix(0, nrow = 0L, ncol = ncol(X_pairs))
					} else {
						X_res = cbind(rep(1, nR), w_r_v, X_r_v)
						if (ncol(X_res) != ncol(X_pairs)) {
							max_cols = max(ncol(X_pairs), ncol(X_res))
							if (ncol(X_pairs) < max_cols) {
								X_pairs = cbind(X_pairs, matrix(0, nrow = nd, ncol = max_cols - ncol(X_pairs)))
							}
							if (ncol(X_res) < max_cols) {
								X_res = cbind(X_res, matrix(0, nrow = nR, ncol = max_cols - ncol(X_res)))
							}
						}
					}
					rbind(X_pairs, X_res)
				}

				X_full = build_combined()
				required = 2L

				if (!is.null(cached_keep) && length(cached_keep) > 0L){
					valid_cached = cached_keep[cached_keep >= 1L & cached_keep <= p]
					if (length(valid_cached) > 0L){
						keep_full = sort(unique(c(required, valid_cached + 2L)))
						keep_full = keep_full[keep_full <= ncol(X_full)]
						if (length(keep_full) > 0L){
							X_try = X_full[, keep_full, drop = FALSE]
							qr_try = qr(X_try)
							if (qr_try$rank == ncol(X_try) && (2L %in% keep_full)){
								private$cached_values$combined_cov_keep = sort(unique(valid_cached))
								return(private$cached_values$combined_cov_keep)
							}
						}
					}
					private$cached_values$combined_cov_keep = NULL
				}

				qr_full = qr(X_full)
				r_full = qr_full$rank
				if (r_full <= 2L){
					private$cached_values$combined_cov_keep = integer(0)
					return(integer(0))
				}
				keep_full = qr_full$pivot[seq_len(r_full)]
				if (!(required %in% keep_full)) keep_full[r_full] = required
				keep_full = sort(unique(keep_full))
				keep_cov = keep_full[keep_full > 2L] - 2L
				keep_cov = keep_cov[keep_cov >= 1L & keep_cov <= p]
				private$cached_values$combined_cov_keep = keep_cov
				keep_cov
			},

		try_combined_fit = function(estimate_only, yT_v, n_k_v, X_diff_v, y_r_v, w_r_v, X_r_v){
			mod = tryCatch(
				fast_cpoisson_combined_with_var_cpp(
					yT_v = as.numeric(yT_v),
					n_k_v = as.numeric(n_k_v),
					X_diff_v = as.matrix(X_diff_v),
					y_r = as.numeric(y_r_v),
					w_r = as.numeric(w_r_v),
					X_r = as.matrix(X_r_v)
				),
				error = function(e) NULL
			)
			if (is.null(mod) || length(mod$b) < 2L || !is.finite(mod$b[2])) return(FALSE)
			if (!estimate_only && (!is.finite(mod$ssq_b_j) || mod$ssq_b_j < 0)) return(FALSE)

			private$cached_values$beta_hat_T = as.numeric(mod$b[2])
			if (!estimate_only) private$cached_values$s_beta_hat_T = sqrt(as.numeric(mod$ssq_b_j))
			private$cached_values$is_z = TRUE
			TRUE
		},

		try_pairs_only = function(estimate_only, yT_v, n_k_v, X_diff_v){
			if (length(yT_v) == 0L) return(FALSE)
			y_prop = yT_v / n_k_v
			Xmm = if (ncol(X_diff_v) > 0L) cbind(1, X_diff_v) else matrix(1, nrow = length(yT_v), ncol = 1L)
			mod = tryCatch({
				res = fast_logistic_regression_weighted_cpp(X = Xmm, y = y_prop, weights = n_k_v)
				if (estimate_only) {
					list(b = res$b, ssq_b_j = NA_real_)
				} else {
					vcov = solve(res$XtWX)
					list(b = res$b, ssq_b_j = vcov[1, 1])
				}
			}, error = function(e) NULL)
			if (is.null(mod) || !is.finite(mod$b[1])) return(FALSE)
			if (!estimate_only && (!is.finite(mod$ssq_b_j) || mod$ssq_b_j < 0)) return(FALSE)

			private$cached_values$beta_hat_T = as.numeric(mod$b[1])
			if (!estimate_only) private$cached_values$s_beta_hat_T = sqrt(as.numeric(mod$ssq_b_j))
			private$cached_values$is_z = TRUE
			TRUE
		},

		try_reservoir_only = function(estimate_only, y_r_v, w_r_v, X_r_v){
			if (length(y_r_v) == 0L) return(FALSE)
			X_full = if (ncol(X_r_v) > 0L) cbind(1, w_r_v, X_r_v) else cbind(1, w_r_v)
			j_treat = 2L
			if (ncol(X_full) > 2L){
				qr_full = qr(X_full)
				r_full = qr_full$rank
				if (r_full < ncol(X_full)){
					keep = qr_full$pivot[seq_len(r_full)]
					if (!(2L %in% keep)) keep[r_full] = 2L
					keep = sort(unique(keep))
					X_full = X_full[, keep, drop = FALSE]
					j_treat = which(keep == 2L)
					if (length(j_treat) != 1L) return(FALSE)
				}
			}
			mod = tryCatch({
				if (estimate_only) {
					list(b = fast_poisson_regression_cpp(X_full, as.numeric(y_r_v))$b, ssq_b_j = NA_real_)
				} else {
					fast_poisson_regression_with_var_cpp(X_full, as.numeric(y_r_v), j = j_treat)
				}
			}, error = function(e) NULL)
			if (is.null(mod) || !is.finite(mod$b[j_treat])) return(FALSE)
			if (!estimate_only && (!is.finite(mod$ssq_b_j) || mod$ssq_b_j < 0)) return(FALSE)

			private$cached_values$beta_hat_T = as.numeric(mod$b[j_treat])
			if (!estimate_only) private$cached_values$s_beta_hat_T = sqrt(as.numeric(mod$ssq_b_j))
			private$cached_values$is_z = TRUE
			TRUE
		},

		# Maximize L_cond_Poisson(pairs) + L_Poisson(reservoir) jointly.
		#
		# Pair contribution (conditional Poisson = weighted binomial logistic):
		#   eta_k = beta_T + X_diff_k' beta_xs
		#   ll_pairs = sum(yT_k * eta_k - n_k * log(1 + exp(eta_k)))
		#
		# Reservoir contribution (marginal Poisson):
		#   eta_i = beta_0 + w_i * beta_T + X_i' beta_xs
		#   ll_res = sum(y_i * eta_i - exp(eta_i))
		#
		# Parameter layout: [beta_0, beta_T, beta_xs (p)] when reservoir present;
		#                   [beta_T, beta_xs (p)]           when pairs-only.
		# The combined case is handled by fast_cpoisson_combined_with_var_cpp
		# (Newton's method with analytic Fisher-information Hessian).
		shared_combined_likelihood = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
				private$cached_values$combined_cov_keep = NULL
			}
			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			p             = if (private$include_covariates()) ncol(private$get_X()) else 0L
			has_reservoir = nRT > 0 && nRC > 0

			# ---- Pair data (conditional Poisson, zero-total pairs discarded) ----
			yT_v     = numeric(0)
			n_k_v    = numeric(0)
			X_diff_v = matrix(nrow = 0L, ncol = p)

			if (m > 0){
				yT    = KKstats$yTs_matched
				yC    = KKstats$yCs_matched
				n_k   = yT + yC
				valid = which(n_k > 0)
				if (length(valid) > 0){
					yT_v  = yT[valid]
					n_k_v = n_k[valid]
					if (p > 0L) {
						# Use the full-width pair-difference matrix here so any later
						# covariate reduction stays aligned with the reservoir matrix.
						X_diff_v = as.matrix(KKstats$X_matched_diffs_full[valid, , drop = FALSE])
					}
				}
			}
			has_pairs = length(yT_v) > 0

			# ---- Reservoir data (marginal Poisson) ----
			y_r_v = numeric(0)
			w_r_v = numeric(0)
			X_r_v = matrix(nrow = 0L, ncol = p)
			if (has_reservoir){
				y_r_v = KKstats$y_reservoir
				w_r_v = KKstats$w_reservoir
				if (p > 0L) X_r_v = as.matrix(KKstats$X_reservoir)
			}

			if (!has_pairs && !has_reservoir){
				private$set_failed_combined_cache()
				return(invisible(NULL))
			}

			if (has_pairs && has_reservoir && p > 0L){
				keep_cov = private$reduce_combined_covariates(X_diff_v, X_r_v, w_r_v)
				if (length(keep_cov) > 0L){
					X_diff_v = X_diff_v[, keep_cov, drop = FALSE]
					X_r_v    = X_r_v[, keep_cov, drop = FALSE]
				} else {
					X_diff_v = matrix(nrow = length(yT_v), ncol = 0L)
					X_r_v    = matrix(nrow = length(y_r_v), ncol = 0L)
				}
			}

			if (has_pairs && has_reservoir){
				if (!private$try_combined_fit(estimate_only, yT_v, n_k_v, X_diff_v, y_r_v, w_r_v, X_r_v)){
					if (private$verbose) message(class(self)[1], ": combined CP/Poisson fit failed; falling back to simpler components.")
					fallback_success = FALSE
					if (has_pairs){
						fallback_success = private$try_pairs_only(estimate_only, yT_v, n_k_v, X_diff_v)
					}
					if (!fallback_success && has_reservoir){
						fallback_success = private$try_reservoir_only(estimate_only, y_r_v, w_r_v, X_r_v)
					}
					if (!fallback_success){
						private$set_failed_combined_cache()
					}
				}
				return(invisible(NULL))
			}

			if (has_pairs){
				if (!private$try_pairs_only(estimate_only, yT_v, n_k_v, X_diff_v)){
					private$set_failed_combined_cache()
				}
			} else if (has_reservoir){
				if (!private$try_reservoir_only(estimate_only, y_r_v, w_r_v, X_r_v)){
					private$set_failed_combined_cache()
				}
			}

			invisible(NULL)
		}
	)
)
