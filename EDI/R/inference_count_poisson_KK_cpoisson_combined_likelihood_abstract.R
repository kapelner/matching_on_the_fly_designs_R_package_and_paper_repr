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
			assertResponseType(des_obj$get_response_type(), "count")
			if (!is(des_obj, "DesignSeqOneByOneKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
			}
			super$initialize(des_obj, verbose)
			assertNoCensoring(private$any_censoring)
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
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared_combined_likelihood(estimate_only = FALSE)
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Compute an asymptotic two-sided p-value for the treatment effect.
		#' @param delta Null treatment effect value.
		#' @return A two-sided p-value.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared_combined_likelihood(estimate_only = FALSE)
			private$assert_finite_se()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(

		# Abstract: subclasses return TRUE (multivariate) or FALSE (univariate).
		include_covariates = function() stop(class(self)[1], " must implement include_covariates()"),

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
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
					if (p > 0L) X_diff_v = as.matrix(KKstats$X_matched_diffs[valid, , drop = FALSE])
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
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			# ---- Case 1: pairs + reservoir — joint BFGS optimization --------
			if (has_pairs && has_reservoir){
				if (p > 0L) {
					# Use the full-width pair differences from the shared C++ preprocessing
					# rather than reconstructing them in R.
					Xd_full = as.matrix(KKstats$X_matched_diffs_full)
					X_diff_v = as.matrix(Xd_full[valid, , drop = FALSE])

					if (!estimate_only) {
						# QR-reduce the combined design to full rank while preserving the
						# intercept and treatment columns required by the C++ parameterization.
						X_rank = rbind(
							cbind(0, 1, X_diff_v),
							cbind(1, w_r_v, X_r_v)
						)
						keep = 1:2
						rank_keep = qr(X_rank[, keep, drop = FALSE])$rank
						if (rank_keep < length(keep)) {
							private$cached_values$beta_hat_T   = NA_real_
							private$cached_values$s_beta_hat_T = NA_real_
							private$cached_values$is_z         = TRUE
							return(invisible(NULL))
						}
						for (j in 3:ncol(X_rank)) {
							trial_keep = c(keep, j)
							trial_rank = qr(X_rank[, trial_keep, drop = FALSE])$rank
							if (trial_rank > rank_keep) {
								keep = trial_keep
								rank_keep = trial_rank
							}
						}
						keep_cov = keep[keep > 2L] - 2L
						X_diff_v = if (length(keep_cov) > 0L) {
							X_diff_v[, keep_cov, drop = FALSE]
						} else {
							matrix(nrow = nrow(X_diff_v), ncol = 0L)
						}
						X_r_v = if (length(keep_cov) > 0L) {
							X_r_v[, keep_cov, drop = FALSE]
						} else {
							matrix(nrow = nrow(X_r_v), ncol = 0L)
						}
					}
				}

				# Combined: Newton's method in C++
				# params layout: [beta_0, beta_T, beta_xs]; beta_T is at 1-based index 2
				mod = tryCatch(
					if (estimate_only) {
						fast_cpoisson_combined_cpp(
							as.double(yT_v), as.double(n_k_v), X_diff_v,
							as.double(y_r_v), as.double(w_r_v), X_r_v
						)
					} else {
						fast_cpoisson_combined_with_var_cpp(
							as.double(yT_v), as.double(n_k_v), X_diff_v,
							as.double(y_r_v), as.double(w_r_v), X_r_v
						)
					},
					error = function(e) NULL
				)
				if (is.null(mod) || !mod$converged){
					private$cached_values$beta_hat_T   = NA_real_
					if (!estimate_only) private$cached_values$s_beta_hat_T = NA_real_
					private$cached_values$is_z         = TRUE
					return(invisible(NULL))
				}
				private$cached_values$beta_hat_T   = as.numeric(mod$b[2L])
				if (!estimate_only) private$cached_values$s_beta_hat_T = sqrt(mod$ssq_b_j)
				private$cached_values$is_z         = TRUE

			# ---- Case 2: pairs only — weighted logistic (C++) -------------------
			} else if (has_pairs){
				y_prop = yT_v / n_k_v
				Xmm    = if (p > 0L) cbind(1, X_diff_v) else matrix(1, nrow = length(yT_v), ncol = 1)
				mod = tryCatch({
					if (estimate_only) {
						res = fast_logistic_regression_weighted_cpp(X = Xmm, y = y_prop, weights = n_k_v)
						list(b = res$b, ssq_b_1 = NA_real_)
					} else {
						res = fast_logistic_regression_weighted_cpp(X = Xmm, y = y_prop, weights = n_k_v)
						# Compute vcov from Fisher information matrix XtWX
						vcov = solve(res$XtWX)
						list(b = res$b, ssq_b_1 = vcov[1, 1])
					}
				}, error = function(e) NULL)
				if (is.null(mod)){
					private$cached_values$beta_hat_T   = NA_real_
					if (!estimate_only) private$cached_values$s_beta_hat_T = NA_real_
					private$cached_values$is_z         = TRUE
					return(invisible(NULL))
				}
				private$cached_values$beta_hat_T   = as.numeric(mod$b[1L])
				if (!estimate_only) private$cached_values$s_beta_hat_T = sqrt(as.numeric(mod$ssq_b_1))
				private$cached_values$is_z         = TRUE

			# ---- Case 3: reservoir only — Poisson regression (C++) --------------
			} else {
				j_beta_T = 2L
				X_full   = if (p > 0L) cbind(1, w_r_v, X_r_v) else cbind(1, w_r_v)
				
				if (!estimate_only) {
					qr_X = qr(X_full)
					if (qr_X$rank < ncol(X_full)){
						keep     = qr_X$pivot[seq_len(qr_X$rank)]
						if (!(j_beta_T %in% keep)) keep[qr_X$rank] = j_beta_T
						keep     = sort(keep)
						X_full   = X_full[, keep, drop = FALSE]
						j_beta_T = which(keep == j_beta_T)
					}
				}
				
				mod = tryCatch({
					if (estimate_only) {
						res = fast_poisson_regression_cpp(X = X_full, y = y_r_v)
						list(b = res$b, ssq_b_j = NA_real_)
					} else {
						res = fast_poisson_regression_cpp(X = X_full, y = y_r_v)
						vcov = solve(res$XtWX)
						list(b = res$b, ssq_b_j = vcov[j_beta_T, j_beta_T])
					}
				}, error = function(e) NULL)
				if (is.null(mod)){
					private$cached_values$beta_hat_T   = NA_real_
					if (!estimate_only) private$cached_values$s_beta_hat_T = NA_real_
					private$cached_values$is_z         = TRUE
					return(invisible(NULL))
				}
				private$cached_values$beta_hat_T   = as.numeric(mod$b[j_beta_T])
				if (!estimate_only) private$cached_values$s_beta_hat_T = sqrt(as.numeric(mod$ssq_b_j))
				private$cached_values$is_z         = TRUE
			}

			invisible(NULL)
		}
	)
)
