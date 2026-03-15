# Abstract Quantile Regression Compound Estimator for KK Matching-on-the-Fly Designs
#
# @description
# An abstract base class providing shared quantile regression logic for KK matching-on-the-fly
# designs. Subclasses override the \code{transform_y_fn} private field to apply a response
# transformation before quantile regression (e.g., \code{identity} for continuous,
# \code{qlogis} for proportion outcomes).
#
# @keywords internal
SeqDesignInferenceAbstractKKQuantileRegrIVWC = R6::R6Class("SeqDesignInferenceAbstractKKQuantileRegrIVWC",
	inherit = SeqDesignInferenceAbstractQuantileRandCI,
	public = list(

		# @param seq_des_obj		A SeqDesign object whose entire n subjects are assigned and response y is recorded within.
		# @param tau				The quantile level for regression, strictly between 0 and 1. The default \code{tau = 0.5}
		# 							estimates the median treatment effect.
		# @param transform_y_fn	A function applied to y values before quantile regression. Subclasses pass
		# 							\code{identity} (continuous) or \code{qlogis} (proportion). Not exposed publicly.
		# @param num_cores			The number of CPU cores to use to parallelize sampling.
		# @param verbose			A flag indicating whether messages should be displayed. Default is \code{FALSE}.
		initialize = function(seq_des_obj, tau = 0.5, transform_y_fn = identity, num_cores = 1, verbose = FALSE){
			assertNumeric(tau, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps)
			if (!requireNamespace("quantreg", quietly = TRUE)) {
				stop("Package 'quantreg' is required. Please install it with install.packages(\"quantreg\").")
			}
			private$tau = tau
			private$transform_y_fn_list = list(fn = transform_y_fn)
			super$initialize(seq_des_obj, num_cores, verbose)
		},

		# @description
		# Computes the appropriate quantile regression compound estimate
		#
		# @return 	The setting-appropriate numeric estimate of the treatment effect
		compute_treatment_estimate = function(){
			if (is.null(private$cached_values$beta_hat_T)){
				private$shared()
			}
			private$cached_values$beta_hat_T
		},

		# @description
		# Computes a 1-alpha level frequentist confidence interval based on asymptotic normality
		# of the quantile regression estimator (z-based).
		#
		# @param alpha					The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		#
		# @return 	A (1 - alpha)-sized frequentist confidence interval for the treatment effect
		compute_mle_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			if (is.null(private$cached_values$s_beta_hat_T)){
				private$shared()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		# @description
		# Computes a 2-sided p-value based on asymptotic normality of the quantile regression estimator.
		#
		# @param delta					The null difference to test against. Default is zero.
		#
		# @return 	The approximate frequentist p-value
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			if (is.null(private$cached_values$s_beta_hat_T)){
				private$shared()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		tau = NULL,
		transform_y_fn_list = NULL,  # list(fn = ...) wrapping avoids R6 treating function as a locked method

		shared = function(){
			KKstats = private$cached_values$KKstats
			if (is.null(KKstats)){
				private$compute_basic_match_data()
				KKstats = private$cached_values$KKstats
			}

			m    = KKstats$m
			nRT  = KKstats$nRT
			nRC  = KKstats$nRC

			beta_m  = NA_real_
			ssq_m   = NA_real_
			beta_r  = NA_real_
			ssq_r   = NA_real_

			if (m > 0){
				res_m = private$quantile_for_matched_pairs()
				beta_m = res_m$beta
				ssq_m  = res_m$ssq
			}
			if (nRT > 0 && nRC > 0){
				res_r = private$quantile_for_reservoir()
				beta_r = res_r$beta
				ssq_r  = res_r$ssq
			}

			m_ok = is.finite(beta_m) && is.finite(ssq_m) && ssq_m > 0
			r_ok = is.finite(beta_r) && is.finite(ssq_r) && ssq_r > 0

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
			invisible(NULL)
		},

		quantile_for_matched_pairs = function(){
			yd = private$transform_y_fn_list$fn(private$cached_values$KKstats$yTs_matched) -
			     private$transform_y_fn_list$fn(private$cached_values$KKstats$yCs_matched)
			Xd = private$cached_values$KKstats$X_matched_diffs
			m  = length(yd)
			tau = private$tau

			# QR-reduce Xd to full rank (same logic as OLS compound)
			if (ncol(Xd) > 0){
				Xd = qr_reduce_full_rank_cpp(Xd)$X_reduced
			}
			p_kept = ncol(Xd)

			# Underdetermined: fall back to intercept-only quantile regression
			if (m <= p_kept + 1){
				fit0 = tryCatch(
					suppressWarnings(quantreg::rq(yd ~ 1, tau = tau)),
					error = function(e) NULL
				)
				if (is.null(fit0) || m < 2){
					return(list(beta = NA_real_, ssq = NA_real_))
				}
				beta = unname(coef(fit0)[1])
				se = private$iqr_se(yd, m)
				if (!is.finite(se) || se <= 0){
					return(list(beta = beta, ssq = NA_real_))
				}
				return(list(beta = beta, ssq = se^2))
			}

			# Main path: rq on differences with reduced covariates
			p = ncol(Xd)
			colnames(Xd) = paste0("xd", seq_len(p))
			dat = as.data.frame(Xd)
			dat$yd__ = yd
			fit = tryCatch(
				suppressWarnings(quantreg::rq(yd__ ~ ., tau = tau, data = dat)),
				error = function(e) NULL
			)
			if (is.null(fit)){
				return(list(beta = NA_real_, ssq = NA_real_))
			}

			se = private$extract_se_from_rq(fit, "(Intercept)")
			if (!is.finite(se) || se <= 0){
				se = private$iqr_se(yd, m)
			}
			beta = tryCatch(coef(fit)[["(Intercept)"]], error = function(e) NA_real_)
			list(beta = beta, ssq = if (is.finite(se) && se > 0) se^2 else NA_real_)
		},

		quantile_for_reservoir = function(){
			y_r = private$transform_y_fn_list$fn(private$cached_values$KKstats$y_reservoir)
			w_r = private$cached_values$KKstats$w_reservoir
			X_r = as.matrix(private$cached_values$KKstats$X_reservoir)
			tau = private$tau

			# Build full design matrix and QR-reduce, always preserving treatment column
			X_full  = cbind(1, w_r, X_r)
			j_treat = 2L
			reduced = qr_reduce_preserve_cols_cpp(X_full, j_treat)
			X_full  = reduced$X_reduced
			j_treat = match(2L, reduced$keep)

			p = ncol(X_full)
			cn = paste0("xr", seq_len(p))
			cn[j_treat] = "trt__"
			colnames(X_full) = cn

			dat = as.data.frame(X_full)
			dat$yr__ = y_r

			fit = tryCatch(
				suppressWarnings(quantreg::rq(yr__ ~ . - 1, tau = tau, data = dat)),
				error = function(e) NULL
			)
			if (is.null(fit)){
				return(list(beta = NA_real_, ssq = NA_real_))
			}

			se   = private$extract_se_from_rq(fit, "trt__")
			if (!is.finite(se) || se <= 0){
				nRT = private$cached_values$KKstats$nRT
				nRC = private$cached_values$KKstats$nRC
				se_T = private$iqr_se(y_r[w_r == 1], nRT)
				se_C = private$iqr_se(y_r[w_r == 0], nRC)
				if (is.finite(se_T) && is.finite(se_C) && se_T > 0 && se_C > 0){
					se = sqrt(se_T^2 + se_C^2)
				}
			}
			beta = tryCatch(coef(fit)[["trt__"]], error = function(e) NA_real_)
			list(beta = beta, ssq = if (is.finite(se) && se > 0) se^2 else NA_real_)
		},

		# IQR-based SE for a sample quantile: CLT approximation SE = IQR/(2*qnorm(0.75)) / sqrt(n)
		iqr_se = function(x, n){
			if (n < 2) return(NA_real_)
			iqr = IQR(x)
			if (!is.finite(iqr) || iqr <= 0) return(NA_real_)
			iqr / (2 * qnorm(0.75)) / sqrt(n)
		},

		# Helper: extract SE from rq fit by coefficient name, trying "nid" then "iid".
		# SEs above 1e6 are treated as invalid (the "nid" sparsity estimator can return
		# astronomically large but finite values when the density at the quantile is near
		# zero, which bypasses the usual !is.finite() || <= 0 guard in callers).
		extract_se_from_rq = function(fit, coef_name){
			.extract_se_from_rq_fit(fit, coef_name)
		},

		# Two-sided sign-flip randomisation test for matched pairs at H0: Q_tau(yd) = delta_0.
		# Within each pair the sign of (transform(yT) - transform(yC)) is equally likely to be
		# + or - under randomisation.  Flipping the assignment also flips X_d, so both are negated.
		compute_rand_pval_matched_pairs = function(delta_0){
			KKstats = private$cached_values$KKstats
			fn  = private$transform_y_fn_list$fn
			tau = private$tau

			yd_adj = fn(KKstats$yTs_matched) - fn(KKstats$yCs_matched) - delta_0
			Xd     = KKstats$X_matched_diffs
			m      = length(yd_adj)

			if (ncol(Xd) > 0){
				Xd = qr_reduce_full_rank_cpp(Xd)$X_reduced
			}

			T_obs = private$qr_intercept_pairs(yd_adj, Xd, tau, m)
			if (!is.finite(T_obs)) return(NA_real_)

			T_rand = replicate(private$nsim_rand, {
				signs = sample(c(-1L, 1L), m, replace = TRUE)
				private$qr_intercept_pairs(
					signs * yd_adj,
					sweep(Xd, 1L, signs, `*`),
					tau, m
				)
			})
			T_rand = T_rand[is.finite(T_rand)]
			if (length(T_rand) == 0L) return(NA_real_)
			mean(abs(T_rand) >= abs(T_obs))
		},

		# Two-sided permutation test for reservoir at H0: Q_tau(transform(y)|w=1,X) - Q_tau(...|w=0,X) = delta_0.
		compute_rand_pval_reservoir = function(delta_0){
			KKstats = private$cached_values$KKstats
			fn  = private$transform_y_fn_list$fn
			tau = private$tau

			w_r = KKstats$w_reservoir
			X_r = as.matrix(KKstats$X_reservoir)
			ty  = fn(KKstats$y_reservoir)

			X_full  = cbind(1, w_r, X_r)
			j_treat = 2L
			reduced = qr_reduce_preserve_cols_cpp(X_full, j_treat)
			X_full  = reduced$X_reduced
			j_treat = match(2L, reduced$keep)
			p  = ncol(X_full)
			cn = paste0("xr", seq_len(p))
			cn[j_treat] = "trt__"
			colnames(X_full) = cn

			y_adj_obs = ty - delta_0 * w_r
			T_obs = private$qr_trt_coef_reservoir(y_adj_obs, X_full, tau)
			if (!is.finite(T_obs)) return(NA_real_)

			T_rand = replicate(private$nsim_rand, {
				w_perm = sample(w_r)
				X_perm = X_full
				X_perm[, j_treat] = w_perm
				private$qr_trt_coef_reservoir(ty - delta_0 * w_perm, X_perm, tau)
			})
			T_rand = T_rand[is.finite(T_rand)]
			if (length(T_rand) == 0L) return(NA_real_)
			mean(abs(T_rand) >= abs(T_obs))
		},

		# Helper: QR intercept from rq(yd ~ Xd) for matched-pair randomisation iterations
		qr_intercept_pairs = function(yd, Xd, tau, m){
			p = ncol(Xd)
			if (p == 0L || m <= p + 1L){
				fit = tryCatch(
					suppressWarnings(quantreg::rq(yd ~ 1, tau = tau)),
					error = function(e) NULL
				)
				if (is.null(fit)) return(NA_real_)
				return(unname(coef(fit)[1L]))
			}
			colnames(Xd) = paste0("xd", seq_len(p))
			dat = as.data.frame(Xd)
			dat$yd__ = yd
			fit = tryCatch(
				suppressWarnings(quantreg::rq(yd__ ~ ., tau = tau, data = dat)),
				error = function(e) NULL
			)
			if (is.null(fit)) return(NA_real_)
			tryCatch(coef(fit)[["(Intercept)"]], error = function(e) NA_real_)
		},

		# Helper: QR treatment coefficient from rq(y_adj ~ X_full - 1) for reservoir randomisation iterations
		qr_trt_coef_reservoir = function(y_adj, X_full, tau){
			dat = as.data.frame(X_full)
			dat$yr__ = y_adj
			fit = tryCatch(
				suppressWarnings(quantreg::rq(yr__ ~ . - 1, tau = tau, data = dat)),
				error = function(e) NULL
			)
			if (is.null(fit)) return(NA_real_)
			tryCatch(coef(fit)[["trt__"]], error = function(e) NA_real_)
		}
	)
)
