# Abstract class for Robust-Regression IVWC Compound Inference for KK Designs
#
# @description
# Fits a variance-weighted compound estimator for KK matching-on-the-fly designs
# with continuous responses using robust linear regression (`MASS::rlm`) for the
# matched-pair and reservoir components separately.
#
# @keywords internal
SeqDesignInferenceAbstractKKRobustRegrIVWC = R6::R6Class("SeqDesignInferenceAbstractKKRobustRegrIVWC",
	inherit = SeqDesignInferenceKKPassThroughCompound,
	public = list(

		# @description
		# Initialize the inference object.
		# @param seq_des_obj		A SeqDesign object (must be a KK design).
		# @param method			Robust-regression fitting method for `MASS::rlm`; one of `"M"` or `"MM"`.
		# @param num_cores			Number of CPU cores for parallel processing.
		# @param verbose			Whether to print progress messages.
		initialize = function(seq_des_obj, method = "MM", num_cores = 1, verbose = FALSE){
			assertResponseType(seq_des_obj$get_response_type(), "continuous")
			assertChoice(method, c("M", "MM"))
			if (!is(seq_des_obj, "SeqDesignKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (SeqDesignKK14 or subclass).")
			}
			super$initialize(seq_des_obj, num_cores, verbose)
			assertNoCensoring(private$any_censoring)
			private$rlm_method = method
		},

		# @description
		# Returns the estimated treatment effect.
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		# @description
		# Computes the approximate confidence interval.
		# @param alpha The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		compute_mle_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		# @description
		# Computes the approximate p-value.
		# @param delta The null difference to test against. For any treatment effect at all this is set to zero (the default).
		compute_mle_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		rlm_method = NULL,

		include_covariates = function() stop(class(self)[1], " must implement include_covariates()"),

		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}

			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			if (m > 0){
				private$robust_for_matched_pairs()
			}
			beta_m = private$cached_values$beta_T_matched
			ssq_m  = private$cached_values$ssq_beta_T_matched
			m_ok   = !is.null(beta_m) && is.finite(beta_m) && !is.null(ssq_m) && is.finite(ssq_m) && ssq_m > 0

			if (nRT > 0 && nRC > 0){
				private$robust_for_reservoir()
			}
			beta_r = private$cached_values$beta_T_reservoir
			ssq_r  = private$cached_values$ssq_beta_T_reservoir
			r_ok   = !is.null(beta_r) && is.finite(beta_r) && !is.null(ssq_r) && is.finite(ssq_r) && ssq_r > 0

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
				stop("KK robust-regression IVWC estimator: could not compute a finite standard error.")
			}
		},

		fit_rlm_with_treatment = function(X, y, j_treat){
			if (nrow(X) <= ncol(X)) return(NULL)

			mod = tryCatch({
				if (identical(private$rlm_method, "M")) {
					MASS::rlm(x = X, y = y, method = "M", psi = MASS::psi.huber)
				} else {
					MASS::rlm(x = X, y = y, method = private$rlm_method)
				}
			}, error = function(e) NULL)
			if (is.null(mod)) return(NULL)

			coef_table = tryCatch(summary(mod)$coefficients, error = function(e) NULL)
			if (is.null(coef_table) || nrow(coef_table) < j_treat) return(NULL)

			beta = as.numeric(coef_table[j_treat, "Value"])
			se   = as.numeric(coef_table[j_treat, "Std. Error"])
			if (!is.finite(beta) || !is.finite(se) || se <= 0) return(NULL)

			list(beta = beta, ssq = se^2)
		},

		robust_for_matched_pairs = function(){
			yd = private$cached_values$KKstats$y_matched_diffs
			Xd = as.matrix(private$cached_values$KKstats$X_matched_diffs)
			m  = length(yd)

			if (private$include_covariates() && ncol(Xd) > 0L){
				qr_Xd = qr(Xd)
				r = qr_Xd$rank
				if (r < ncol(Xd)) {
					Xd = Xd[, qr_Xd$pivot[seq_len(r)], drop = FALSE]
				}
			}

			X = if (private$include_covariates() && ncol(Xd) > 0L) {
				cbind(1, Xd)
			} else {
				matrix(1, nrow = m, ncol = 1L)
			}

			fit = private$fit_rlm_with_treatment(X, yd, 1L)
			if (is.null(fit)) {
				private$cached_values$beta_T_matched     = if (m >= 1) mean(yd) else NA_real_
				private$cached_values$ssq_beta_T_matched = if (m >= 2) var(yd) / m else NA_real_
			} else {
				private$cached_values$beta_T_matched     = fit$beta
				private$cached_values$ssq_beta_T_matched = fit$ssq
			}
		},

		robust_for_reservoir = function(){
			y_r = private$cached_values$KKstats$y_reservoir
			w_r = private$cached_values$KKstats$w_reservoir
			X_r = as.matrix(private$cached_values$KKstats$X_reservoir)
			j_treat = 2L

			if (private$include_covariates()){
				X_full = cbind(1, w_r, X_r)
				qr_full = qr(X_full)
				r_full = qr_full$rank
				if (r_full < ncol(X_full)){
					keep = qr_full$pivot[seq_len(r_full)]
					if (!(2L %in% keep)) keep[r_full] = 2L
					keep = sort(keep)
					X_full = X_full[, keep, drop = FALSE]
					j_treat = which(keep == 2L)
				}
			} else {
				X_full = cbind(1, w_r)
			}

			fit = private$fit_rlm_with_treatment(X_full, y_r, j_treat)
			if (is.null(fit)) {
				private$cached_values$beta_T_reservoir     = NA_real_
				private$cached_values$ssq_beta_T_reservoir = NA_real_
			} else {
				private$cached_values$beta_T_reservoir     = fit$beta
				private$cached_values$ssq_beta_T_reservoir = fit$ssq
			}
		}
	)
)
