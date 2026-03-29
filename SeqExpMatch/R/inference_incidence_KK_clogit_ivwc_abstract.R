#' Abstract class for Conditional Logistic Compound Inference
#'
#' @keywords internal
InferenceAbstractKKClogitIVWC = R6::R6Class("InferenceAbstractKKClogitIVWC",
	inherit = InferenceKKPassThrough,
	public = list(

		# @description
		# Initialize the inference object.
		# @param des_obj		A DesignSeqOneByOne object (must be a KK design).
		# @param num_cores			Number of CPU cores for parallel processing.
		# @param verbose			Whether to print progress messages.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE, make_fork_cluster = NULL){
			assertResponseType(des_obj$get_response_type(), "incidence")
			if (!is(des_obj, "DesignSeqOneByOneKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
			}
			super$initialize(des_obj, num_cores, verbose, make_fork_cluster = make_fork_cluster)
			assertNoCensoring(private$any_censoring)
		},

		# @description
		# Returns the estimated treatment effect.
		compute_treatment_estimate = function(){
			private$shared()
			private$cached_values$beta_hat_T
		},

		# @description
		# Computes the asymptotic confidence interval.
		# @param alpha					The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		# @description
		# Computes the asymptotic p-value.
		# @param delta					The null difference to test against. For any treatment effect at all this is set to zero (the default).
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$assert_finite_se()
			if (delta == 0){
				private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			} else {
				stop("TO-DO")
			}
		}
	),

	private = list(

		# Abstract: subclasses return TRUE (multivariate) or FALSE (univariate).
		include_covariates = function() stop(class(self)[1], " must implement include_covariates()"),

		shared = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			# Recompute KKstats if cache was cleared (e.g., after y transformation for rand CI)
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}

			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			# --- Matched pairs: clogit ---
			if (m > 0){
				private$clogit_for_matched_pairs()
			}
			beta_m   = private$cached_values$beta_T_matched
			ssq_m    = private$cached_values$ssq_beta_T_matched
			m_ok     = !is.null(beta_m) && is.finite(beta_m) &&
			           !is.null(ssq_m)  && is.finite(ssq_m) && ssq_m > 0

			# --- Reservoir: logistic regression ---
			if (nRT > 0 && nRC > 0){
				private$logistic_for_reservoir()
			}
			beta_r   = private$cached_values$beta_T_reservoir
			ssq_r    = private$cached_values$ssq_beta_T_reservoir
			r_ok     = !is.null(beta_r) && is.finite(beta_r) &&
			           !is.null(ssq_r)  && is.finite(ssq_r) && ssq_r > 0

			# --- Variance-weighted combination (mirrors InferenceContinMultOLSKK) ---
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
				stop("Clogit/logistic compound estimator: could not compute a finite standard error (possible perfect separation or insufficient data).")
			}
		},

		clogit_for_matched_pairs = function(){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(0L, private$n)
			m_vec[is.na(m_vec)] = 0L

			i_matched = which(m_vec > 0)
			y_m       = private$y[i_matched]
			w_m       = private$w[i_matched]
			strata_m  = m_vec[i_matched]
			X_m       = if (private$include_covariates()) as.data.frame(private$X[i_matched, , drop = FALSE]) else data.frame()

			mod = clogit_helper(y_m, X_m, w_m, strata_m)
			if (is.null(mod)) return(invisible(NULL))

			beta = as.numeric(mod$b[1])
			ssq  = as.numeric(mod$ssq_b_j)
			private$cached_values$beta_T_matched     = if (is.finite(beta)) beta else NA_real_
			private$cached_values$ssq_beta_T_matched = if (is.finite(ssq) && ssq > 0) ssq else NA_real_
		},

		logistic_for_reservoir = function(){
			y_r    = private$cached_values$KKstats$y_reservoir
			w_r    = private$cached_values$KKstats$w_reservoir
			X_r    = as.matrix(private$cached_values$KKstats$X_reservoir)
			j_treat = 2L

			if (private$include_covariates()){
				X_full = cbind(1, w_r, X_r)
				# QR-reduce to full rank while always preserving the treatment column,
				# mirroring the approach used in ols_for_reservoir().
				qr_full = qr(X_full)
				r_full  = qr_full$rank
				if (r_full < ncol(X_full)){
					keep = qr_full$pivot[seq_len(r_full)]
					if (!(2L %in% keep)) keep[r_full] = 2L
					keep    = sort(keep)
					X_full  = X_full[, keep, drop = FALSE]
					j_treat = which(keep == 2L)
				}
			} else {
				X_full = cbind(1, w_r)
			}

			mod = tryCatch(
				fast_logistic_regression_with_var(X_full, y_r, j = j_treat),
				error = function(e) NULL
			)
			# Fallback: if the covariates model failed, retry with just intercept + treatment
			if (is.null(mod) && private$include_covariates()){
				mod = tryCatch(
					fast_logistic_regression_with_var(cbind(1, w_r), y_r, j = 2L),
					error = function(e) NULL
				)
				j_treat = 2L
			}
			if (is.null(mod)) return(invisible(NULL))

			beta = as.numeric(mod$b[j_treat])
			ssq  = as.numeric(mod$ssq_b_j)
			private$cached_values$beta_T_reservoir     = if (is.finite(beta)) beta else NA_real_
			private$cached_values$ssq_beta_T_reservoir = if (is.finite(ssq) && ssq > 0) ssq else NA_real_
		}
	)
)
