#' Abstract class for Conditional Logistic Compound Inference
#'
#' Under \code{harden = TRUE}, the matched-pair clogit and reservoir logistic
#' components preserve the treatment column and retry reduced covariate sets after
#' QR-based rank reduction. Extreme finite coefficients / standard errors are
#' rejected and treated as non-estimable.
#'
#' @keywords internal
InferenceAbstractKKClogitIVWC = R6::R6Class("InferenceAbstractKKClogitIVWC",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj		A DesignSeqOneByOne object (must be a KK design).
		#' @param verbose			Whether to print progress messages.
		initialize = function(des_obj,  verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "incidence")
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
					stop("TO-DO")
				}
				NA_real_
			}
		}
	),

	private = list(
		max_abs_reasonable_coef = 1e4,
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
				i_matched = which(private$m > 0)
				y_m       = private$y[i_matched]
				w_m       = private$w[i_matched]
				strata_m  = private$m[i_matched]
				
				if (private$include_covariates() && !is.null(private$best_Xmm_colnames_matched)){
					X_m = private$get_X()[i_matched, intersect(private$best_Xmm_colnames_matched, colnames(private$get_X())), drop = FALSE]
					mod_m = clogit_helper(y_m, as.data.frame(X_m), w_m, strata_m)
				} else {
					mod_m = clogit_helper(y_m, data.frame(), w_m, strata_m)
				}
				if (!is.null(mod_m) && is.finite(mod_m$b[1])) beta_m = as.numeric(mod_m$b[1])
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
				
				mod_r = tryCatch(fast_logistic_regression_cpp(X = Xmm, y = as.numeric(y_r)), error = function(e) NULL)
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

		# Abstract: subclasses return TRUE (multivariate) or FALSE (univariate).
		include_covariates = function() stop(class(self)[1], " must implement include_covariates()"),

		component_is_usable = function(beta, ssq){
			!is.null(beta) && is.finite(beta) &&
				abs(beta) <= private$max_abs_reasonable_coef &&
				!is.null(ssq) && is.finite(ssq) && ssq > 0 &&
				sqrt(ssq) <= private$max_abs_reasonable_coef
		},

		cache_failed_component = function(which_component){
			private$cached_values[[paste0("beta_T_", which_component)]] = NA_real_
			private$cached_values[[paste0("ssq_beta_T_", which_component)]] = NA_real_
			invisible(NULL)
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			private$clear_nonestimable_state()

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
			m_ok     = private$component_is_usable(beta_m, ssq_m)

			# --- Reservoir: logistic regression ---
			if (nRT > 0 && nRC > 0){
				private$logistic_for_reservoir()
			}
			beta_r   = private$cached_values$beta_T_reservoir
			ssq_r    = private$cached_values$ssq_beta_T_reservoir
			r_ok     = private$component_is_usable(beta_r, ssq_r)

			# --- Variance-weighted combination (mirrors InferenceContinMultOLSKK) ---
			if (m_ok && r_ok){
				if (estimate_only) {
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
				private$cache_nonestimable_estimate("kk_clogit_ivwc_no_usable_component")
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}
			if (is.finite(private$cached_values$beta_hat_T) &&
			    abs(private$cached_values$beta_hat_T) > private$max_abs_reasonable_coef){
				private$cache_nonestimable_estimate("kk_clogit_ivwc_extreme_estimate")
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}
			if (!estimate_only &&
			    is.finite(private$cached_values$s_beta_hat_T) &&
			    private$cached_values$s_beta_hat_T > private$max_abs_reasonable_coef){
				private$cache_nonestimable_se("kk_clogit_ivwc_extreme_standard_error")
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}
			private$clear_nonestimable_state()
			private$cached_values$is_z = TRUE
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},

		clogit_for_matched_pairs = function(){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			i_matched = which(m_vec > 0)
			y_m       = private$y[i_matched]
			w_m       = private$w[i_matched]
			strata_m  = m_vec[i_matched]
			if (!private$include_covariates()){
				mod = clogit_helper(y_m, data.frame(), w_m, strata_m)
				if (is.null(mod)) {
					private$cache_failed_component("matched")
					return(invisible(NULL))
				}
				beta = as.numeric(mod$b[1])
				ssq  = as.numeric(mod$ssq_b_j)
				private$cached_values$beta_T_matched     = if (is.finite(beta) && abs(beta) <= private$max_abs_reasonable_coef) beta else NA_real_
				private$cached_values$ssq_beta_T_matched = if (is.finite(ssq) && ssq > 0 && sqrt(ssq) <= private$max_abs_reasonable_coef) ssq else NA_real_
				return(invisible(NULL))
			}

			X_m = as.matrix(private$get_X()[i_matched, , drop = FALSE])
			X_full = cbind(w = w_m, X_m)
			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 1L,
				fit_fun = function(X_fit){
					clogit_helper(
						y_m,
						as.data.frame(X_fit[, -1, drop = FALSE]),
						X_fit[, 1],
						strata_m
					)
				},
				fit_ok = function(mod, X_fit, keep){
					!is.null(mod) &&
						is.finite(mod$b[1]) &&
						abs(as.numeric(mod$b[1])) <= private$max_abs_reasonable_coef &&
						is.finite(mod$ssq_b_j) &&
						mod$ssq_b_j > 0 &&
						sqrt(as.numeric(mod$ssq_b_j)) <= private$max_abs_reasonable_coef
				}
			)
			mod = attempt$fit
			if (is.null(mod)) {
				private$cache_failed_component("matched")
				return(invisible(NULL))
			}

			if (private$include_covariates()){
				private$best_Xmm_colnames_matched = setdiff(colnames(attempt$X_fit), "w")
			}

			beta = as.numeric(mod$b[1])
			ssq  = as.numeric(mod$ssq_b_j)
			private$cached_values$beta_T_matched     = if (is.finite(beta) && abs(beta) <= private$max_abs_reasonable_coef) beta else NA_real_
			private$cached_values$ssq_beta_T_matched = if (is.finite(ssq) && ssq > 0 && sqrt(ssq) <= private$max_abs_reasonable_coef) ssq else NA_real_
		},

		logistic_for_reservoir = function(){
			y_r    = private$cached_values$KKstats$y_reservoir
			w_r    = private$cached_values$KKstats$w_reservoir
			X_r    = as.matrix(private$cached_values$KKstats$X_reservoir)
			j_treat = 2L

			if (private$include_covariates()){
				X_full = cbind(Intercept = 1, w = w_r, X_r)
				attempt = private$fit_with_hardened_qr_column_dropping(
					X_full = X_full,
					required_cols = c(1L, 2L),
					fit_fun = function(X_fit){
						j_treat_fit = match("w", colnames(X_fit))
						fast_logistic_regression_with_var(X_fit, y_r, j = j_treat_fit)
					},
					fit_ok = function(mod, X_fit, keep){
						if (is.null(mod)) return(FALSE)
						j_treat_fit = match("w", colnames(X_fit))
						if (!is.finite(j_treat_fit) || is.na(j_treat_fit)) return(FALSE)
						beta = suppressWarnings(as.numeric(mod$b[j_treat_fit]))
						ssq = suppressWarnings(as.numeric(mod$ssq_b_j))
						is.finite(beta) &&
							abs(beta) <= private$max_abs_reasonable_coef &&
							is.finite(ssq) &&
							ssq > 0 &&
							sqrt(ssq) <= private$max_abs_reasonable_coef
					}
				)
				mod = attempt$fit
				j_treat = if (!is.null(attempt$X_fit)) match("w", colnames(attempt$X_fit)) else NA_integer_
				if (!is.null(mod) && private$include_covariates()){
					private$best_Xmm_colnames_reservoir = setdiff(colnames(attempt$X_fit), c("Intercept", "w"))
					private$best_Xmm_j_treat_reservoir = j_treat
				}
			} else {
				X_full = cbind(Intercept = 1, w = w_r)
				mod = tryCatch(
					fast_logistic_regression_with_var(X_full, y_r, j = 2L),
					error = function(e) NULL
				)
				j_treat = 2L
			}
			if (is.null(mod) || !is.finite(j_treat) || is.na(j_treat)) {
				private$cache_failed_component("reservoir")
				return(invisible(NULL))
			}

			beta = as.numeric(mod$b[j_treat])
			ssq  = as.numeric(mod$ssq_b_j)
			private$cached_values$beta_T_reservoir     = if (is.finite(beta) && abs(beta) <= private$max_abs_reasonable_coef) beta else NA_real_
			private$cached_values$ssq_beta_T_reservoir = if (is.finite(ssq) && ssq > 0 && sqrt(ssq) <= private$max_abs_reasonable_coef) ssq else NA_real_
		}
	)
)
