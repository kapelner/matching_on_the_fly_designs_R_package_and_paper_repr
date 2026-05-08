#' Abstract class for Stratified Cox / Standard Cox Compound Inference
#'
#' This class implements a compound estimator for KK matching-on-the-fly designs with
#' survival responses. For matched pairs, it uses stratified Cox proportional hazards
#' regression (each pair is a stratum). For reservoir subjects, it uses standard Cox
#' regression. The two estimates (both log-hazard ratios) are combined via a
#' variance-weighted linear combination.
#'
#' Under \code{harden = TRUE}, multivariate fits preserve the treatment column and
#' progressively retry reduced covariate sets after QR-based rank reduction and
#' correlation-based pruning. Extreme finite coefficients / standard errors are
#' rejected and treated as non-estimable.
#'
#' @keywords internal
InferenceAbstractKKStratCoxIVWC = R6::R6Class("InferenceAbstractKKStratCoxIVWC",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj  	A DesignSeqOneByOne object (must be a KK design).
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose  		Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL,  verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
			}
			if (should_run_asserts()) {
				if (!inherits(des_obj, "DesignSeqOneByOneKK14") && !inherits(des_obj, "FixedDesignBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
				}
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
		},

		#' @description
		#' Returns the estimated treatment effect (log-hazard ratio).
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
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
		compute_asymp_two_sided_pval = function(delta = 0){
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
		max_abs_reasonable_coef = 1e4,

		# Abstract: subclasses return TRUE (multivariate) or FALSE (univariate).
		cox_design_candidates = function(w, X){
			X_full = matrix(w, ncol = 1)
			colnames(X_full) = "w"
			X_covs = as.matrix(X)
			if (ncol(X_covs) > 0L){
				X_full = cbind(X_full, X_covs)
			}
			if (!private$harden || ncol(X_full) <= 1L){
				return(list(X_full))
			}

			attempt = private$fit_with_hardened_qr_column_dropping(
				X_full = X_full,
				required_cols = 1L,
				fit_fun = function(X_fit) X_fit,
				fit_ok = function(mod, X_fit, keep) TRUE
			)
			candidates = list(attempt$X)
			keys = paste(colnames(candidates[[1L]]), collapse = "|")

			thresholds = c(0.99, 0.95, 0.90, 0.85, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20, 0.10)
			X_cov_orig = X_full[, -1, drop = FALSE]
			for (thresh in thresholds){
				X_cov = drop_highly_correlated_cols(X_cov_orig, threshold = thresh)$M
				X_try = matrix(w, ncol = 1)
				colnames(X_try) = "w"
				if (ncol(X_cov) > 0){
					X_try = cbind(X_try, X_cov)
				}
				attempt_try = private$fit_with_hardened_qr_column_dropping(
					X_full = X_try,
					required_cols = 1L,
					fit_fun = function(X_fit) X_fit,
					fit_ok = function(mod, X_fit, keep) TRUE
				)
				key = paste(colnames(attempt_try$X), collapse = "|")
				if (!(key %in% keys)){
					candidates[[length(candidates) + 1L]] = attempt_try$X
					keys = c(keys, key)
				}
			}
			candidates
		},

		rcpp_cox_fit_is_usable = function(res, estimate_only = FALSE){
			if (is.null(res) || !isTRUE(res$converged)) return(FALSE)
			beta = res$coefficients[1L]
			if (!is.finite(beta) || abs(beta) > private$max_abs_reasonable_coef) return(FALSE)
			if (estimate_only) return(TRUE)
			se = tryCatch(sqrt(res$vcov[1L, 1L]), error = function(e) NA_real_)
			is.finite(se) && se > 0 && se <= private$max_abs_reasonable_coef
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

			# --- Matched pairs: Stratified Cox ---
			if (m > 0){
				private$strat_cox_for_matched_pairs(estimate_only = estimate_only)
			}
			beta_m   = private$cached_values$beta_T_matched
			ssq_m    = private$cached_values$ssq_beta_T_matched
			m_ok     = !is.null(beta_m) && is.finite(beta_m) &&
			           (!estimate_only && !is.null(ssq_m) && is.finite(ssq_m) && ssq_m > 0 || estimate_only)

			# --- Reservoir: Standard Cox ---
			if (nRT > 0 && nRC > 0){
				private$cox_for_reservoir(estimate_only = estimate_only)
			}
			beta_r   = private$cached_values$beta_T_reservoir
			ssq_r    = private$cached_values$ssq_beta_T_reservoir
			r_ok     = !is.null(beta_r) && is.finite(beta_r) &&
			           (!estimate_only && !is.null(ssq_r) && is.finite(ssq_r) && ssq_r > 0 || estimate_only)

			# --- Variance-weighted combination ---
			if (m_ok && r_ok){
				w_star = ssq_r / (ssq_r + ssq_m)
				private$cached_values$beta_hat_T   = w_star * beta_m + (1 - w_star) * beta_r
				if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = sqrt(ssq_m * ssq_r / (ssq_m + ssq_r))
			} else if (m_ok){
				private$cached_values$beta_hat_T   = beta_m
				private$cached_values$s_beta_hat_T = if (estimate_only) NA_real_ else sqrt(ssq_m)
			} else if (r_ok){
				private$cached_values$beta_hat_T   = beta_r
				private$cached_values$s_beta_hat_T = if (estimate_only) NA_real_ else sqrt(ssq_r)
			} else {
				private$cache_nonestimable_estimate("kk_strat_cox_ivwc_no_usable_component")
				return(invisible(NULL))
			}
			if (is.finite(private$cached_values$beta_hat_T) &&
			    abs(private$cached_values$beta_hat_T) > private$max_abs_reasonable_coef){
				private$cache_nonestimable_estimate("kk_strat_cox_ivwc_extreme_estimate")
				return(invisible(NULL))
			}
			if (!estimate_only &&
			    (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0 ||
			     private$cached_values$s_beta_hat_T > private$max_abs_reasonable_coef)){
				private$cache_nonestimable_se("kk_strat_cox_ivwc_standard_error_unavailable")
				return(invisible(NULL))
			}
			private$clear_nonestimable_state()
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},

		strat_cox_for_matched_pairs = function(estimate_only = FALSE){
			private$cached_values$beta_T_matched = NULL
			private$cached_values$ssq_beta_T_matched = NULL

			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			i_matched = which(m_vec > 0)
			y_m       = private$y[i_matched]
			dead_m    = private$dead[i_matched]
			w_m       = private$w[i_matched]
			strata_m  = m_vec[i_matched]

			# Filter strata that have no events (provide no information for Cox partial likelihood)
			strata_with_events = unique(strata_m[dead_m == 1])
			if (length(strata_with_events) == 0) return(invisible(NULL))

			i_valid = which(strata_m %in% strata_with_events)
			y_v     = y_m[i_valid]
			dead_v  = dead_m[i_valid]
			w_v     = w_m[i_valid]
			strata_v = as.integer(strata_m[i_valid])

			res = NULL
			if (ncol(as.matrix(private$X)) > 0){
				X_m = as.matrix(private$get_X()[i_matched[i_valid], drop = FALSE])
				for (X_candidate in private$cox_design_candidates(w_v, X_m)){
					res_try = tryCatch(
						fast_stratified_coxph_regression_cpp(X_candidate, y_v, dead_v, strata_v, estimate_only = estimate_only),
						error = function(e) NULL
					)
					if (private$rcpp_cox_fit_is_usable(res_try, estimate_only = estimate_only)){
						res = res_try
						break
					}
				}
			} else {
				X_mat = matrix(w_v, ncol = 1L)
				colnames(X_mat) = "w"
				res = tryCatch(
					fast_stratified_coxph_regression_cpp(X_mat, y_v, dead_v, strata_v, estimate_only = estimate_only),
					error = function(e) NULL
				)
			}
			if (is.null(res)) return(invisible(NULL))

			beta = res$coefficients[1L]
			private$cached_values$beta_T_matched = beta
			
			if (!estimate_only) {
				se   = sqrt(res$vcov[1L, 1L])
				if (!is.finite(beta) || !is.finite(se) || se <= 0 || se > private$max_abs_reasonable_coef || abs(beta) > private$max_abs_reasonable_coef) return(invisible(NULL))
				private$cached_values$ssq_beta_T_matched = se^2
			}
		},

		cox_for_reservoir = function(estimate_only = FALSE){
			private$cached_values$beta_T_reservoir = NULL
			private$cached_values$ssq_beta_T_reservoir = NULL

			y_r    = private$cached_values$KKstats$y_reservoir
			w_r    = private$cached_values$KKstats$w_reservoir
			m_vec_safe = private$m
			if (is.null(m_vec_safe)) m_vec_safe = rep(0L, private$n)
			m_vec_safe[is.na(m_vec_safe)] = 0L
			dead_r = private$dead[m_vec_safe == 0]
			X_r    = as.matrix(private$cached_values$KKstats$X_reservoir)

			res = NULL
			if (ncol(as.matrix(private$X)) > 0){
				for (X_candidate in private$cox_design_candidates(w_r, X_r)){
					res_try = tryCatch(
						fast_coxph_regression_cpp(X_candidate, y_r, dead_r, estimate_only = estimate_only),
						error = function(e) NULL
					)
					if (private$rcpp_cox_fit_is_usable(res_try, estimate_only = estimate_only)){
						res = res_try
						break
					}
				}
			}
			if (is.null(res)){
				X_mat = matrix(w_r, ncol = 1L)
				colnames(X_mat) = "w"
				res = tryCatch(
					fast_coxph_regression_cpp(X_mat, y_r, dead_r, estimate_only = estimate_only),
					error = function(e) NULL
				)
			}

			if (is.null(res)) return(invisible(NULL))

			beta = res$coefficients[1L]
			private$cached_values$beta_T_reservoir = beta
			
			if (!estimate_only) {
				se   = sqrt(res$vcov[1L, 1L])
				if (!is.finite(beta) || !is.finite(se) || se <= 0 || se > private$max_abs_reasonable_coef || abs(beta) > private$max_abs_reasonable_coef) return(invisible(NULL))
				private$cached_values$ssq_beta_T_reservoir = se^2
			}
		}
	)
)
