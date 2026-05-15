#' Abstract class for LWA-style Marginal Cox / Standard Cox Compound Inference
#'
#' This class implements a compound estimator for KK matching-on-the-fly designs with
#' survival responses. For matched pairs, it uses a marginal Cox proportional hazards
#' model with Lee-Wei-Amato style cluster-robust variance, treating each pair as a
#' cluster of size two. For reservoir subjects, it uses standard Cox regression. The
#' two estimates (both log-hazard ratios) are combined via a variance-weighted linear
#' combination.
#'
#' Under \code{harden = TRUE}, multivariate component fits preserve the treatment
#' column and retry reduced covariate sets after QR-based rank reduction and
#' correlation-based pruning. Extreme finite coefficients / standard errors are
#' rejected and treated as non-estimable.
#'
#' @keywords internal
InferenceAbstractKKLWACoxIVWC = R6::R6Class("InferenceAbstractKKLWACoxIVWC",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = c(InferenceMixinKKPassThrough$public, list(
		#' @description Initialize the inference object.
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
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			private$init_kk_passthrough(des_obj)
		},
		#' @description Returns the estimated treatment effect (log-hazard ratio).
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Computes the asymptotic confidence interval.
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
		#' @description Computes the asymptotic p-value.
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
	)),
	private = c(InferenceMixinKKPassThrough$private, list(
		compute_basic_match_data = function() private$compute_basic_kk_match_data_impl(),
		supports_likelihood_tests = function() FALSE,
		max_abs_reasonable_coef = 1e4,
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			if (m > 0){
				private$lwa_cox_for_matched_pairs()
			}
			beta_m = private$cached_values$beta_T_matched
			ssq_m  = private$cached_values$ssq_beta_T_matched
			m_ok   = !is.null(beta_m) && is.finite(beta_m) &&
			         !is.null(ssq_m) && is.finite(ssq_m) && ssq_m > 0
			if (nRT > 0 && nRC > 0){
				private$cox_for_reservoir()
			}
			beta_r = private$cached_values$beta_T_reservoir
			ssq_r  = private$cached_values$ssq_beta_T_reservoir
			r_ok   = !is.null(beta_r) && is.finite(beta_r) &&
			         !is.null(ssq_r) && is.finite(ssq_r) && ssq_r > 0
			if (m_ok && r_ok){
				w_star = ssq_r / (ssq_r + ssq_m)
				private$cached_values$beta_hat_T   = w_star * beta_m + (1 - w_star) * beta_r
				if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = sqrt(ssq_m * ssq_r / (ssq_m + ssq_r))
			} else if (m_ok){
				private$cached_values$beta_hat_T   = beta_m
				private$cached_values$s_beta_hat_T = sqrt(ssq_m)
			} else if (r_ok){
				private$cached_values$beta_hat_T   = beta_r
				private$cached_values$s_beta_hat_T = sqrt(ssq_r)
			} else {
				private$cache_nonestimable_estimate("kk_lwa_cox_ivwc_no_usable_component")
				return(invisible(NULL))
			}
			if (is.finite(private$cached_values$beta_hat_T) &&
			    abs(private$cached_values$beta_hat_T) > private$max_abs_reasonable_coef){
				private$cache_nonestimable_estimate("kk_lwa_cox_ivwc_extreme_estimate")
				return(invisible(NULL))
			}
			if (!estimate_only &&
			    (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0 ||
			     private$cached_values$s_beta_hat_T > private$max_abs_reasonable_coef)){
				private$cache_nonestimable_se("kk_lwa_cox_ivwc_standard_error_unavailable")
				return(invisible(NULL))
			}
			private$clear_nonestimable_state()
		},
		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},
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
		# Fit Cox model with optional cluster-robust variance.
		# cluster: integer vector of cluster IDs (length = length(y)), or NULL for standard Cox.
		fit_cox_model = function(y, dead, w, X, cluster = NULL){
			if (length(y) == 0L || sum(dead) == 0L) return(NULL)
			cl_int = if (!is.null(cluster)) as.integer(cluster) else NULL
			for (X_full in private$cox_design_candidates(w, X)){
				res = tryCatch(
					fast_coxph_regression_cpp(X_full, y, dead, cluster = cl_int),
					error = function(e) NULL
				)
				if (is.null(res)) next
				beta = tryCatch(res$coefficients[1L], error = function(e) NA_real_)
				se   = tryCatch(sqrt(res$vcov[1L, 1L]), error = function(e) NA_real_)
				if (!is.finite(beta) || abs(beta) > private$max_abs_reasonable_coef ||
				    !is.finite(se) || se <= 0 || se > private$max_abs_reasonable_coef) next
				return(list(beta = beta, ssq = se^2))
			}
			NULL
		},
		lwa_cox_for_matched_pairs = function(){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			i_matched = which(m_vec > 0)
			if (length(i_matched) == 0L) return(invisible(NULL))
			fit = private$fit_cox_model(
				y = private$y[i_matched],
				dead = private$dead[i_matched],
				w = private$w[i_matched],
				X = private$get_X()[i_matched, drop = FALSE],
				cluster = m_vec[i_matched]
			)
			if (is.null(fit)) return(invisible(NULL))
			private$cached_values$beta_T_matched = fit$beta
			private$cached_values$ssq_beta_T_matched = fit$ssq
		},
		cox_for_reservoir = function(){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L
			i_reservoir = which(m_vec == 0L)
			if (length(i_reservoir) == 0L) return(invisible(NULL))
			fit = private$fit_cox_model(
				y = private$y[i_reservoir],
				dead = private$dead[i_reservoir],
				w = private$w[i_reservoir],
				X = private$get_X()[i_reservoir, drop = FALSE]
			)
			if (is.null(fit)) return(invisible(NULL))
			private$cached_values$beta_T_reservoir = fit$beta
			private$cached_values$ssq_beta_T_reservoir = fit$ssq
		}
	))
)
