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
	inherit = InferenceAsympLik,
	public = as.list(modifyList(as.list(InferenceMixinKKPassThrough$public), list(
		#' @description Initialize the inference object.
		#' @param des_obj  	A DesignSeqOneByOne object (must be a KK design).
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose  		Whether to print progress messages.
		#' @param smart_cold_start_default   Whether to use smart cold start values.
		initialize = function(des_obj, model_formula = NULL,  verbose = FALSE, smart_cold_start_default = TRUE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
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
		},
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect.
		#' @param B  					Number of bootstrap samples.
		#' @param show_progress Whether to show a progress bar.
		#' @param debug         Whether to return diagnostics.
		#' @param bootstrap_type Optional resampling scheme.
		#' @return A numeric vector of bootstrap estimates.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			InferenceMixinKKPassThrough$public$approximate_bootstrap_distribution_beta_hat_T(B, show_progress, debug, bootstrap_type)
		},
		#' @description Gated off for IVWC.
		#' @param subject_or_block_weights weights.
		#' @param estimate_only flag.
		compute_estimate_with_bootstrap_weights = function(subject_or_block_weights, estimate_only = FALSE) {
			stop_bayesian_bootstrap_for_ivwc(self)
		},
		#' @description Gated off for IVWC.
		#' @param B replicates.
		#' @param show_progress flag.
		#' @param debug flag.
		#' @param weighting_unit_type type.
		approximate_bayesian_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, weighting_unit_type = NULL) {
			stop_bayesian_bootstrap_for_ivwc(self)
		},
		#' @description Gated off for IVWC.
		#' @param delta null.
		#' @param B replicates.
		#' @param type type.
		#' @param na.rm flag.
		#' @param show_progress flag.
		#' @param min_number_usable_samples count.
		#' @param weighting_unit_type type.
		compute_bayesian_bootstrap_two_sided_pval = function(delta = 0, B = 501, type = NULL, na.rm = FALSE, show_progress = TRUE, min_number_usable_samples = 5L, weighting_unit_type = NULL) {
			stop_bayesian_bootstrap_for_ivwc(self)
		},
		#' @description Gated off for IVWC.
		#' @param alpha level.
		#' @param B replicates.
		#' @param type type.
		#' @param na.rm flag.
		#' @param show_progress flag.
		#' @param min_number_usable_samples count.
		#' @param weighting_unit_type type.
		compute_bayesian_bootstrap_confidence_interval = function(alpha = 0.05, B = 501, type = NULL, na.rm = TRUE, show_progress = TRUE, min_number_usable_samples = 5L, weighting_unit_type = NULL) {
			stop_bayesian_bootstrap_for_ivwc(self)
		}
	))),
	private = as.list(modifyList(as.list(InferenceMixinKKPassThrough$private), list(
		compute_basic_match_data = function() private$compute_basic_kk_match_data_impl(),
		supports_likelihood_tests = function() FALSE,
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
				X_try_df = as.data.frame(X_try)
				key = paste(colnames(X_try_df), collapse = "|")
				if (key %in% keys) next
				keys = c(keys, key)
				attempt_try = private$fit_with_hardened_qr_column_dropping(
					X_full = X_try,
					required_cols = 1L,
					fit_fun = function(X_fit) X_fit,
					fit_ok = function(mod, X_fit, keep) TRUE
				)
				candidates[[length(candidates) + 1L]] = attempt_try$X
			}
			candidates
		},
		rcpp_cox_fit_is_usable = function(fit, estimate_only = FALSE){
			if (is.null(fit) || !isTRUE(fit$converged)) return(FALSE)
			coef1 = fit$coefficients[1L]
			if (!is.finite(coef1) || abs(coef1) > private$max_abs_reasonable_coef) return(FALSE)
			if (estimate_only) return(TRUE)
			se = sqrt(fit$vcov[1L, 1L])
			is.finite(se) && se > 0 && se <= private$max_abs_reasonable_coef
		},
		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))
			private$clear_nonestimable_state()
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}
			KKstats = private$cached_values$KKstats
			if (is.null(KKstats)) return(invisible(NULL))
			m = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			if (sum(private$dead) == 0L){
				private$cache_nonestimable_estimate("kk_strat_cox_ivwc_no_events")
				return(invisible(NULL))
			}
			# ── Matched pairs ───────────────────────────────────────────────
			beta_m = NA_real_; ssq_m = NA_real_
			if (m > 0){
				y_m = KKstats$y_matched_long
				dead_m = KKstats$dead_matched_long
				w_m = KKstats$w_matched_long
				strata_m = KKstats$m_matched_long
				X_cov_m = KKstats$X_matched_long
				X_m = cbind(w = w_m, X_cov_m)
				
				res = tryCatch(
					fast_stratified_coxph_regression_cpp(X_m, y_m, dead_m, as.integer(strata_m), estimate_only = estimate_only),
					error = function(e) NULL
				)
				if (private$rcpp_cox_fit_is_usable(res, estimate_only = estimate_only)){
					beta_m = res$coefficients[1L]
					if (!estimate_only) ssq_m = res$vcov[1L, 1L]
				}
			}
			# ── Reservoir ───────────────────────────────────────────────────
			beta_r = NA_real_; ssq_r = NA_real_
			res = NULL
			if (nRT > 0 && nRC > 0){
				y_r = KKstats$y_reservoir
				dead_r = KKstats$dead_reservoir
				w_r = KKstats$w_reservoir
				X_cov_r = as.matrix(KKstats$X_reservoir)
				candidates = private$cox_design_candidates(w_r, X_cov_r)
				for (X_candidate in candidates){
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
				X_mat = matrix(private$w[private$m == 0], ncol = 1L)
				colnames(X_mat) = "w"
				res = tryCatch(
					fast_coxph_regression_cpp(X_mat, y_r, dead_r, estimate_only = estimate_only),
					error = function(e) NULL
				)
			}
			if (!is.null(res)){
				beta_r = res$coefficients[1L]
				if (!estimate_only) ssq_r = res$vcov[1L, 1L]
			}
			# ── Combine ─────────────────────────────────────────────────────
			m_ok = is.finite(beta_m) && (estimate_only || (is.finite(ssq_m) && ssq_m > 0))
			r_ok = is.finite(beta_r) && (estimate_only || (is.finite(ssq_r) && ssq_r > 0))
			if (m_ok && r_ok){
				if (estimate_only){
					# Use simple sample size weighting as fallback if SEs not requested
					w_star = m / (m + (nRT+nRC)/2)
					private$cached_values$beta_hat_T = w_star * beta_m + (1 - w_star) * beta_r
				} else {
					w_star = ssq_r / (ssq_r + ssq_m)
					private$cached_values$beta_hat_T = w_star * beta_m + (1 - w_star) * beta_r
					private$cached_values$s_beta_hat_T = sqrt(ssq_m * ssq_r / (ssq_m + ssq_r))
				}
			} else if (m_ok){
				private$cached_values$beta_hat_T = beta_m
				if (!estimate_only) private$cached_values$s_beta_hat_T = sqrt(ssq_m)
			} else if (r_ok){
				private$cached_values$beta_hat_T = beta_r
				if (!estimate_only) private$cached_values$s_beta_hat_T = sqrt(ssq_r)
			} else {
				private$cache_nonestimable_estimate("kk_strat_cox_ivwc_both_failed")
			}
		}
	)))
)

#' Abstract class for Stratified Cox Combined-Likelihood Inference
#'
#' @keywords internal
#' @noRd
InferenceAbstractKKStratCoxOneLik = R6::R6Class("InferenceAbstractKKStratCoxOneLik",
	lock_objects = FALSE,
	inherit = InferenceAsympLik,
	public = as.list(modifyList(as.list(InferenceMixinKKPassThrough$public), list(
		#' @description Initialize the inference object.
		initialize = function(des_obj, model_formula = NULL,  verbose = FALSE, smart_cold_start_default = TRUE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula, smart_cold_start_default = smart_cold_start_default)
			private$init_kk_passthrough(des_obj)
		},
		#' @description Returns the combined-likelihood estimate of the treatment effect.
		compute_estimate = function(estimate_only = FALSE){
			private$shared_combined_likelihood(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},
		#' @description Computes an asymptotic confidence interval for the treatment effect.
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
		#' @description Returns a 2-sided p-value for H0: beta_T = delta.
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared_combined_likelihood(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},
		#' @description Creates the bootstrap distribution of the estimate for the treatment effect.
		approximate_bootstrap_distribution_beta_hat_T = function(B = 501, show_progress = TRUE, debug = FALSE, bootstrap_type = NULL){
			InferenceMixinKKPassThrough$public$approximate_bootstrap_distribution_beta_hat_T(B, show_progress, debug, bootstrap_type)
		}
	))),
	private = as.list(modifyList(as.list(InferenceMixinKKPassThrough$private), list(
		compute_basic_match_data = function() private$compute_basic_kk_match_data_impl(),
		max_abs_reasonable_coef = 1e4,
		best_X_colnames = NULL,
		optimization_alg = "lbfgs",
		shared_combined_likelihood = function(estimate_only = FALSE){
			# Placeholder for combined strat Cox logic
		}
	)))
)

#' Stratified Cox / Standard Cox Compound Inference for KK Designs
#' @export
InferenceSurvivalKKStratCoxIVWC = R6::R6Class("InferenceSurvivalKKStratCoxIVWC",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKStratCoxIVWC,
	public = list()
)

#' Stratified Cox Combined-Likelihood Compound Inference for KK Designs
#' @export
InferenceSurvivalKKStratCoxOneLik = R6::R6Class("InferenceSurvivalKKStratCoxOneLik",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKStratCoxOneLik,
	public = list()
)
