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
		#' @param des_obj		A DesignSeqOneByOne object (must be a KK design).
		#' @param verbose			Whether to print progress messages.
		initialize = function(des_obj,  verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "survival")
			}
			if (should_run_asserts()) {
				if (!is(des_obj, "DesignSeqOneByOneKK14") && !is(des_obj, "FixedDesignBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
				}
			}
			super$initialize(des_obj, verbose)
		},

		#' @description
		#' Returns the estimated treatment effect (log-hazard ratio).
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
			if (should_run_asserts()) {
				if (delta == 0){
					private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
				} else {
					stop("Testing non-zero delta is not yet implemented for this class.")
				}
			}
		}
	),

	private = list(
		max_abs_reasonable_coef = 1e4,

		# Abstract: subclasses return TRUE (multivariate) or FALSE (univariate).
		include_covariates = function() stop(class(self)[1], " must implement include_covariates()"),

		cox_design_candidates = function(w, X){
			X_full = cbind(w = w, as.matrix(X))
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
				X_try = cbind(w = w, X_cov)
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

		cox_fit_is_usable = function(mod){
			if (is.null(mod)) return(FALSE)
			beta = tryCatch(as.numeric(coef(mod)["w"]), error = function(e) NA_real_)
			se = tryCatch(sqrt(as.numeric(vcov(mod)["w", "w"])), error = function(e) NA_real_)
			is.finite(beta) && abs(beta) <= private$max_abs_reasonable_coef &&
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
				private$strat_cox_for_matched_pairs()
			}
			beta_m   = private$cached_values$beta_T_matched
			ssq_m    = private$cached_values$ssq_beta_T_matched
			m_ok     = !is.null(beta_m) && is.finite(beta_m) &&
			           !is.null(ssq_m)  && is.finite(ssq_m) && ssq_m > 0

			# --- Reservoir: Standard Cox ---
			if (nRT > 0 && nRC > 0){
				private$cox_for_reservoir()
			}
			beta_r   = private$cached_values$beta_T_reservoir
			ssq_r    = private$cached_values$ssq_beta_T_reservoir
			r_ok     = !is.null(beta_r) && is.finite(beta_r) &&
			           !is.null(ssq_r)  && is.finite(ssq_r) && ssq_r > 0

			# --- Variance-weighted combination ---
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
				private$cache_nonestimable_estimate("kk_strat_cox_ivwc_no_usable_component")
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}
			if (is.finite(private$cached_values$beta_hat_T) &&
			    abs(private$cached_values$beta_hat_T) > private$max_abs_reasonable_coef){
				private$cache_nonestimable_estimate("kk_strat_cox_ivwc_extreme_estimate")
				private$cached_values$is_z = TRUE
				return(invisible(NULL))
			}
			if (!estimate_only &&
			    (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0 ||
			     private$cached_values$s_beta_hat_T > private$max_abs_reasonable_coef)){
				private$cache_nonestimable_se("kk_strat_cox_ivwc_standard_error_unavailable")
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

		strat_cox_for_matched_pairs = function(){
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

			dat = data.frame(y = y_m[i_valid], dead = dead_m[i_valid], w = w_m[i_valid], strata = strata_m[i_valid])
			formula_str = "survival::Surv(y, dead) ~ w + strata(strata)"

			mod = NULL
			if (private$include_covariates()){
				X_m = as.matrix(private$get_X()[i_matched[i_valid], , drop = FALSE])
				for (X_candidate in private$cox_design_candidates(dat$w, X_m)){
					dat_try = dat
					formula_try = formula_str
					X_covs = X_candidate[, colnames(X_candidate) != "w", drop = FALSE]
					if (ncol(X_covs) > 0){
						colnames(X_covs) = paste0("x", seq_len(ncol(X_covs)))
						dat_try = cbind(dat_try, X_covs)
						formula_try = paste(formula_try, "+", paste(colnames(X_covs), collapse = " + "))
					}
					mod_try = tryCatch(
						survival::coxph(as.formula(formula_try), data = dat_try),
						error = function(e) NULL,
						warning = function(w) {
							if (grepl("converge|infinite", w$message)) return(NULL)
							invokeRestart("muffleWarning")
						}
					)
					if (private$cox_fit_is_usable(mod_try)){
						mod = mod_try
						break
					}
				}
			} else {
				mod = tryCatch(
					survival::coxph(as.formula(formula_str), data = dat),
					error = function(e) NULL,
					warning = function(w) {
						if (grepl("converge|infinite", w$message)) return(NULL)
						invokeRestart("muffleWarning")
					}
				)
			}
			if (is.null(mod)) return(invisible(NULL))

			beta = as.numeric(coef(mod)["w"])
			se   = sqrt(as.numeric(vcov(mod)["w", "w"]))
			
			# If estimate is clearly degenerate, treat as failure
			if (!is.finite(beta) || !is.finite(se) || se <= 0 || se > private$max_abs_reasonable_coef || abs(beta) > private$max_abs_reasonable_coef) return(invisible(NULL))

			private$cached_values$beta_T_matched     = beta
			private$cached_values$ssq_beta_T_matched = se^2
		},

		cox_for_reservoir = function(){
			private$cached_values$beta_T_reservoir = NULL
			private$cached_values$ssq_beta_T_reservoir = NULL

			y_r    = private$cached_values$KKstats$y_reservoir
			w_r    = private$cached_values$KKstats$w_reservoir
			m_vec_safe = private$m
			if (is.null(m_vec_safe)) m_vec_safe = rep(0L, private$n)
			m_vec_safe[is.na(m_vec_safe)] = 0L
			dead_r = private$dead[m_vec_safe == 0]
			X_r    = as.matrix(private$cached_values$KKstats$X_reservoir)

			dat = data.frame(y = y_r, dead = dead_r, w = w_r)
			formula_str = "survival::Surv(y, dead) ~ w"

			mod = NULL
			if (private$include_covariates()){
				for (X_candidate in private$cox_design_candidates(dat$w, X_r)){
					dat_try = dat
					formula_try = formula_str
					X_covs = X_candidate[, colnames(X_candidate) != "w", drop = FALSE]
					if (ncol(X_covs) > 0){
						colnames(X_covs) = paste0("x", seq_len(ncol(X_covs)))
						dat_try = cbind(dat_try, X_covs)
						formula_try = paste(formula_try, "+", paste(colnames(X_covs), collapse = " + "))
					}
					mod_try = tryCatch(
						survival::coxph(as.formula(formula_try), data = dat_try),
						error = function(e) NULL,
						warning = function(w) {
							if (grepl("converge|infinite", w$message)) return(NULL)
							invokeRestart("muffleWarning")
						}
					)
					if (private$cox_fit_is_usable(mod_try)){
						mod = mod_try
						break
					}
				}
			}
			if (is.null(mod)){
				mod = tryCatch(
					survival::coxph(survival::Surv(y, dead) ~ w, data = data.frame(y = y_r, dead = dead_r, w = w_r)),
					error = function(e) NULL,
					warning = function(w) {
						if (grepl("converge|infinite", w$message)) return(NULL)
						invokeRestart("muffleWarning")
					}
				)
			}

			if (is.null(mod)) return(invisible(NULL))

			beta = as.numeric(coef(mod)["w"])
			se   = sqrt(as.numeric(vcov(mod)["w", "w"]))
			
			if (!is.finite(beta) || !is.finite(se) || se <= 0 || se > private$max_abs_reasonable_coef || abs(beta) > private$max_abs_reasonable_coef) return(invisible(NULL))

			private$cached_values$beta_T_reservoir     = beta
			private$cached_values$ssq_beta_T_reservoir = se^2
		}
	)
)
