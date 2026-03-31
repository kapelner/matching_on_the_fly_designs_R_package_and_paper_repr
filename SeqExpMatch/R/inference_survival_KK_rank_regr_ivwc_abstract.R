#' Abstract class for Survival Rank-based Regression (AFT) Compound Inference
#'
#' This class implements a robust compound estimator for KK matching-on-the-fly
#' designs with survival responses using rank-based estimating equations via the
#' \pkg{aftgee} package. For matched pairs, it fits a rank-based AFT model with
#' clustering. For reservoir subjects, it fits a standard rank-based AFT model.
#' The two estimates (both log-time ratios) are combined via a variance-weighted
#' linear combination.
#'
#' @details
#' This class requires the \pkg{aftgee} package.
#'
#' @keywords internal
InferenceAbstractKKSurvivalRankRegrIVWC = R6::R6Class("InferenceAbstractKKSurvivalRankRegrIVWC",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj		A DesignSeqOneByOne object (must be a KK design).
		#' @param num_cores			Number of CPU cores for parallel processing.
		#' @param verbose			Whether to print progress messages.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE){
			res_type = des_obj$get_response_type()
			if (res_type == "incidence"){
				stop("Rank-based regression is not recommended for incidence data; clogit and compound mean diff is recommended.")
			}
			assertResponseType(res_type, "survival")
			if (!is(des_obj, "DesignSeqOneByOneKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
			}
			super$initialize(des_obj, num_cores, verbose)
			if (!requireNamespace("aftgee", quietly = TRUE)) {
				stop("Package 'aftgee' is required for ", class(self)[1], ". Please install it.")
			}
		},

		#' @description
		#' Returns the estimated treatment effect (log-time ratio).
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
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes the asymptotic p-value.
		#' @param delta                                   The null difference to test against. For
		#'   any treatment effect at all this is set to zero (the default).
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared()
			private$assert_finite_se()
			if (delta == 0){
				private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
			} else {
				stop("Testing non-zero delta is not yet implemented for the combined rank-regression estimator.")
			}
		}
	),

	private = list(

		# Abstract: subclasses return TRUE (multivariate) or FALSE (univariate).
		include_covariates = function() stop(class(self)[1], " must implement include_covariates()"),

		extract_term_estimate = function(mod, term_name = "w"){
			coefs = tryCatch(stats::coef(mod), error = function(e) NULL)
			if (is.null(coefs) || is.null(names(coefs)) || !(term_name %in% names(coefs))){
				return(NA_real_)
			}
			as.numeric(coefs[[term_name]])
		},

		extract_term_se = function(mod, term_name = "w"){
			coef_table = tryCatch(summary(mod)$coefficients, error = function(e) NULL)
			if (is.null(coef_table) || is.null(dim(coef_table))){
				return(NA_real_)
			}
			if (is.null(rownames(coef_table)) || !(term_name %in% rownames(coef_table))){
				return(NA_real_)
			}
			se_col = intersect(colnames(coef_table), c("StdErr", "Std.Err", "Std.err", "Std.Error"))[1]
			if (is.na(se_col) || length(se_col) == 0L){
				return(NA_real_)
			}
			as.numeric(coef_table[term_name, se_col])
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			# --- Matched pairs: aftgee with clustering ---
			if (m > 0){
				private$aftgee_for_matched_pairs()
			}
			beta_m   = private$cached_values$beta_T_matched
			ssq_m    = private$cached_values$ssq_beta_T_matched
			m_ok     = !is.null(beta_m) && is.finite(beta_m) &&
			           !is.null(ssq_m)  && is.finite(ssq_m) && ssq_m > 0

			# --- Reservoir: aftgee (independent) ---
			if (nRT > 0 && nRC > 0){
				private$aftgee_for_reservoir()
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
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
			}
			private$cached_values$is_z = TRUE
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},

		aftgee_for_matched_pairs = function(){
			m_vec = private$m
			if (is.null(m_vec)) m_vec = rep(NA_integer_, private$n)
			m_vec[is.na(m_vec)] = 0L

			i_matched = which(m_vec > 0)
			y_m       = private$y[i_matched]
			dead_m    = private$dead[i_matched]
			w_m       = private$w[i_matched]
			strata_m  = m_vec[i_matched]

			# Filter strata that have no events (aftgee needs at least some events)
			if (sum(dead_m) < 2) return(invisible(NULL))

			dat = data.frame(y = y_m, dead = dead_m, w = w_m, id = strata_m)
			formula_str = "survival::Surv(y, dead) ~ w"

			if (private$include_covariates()){
				X_m = as.matrix(private$X[i_matched, , drop = FALSE])
				colnames(X_m) = paste0("x", 1:ncol(X_m))
				dat = cbind(dat, X_m)
				formula_str = paste(formula_str, "+", paste(colnames(X_m), collapse = " + "))
			}

			mod = tryCatch({
				# Use Gehan-type rank estimator (default)
				suppressMessages(aftgee::aftgee(as.formula(formula_str), id = id, data = dat, corstr = "exchangeable"))
			}, error = function(e) NULL)

			if (is.null(mod)) return(invisible(NULL))

			beta = private$extract_term_estimate(mod, "w")
			se   = private$extract_term_se(mod, "w")

			private$cached_values$beta_T_matched     = if (is.finite(beta)) beta else NA_real_
			private$cached_values$ssq_beta_T_matched = if (is.finite(se) && se > 0) se^2 else NA_real_
		},

		aftgee_for_reservoir = function(){
			y_r    = private$cached_values$KKstats$y_reservoir
			w_r    = private$cached_values$KKstats$w_reservoir
			dead_r = private$dead[private$m == 0]
			X_r    = as.matrix(private$cached_values$KKstats$X_reservoir)

			if (sum(dead_r) < 2) return(invisible(NULL))

			dat = data.frame(y = y_r, dead = dead_r, w = w_r)
			formula_str = "survival::Surv(y, dead) ~ w"

			if (private$include_covariates()){
				X_covs = as.data.frame(X_r)
				colnames(X_covs) = paste0("x", 1:ncol(X_covs))
				dat = cbind(dat, X_covs)
				formula_str = paste(formula_str, "+", paste(colnames(X_covs), collapse = " + "))
			}

			mod = tryCatch({
				suppressMessages(aftgee::aftgee(as.formula(formula_str), data = dat))
			}, error = function(e) NULL)

			if (is.null(mod)) return(invisible(NULL))

			beta = private$extract_term_estimate(mod, "w")
			se   = private$extract_term_se(mod, "w")

			private$cached_values$beta_T_reservoir     = if (is.finite(beta)) beta else NA_real_
			private$cached_values$ssq_beta_T_reservoir = if (is.finite(se) && se > 0) se^2 else NA_real_
		}
	)
)
