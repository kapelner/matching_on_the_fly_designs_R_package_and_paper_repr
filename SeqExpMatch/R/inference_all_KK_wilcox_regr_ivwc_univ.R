#' Univariate Wilcox Rank-based Regression Compound Inference for KK Designs
#'
#' @description
#' Fits a robust compound estimator for KK matching-on-the-fly designs using rank-based
#' regression (R-estimation) with only the treatment indicator (no additional covariates).
#' For matched pairs, it uses R-estimation on the within-pair response differences.
#' For reservoir subjects, it uses a standard rank-based linear model with the treatment
#' assignment as the sole predictor. This is the univariate (covariate-free) variant.
#'
#' @export
InferenceAllKKWilcoxRegrUnivIVWC = R6::R6Class("InferenceAllKKWilcoxRegrUnivIVWC",
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj A DesignSeqOneByOne object (must be a KK design).
		#' @param num_cores Number of CPU cores for parallel processing.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, num_cores = 1, verbose = FALSE, make_fork_cluster = NULL){
			res_type = des_obj$get_response_type()
			if (res_type == "incidence"){
				stop("Rank-based regression is not recommended for incidence data; clogit and compound mean diff is recommended.")
			}
			assertResponseType(res_type, c("continuous", "count", "proportion", "survival", "ordinal"))
			if (!is(des_obj, "DesignSeqOneByOneKK14")){
				stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
			}
			super$initialize(des_obj, num_cores, verbose, make_fork_cluster = make_fork_cluster)
			if (private$any_censoring){
				stop(class(self)[1], " does not support censored survival data.")
			}
			if (!requireNamespace("Rfit", quietly = TRUE)) {
				stop("Package 'Rfit' is required for ", class(self)[1], ". Please install it.")
			}
		},

		#' @description
		#' Returns the estimated treatment effect.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes the asymptotic confidence interval.
		#' @param alpha The confidence level is 1 - \code{alpha}. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes the asymptotic two-sided p-value.
		#' @param delta The null difference to test against. Default is zero.
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

		include_covariates = function(){
			FALSE
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			# Recompute KKstats if cache was cleared (e.g., after y transformation for rand CI)
			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}

			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			# --- Matched pairs: Rfit on differences ---
			if (m > 0){
				private$rfit_for_matched_pairs()
			}
			beta_m   = private$cached_values$beta_T_matched
			ssq_m    = private$cached_values$ssq_beta_T_matched
			m_ok     = !is.null(beta_m) && is.finite(beta_m) &&
			           !is.null(ssq_m)  && is.finite(ssq_m) && ssq_m > 0

			# --- Reservoir: Rfit on level data ---
			if (nRT > 0 && nRC > 0){
				private$rfit_for_reservoir()
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
				stop("Rank-regression compound estimator: could not compute a finite standard error.")
			}
		},

		rfit_for_matched_pairs = function(){
			KKstats = private$cached_values$KKstats
			dy = KKstats$y_matched_diffs

			if (private$include_covariates()){
				dX = KKstats$X_matched_diffs
				# Formula: dy ~ 1 + dX. Intercept is the treatment effect.
				dat = as.data.frame(dX)
				colnames(dat) = paste0("x", 1:ncol(dat))
				dat$dy = dy
				mod = tryCatch({
					Rfit::rfit(dy ~ ., data = dat)
				}, error = function(e) NULL)

				if (is.null(mod)) return(invisible(NULL))

				summ = tryCatch(summary(mod), error = function(e) NULL)
				if (is.null(summ) || is.null(summ$coefficients)) return(invisible(NULL))
				if (!"(Intercept)" %in% rownames(summ$coefficients)) return(invisible(NULL))

				beta = as.numeric(summ$coefficients["(Intercept)", "Estimate"])
				se   = as.numeric(summ$coefficients["(Intercept)", "Std. Error"])

				private$cached_values$beta_T_matched     = if (length(beta) == 1L && is.finite(beta)) beta else NA_real_
				private$cached_values$ssq_beta_T_matched = if (length(se) == 1L && is.finite(se) && se > 0) se^2 else NA_real_
			} else {
				# Univariate case: Rfit fails on dy ~ 1. Use Wilcoxon HL estimate.
				# This is mathematically consistent with the RankRegr goal.
				if (all(dy == 0)) return(invisible(NULL))
				mod = tryCatch({
					stats::wilcox.test(dy, conf.int = TRUE)
				}, error = function(e) NULL)

				if (is.null(mod)) return(invisible(NULL))

				beta = as.numeric(mod$estimate)
				se   = (as.numeric(mod$conf.int[2]) - as.numeric(mod$conf.int[1])) / (2 * 1.96)

				private$cached_values$beta_T_matched     = if (is.finite(beta)) beta else NA_real_
				private$cached_values$ssq_beta_T_matched = if (is.finite(se) && se > 0) se^2 else NA_real_
			}
		},

		rfit_for_reservoir = function(){
			y_r = private$cached_values$KKstats$y_reservoir
			w_r = private$cached_values$KKstats$w_reservoir
			X_r = as.matrix(private$cached_values$KKstats$X_reservoir)

			if (private$include_covariates()){
				dat = data.frame(y = y_r, w = w_r)
				X_covs = as.data.frame(X_r)
				colnames(X_covs) = paste0("x", 1:ncol(X_covs))
				dat = cbind(dat, X_covs)

				mod = tryCatch({
					Rfit::rfit(y ~ ., data = dat)
				}, error = function(e) NULL)

				if (is.null(mod)) return(invisible(NULL))

				summ = tryCatch(summary(mod), error = function(e) NULL)
				if (is.null(summ) || is.null(summ$coefficients)) return(invisible(NULL))
				if (!"w" %in% rownames(summ$coefficients)) return(invisible(NULL))

				beta = as.numeric(summ$coefficients["w", "Estimate"])
				se   = as.numeric(summ$coefficients["w", "Std. Error"])

				private$cached_values$beta_T_reservoir     = if (length(beta) == 1L && is.finite(beta)) beta else NA_real_
				private$cached_values$ssq_beta_T_reservoir = if (length(se) == 1L && is.finite(se) && se > 0) se^2 else NA_real_
			} else {
				# Univariate case: use Wilcoxon rank-sum HL estimate
				yT = y_r[w_r == 1]
				yC = y_r[w_r == 0]
				mod = tryCatch({
					stats::wilcox.test(yT, yC, conf.int = TRUE)
				}, error = function(e) NULL)

				if (is.null(mod)) return(invisible(NULL))

				beta = as.numeric(mod$estimate)
				se   = (as.numeric(mod$conf.int[2]) - as.numeric(mod$conf.int[1])) / (2 * 1.96)

				private$cached_values$beta_T_reservoir     = if (is.finite(beta)) beta else NA_real_
				private$cached_values$ssq_beta_T_reservoir = if (is.finite(se) && se > 0) se^2 else NA_real_
			}
		}
	)
)
