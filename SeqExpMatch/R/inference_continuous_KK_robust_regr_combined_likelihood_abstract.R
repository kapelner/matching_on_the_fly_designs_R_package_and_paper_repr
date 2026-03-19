# Abstract class for Robust-Regression Combined-Likelihood Inference for KK Designs
#
# @description
# Fits a single stacked robust regression over matched-pair differences and reservoir
# observations for KK matching-on-the-fly designs with continuous responses.
#
# @keywords internal
SeqDesignInferenceAbstractKKRobustRegrCombinedLikelihood = R6::R6Class("SeqDesignInferenceAbstractKKRobustRegrCombinedLikelihood",
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
		# Returns the combined robust-regression estimate of the treatment effect.
		compute_treatment_estimate = function(){
			private$fit_combined()
			private$cached_values$beta_hat_T
		},

		# @description
		# Computes the approximate confidence interval.
		# @param alpha The confidence level in the computed confidence interval is 1 - \code{alpha}. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$fit_combined()
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		# @description
		# Computes the approximate p-value.
		# @param delta The null difference to test against. For any treatment effect at all this is set to zero (the default).
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$fit_combined()
			private$assert_finite_se()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},

		duplicate = function(verbose = FALSE){
			i = super$duplicate(verbose = verbose)
			i
		}
	),

	private = list(
		rlm_method = NULL,
		rlm_force_M = FALSE,

		include_covariates = function() stop(class(self)[1], " must implement include_covariates()"),

		reduce_design_matrix_once = function(X, j_treat, cache_key){
			cached = private$cached_values[[cache_key]]
			if (!is.null(cached)) return(cached)

			qr_X = qr(X)
			if (qr_X$rank < ncol(X)){
				keep = qr_X$pivot[seq_len(qr_X$rank)]
				if (!(j_treat %in% keep)) keep[qr_X$rank] = j_treat
				keep = sort(unique(keep))
				X = X[, keep, drop = FALSE]
				j_treat = which(keep == j_treat)
			}

			cached = list(X = X, j_treat = j_treat)
			private$cached_values[[cache_key]] = cached
			cached
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				stop("KK robust-regression combined-likelihood estimator: could not compute a finite standard error.")
			}
		},

		fit_rlm = function(X, y, j_treat){
			if (nrow(X) <= ncol(X)) return(NULL)

			run_rlm = function(method){
				tryCatch({
					if (identical(method, "M")) {
						MASS::rlm(x = X, y = y, method = "M", psi = MASS::psi.huber)
					} else {
						MASS::rlm(x = X, y = y, method = method)
					}
				}, error = function(e) e)
			}

			method_to_try = if (isTRUE(private$rlm_force_M)) "M" else private$rlm_method
			mod = run_rlm(method_to_try)
			if (inherits(mod, "error") && identical(method_to_try, "MM")){
				msg = if (length(mod$message) == 0L) "" else mod$message
				if (grepl("'lqs' failed", msg, fixed = TRUE) || grepl("singular", msg, ignore.case = TRUE)) {
					private$rlm_force_M = TRUE
					mod = run_rlm("M")
				}
			}
			if (inherits(mod, "error")) return(NULL)
			if (is.null(mod)) return(NULL)

			coef_table = tryCatch(summary(mod)$coefficients, error = function(e) NULL)
			if ((is.null(coef_table) || nrow(coef_table) < j_treat) && identical(method_to_try, "MM")){
				private$rlm_force_M = TRUE
				mod = run_rlm("M")
				if (inherits(mod, "error") || is.null(mod)) return(NULL)
				coef_table = tryCatch(summary(mod)$coefficients, error = function(e) NULL)
			}
			if (is.null(coef_table) || nrow(coef_table) < j_treat) return(NULL)
			beta = as.numeric(coef_table[j_treat, "Value"])
			se   = as.numeric(coef_table[j_treat, "Std. Error"])
			if ((!is.finite(beta) || !is.finite(se) || se <= 0) && identical(method_to_try, "MM")){
				private$rlm_force_M = TRUE
				mod = run_rlm("M")
				if (inherits(mod, "error") || is.null(mod)) return(NULL)
				coef_table = tryCatch(summary(mod)$coefficients, error = function(e) NULL)
				if (is.null(coef_table) || nrow(coef_table) < j_treat) return(NULL)
				beta = as.numeric(coef_table[j_treat, "Value"])
				se   = as.numeric(coef_table[j_treat, "Std. Error"])
			}
			if (!is.finite(beta) || !is.finite(se) || se <= 0) return(NULL)

			list(beta = beta, se = se, mod = mod)
		},

		fit_combined = function(){
			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))

			KKstats = private$cached_values$KKstats
			if (is.null(KKstats)){
				private$compute_basic_match_data()
				KKstats = private$cached_values$KKstats
			}

			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC
			nR  = nRT + nRC

			if (m > 0 && nRT > 0 && nRC > 0){
				if (private$include_covariates()){
					Xd_full = as.matrix(KKstats$X_matched_diffs_full)
					p = ncol(Xd_full)
				} else {
					Xd_full = matrix(nrow = m, ncol = 0L)
					p = 0L
				}
				X_comb = rbind(
					if (p > 0L) cbind(0, 1, Xd_full) else matrix(c(0, 1), nrow = m, ncol = 2L, byrow = TRUE),
					if (p > 0L) cbind(rep(1, nR), KKstats$w_reservoir, as.matrix(KKstats$X_reservoir)) else cbind(rep(1, nR), KKstats$w_reservoir)
				)
				y_comb = c(KKstats$y_matched_diffs, KKstats$y_reservoir)
				j_treat = 2L
			} else if (m > 0){
				Xd = if (private$include_covariates()) as.matrix(KKstats$X_matched_diffs) else matrix(nrow = m, ncol = 0L)
				X_comb = if (ncol(Xd) > 0L) cbind(1, Xd) else matrix(1, nrow = m, ncol = 1L)
				y_comb = KKstats$y_matched_diffs
				j_treat = 1L
			} else if (nRT > 0 && nRC > 0){
				X_r = if (private$include_covariates()) as.matrix(KKstats$X_reservoir) else matrix(nrow = nR, ncol = 0L)
				X_comb = if (ncol(X_r) > 0L) cbind(rep(1, nR), KKstats$w_reservoir, X_r) else cbind(rep(1, nR), KKstats$w_reservoir)
				y_comb = KKstats$y_reservoir
				j_treat = 2L
			} else {
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			reduced = private$reduce_design_matrix_once(
				X_comb,
				j_treat,
				cache_key = "kk_robust_combined_reduced_design"
			)
			X_comb = reduced$X
			j_treat = reduced$j_treat

			fit = private$fit_rlm(X_comb, y_comb, j_treat)
			if (is.null(fit)){
				private$cached_values$beta_hat_T   = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
				private$cached_values$is_z         = TRUE
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T   = fit$beta
			private$cached_values$s_beta_hat_T = fit$se
			private$cached_values$is_z         = TRUE
			private$cached_values$full_coefficients = stats::coef(fit$mod)
			private$cached_values$full_vcov = tryCatch(stats::vcov(fit$mod), error = function(e) NULL)
			invisible(NULL)
		}
	)
)
