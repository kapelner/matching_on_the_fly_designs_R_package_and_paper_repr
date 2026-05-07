#' Abstract class for Conditional Poisson / Negative Binomial Compound Inference
#'
#' This class implements a compound estimator for KK matching-on-the-fly designs with
#' count responses. For matched pairs, it uses conditional Poisson regression
#' (equivalent to weighted logistic regression) and for the reservoir, it uses
#' standard negative binomial regression. The two estimates are combined via
#' inverse-variance weighted combination (IVWC).
#'
#' @keywords internal
InferenceAbstractKKPoissonCPoissonIVWC = R6::R6Class("InferenceAbstractKKPoissonCPoissonIVWC",
	lock_objects = FALSE,
	inherit = InferenceKKPassThrough,
	public = list(

		#' @description
		#' Initialize the inference object.
		#' @param des_obj		A DesignSeqOneByOne object (must be a KK design).
		#' @param model_formula   Optional formula for covariate adjustment. If \code{NULL} (default),
		#'   the formula from the design object is used and its pre-computed design matrix is
		#'   reused. If a formula is provided, a new design matrix is constructed from the
		#'   design's imputed covariates.
		#' @param verbose			Whether to print progress messages.
		initialize = function(des_obj, model_formula = NULL,  verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
			}
			if (should_run_asserts()) {
				if (!inherits(des_obj, "DesignSeqOneByOneKK14") && !inherits(des_obj, "FixedDesignBinaryMatch")){
					stop(class(self)[1], " requires a KK matching-on-the-fly design (DesignSeqOneByOneKK14 or subclass).")
				}
			}
			super$initialize(des_obj, verbose = verbose, model_formula = model_formula)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
		},

		#' @description
		#' Returns the estimated treatment effect (log-rate ratio).
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
		#' @param delta                                   The null difference to test against. Default is 0.
		compute_asymp_two_sided_pval = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared()
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		},

		#' @description
		#' Duplicates the object while preserving caches.
		#' @param verbose Whether the duplicate should be verbose.
		#' @param make_fork_cluster Whether the duplicate should be allowed to create a fork cluster.
		duplicate = function(verbose = FALSE, make_fork_cluster = FALSE){
			inf_obj = super$duplicate(verbose = verbose, make_fork_cluster = make_fork_cluster)
			inf_obj
		}
	),

	private = list(

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Use the same joint-likelihood logic for the point estimate
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		compute_treatment_estimate_during_power_simulation = function(y, permutations, delta, transform_responses, zero_one_logit_clamp = 1e-4){
			private$compute_treatment_estimate_during_power_simulation_worker(y, permutations, delta, transform_responses, zero_one_logit_clamp = zero_one_logit_clamp)
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T)){
				return(invisible(NULL))
			}
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (is.null(private$cached_values$KKstats)){
				private$compute_basic_match_data()
			}

			KKstats = private$cached_values$KKstats
			m   = KKstats$m
			nRT = KKstats$nRT
			nRC = KKstats$nRC

			if (m > 0){
				private$cpoisson_for_matched_pairs(estimate_only = estimate_only)
			}
			beta_m = private$cached_values$beta_T_matched
			ssq_m = private$cached_values$ssq_beta_T_matched
			m_ok = !is.null(beta_m) && is.finite(beta_m) && 
			       (!estimate_only && !is.null(ssq_m) && is.finite(ssq_m) && ssq_m > 0 || estimate_only)

			if (nRT > 0 && nRC > 0){
				private$negbin_for_reservoir(estimate_only = estimate_only)
			}
			beta_r = private$cached_values$beta_T_reservoir
			ssq_r = private$cached_values$ssq_beta_T_reservoir
			r_ok = !is.null(beta_r) && is.finite(beta_r) && 
			       (!estimate_only && !is.null(ssq_r) && is.finite(ssq_r) && ssq_r > 0 || estimate_only)

			if (m_ok && r_ok){
				w_star = ssq_r / (ssq_r + ssq_m)
				private$cached_values$beta_hat_T = w_star * beta_m + (1 - w_star) * beta_r
				if (estimate_only) return(invisible(NULL))
				private$cached_values$s_beta_hat_T = sqrt(ssq_m * ssq_r / (ssq_m + ssq_r))
			} else if (m_ok){
				private$cached_values$beta_hat_T = beta_m
				private$cached_values$s_beta_hat_T = if (estimate_only) NA_real_ else sqrt(ssq_m)
			} else if (r_ok){
				private$cached_values$beta_hat_T = beta_r
				private$cached_values$s_beta_hat_T = if (estimate_only) NA_real_ else sqrt(ssq_r)
			} else {
				private$cached_values$beta_hat_T = NA_real_
				private$cached_values$s_beta_hat_T = NA_real_
			}
		},

		cpoisson_for_matched_pairs = function(estimate_only = FALSE){
			KKstats = private$cached_values$KKstats
			yT = KKstats$yTs_matched
			yC = KKstats$yCs_matched

			# Filter pairs where total count is zero (provide no information for the conditional likelihood)
			y_total = yT + yC
			valid_idx = which(y_total > 0)
			if (length(valid_idx) == 0) return(invisible(NULL))

			y_prop = yT[valid_idx] / y_total[valid_idx]
			weights = y_total[valid_idx]

			# If no covariates, use weighted intercept-only logistic
			X = matrix(1, nrow = length(valid_idx), ncol = 1)
			colnames(X) = "(Intercept)"
			if (ncol(as.matrix(private$X)) > 0 && !is.null(private$best_X_colnames_matched)){
				X = cbind(1, KKstats$X_matched_diffs[valid_idx, intersect(private$best_X_colnames_matched, colnames(KKstats$X_matched_diffs)), drop = FALSE])
			}

			mod = tryCatch({
				res = fast_logistic_regression_weighted_cpp(X = X, y = y_prop, weights = weights)
				ssq_b_1 = if (estimate_only) NA_real_ else {
					vcov_mat = tryCatch(solve(res$XtWX), error = function(e) NULL)
					if (!is.null(vcov_mat)) vcov_mat[1L, 1L] else NA_real_
				}
				list(b = res$b, ssq_b_1 = ssq_b_1, X_fit = X)
			}, error = function(e) NULL)

			if (!is.null(mod) && is.finite(mod$b[1])){
				private$cached_values$beta_T_matched = as.numeric(mod$b[1])
				private$best_X_colnames_matched = setdiff(colnames(mod$X_fit), "(Intercept)")
				if (!estimate_only) private$cached_values$ssq_beta_T_matched = as.numeric(mod$ssq_b_1)
			}
		},

		negbin_for_reservoir = function(estimate_only = FALSE){
			y_r    = private$cached_values$KKstats$y_reservoir
			w_r    = private$cached_values$KKstats$w_reservoir
			X_r    = as.matrix(private$cached_values$KKstats$X_reservoir)
			j_treat = 2L

			if (ncol(as.matrix(private$X)) > 0){
				X_full = cbind(1, w_r, X_r)
				attempt = private$fit_with_hardened_qr_column_dropping(
					X_full = X_full,
					required_cols = 2L, # intercept and treatment
					fit_fun = function(X_fit){
						fast_neg_bin_with_var_cpp(X = X_fit, y = as.integer(y_r))
					},
					fit_ok = function(mod, X_fit, keep){
						if (is.null(mod) || !is.finite(mod$b[2L])) return(FALSE)
						if (estimate_only) return(TRUE)
						j_col = which(colnames(X_fit) == "w_r")
						if (length(j_col) == 0L) j_col = 2L
						ssq = mod$vcov[j_col, j_col]
						is.finite(ssq) && ssq > 0
					}
				)
				mod = attempt$fit
				if (!is.null(mod)){
					j_treat = which(colnames(attempt$X) == "w_r")
					if (length(j_treat) == 0) j_treat = 2L
				}
			} else {
				X = cbind(1, w_r)
				mod = tryCatch(fast_neg_bin_with_var_cpp(X = X, y = as.integer(y_r)), error = function(e) NULL)
			}

			if (!is.null(mod) && is.finite(mod$b[j_treat])){
				private$cached_values$beta_T_reservoir = as.numeric(mod$b[j_treat])
				if (ncol(as.matrix(private$X)) > 0 && !is.null(attempt$fit)){
					private$best_X_colnames_reservoir = setdiff(colnames(attempt$X), c("(Intercept)", "w_r"))
					private$best_X_j_treat_reservoir = j_treat
				}
			}
			if (!estimate_only && !is.null(mod)) private$cached_values$ssq_beta_T_reservoir = as.numeric(mod$vcov[j_treat, j_treat])
		},

		best_X_colnames_matched = NULL,
		best_X_colnames_reservoir = NULL,
		best_X_j_treat_reservoir = 2L
	)
)
