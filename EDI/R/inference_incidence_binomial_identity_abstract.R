#' Binomial Identity-Link Inference for Binary Responses
#'
#' Internal base class for binomial identity-link regression inference on
#' incidence outcomes. The treatment effect is reported on the risk-difference
#' scale.
#'
#' @keywords internal
#' @noRd
InferenceIncidBinomialIdentityAbstract = R6::R6Class("InferenceIncidBinomialIdentityAbstract",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(
		#' @description
		#' Returns the treatment effect estimate (risk difference).
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes an asymptotic confidence interval for the risk difference.
		#' @param alpha The confidence level is 1 - alpha. Default 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes an asymptotic two-sided p-value for the treatment effect.
		#' @param delta Null treatment effect. Default 0.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),
	private = list(
		best_Xmm_colnames = NULL,

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			# Ensure we have the best design from the original data
			if (is.null(private$best_Xmm_colnames)){
				private$shared(estimate_only = TRUE)
			}
			# Fallback if initial fit failed
			if (is.null(private$best_Xmm_colnames)){
				return(self$compute_treatment_estimate(estimate_only = estimate_only))
			}

			# Use the same design matrix structure as the original fit
			Xmm_cols = private$best_Xmm_colnames
			X_data = private$get_X()
			# Filter columns that are in Xmm_cols (excluding treatment which we provide separately)
			X_cov = X_data[, intersect(Xmm_cols, colnames(X_data)), drop = FALSE]
			Xmm = cbind(treatment = private$w, X_cov)
			# Add intercept
			Xmm = cbind("(Intercept)" = 1, Xmm)

			mod = tryCatch(private$fit_constrained_binomial_maybe_fast(Xmm, j_treat = 2L, estimate_only = estimate_only), error = function(e) NULL)
			if (is.null(mod) || !is.finite(mod$b[2])){
				return(NA_real_)
			}
			as.numeric(mod$b[2])
		},

		supports_reusable_bootstrap_worker = function(){
			TRUE
		},

		create_bootstrap_worker_state = function(){
			private$create_design_backed_bootstrap_worker_state()
		},

		load_bootstrap_sample_into_worker = function(worker_state, indices){
			private$load_bootstrap_sample_into_design_backed_worker(worker_state, indices)
		},

		compute_bootstrap_worker_estimate = function(worker_state){
			private$compute_bootstrap_worker_estimate_via_compute_treatment_estimate(worker_state)
		},

		build_design_matrix = function() stop(class(self)[1], " must implement build_design_matrix()."),

		fit_constrained_binomial_maybe_fast = function(X_fit, j_treat, estimate_only = FALSE){
			if (estimate_only) {
				fast_identity_binomial_regression_cpp(X_fit, as.numeric(private$y))
			} else {
				fast_identity_binomial_regression_with_var_cpp(X_fit, as.numeric(private$y), j = j_treat)
			}
		},

		method_label = function() "Binomial identity-link regression",

		get_covariate_names = function(){
			X = private$get_X()
			p = ncol(X)
			x_names = colnames(X)
			if (is.null(x_names)) x_names = paste0("x", seq_len(p))
			x_names
		},

		set_failed_fit_cache = function(){
			private$cache_nonestimable_estimate("binomial_identity_fit_unavailable")
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				return(invisible(NULL))
			}
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			X_full = private$build_design_matrix()
			if (is.null(dim(X_full))){
				X_full = matrix(X_full, ncol = 2L)
			}
			if (is.null(colnames(X_full))) {
				colnames(X_full) = c("(Intercept)", "treatment", if (ncol(X_full) > 2L) private$get_covariate_names() else NULL)
			}
			fit_candidate = function(X_cand){
				reduced = private$reduce_design_matrix_preserving_treatment_fixed_covariates(X_cand)
				X_fit = reduced$X
				j_treat = reduced$j_treat
				if (is.null(X_fit) || !is.finite(j_treat) || nrow(X_fit) <= ncol(X_fit)) return(NULL)
				mod = tryCatch(private$fit_constrained_binomial_maybe_fast(X_fit, j_treat, estimate_only = estimate_only), error = function(e) NULL)
				if (is.null(mod) || !isTRUE(mod$converged) || (!estimate_only && !is.finite(mod$ssq_b_j))) return(NULL)
				list(mod = mod, j_treat = j_treat, X_fit = X_fit)
			}

			result = fit_candidate(X_full)

			if (private$harden && is.null(result) && ncol(X_full) > 2L){
				X_cov_orig = as.matrix(X_full[, -(1:2), drop = FALSE])
				for (thresh in c(0.95, 0.90, 0.80, 0.70)){
					X_cov_red = drop_highly_correlated_cols(X_cov_orig, threshold = thresh)$M
					if (ncol(X_cov_red) == 0L) next
					X_cand = cbind(X_full[, 1:2, drop = FALSE], X_cov_red)
					private$fixed_covariate_keep_cache = NULL
					private$reduced_design_keep_cache = NULL
					result = fit_candidate(X_cand)
					if (!is.null(result)) break
				}
				if (is.null(result)){
					private$fixed_covariate_keep_cache = NULL
					private$reduced_design_keep_cache = NULL
					result = fit_candidate(X_full[, 1:2, drop = FALSE])
				}
			}

			if (is.null(result)){
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}

			private$best_Xmm_colnames = setdiff(colnames(result$X_fit), c("(Intercept)", "treatment"))
			private$cached_values$beta_hat_T = as.numeric(result$mod$b[result$j_treat])
			if (estimate_only) return(invisible(NULL))
			private$cached_values$s_beta_hat_T = sqrt(as.numeric(result$mod$ssq_b_j))
			private$cached_values$is_z = TRUE
			private$cached_values$df = nrow(result$X_fit) - ncol(result$X_fit)
		}
	)
)
