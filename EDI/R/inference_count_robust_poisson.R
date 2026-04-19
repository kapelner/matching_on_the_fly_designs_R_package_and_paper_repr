#' Robust Poisson Regression Inference for Count Responses
#'
#' Fits a Poisson log-link regression for count responses using the treatment
#' indicator and, optionally, all recorded covariates as predictors. Inference
#' uses a Huber-White sandwich variance.
#'
#' @export
InferenceCountRobustPoisson = R6::R6Class("InferenceCountRobustPoisson",
	lock_objects = FALSE,
	inherit = InferenceMLEorKMforGLMs,
	public = list(

		#' @description
		#' Initialize a robust Poisson regression inference object.
		#' @param des_obj A completed \code{Design} object with a count response.
		#' @param include_covariates Logical. If \code{TRUE}, all covariates in the design
		#'   are included as predictors. If \code{FALSE}, only the treatment indicator
		#'   is used. If \code{NULL} (default), it is set to \code{TRUE} if the design
		#'   contains covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, include_covariates = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "count")
				assertFlag(include_covariates, null.ok = TRUE)
			}
			super$initialize(des_obj, verbose)
			if (should_run_asserts()) {
				assertNoCensoring(private$any_censoring)
			}
			
			if (is.null(include_covariates)) {
				include_covariates = des_obj$has_covariates()
			}
			private$include_covariates = include_covariates
		},

		#' @description
		#' Computes the robust Poisson estimate of the treatment effect.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
		}
	),

	private = list(
		include_covariates = NULL,
		best_Xmm_colnames = NULL,

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			if (is.null(private$best_Xmm_colnames)){
				private$shared(estimate_only = TRUE)
			}
			if (is.null(private$best_Xmm_colnames)){
				return(self$compute_treatment_estimate(estimate_only = estimate_only))
			}

			Xmm_cols = private$best_Xmm_colnames
			X_data = private$get_X()
			
			if (length(Xmm_cols) == 0L){
				Xmm = cbind(1, private$w)
			} else {
				X_cov = X_data[, intersect(Xmm_cols, colnames(X_data)), drop = FALSE]
				Xmm = cbind(1, treatment = private$w, X_cov)
			}

			res = fast_poisson_regression_cpp(X = Xmm, y = as.numeric(private$y))
			if (is.null(res) || !is.finite(res$b[2])){
				return(NA_real_)
			}
			as.numeric(res$b[2])
		},

		supports_reusable_bootstrap_worker = function(){
			TRUE
		},

		build_design_matrix = function(){
			if (private$include_covariates) {
				private$create_design_matrix()
			} else {
				X = cbind(1, private$w)
				colnames(X) = c("(Intercept)", "treatment")
				X
			}
		},

		fit_count_model_with_var = function(Xmm, estimate_only = FALSE){
			reduced = private$reduce_design_matrix_preserving_treatment_fixed_covariates(Xmm)
			X_fit = reduced$X
			j_treat = reduced$j_treat
			if (is.null(X_fit) || !is.finite(j_treat) || nrow(X_fit) <= ncol(X_fit)){
				return(list(b = rep(NA_real_, ncol(Xmm)), ssq_b_2 = NA_real_))
			}

			mod = tryCatch(
				fast_poisson_regression_cpp(X = X_fit, y = as.numeric(private$y)),
				error = function(e) NULL
			)
			if (is.null(mod)){
				return(list(b = rep(NA_real_, ncol(Xmm)), ssq_b_2 = NA_real_))
			}

			coef_hat = as.numeric(mod$b)
			if (length(coef_hat) != ncol(X_fit) || any(!is.finite(coef_hat))){
				return(list(b = rep(NA_real_, ncol(Xmm)), ssq_b_2 = NA_real_))
			}

			if (estimate_only){
				b_full = rep(NA_real_, ncol(Xmm))
				b_full[reduced$keep] = coef_hat
				names(b_full) = colnames(Xmm)
				return(list(b = b_full, ssq_b_2 = NA_real_))
			}

			mu_hat = as.numeric(exp(X_fit %*% coef_hat))
			if (length(mu_hat) != nrow(X_fit) || any(!is.finite(mu_hat)) || any(mu_hat <= 0)){
				return(list(b = rep(NA_real_, ncol(Xmm)), ssq_b_2 = NA_real_))
			}

			bread = NULL
			var_keep = seq_len(ncol(X_fit))
			X_var = X_fit
			cross_mat = tryCatch(
				crossprod(X_fit, X_fit * mu_hat),
				error = function(e) NULL
			)
			if (!is.null(cross_mat)) {
				bread = tryCatch(solve(cross_mat), error = function(e) NULL)
			}
			if (is.null(bread)) {
				sqrt_mu = sqrt(mu_hat)
				X_weighted = X_fit * sqrt_mu
				reduced_weighted = private$reduce_design_matrix_preserving_treatment(X_weighted)
				keep_sub = as.integer(reduced_weighted$keep)
				if (length(keep_sub) == 0L) {
					return(list(b = rep(NA_real_, ncol(Xmm)), ssq_b_2 = NA_real_))
				}
				X_var_candidate = as.matrix(X_fit[, keep_sub, drop = FALSE])
				cross_mat = tryCatch(
					crossprod(X_var_candidate, X_var_candidate * mu_hat),
					error = function(e) NULL
				)
				if (!is.null(cross_mat)) {
					bread = tryCatch(solve(cross_mat), error = function(e) NULL)
					if (!is.null(bread)) {
						var_keep = keep_sub
						X_var = X_var_candidate
					}
				}
			}
			if (is.null(bread)){
				return(list(b = rep(NA_real_, ncol(Xmm)), ssq_b_2 = NA_real_))
			}

			resid = as.numeric(private$y) - mu_hat
			meat = crossprod(X_var, X_var * (resid^2))
			vcov_trim = bread %*% meat %*% bread
			vcov_full = matrix(NA_real_, ncol(X_fit), ncol(X_fit))
			if (length(var_keep) == ncol(X_fit)) {
				vcov_full = vcov_trim
			} else {
				for (i in seq_along(var_keep)){
					for (j in seq_along(var_keep)){
						vcov_full[var_keep[i], var_keep[j]] = vcov_trim[i, j]
					}
				}
			}
			vcov_robust = (vcov_full + t(vcov_full)) / 2

			ssq_b_j = as.numeric(vcov_robust[j_treat, j_treat])
			if (!is.finite(ssq_b_j) || ssq_b_j < 0){
				ssq_b_j = NA_real_
			}

			b_full = rep(NA_real_, ncol(Xmm))
			b_full[reduced$keep] = coef_hat
			names(b_full) = colnames(Xmm)
			list(b = b_full, ssq_b_2 = ssq_b_j)
		},

		generate_mod = function(estimate_only = FALSE){
			private$fit_count_model_with_var(private$build_design_matrix(), estimate_only = estimate_only)
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			has_cached_se = !is.null(private$cached_values$s_beta_hat_T) &&
				length(private$cached_values$s_beta_hat_T) == 1L &&
				is.finite(private$cached_values$s_beta_hat_T)
			if (!is.null(private$cached_values$is_z) && (estimate_only || has_cached_se)) return(invisible(NULL))
			model_output = private$generate_mod(estimate_only = estimate_only)
			private$cached_values$beta_hat_T = model_output$b[2]
			if (estimate_only) return(invisible(NULL))

			ssq = model_output$ssq_b_2
			if (!is.null(ssq) && !is.na(ssq) && ssq > 0) {
				private$cached_values$s_beta_hat_T = sqrt(ssq)
			} else {
				private$cached_values$s_beta_hat_T = NA_real_
			}
			private$cached_values$is_z = TRUE
			private$cached_values$df = Inf
		}
	)
)
