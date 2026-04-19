#' Fractional Logit Inference for Proportion Responses
#'
#' Fits a fractional logit model for proportion responses using a binomial logit
#' quasi-likelihood with sandwich-robust variance. The treatment indicator and,
#' optionally, all recorded covariates are used as predictors. The treatment
#' effect is reported on the log-odds scale.
#'
#' @export
InferencePropFractionalLogit = R6::R6Class("InferencePropFractionalLogit",
	lock_objects = FALSE,
	inherit = InferenceAsymp,
	public = list(

		#' @description
		#' Initialize a fractional-logit inference object.
		#' @param des_obj A completed \code{Design} object with a proportion response.
		#' @param include_covariates Logical. If \code{TRUE}, all covariates in the design
		#'   are included as predictors. If \code{FALSE}, only the treatment indicator
		#'   is used. If \code{NULL} (default), it is set to \code{TRUE} if the design
		#'   contains covariates.
		#' @param verbose Whether to print progress messages.
		initialize = function(des_obj, include_covariates = NULL, verbose = FALSE){
			if (should_run_asserts()) {
				assertResponseType(des_obj$get_response_type(), "proportion")
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
		#' Computes the fractional-logit estimate of the treatment effect on the
		#' log-odds scale.
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = TRUE)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Computes a 1 - \code{alpha} confidence interval using the sandwich
		#' standard error.
		#' @param alpha The confidence level. The default is 0.05.
		compute_asymp_confidence_interval = function(alpha = 0.05){
			if (should_run_asserts()) {
				assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			}
			private$shared(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Computes a two-sided p-value for the treatment effect.
		#' @param delta The null treatment effect on the log-odds scale.
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			if (should_run_asserts()) {
				assertNumeric(delta)
			}
			private$shared(estimate_only = FALSE)
			if (should_run_asserts()) {
				private$assert_finite_se()
			}
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		include_covariates = NULL,
		best_Xmm_colnames = NULL,

		compute_treatment_estimate_during_randomization_inference = function(estimate_only = TRUE){
			private$shared(estimate_only = estimate_only)
			private$cached_values$beta_hat_T
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

		get_covariate_names = function(){
			X = private$get_X()
			p = ncol(X)
			x_names = colnames(X)
			if (is.null(x_names)){
				x_names = paste0("x", seq_len(p))
			}
			x_names
		},

		select_covariate_to_drop = function(X_curr, coef_hat){
			covariate_cols = seq.int(3L, ncol(X_curr))
			if (length(covariate_cols) == 0L) return(NA_integer_)
			coef_mags = abs(coef_hat[covariate_cols])
			if (length(coef_mags) == 0L || all(!is.finite(coef_mags))){
				return(tail(covariate_cols, 1L))
			}
			covariate_cols[which.max(replace(coef_mags, !is.finite(coef_mags), -Inf))]
		},

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				return(invisible(NULL))
			}
		},

		set_failed_fit_cache = function(){
			private$cache_nonestimable_estimate("fractional_logit_fit_unavailable")
			private$cached_values$full_coefficients = NULL
			private$cached_values$full_vcov = NULL
			private$cached_values$summary_table = NULL
		},

		fit_fractional_logit = function(X_full, estimate_only = FALSE){
			X_curr = X_full

			repeat {
				reduced = private$reduce_design_matrix_preserving_treatment(X_curr)
				X_fit = reduced$X
				j_treat = reduced$j_treat
				if (is.null(X_fit) || !is.finite(j_treat) || nrow(X_fit) <= ncol(X_fit)){
					return(NULL)
				}

				mod = tryCatch(
					fast_logistic_regression_cpp(X = X_fit, y = as.numeric(private$y)),
					error = function(e) NULL
				)
				if (is.null(mod)){
					return(NULL)
				}

				coef_hat = as.numeric(mod$b)
				converged = all(is.finite(coef_hat))
				if (!converged){
					if (ncol(X_curr) <= 2L) return(NULL)
					drop_col = private$select_covariate_to_drop(X_curr, coef_hat)
					if (!is.finite(drop_col)) return(NULL)
					X_curr = X_curr[, -drop_col, drop = FALSE]
					next
				}

				if (estimate_only){
					return(list(
						beta_hat = coef_hat[j_treat],
						se = NA_real_,
						coefficients = coef_hat,
						vcov = NULL,
						summary_table = NULL,
						X_fit = X_fit
					))
				}

				mu_hat = inv_logit(X_fit %*% coef_hat)
				mu_hat = pmin(pmax(as.numeric(mu_hat), .Machine$double.eps), 1 - .Machine$double.eps)
				W = mu_hat * (1 - mu_hat)
				if (any(!is.finite(W)) || any(W <= 0)){
					if (ncol(X_curr) <= 2L) return(NULL)
					drop_col = private$select_covariate_to_drop(X_curr, coef_hat)
					if (!is.finite(drop_col)) return(NULL)
					X_curr = X_curr[, -drop_col, drop = FALSE]
					next
				}

				post_fit = tryCatch(
					glm_sandwich_post_fit_cpp(
						X_fit = X_fit,
						y = as.numeric(private$y),
						coef_hat = coef_hat,
						mu_hat = mu_hat,
						working_weights = W,
						j_treat = j_treat
					),
					error = function(e) NULL
				)
				if (is.null(post_fit)){
					if (ncol(X_curr) <= 2L) return(NULL)
					drop_col = private$select_covariate_to_drop(X_curr, coef_hat)
					if (!is.finite(drop_col)) return(NULL)
					X_curr = X_curr[, -drop_col, drop = FALSE]
					next
				}

				coef_names = colnames(X_fit)
				names(coef_hat) = coef_names
				vcov_robust = post_fit$vcov
				colnames(vcov_robust) = rownames(vcov_robust) = coef_names
				std_err = post_fit$std_err
				names(std_err) = coef_names
				z_vals = post_fit$z_vals
				names(z_vals) = coef_names

				return(list(
					beta_hat = post_fit$beta_hat,
					se = post_fit$se,
					coefficients = coef_hat,
					vcov = vcov_robust,
					summary_table = cbind(
						Value = coef_hat,
						`Std. Error` = std_err,
						`z value` = z_vals,
						`Pr(>|z|)` = 2 * stats::pnorm(-abs(z_vals))
					),
					X_fit = X_fit
				))
			}
		},

		shared = function(estimate_only = FALSE){
			if (estimate_only && !is.null(private$cached_values$beta_hat_T)) return(invisible(NULL))
			if (!estimate_only && !is.null(private$cached_values$s_beta_hat_T)) return(invisible(NULL))

			if (!is.null(private$cached_values$beta_hat_T) && (estimate_only || !is.null(private$cached_values$summary_table))) return(invisible(NULL))

			X_full = private$build_design_matrix()
			if (is.null(dim(X_full))){
				X_full = matrix(X_full, ncol = 2L)
			}
			if (is.null(colnames(X_full))) {
				full_names = c("(Intercept)", "treatment", if (ncol(X_full) > 2L) private$get_covariate_names() else NULL)
				colnames(X_full) = full_names[seq_len(ncol(X_full))]
			}

			fit = private$fit_fractional_logit(X_full, estimate_only = estimate_only)
			if (is.null(fit)){
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T = fit$beta_hat
			private$best_Xmm_colnames = setdiff(colnames(fit$X_fit), c("(Intercept)", "treatment"))
			
			if (estimate_only) return(invisible(NULL))
			private$cached_values$s_beta_hat_T = fit$se
			private$cached_values$is_z = TRUE
			private$cached_values$df = nrow(X_full) - (if (is.null(fit$vcov)) 2L else ncol(fit$vcov))
			private$cached_values$full_coefficients = fit$coefficients
			private$cached_values$full_vcov = fit$vcov
			private$cached_values$summary_table = fit$summary_table
		}
	)
)
