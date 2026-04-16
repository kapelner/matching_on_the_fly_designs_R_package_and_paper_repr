#' Abstract class for all-subject modified-Poisson inference in KK designs
#'
#' @keywords internal
InferenceAbstractKKModifiedPoisson = R6::R6Class("InferenceAbstractKKModifiedPoisson",
	lock_objects = FALSE,
	inherit = InferenceAbstractKKMarginalIncid,
	public = list(

		#' @description
		#' Compute treatment estimate
		#' @param estimate_only If TRUE, skip variance component calculations.
		compute_treatment_estimate = function(estimate_only = FALSE){
			private$shared(estimate_only = TRUE)
			private$cached_values$beta_hat_T
		},

		#' @description
		#' Compute asymp confidence interval
		#' @param alpha Description for alpha
		compute_asymp_confidence_interval = function(alpha = 0.05){
			assertNumeric(alpha, lower = .Machine$double.xmin, upper = 1 - .Machine$double.xmin)
			private$shared(estimate_only = FALSE)
			private$assert_finite_se()
			private$compute_z_or_t_ci_from_s_and_df(alpha)
		},

		#' @description
		#' Compute asymp two sided pval for treatment effect
		#' @param delta Description for delta
		compute_asymp_two_sided_pval_for_treatment_effect = function(delta = 0){
			assertNumeric(delta)
			private$shared(estimate_only = FALSE)
			private$assert_finite_se()
			private$compute_z_or_t_two_sided_pval_from_s_and_df(delta)
		}
	),

	private = list(
		max_abs_reasonable_coef = 1e4,

		build_design_matrix = function() stop(class(self)[1], " must implement build_design_matrix()."),

		assert_finite_se = function(){
			if (!is.finite(private$cached_values$s_beta_hat_T) || private$cached_values$s_beta_hat_T <= 0){
				return(invisible(NULL))
			}
		},

		set_failed_fit_cache = function(){
			private$cache_nonestimable_estimate("kk_modified_poisson_fit_unavailable")
			private$cached_values$full_coefficients = NULL
			private$cached_values$full_vcov = NULL
			private$cached_values$summary_table = NULL
		},

		coefficients_are_usable = function(coef_hat){
			length(coef_hat) > 0L &&
				all(is.finite(coef_hat)) &&
				max(abs(coef_hat), na.rm = TRUE) <= private$max_abs_reasonable_coef
		},

		fit_modified_poisson = function(X_fit, j_treat, estimate_only = FALSE){
			mod = tryCatch(
				fast_poisson_regression_cpp(X = X_fit, y = as.numeric(private$y)),
				error = function(e) NULL
			)
			if (is.null(mod)){
				return(NULL)
			}

			coef_hat = as.numeric(mod$b)
			if (length(coef_hat) != ncol(X_fit) || !private$coefficients_are_usable(coef_hat)){
				return(NULL)
			}

			if (estimate_only){
				return(list(
					beta_hat = coef_hat[j_treat],
					se = NA_real_,
					coefficients = coef_hat,
					vcov = NULL,
					summary_table = NULL
				))
			}

			mu_hat = as.numeric(exp(X_fit %*% coef_hat))
			if (length(mu_hat) != nrow(X_fit) || any(!is.finite(mu_hat)) || any(mu_hat <= 0)){
				return(NULL)
			}

			post_fit = tryCatch(
				glm_cluster_sandwich_post_fit_cpp(
					X_fit = X_fit,
					y = as.numeric(private$y),
					coef_hat = coef_hat,
					mu_hat = mu_hat,
					working_weights = mu_hat,
					cluster_id = private$get_cluster_ids(),
					j_treat = j_treat
				),
				error = function(e) NULL
			)
			if (is.null(post_fit)){
				return(NULL)
			}

			coef_names = colnames(X_fit)
			names(coef_hat) = coef_names
			vcov_robust = post_fit$vcov
			colnames(vcov_robust) = rownames(vcov_robust) = coef_names
			std_err = post_fit$std_err
			names(std_err) = coef_names
			z_vals = post_fit$z_vals
			names(z_vals) = coef_names

			list(
				beta_hat = post_fit$beta_hat,
				se = post_fit$se,
				coefficients = coef_hat,
				vcov = vcov_robust,
				summary_table = cbind(
					Value = coef_hat,
					`Std. Error` = std_err,
					`z value` = z_vals,
					`Pr(>|z|)` = 2 * stats::pnorm(-abs(z_vals))
				)
			)
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
				colnames(X_full) = c("(Intercept)", "treatment", if (ncol(X_full) > 2L) private$get_covariate_names() else NULL)
			}

			reduced = private$reduce_design_matrix_preserving_treatment(X_full)
			X_fit = reduced$X
			j_treat = reduced$j_treat
			if (is.null(X_fit) || !is.finite(j_treat) || nrow(X_fit) <= ncol(X_fit)){
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}

			fit = private$fit_modified_poisson(X_fit, j_treat, estimate_only = estimate_only)
			if (private$harden && is.null(fit) && ncol(X_full) > 2L){
				reduced = private$reduce_design_matrix_preserving_treatment(X_full[, 1:2, drop = FALSE])
				X_fit = reduced$X
				j_treat = reduced$j_treat
				if (!is.null(X_fit) && is.finite(j_treat) && nrow(X_fit) > ncol(X_fit)){
					fit = private$fit_modified_poisson(X_fit, j_treat, estimate_only = estimate_only)
				}
			}
			if (is.null(fit)){
				private$set_failed_fit_cache()
				return(invisible(NULL))
			}

			private$cached_values$beta_hat_T = fit$beta_hat
			if (estimate_only) return(invisible(NULL))
			private$cached_values$s_beta_hat_T = fit$se
			private$cached_values$is_z = TRUE
			private$cached_values$df = nrow(X_fit) - ncol(X_fit)
			private$cached_values$full_coefficients = fit$coefficients
			private$cached_values$full_vcov = fit$vcov
			private$cached_values$summary_table = fit$summary_table
		}
	)
)
